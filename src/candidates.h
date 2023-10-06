#ifndef STOCKFISH_CANDIDATES_H
#define STOCKFISH_CANDIDATES_H

#include "main.h"
#include "matrix.h"
#include "poly.h"
#include "time.h"

inline int gray_to_lex(int gray) {
    int result = gray ^ (gray >> 8);
    result ^= (result >> 4);
    result ^= (result >> 2);
    result ^= (result >> 1);
    return result;
}

inline int lex_to_gray(int lex) {
    return (lex >> 1) ^ lex;
}

template<typename Params>
struct Candidates {
    static Params params;
    MatrixGF2 candidates;
    // We keep track of the candidates (alpha, g(alpha)) where g(alpha) != 0.
    int nr_candidates = params.q * (params.q - 1);

    std::vector<bool> deleted;
    std::vector<typename decltype(params.gf2m)::Element> _idx2a;
    std::vector<typename decltype(params.gf2m)::Element> _idx2ga;

    // static int aga2idx(int a, int ga) { return (ga - 1) * params.q + a; }
    typename decltype(params.gf2m)::Element idx2a(int idx) { return _idx2a[idx]; }

    typename decltype(params.gf2m)::Element idx2ga(int idx) { return _idx2ga[idx]; }

    void generate_candidate_matrix() {
        auto loop = [this](int j_start, int j_end) {
            for (int j = j_start; j < j_end; j++) {
                MatrixGF2::pack_t *packptr = &candidates.packs[j * candidates.packs_per_column];
                if (j % 100000 == 0)
                    std::cout << j << "/" << nr_candidates << "\r" << std::flush;
                typename decltype(params.gf2m)::Element el = idx2ga(j).inv();
                for (int i = 0; i < params.m * params.t; i += params.m) {
                    MatrixGF2::pack_t val = el.val;

                    int value_shift = i % candidates.rows_per_pack;
                    MatrixGF2::pack_t value_mask = ((1ULL << params.m) - 1) << value_shift;

                    // write val into pack.
                    *packptr &= ~value_mask;
                    *packptr |= val << value_shift;

                    if (value_shift + params.m == candidates.rows_per_pack) {
                        packptr++;
                    } else if (value_shift + params.m > candidates.rows_per_pack) {
                        // We split into two.

                        // Overflow into the next pack.
                        int overflow_bits = value_shift + params.m - candidates.rows_per_pack;
                        MatrixGF2::pack_t overflow_mask = (1ULL << overflow_bits) - 1;
                        packptr++;
                        *packptr &= ~overflow_mask;
                        *packptr |= val >> (params.m - overflow_bits);
                    }
                    el = el * idx2a(j);
                }
            }
        };
        pool.parallelize_loop(0, nr_candidates, loop).wait();
    }

    Candidates(bool botan=false) {
        resize();
        for (int j = 0; j < nr_candidates; j++) {
            deleted[j] = false;
            _idx2a[j] = j % params.q;
            _idx2ga[j] = (j / params.q) + 1;
        }

        std::string fname = std::string("cache-") + params.name + ".bin";
        if (std::filesystem::exists(fname)) {
            spdlog::info("Loading precomputed candidate matrix.");
            load(fname);
        } else {
            spdlog::info("Regenerating candidate matrix, this may take a while...");
            generate_candidate_matrix();

            if (botan) {
                // With Botan, we can actually filter the initial code book, as all possible values for alpha are known.
                if (params.n < params.q) {
                    spdlog::info("Filtering initial codebook for alpha candidates");
                    int nr_filtered = 0;
                    for (int i = params.n; i < params.q; i++) {
                        int val = lex_to_gray(i);
                        nr_filtered += filter_candidates(params.gf2m.from_int(val));
                    }
                    defragmentation();
                    spdlog::info("Removed {0} candidates based on initial alpha candidates", nr_filtered);
                }
            }
            spdlog::info("Saving candidate matrix for later, next time it should be quick.");
            save(fname);
        }
    };

    int parallel_blocks() {
        if (nr_candidates <= params.q)
            return std::min(params.m,
                            (int) pool.get_thread_count()); // any small number will do fine here, so we scale a bit with problem size.
        else
            return 0; // default, nr of threads in pool
    }

    // This adds all rows in row_mask to row)
    void add_rows_to_row(int r, MatrixGF2 &row_mask) {
        // Code is thread-safe because it only operates on individual columns.
        auto loop = [this, r, row_mask](int a, int b) {
            for (int j = a; j < b; j++) {
                if (deleted[j])
                    continue;
                int pack_idx = j * candidates.packs_per_column;

                // First we determine the weight (modulo 2) of the inner product of the column with row_mask
                int weight = 0;
                for (int i = 0; i < params.m * params.t; i += MatrixGF2::rows_per_pack) {
                    MatrixGF2::pack_t pack = candidates.packs[pack_idx + i / MatrixGF2::rows_per_pack];
                    MatrixGF2::pack_t mask = row_mask.packs[i / MatrixGF2::rows_per_pack];
                    weight ^= __builtin_popcountll(pack & mask);
                }
                // set_entry takes care of only looking at the lowest bit.
                candidates.set_entry(r, j, candidates.get_entry(r, j) + weight);
            }
        };
        pool.parallelize_loop(0, nr_candidates, loop, parallel_blocks()).wait();
    }

    void add_row_to_row(int dest_r, int source_r) {
        auto loop = [this, dest_r, source_r](int j_start, int j_end) {
            for (int j = j_start; j < j_end; j++) {
                if (deleted[j])
                    continue;
                int old_val = candidates.get_entry(dest_r, j);
                int other_val = candidates.get_entry(source_r, j);
                candidates.set_entry(dest_r, j, old_val ^ other_val);
            }
        };
        pool.parallelize_loop(0, nr_candidates, loop, parallel_blocks()).wait();
    }

    // This adds row to all rows in row_mask.
    void add_row_to_rows(MatrixGF2 &row_mask, int row) {
        auto loop = [this, row_mask, row](int j_start, int j_end) {
            for (int j = j_start; j < j_end; j++) {
                if (deleted[j])
                    continue;

                if (candidates.get_entry(row, j)) {
                    candidates.add_column(j, row_mask);
                }
            }
        };
        pool.parallelize_loop(0, nr_candidates, loop, parallel_blocks()).wait();
    }

    struct FindResult {
        int nr_candidates;
        uint64_t elapsed;
        int best_distance;
        int best_idx;
        typename decltype(params.gf2m)::Element best_a;
        typename decltype(params.gf2m)::Element best_ga;

        FindResult() {
            best_distance = INT_MAX;
            best_idx = 0;
            best_a = 0;
            best_ga = 0;
        }

        std::string to_string() {
            std::stringstream os;
            if (best_distance == INT_MAX) {
                os << "[FindResult: no result]";
            } else {
                os << "[FindResult: alpha=" << best_a << " and g(alpha)=" << best_ga << " with distance "
                   << best_distance << " at index " << best_idx << "]";
            }
            return os.str();
        }

        friend std::ostream &operator<<(std::ostream &os, const FindResult &r) {
            return os << r.to_string();
        }
    };

    FindResult find_target(const MatrixGF2 &target) {
        auto begin_find_target = std::chrono::high_resolution_clock::now();
        FindResult best_result;

        auto loop = [this, target](const int a, const int b) {
            int nr_candidates_not_deleted = 0;
            FindResult result;
            for (int j = a; j < b; j++) {
                if (deleted[j])
                    continue;
                nr_candidates_not_deleted++;
                int distance = 0;

                // naive:
                // for (int i = 0; i < params.m * params.t; i++)
                //     distance = distance + (candidates.get_entry(j, i) ^ target.get_entry(i, 0));
                for (int i = 0; i < params.m * params.t; i += candidates.rows_per_pack) {
                    int pack_idx = j * candidates.packs_per_column + i / candidates.rows_per_pack;
                    MatrixGF2::pack_t pack = candidates.packs[pack_idx];
                    MatrixGF2::pack_t sum = pack ^ target.packs[i / candidates.rows_per_pack];
                    distance = distance + __builtin_popcountll(sum);
                }
                if (distance < result.best_distance) {
                    result.best_distance = distance;
                    result.best_idx = j;
                    result.best_a = idx2a(j);
                    result.best_ga = idx2ga(j);
                }
            }
            result.nr_candidates = nr_candidates_not_deleted;
            return result;
        };
        BS::multi_future<FindResult> mf = pool.parallelize_loop(0, nr_candidates, loop, parallel_blocks());
        std::vector<FindResult> results = mf.get();
        int nr_candidates_not_deleted = 0;
        for (auto &result: results) {
            nr_candidates_not_deleted += result.nr_candidates;
            if (result.best_distance < best_result.best_distance)
                best_result = result;
        }
        auto end_find_target = std::chrono::high_resolution_clock::now();
        best_result.nr_candidates = nr_candidates_not_deleted;
        best_result.elapsed = get_duration(begin_find_target, end_find_target);
        return best_result;
    }

    int filter_candidates(typename decltype(params.gf2m)::Element a) {
        auto loop = [this, a](int j_start, int j_end) {
            int count = 0;
            for (int j = j_start; j < j_end; j++) {
                if (deleted[j])
                    continue;
                // Remove all columns where a has the given value.
                if (idx2a(j) == a) {
                    deleted[j] = true;
                    count++;
                }
            }
            return count;
        };

        int total = 0;
        BS::multi_future<int> mf = pool.parallelize_loop(0, nr_candidates, loop, parallel_blocks());
        std::vector<int> counts = mf.get();
        for (auto &count: counts) {
            total += count;
        }
        return total;
    }

    int filter_candidates(Poly<decltype(params.gf2m)> g) {
        typename decltype(params.gf2m)::Element g_eval[params.q];
        for (int a = 0; a < params.q; a++) {
            g_eval[a] = g.evaluate(params.gf2m.from_int(a));
        }

        auto loop = [this, g_eval](int j_start, int j_end) {
            int count = 0;
            for (int j = j_start; j < j_end; j++) {
                if (deleted[j])
                    continue;
                // Remove all columns where a has the given value.
                auto a = idx2a(j);
                auto ga = idx2ga(j);
                // Remove all columns where a and ga do not fit the polynomial.
                if (g_eval[a.val] != ga) {
                    deleted[j] = true;
                    count++;
                }
            }
            return count;
        };

        int total = 0;
        BS::multi_future<int> mf = pool.parallelize_loop(0, nr_candidates, loop, parallel_blocks());
        std::vector<int> counts = mf.get();
        for (auto &count: counts) {
            total += count;
        }
        return total;
    }

    // not embarrasingly parallelizable
    int defragmentation() {
        int nr_removed = 0;
        int dst = 0;
        for (int src = 0; src < nr_candidates; src++) {
            if (deleted[src]) {
                nr_removed++;
                continue;
            }

            if (src != dst) {
                auto src_pack = candidates.packs.begin() + src * candidates.packs_per_column;
                auto dst_pack = candidates.packs.begin() + dst * candidates.packs_per_column;
                // Copy src to dst.
                std::copy(src_pack, src_pack + candidates.packs_per_column, dst_pack);
                deleted[dst] = false;
                _idx2a[dst] = _idx2a[src];
                _idx2ga[dst] = _idx2ga[src];
            }
            dst++;
        }
        nr_candidates -= nr_removed;
        resize();
        return nr_removed;
    }

    void resize() {
        candidates.resize(params.m * params.t, nr_candidates);
        deleted.resize(nr_candidates);
        _idx2a.resize(nr_candidates);
        _idx2ga.resize(nr_candidates);
    }

    int get_distance(int j1, int j2) {
        int distance = 0;
        MatrixGF2::pack_t *ipacks = &candidates.packs[j1 * candidates.packs_per_column];
        MatrixGF2::pack_t *jpacks = &candidates.packs[j2 * candidates.packs_per_column];
        for (int k = 0; k < candidates.packs_per_column; k++)
            distance = distance + __builtin_popcountll(ipacks[k] ^ jpacks[k]);
        return distance;
    }

    std::vector<uint64_t> get_distances_report() {
        auto loop = [this](int j_start, int j_end) {
            std::vector<uint64_t> distances(params.t * params.m);

            for (int j1 = j_start; j1 < j_end; j1++) {
                // We use the trick of Gauss as a pupil to get a constant workload in each thread.
                for (int j2 = j1 + 1; j2 < nr_candidates; j2++) {
                    int distance = get_distance(j1, j2);
                    distances[distance]++;
                }
                int j3 = nr_candidates - 1 - j1;

                // We have to avoid double counting around the center for odd nr_candidates.
                if (j3 == j1)
                    continue;

                for (int j2 = j3 + 1; j2 < nr_candidates; j2++) {
                    int distance = get_distance(j3, j2);
                    distances[distance]++;
                }
            }
            return distances;
        };

        // We remove all duplicates here, so we don't need to check deleted flags.
        defragmentation();

        std::vector<uint64_t> total_distances(params.t * params.m);
        auto mf = pool.parallelize_loop(0, (nr_candidates + 1) / 2, loop);
        auto distances = mf.get();
        for (auto &distance: distances) {
            for (int k = 0; k < params.t * params.m; k++)
                total_distances[k] += distance[k];
        }
        return total_distances;
    }

    void load(const std::string& fname) {
        FILE *infp;
        infp = fopen(fname.c_str(), "rb");
        if (!infp)
            throw std::runtime_error("Error opening file: " + fname);
        fread(&nr_candidates, sizeof(nr_candidates), 1, infp);
        resize();
        fread(candidates.packs.data(), sizeof(MatrixGF2::pack_t), candidates.packs.size(), infp);
        fread(_idx2a.data(), sizeof(typename decltype(params.gf2m)::Element), _idx2a.size(), infp);
        fread(_idx2ga.data(), sizeof(typename decltype(params.gf2m)::Element), _idx2ga.size(), infp);
        fclose(infp);
    }

    void save(const std::string& fname) {
        FILE *outfp;
        outfp = fopen(fname.c_str(), "wb");
        fwrite(&nr_candidates, sizeof(nr_candidates), 1, outfp);
        resize();
        fwrite(candidates.packs.data(), sizeof(MatrixGF2::pack_t), candidates.packs.size(), outfp);
        fwrite(_idx2a.data(), sizeof(typename decltype(params.gf2m)::Element), _idx2a.size(), outfp);
        fwrite(_idx2ga.data(), sizeof(typename decltype(params.gf2m)::Element), _idx2ga.size(), outfp);
        fclose(outfp);
    }

};

template<typename Params>
Params Candidates<Params>::params;

#endif //STOCKFISH_CANDIDATES_H
