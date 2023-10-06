#ifndef STOCKFISH_BOTAN_H
#define STOCKFISH_BOTAN_H

#include <iostream>

#include "3rdparty/spdlog/spdlog.h"
#include "3rdparty/nlohmann/json.hpp"

#include "leakfile.h"
#include "randomness.h"
#include "candidates.h"
#include "poly.h"
#include "time.h"


template<typename Params>
struct BotanSecret {
    static Params params;

    std::vector<uint8_t> sk;
    Poly<decltype(params.gf2m)> goppa;
    std::vector<typename decltype(params.gf2m)::Element> alpha;

    BotanSecret(const std::string& aSk) {
        sk.resize(aSk.length()/2);
        for (int i = 0; i < aSk.length()/2; i++) {
            sk[i] = (hexdigit(aSk[2*i]) << 4) + hexdigit(aSk[2*i+1]);
        }

        const int element_len = 2;
        const int goppa_len = (params.t + 1) * element_len;
        // Botan only uses a permutation of (0, 1, ..., n) for alpha.
        const int Linv_len = params.n * element_len;

        for (int i = 0; i < params.t + 1; i++) {
            int idx = element_len * i;
            int coef = (sk[idx] << 8) + sk[idx + 1];
            goppa.set(i, coef);
        }

        const int Linv_offset = goppa_len;
        std::vector<int> Linv(params.n);
        for (int i = 0; i < params.n; i++) {
            int idx = Linv_offset + element_len * i;
            int inv = (sk[idx] << 8) + sk[idx + 1];
            Linv[i] = inv;
        }

        // Botan does Gauss eliminiation from the right. This means that the alphas we are looking for
        // are at position n-tm to n-1. Then the public key follows at position 0 to n-tm-1.
        // But, also, Botan does not use L as values for alpha. Instead, it converts all values in L
        // with lex_to_gray into a gray code number in 2^m. This does not affect the security (as support splitting
        // is still possible), but removes some biases in the bits I suppose.

        alpha.resize(params.n);
        std::vector<int> L(params.n);
        for (int i = 0; i < params.n; i++) {
            L[Linv[i]] = i;
        }
        for (int i = 0; i < params.n; i++) {
            int idx = (i + params.k) % params.n;
            int val = lex_to_gray(L[idx]);
            alpha[i] = params.gf2m.from_int(val);
        }
    }

    bool check_alpha(int column, typename decltype(params.gf2m)::Element alpha_guess) {
        return alpha[column] == alpha_guess;
    }

    bool check_goppa(Poly<decltype(params.gf2m)>& goppa_guess) {
        return goppa == goppa_guess;
    }

};

template<typename Params>
Params BotanSecret<Params>::params;


template<typename Params>
struct Botan {
    static Params params;
    bool verbose;

    Botan() {
    }

    void parse_leak(MatrixGF2& leak, Leakfile leakfile, Randomness& rng) {
        std::cout << "Parsing leak in " << leakfile.filename << std::endl;

        // Leakage model:
        // We leak all bits of the first m*t columns of H, except for the diagonal.
        int leak_size = (params.m * params.t - 1) * params.m * params.t;

        std::string leakdata = leakfile.leak.substr(leakfile.leak.length() - leak_size * 2);
        int leakdatapos = 1; // We only look at the second digit.
        for (int j = 0; j < params.m * params.t; j++) {
            // Botan has two loops, one from j+1 to mt-1, and one from j-1 to 0! Order matters!
            leak.set_entry(j, j, rng.sample_uniform_bit());
            for (int i = j + 1; i < params.m * params.t; i++) {
                char leakbyte = leakdata[leakdatapos];
                leak.set_entry(i, j, leakbyte == '1' ? 1 : 0);
                leakdatapos = leakdatapos + 2;
            }
            for (int i = j - 1; i >= 0; i--) {
                char leakbyte = leakdata[leakdatapos];
                leak.set_entry(i, j, leakbyte == '1' ? 1 : 0);
                leakdatapos = leakdatapos + 2;
            }
        }
    }

    // Public key is stored in transposed form, n-mt rows and mt columns.
    // The byte data is by row, and each row is 32-bit aligned.
    void parse_pk(MatrixGF2& pk, Leakfile leakfile) {
        spdlog::info("Parsing public key in {}", leakfile.filename);

        auto &aPk = leakfile.pk;

        // rowlen in bytes
        int rowlen = ((params.m * params.t + 31) / 32) * 4;
        for (int idx = 0; idx < aPk.size()/2; idx++) {
            auto data = (hexdigit(aPk[2*idx]) << 4) + hexdigit(aPk[2*idx+1]);
            for (int b = 0; b < 8; b++) {
                int val = (data >> b) & 1;
                int j = idx / rowlen;
                int i = (idx % rowlen) * 8 + b;
                if (j < pk.cols)
                    pk.set_entry(i, j, val);
            }
        }
    }

    int add_noise(MatrixGF2& leak, double readout_error, Randomness &rng) {
        if (readout_error == 0)
            return 0;

        int bits_flipped = 0;
        spdlog::info("Adding noise with readout error rate {}", readout_error);

        for (int j = 0; j < params.m * params.t; j++) {
            for (int i = 0; i < params.m * params.t; i++) {
                if (i == j)
                    // In Classic McEliece, the diagonal is already a systematic error (and chosen randomly by parse_leak).
                    continue;

                double draw = rng.sample_uniform_unit();
                // std::cout << "Draw " << draw << " compare with " << readout_error << ": " << (draw < readout_error) << std::endl;
                if (draw < readout_error) {
                    bits_flipped++;
                    leak.set_entry(i, j, 1 - leak.get_entry(i, j));
                }
            }
        }
        return bits_flipped;
    }

    // readout_error rate
    bool recover_from_leak(Leakfile leakfile, unsigned int seed = 0, double readout_error = 0,
                           bool with_report = false,
                           std::string report_fname="result.json") {
        auto begin_recover_from_leak = std::chrono::high_resolution_clock::now();
        Randomness rng(seed);

        // The candidates for columns of H during Gauss elim.
        Candidates<Params> candidates(true);
        std::vector<typename Candidates<Params>::FindResult> results;
        // The known points of the goppa polynomial (until the polynomial is fully known).
        std::vector<std::pair<typename decltype(params.gf2m)::Element, typename decltype(params.gf2m)::Element>> goppa_points;
        // The goppa polynomial once it is known.
        Poly<decltype(params.gf2m)> goppa;

        // We access the secret key only for verification, so we can abort early if we do not hit.
        BotanSecret<Params> sk(leakfile.sk);

        // We need the public key to complete the attack.
        MatrixGF2 pk = MatrixGF2(params.m * params.t, params.k);
        parse_pk(pk, leakfile);

        // Prepare the leak file, adding the desired readout error.
        MatrixGF2 leak = MatrixGF2(params.m * params.t, params.m * params.t);
        parse_leak(leak, leakfile, rng);
        int bits_flipped = add_noise(leak, readout_error, rng);
        // The estimated error is the readout error applied to all bits not on the diagonal, plus half the bits on the diagonal.
        double estimated_error = 0;
        estimated_error += readout_error * (params.m * params.t) * (params.m * params.t - 1);
        estimated_error += ((double)params.m * params.t) / 2;
        estimated_error = estimated_error / (params.m * params.t * params.m * params.t);
        spdlog::info("Estimated error rate is {}", estimated_error);
        // We keep track of the actual number of bits we correct.
        int errors_corrected = 0;

        for (int current_column = 0; current_column < params.n; current_column++) {
            spdlog::info("Recovering column {0}...", current_column);

            // The target vector we are looking for (either from the leak or from the public key).
            MatrixGF2 target;
            if (current_column < params.m * params.t) {
                spdlog::info("Target at column {} taken from leak data.", current_column);
                target = leak.get_column(current_column);
            } else {
                int pk_column = current_column - params.m * params.t;
                spdlog::info("Target at column {} taken from column {} in public key.", current_column, pk_column);
                target = pk.get_column(pk_column);
            }

#if 0
            // Optionally output the desired target vector for debugging.
            if (1) {
                auto target_a = sk.alpha[current_column];
                auto target_ga = sk.goppa.evaluate(target_a);
                int idx;
                std::cout << candidates.nr_candidates << std::endl;
                std::cout << (target_ga.inv()).val << std::endl;
                std::cout << "Target: ";
                for (int i = 0; i < params.m * params.t; i++)
                    std::cout << target.get_entry(i, 0);
                std::cout << std::endl;
                bool found = false;
                for (idx = 0; idx < candidates.nr_candidates; idx++) {
                    if (candidates.idx2a(idx) == target_a && candidates.idx2ga(idx) == target_ga) {
                        std::cout << target_a << ", " << target_ga << ", " << idx << std::endl;
                        std::cout << "Candi.: ";
                        auto tc = candidates.candidates.get_column(idx);
                        for (int i = 0; i < params.m * params.t; i++)
                            //std::cout << candidates.candidates.get_entry(i, idx);
                            std::cout << tc.get_entry(i, 0);
                        std::cout << std::endl;
                        found = true;
                        break;
                    }
                }
                if (! found) {
                    std::cout << "Candi.: NOT FOUND" << std::endl;
                }
            }
#endif

            auto result = candidates.find_target(target);
            results.push_back(result);
            auto best_a = result.best_a;
            auto best_ga = result.best_ga;
            spdlog::info("Found {0}: {1}", current_column, result.to_string());

            // Sanity check.
            if (readout_error == 0 && result.best_distance > 1) {
                spdlog::critical("Distance is too large in error-free case!");
            }

            if (current_column >= params.m * params.t && result.best_distance > 0) {
                spdlog::critical("Distance is too large for decoding public-key, we must have failed!");
            }

            // Find goppa polynomial.
            if (current_column <= params.t) {
                goppa_points.push_back(std::make_pair(best_a, best_ga));
                spdlog::debug("Remembering goppa point.");
            }
            // If current_column is t, we have t+1 points (starting from 0).
            // Verify goppa polynomial.
            if (current_column == params.t) {
                // Only do this once.
                goppa = goppa.lagrange(goppa_points);
                spdlog::info("Goppa polynomial: {}", goppa.to_str());
                bool goppa_correct = sk.check_goppa(goppa);
                if (!goppa_correct) {
                    spdlog::critical("Secret key goppa differs! {}", sk.goppa.to_str());
                    break;
                }
            }
            if (current_column >= params.t) {
                auto ga_based_on_goppa = goppa.evaluate(best_a);
                if (ga_based_on_goppa == best_ga) {
                    spdlog::info("Verified the solution also based on the determined goppa polynomial!");
                } else {
                    spdlog::critical ("Goppa polynomial evaluation differs! best_ga = {}, goppa(best_a) = {}", best_ga.to_str(),
                                      ga_based_on_goppa.to_str());
                    break;
                }
            }


            // Verify the solution based on secret key. This allows us to abort failed trials early.
            bool alpha_correct = sk.check_alpha(current_column, best_a);
            if  (alpha_correct) {
                spdlog::info("Verified the solution also based on the secret key alpha!");
            } else {
                spdlog::critical("Secret key alpha differs! best_a = {}, a = {}", best_a.to_str(), sk.alpha[current_column].to_str());
                break;
            }

            {
                auto ga_based_on_sk_goppa = sk.goppa.evaluate(best_a);
                if (ga_based_on_sk_goppa == best_ga) {
                    spdlog::info("Verified the solution also based on the secret key goppa polynomial!");
                } else {
                    spdlog::critical ("Secret key Goppa polynomial evaluation differs! best_ga = {}, sk.goppa(best_a) = {}", best_ga.to_str(),
                                      ga_based_on_sk_goppa.to_str());
                    break;
                }
            }

            // We have to get the real column before defragmentation! Because result.best_idx is not stable.
            auto real_column = candidates.candidates.get_column(result.best_idx);
            for (int i = 0; i < params.m * params.t; i++) {
                if (target.get_entry(i, 0) != real_column.get_entry(i, 0))
                    errors_corrected++;
            }

            // OPTIMIZATION: We now invalidate best_idx by filtering the candidate matrix to remove all columns that have the same a.
            //std::cout << "Filtering based on alpha..." << std::endl;
            int nr_filtered = candidates.filter_candidates(best_a);
            spdlog::info("Removed {0} candidates for the same alpha.", nr_filtered);

            if (current_column == params.t) {
                // Only do this once.
                spdlog::info("Filtering based on goppa polynomial...");
                int nr_filtered_goppa = candidates.filter_candidates(goppa);
                spdlog::info("Removed {0} candidates based on goppa polynomial.", nr_filtered_goppa);

                // std::cout << "Defragmentation..." << std::endl;
                int defragmented = candidates.defragmentation();
                spdlog::info("Defragmentation removed a total of {0} columns", defragmented);
            }

            if (current_column < params.m * params.t) {
                spdlog::info("Applying Gauss elimination to candidate matrix.");
                // Based on this actual real first column, we can now forward-propagate the Gaussian
                // elimination steps done by Botan McEliece implementation.
                // From pk_gen.candidates
                int current_row = current_column;
                // Find all rows that are added to the current row. If the current column
                // is 0 at the current row, this is the first one with a 1 in that column,
                // followed by all rows with a 0 in that column. If the current column is
                // already 1, it's simply all rows with a 0.
                MatrixGF2 row_mask = MatrixGF2(params.m * params.t, 1);
                int i = current_row;
                bool encountered_one = real_column.get_entry(current_row, 0);
                for (i = current_row; i < params.m * params.t; i++) {
                    if (real_column.get_entry(i, 0))
                        break;
                }
                if (i == params.m * params.t) {
                    spdlog::critical("Did not find a 1 in the real column.");
                    break;
                }
                if (i != current_row)
                    candidates.add_row_to_row(current_row, i);

                // real_column already is almost the correct mask to add
                // to all columns where the current row is 1. We only need to clear the current row bit.
                auto what_we_add = real_column.get_column(0);
                what_we_add.set_entry(current_row, 0, 0);
                candidates.add_row_to_rows(what_we_add, current_row);
            }
        }
        auto end_recover_from_leak = std::chrono::high_resolution_clock::now();

        bool success = results.size() == params.n;
        double actual_error = ((double) errors_corrected)/(params.m*params.t * params.m*params.t);
        if (success) {
            spdlog::info("Corrected {} bits, yielding an actual error rate of {} (estimated: {})", errors_corrected,
                         actual_error, estimated_error);
        }

        if (with_report) {
            nlohmann::ordered_json columns = nlohmann::ordered_json::array();
            for (int j = 0; j < results.size(); j++) {
                auto result = results[j];
                nlohmann::ordered_json column = {
                        { "j", j },
                        { "nr_of_candidates", result.nr_candidates },
                        { "elapsed", result.elapsed },
                        { "distance", result.best_distance },
                        { "alpha", result.best_a.val },
                        { "g_of_alpha", result.best_ga.val },
                        { "correct", success || (j < results.size() - 1) }
                };
                columns.push_back(column);
            }
            nlohmann::ordered_json summary = {
                    { "leakfile", leakfile.filename },
                    { "kem", leakfile.kem },
                    { "params", {
                            { "m", params.m },
                            { "t", params.t },
                            { "n", params.n },
                            { "q", params.q },
                            { "k", params.k }
                    } },
                    { "error", readout_error },
                    { "estimated_error", estimated_error },
                    { "seed", rng.seed },
                    { "success", success },
                    { "elapsed", get_duration(begin_recover_from_leak, end_recover_from_leak) },
                    { "recovered_columns", success ? results.size() : results.size() - 1 },
            };
            if (success) {
                summary.push_back({"errors_corrected", errors_corrected });
                summary.push_back({"actual_error", actual_error });
            }
            summary.push_back({"columns", columns });
            std::ofstream file(report_fname);
            if (!file)
                throw std::runtime_error("Error opening result file: " + report_fname);
            file << summary.dump(2);
        }

        return success;
    }
};

template<typename Params>
Params Botan<Params>::params;

#endif // STOCKFISH_BOTAN_H
