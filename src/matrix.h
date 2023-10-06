#ifndef STOCKFISH_MATRIX_H
#define STOCKFISH_MATRIX_H

// This is an implementation of matrices over F2, optimized to access all rows of a column
// in packs of 64 bit each.
struct MatrixGF2 {
    int rows;
    int cols;

    // We pack the rows of each column together.
    typedef uint64_t pack_t;
    constexpr static int rows_per_pack = sizeof(pack_t) * 8;

    int packs_per_column;
    std::vector<pack_t> packs;

    MatrixGF2(int aRows = 0, int aColumns = 0) {
        resize(aRows, aColumns);
    }

    void resize(int aRows, int aColumns) {
        rows = aRows;
        cols = aColumns;

        packs_per_column = (rows + rows_per_pack - 1) / rows_per_pack;
        packs.resize(packs_per_column * cols);
    }

    // Get a single bit. This is somewhat slow.
    int get_entry(int i, int j) const {
        int pack = j * packs_per_column + i / rows_per_pack;
        int bit = i % rows_per_pack;
        return (packs[pack] >> bit) & 1;
    }

    // Set a single bit to val (which is 0 or 1). This is somewhat slow.
    void set_entry(int i, int j, int val) {
        int pack = j * packs_per_column + i / rows_per_pack;
        int bit = i % rows_per_pack;
        pack_t val_mask = ((pack_t) 1) << bit;
        pack_t val_bit = (((pack_t) val) & 1) << bit;
        packs[pack] = (packs[pack] & ~val_mask) | val_bit;
    }

    void randomize() {
        uint64_t mask_last_pack_in_column;
        if (rows % rows_per_pack == 0)
            mask_last_pack_in_column = ~0;
        else {
            mask_last_pack_in_column = 1;
            mask_last_pack_in_column <<= rows % rows_per_pack;
            mask_last_pack_in_column -= 1;
        }
        for (int i = 0; i < packs.size(); i++) {
            int row = i / packs_per_column;
            int packcol = i % packs_per_column;
            uint64_t r1 = std::rand();
            uint64_t r2 = std::rand();
            uint64_t r3 = std::rand();
            uint64_t r4 = std::rand();

            pack_t pack = (r1 << 48) ^ (r2 << 32) ^ (r3 << 16) ^ r4;
            if (packcol == packs_per_column - 1) {
                // clamp last value in row.
                pack = pack & mask_last_pack_in_column;
            }
            packs[i] = pack;
        }
    }

    void debug() {
        std::cout << "MatrixGF2(" << rows << ", " << cols << ")" << std::endl;
#if 0
        for (int i = 0; i < 20 && i < rows; i++) {
            for (int j = 0; j < 78 && j < cols; j++) {
                std::cout << get_entry(i, j);
            }
            std::cout << std::endl;
        }
#endif
        std::cout << "packs: " << packs.size() << " (" << packs_per_column << " per row)" << std::endl;
        for (int i = 0; i < 40 && i < packs.size(); i++) {
            printf("%016llx ", (long long unsigned int) packs[i]);
            if ((i+1) % 4 == 0)
                printf("\n");
        }
        printf("\n");
    }

    MatrixGF2 get_column(int j) {
        MatrixGF2 result(rows, 1);
        pack_t* packptr = &packs[j * packs_per_column];
        for (int i = 0; i < rows; i += rows_per_pack) {
            result.packs[i/rows_per_pack] = packptr[i/rows_per_pack];
        }
        return result;
    }

    void add_column(int j, const MatrixGF2& c) {
        pack_t* packptr = &packs[j * packs_per_column];
        for (int i = 0; i < rows; i += rows_per_pack) {
            packptr[i/rows_per_pack] ^= c.packs[i/rows_per_pack];
        }
    }

    // Add row src to row dst.
    void add_row(int src, int dst) {
        for (int j = 0; j < cols; j++) {
            int src_bit = get_entry(src, j);
            int dst_bit = get_entry(dst, j);
            set_entry(dst, j, dst_bit ^ src_bit);
        }
    }
};

#endif //STOCKFISH_MATRIX_H
