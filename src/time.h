#ifndef STOCKFISH_TIME_H
#define STOCKFISH_TIME_H

#include <chrono>

static uint64_t get_duration(std::chrono::high_resolution_clock::time_point start, std::chrono::high_resolution_clock::time_point end) {
    return std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
}

#endif //STOCKFISH_TIME_H
