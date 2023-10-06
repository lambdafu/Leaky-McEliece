#ifndef STOCKFISH_RANDOMNESS_H
#define STOCKFISH_RANDOMNESS_H

#include <cstdint>
#include <chrono>
#include <random>

struct Randomness {
    unsigned int seed;
    std::mt19937 engine;

    Randomness(unsigned int aSeed = 0) {
        if (aSeed == 0) {
            std::random_device os_seed;
            seed = os_seed();
        } else {
            seed = aSeed;
        }
        std::cout << "Using seed: " << seed << std::endl;
        engine.seed(seed);
    }

    int sample_uniform_bit() {
        auto distribution = std::uniform_int_distribution<int>(0, 1);
        return distribution(engine);
    }

    // Number between 0 and up_to-1.
    int sample_uniform_number(int up_to) {
        auto distribution = std::uniform_int_distribution<int>(0, up_to - 1);
        return distribution(engine);
    }

    double sample_uniform_unit() {
        auto distribution = std::uniform_real_distribution<double>(0, 1);
        return distribution(engine);
    }

    // FIXME: Add test for randomness.
};

#endif //STOCKFISH_RANDOMNESS_H
