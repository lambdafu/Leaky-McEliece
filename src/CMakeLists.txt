cmake_minimum_required(VERSION 3.22)
project(Stockfish)

set(CMAKE_CXX_STANDARD 17)

include_directories(3rdparty)
add_executable(Stockfish main.cpp gf2n.h matrix.h candidates.h leakfile.h poly.h classic.h randomness.h main.h time.h botan.h)
