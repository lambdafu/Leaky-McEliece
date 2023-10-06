#ifndef STOCKFISH_MCELIECE_H
#define STOCKFISH_MCELIECE_H

#include <vector>
#include <sstream>
#include <cassert>
#include <cstdint>
#include <iostream>
#include "gf2n.h"

enum McElieceParameters {
    // Testing
    Classic51220,
    Classic102450,

    // NIST submission
    Classic348864,
    Classic460896,
    Classic6688128,
    Classic6960119,
    Classic8192128,

    Botan51220,
    Botan102450,

    // NIST submission
    Botan348864,
    Botan460896,
    Botan6688128,
    Botan6960119,
    Botan8192128,


};

template<McElieceParameters p>
struct McEliece;

template<>
struct McEliece<Classic51220> {
    McElieceParameters type = Classic51220;
    static constexpr const char*const name = "Classic51220";
    static constexpr int m = 9;
    static constexpr int q = 1 << m;
    static constexpr int n = 512;
    static constexpr int t = 20;
    static constexpr int k = n - t * m;
    static constexpr int irr = 0b10'0000'0011; // z^9 + z^1 + 1
    // This also generates the static multiplication and inverse tables.
    GF2N<m, irr> gf2m;
    typedef GF2N<m, irr>::Element element_t;
};

template<>
struct McEliece<Classic102450> {
    McElieceParameters type = Classic102450;
    static constexpr const char*const name = "Classic102450";
    static constexpr int m = 10;
    static constexpr int q = 1 << m;
    static constexpr int n = 1024;
    static constexpr int t = 50;
    static constexpr int k = n - t * m;
    static constexpr int irr = 0b100'0000'1001; // z^10 + z^3 + 1
    // This also generates the static multiplication and inverse tables.
    GF2N<m, irr> gf2m;
    typedef GF2N<m, irr>::Element element_t;
};

template<>
struct McEliece<Classic348864> {
    McElieceParameters type = Classic348864;
    static constexpr const char*const name = "Classic348864";
    static constexpr int m = 12;
    static constexpr int q = 1 << m;
    static constexpr int n = 3488;
    static constexpr int t = 64;
    static constexpr int k = n - t * m;
    static constexpr int irr = 0b1'0000'0000'1001; // z^12 + z^3 + 1
    // This also generates the static multiplication and inverse tables.
    GF2N<m, irr> gf2m;
    typedef GF2N<m, irr>::Element element_t;
};

template<>
struct McEliece<Classic460896> {
    McElieceParameters type = Classic460896;
    static constexpr const char*const name = "Classic460896";
    static constexpr int m = 13;
    static constexpr int q = 1 << m;
    static constexpr int n = 4608;
    static constexpr int t = 96;
    static constexpr int k = n - t * m;
    static constexpr int irr = 0b10'0000'0001'1011; // z^13+z^4+z^3+z+1
    // This also generates the static multiplication and inverse tables.
    GF2N<m, irr> gf2m;
    typedef GF2N<m, irr>::Element element_t;
};

template<>
struct McEliece<Classic6688128> {
    McElieceParameters type = Classic6688128;
    static constexpr const char*const name = "Classic6688128";
    static constexpr int m = 13;
    static constexpr int q = 1 << m;
    static constexpr int n = 6688;
    static constexpr int t = 128;
    static constexpr int k = n - t * m;
    static constexpr int irr = 0b10'0000'0001'1011; // z^13+z^4+z^3+z+1
    // This also generates the static multiplication and inverse tables.
    GF2N<m, irr> gf2m;
    typedef GF2N<m, irr>::Element element_t;
};

template<>
struct McEliece<Classic6960119> {
    McElieceParameters type = Classic6960119;
    static constexpr const char*const name = "Classic6960119";
    static constexpr int m = 13;
    static constexpr int q = 1 << m;
    static constexpr int n = 6960;
    static constexpr int t = 119;
    static constexpr int k = n - t * m;
    static constexpr int irr = 0b10'0000'0001'1011; // z^13+z^4+z^3+z+1
    // This also generates the static multiplication and inverse tables.
    GF2N<m, irr> gf2m;
    typedef GF2N<m, irr>::Element element_t;
};

template<>
struct McEliece<Classic8192128> {
    McElieceParameters type = Classic8192128;
    static constexpr const char*const name = "Classic8192128";
    static constexpr int m = 13;
    static constexpr int q = 1 << m;
    static constexpr int n = 8192;
    static constexpr int t = 128;
    static constexpr int k = n - t * m;
    static constexpr int irr = 0b10'0000'0001'1011; // z^13+z^4+z^3+z+1
    // This also generates the static multiplication and inverse tables.
    GF2N<m, irr> gf2m;
    typedef GF2N<m, irr>::Element element_t;
};

template<>
struct McEliece<Botan51220> {
    McElieceParameters type = Botan51220;
    static constexpr const char*const name = "Botan51220";
    static constexpr int m = 9;
    static constexpr int q = 1 << m;
    static constexpr int n = 512;
    static constexpr int t = 20;
    static constexpr int k = n - t * m;
    static constexpr int irr = 01041; // see botan/lib/pubkey/mce/gf2m_small_m.cpp
    // This also generates the static multiplication and inverse tables.
    GF2N<m, irr> gf2m;
    typedef GF2N<m, irr>::Element element_t;
};

template<>
struct McEliece<Botan102450> {
    McElieceParameters type = Botan102450;
    static constexpr const char*const name = "Botan102450";
    static constexpr int m = 10;
    static constexpr int q = 1 << m;
    static constexpr int n = 1024;
    static constexpr int t = 50;
    static constexpr int k = n - t * m;
    static constexpr int irr = 02011; // see botan/lib/pubkey/mce/gf2m_small_m.cpp
    // This also generates the static multiplication and inverse tables.
    GF2N<m, irr> gf2m;
    typedef GF2N<m, irr>::Element element_t;
};

template<>
struct McEliece<Botan348864> {
    McElieceParameters type = Botan348864;
    static constexpr const char*const name = "Botan348864";
    static constexpr int m = 12;
    static constexpr int q = 1 << m;
    static constexpr int n = 3488;
    static constexpr int t = 64;
    static constexpr int k = n - t * m;
    static constexpr int irr = 010123; // see botan/lib/pubkey/mce/gf2m_small_m.cpp
    // This also generates the static multiplication and inverse tables.
    GF2N<m, irr> gf2m;
    typedef GF2N<m, irr>::Element element_t;
};

template<>
struct McEliece<Botan460896> {
    McElieceParameters type = Botan460896;
    static constexpr const char*const name = "Botan460896";
    static constexpr int m = 13;
    static constexpr int q = 1 << m;
    static constexpr int n = 4608;
    static constexpr int t = 96;
    static constexpr int k = n - t * m;
    static constexpr int irr = 020033; // see botan/lib/pubkey/mce/gf2m_small_m.cpp
    // This also generates the static multiplication and inverse tables.
    GF2N<m, irr> gf2m;
    typedef GF2N<m, irr>::Element element_t;
};

template<>
struct McEliece<Botan6688128> {
    McElieceParameters type = Botan6688128;
    static constexpr const char*const name = "Botan6688128";
    static constexpr int m = 13;
    static constexpr int q = 1 << m;
    static constexpr int n = 6688;
    static constexpr int t = 128;
    static constexpr int k = n - t * m;
    static constexpr int irr = 020033; // see botan/lib/pubkey/mce/gf2m_small_m.cpp
    // This also generates the static multiplication and inverse tables.
    GF2N<m, irr> gf2m;
    typedef GF2N<m, irr>::Element element_t;
};

template<>
struct McEliece<Botan6960119> {
    McElieceParameters type = Botan6960119;
    static constexpr const char*const name = "Botan6960119";
    static constexpr int m = 13;
    static constexpr int q = 1 << m;
    static constexpr int n = 6960;
    static constexpr int t = 119;
    static constexpr int k = n - t * m;
    static constexpr int irr = 020033; // see botan/lib/pubkey/mce/gf2m_small_m.cpp
    // This also generates the static multiplication and inverse tables.
    GF2N<m, irr> gf2m;
    typedef GF2N<m, irr>::Element element_t;
};

template<>
struct McEliece<Botan8192128> {
    McElieceParameters type = Botan8192128;
    static constexpr const char*const name = "Botan8192128";
    static constexpr int m = 13;
    static constexpr int q = 1 << m;
    static constexpr int n = 8192;
    static constexpr int t = 128;
    static constexpr int k = n - t * m;
    static constexpr int irr = 020033; // see botan/lib/pubkey/mce/gf2m_small_m.cpp
    // This also generates the static multiplication and inverse tables.
    GF2N<m, irr> gf2m;
    typedef GF2N<m, irr>::Element element_t;
};


//010123,  /* extension degree 12 */
//020033,  /* extension degree 13 */


#endif //STOCKFISH_MCELIECE_H
