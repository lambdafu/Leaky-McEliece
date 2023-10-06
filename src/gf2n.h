#ifndef STOCKFISH_GF2N_H
#define STOCKFISH_GF2N_H

#include <vector>
#include <sstream>
#include <cassert>
#include <cstdint>
#include <iostream>

// This class defines the general field with 2^m elements, represented by the polynomial ring
// over GF(2) modulo the irreducible polynomial irr (given as a bitmask, i.e. 0x1011 == x^3 + x + 1).
template<int m, int irr>
struct GF2N {
    // q is the number of elements in the field.
    const static int q = 1 << m;

    // The class for an element.
    struct Element {
        // The integer representation of that element.
        // Low order bits represent low order monomials.
        uint16_t val;

        // Create a new element with the given integer represenation (or 0).
        Element(uint16_t aVal = 0) {
            val = aVal;
        }

        // Test if the monomial POWER is contained in the polynomial
        // represenation of this element.
        int test(int power) const {
            return (val >> power) & 1;
        }

        // Get the inverse of this element.
        Element inv() const {
            return Element(GF2N<m, irr>::tables.inverse(val));
        }

        // Add this element to another element.
        Element operator+(const Element& rval) {
            return Element(val ^ rval.val);
        }

        // Subtract this element from another element.
        Element operator-(const Element& rval) {
            return Element(val ^ rval.val);
        }

        // Multiply this element with another element.
        Element operator*(const Element& rval) {
            return Element(GF2N<m, irr>::tables.multiply(val, rval.val));
        }

        // Divide this element by another element.
        Element operator/(const Element& rval) {
            return this->operator*(rval.inv());
        }

        // Test if this element is equal to another.
        bool operator==(const Element &rhs) const {
            return val == rhs.val;
        }

        // Test if this element is not equal to another.
        bool operator!=(const Element &rhs) const {
            return !(rhs == *this);
        }

        // Return the polynomial representation of this element as a string.
        const std::string to_str() const {
            std::ostringstream out;
            bool any = false;
            for (int i = m - 1; i >= 0 ; i--) {
                if (test(i)) {
                    if (any) {
                        out << " + ";
                    } else {
                        any = true;
                    }
                    if (i == 0)
                        out << "1";
                    else if (i == 1)
                        out << "z";
                    else
                        out << "z^" << i;
                }
            }
            if (any == false)
                out << "0";
            return out.str();
        }

        // Output the polynomial representation of an element to the ostream.
        friend std::ostream& operator<< ( std::ostream &os, const Element &x )
        {
            return os << x.to_str();
        }

    };

    // This type is defined just for automatic initialization of the inverse
    // table and the multiplication table.
    struct static_tables {
        // We store the product a*b at index (a*(a+1)/2+b), where a is >= b.
        // The order is 0*0, 1*0, 1*1, 2*0, 2*1, 2*2, 3*0, etc.
        std::vector<int> mult;
        // Well, inverse of 0 is set to 0 (but should not be used).
        std::vector<int> inv;

        // This function can calculate the product of two elements based on their
        // integer representation without using a lookup table.
        static int _mult(int a, int b) {
            int p = 0;
            while (a != 0 && b != 0) {
                if (b & 1)
                    p ^= a;
                a <<= 1;
                if (a & q)
                    a ^= irr;
                b >>= 1;
            }
            return p;
        }

        inline int multiply(int a, int b) {
            if (a < b)
                std::swap(a,b);
            return mult[a * (a + 1) / 2 + b];
        }

        inline int inverse(int a) {
            return inv[a];
        }

        // The constructor for this class generates the inverse table
        // and the multiplication table.
        static_tables() {
            mult.resize(q * (q+1) / 2 + q);
            inv.resize(q);

            inv[0] = 0;
            for (int a = 0; a < q; a++) {
                for (int b = 0; b <= a; b++) {
                    int p = _mult(a, b);
                    // As a side effect of calculating all products, we also collect the inverses.
                    if (p == 1) {
                        // In one case this will be redundant (a = b = 1), but it doesn't matter.
                        inv[a] = b;
                        inv[b] = a;
                    }
                    mult[a * (a + 1) / 2 + b] = p;
                }
            }
        }
    };
    static inline static_tables tables;

    GF2N() {
    }

    static Element from_int(int val) {
        return Element(val);
    }

    const std::string to_str() const {
        std::ostringstream out;
        out << "F_2^" << m;
        return out.str();
    }

    friend std::ostream &operator<<(std::ostream &os, GF2N<m, irr> const &f) {
        return os << f.to_str();
    }
};

#endif //STOCKFISH_GF2N_H
