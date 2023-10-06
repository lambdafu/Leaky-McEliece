#ifndef STOCKFISH_POLY_H
#define STOCKFISH_POLY_H

#include <vector>
#include <sstream>
#include <cassert>
#include <cstdint>
#include <iostream>

#include <spdlog/spdlog.h>

// This class defines a polynomial over a field. It is mainly used to calculate and
// evaluate the Lagrange polynomial.
template<typename F>
struct Poly {
    std::vector<typename F::Element> coefficients;

    Poly() {
    }

    // Return the degree of the polynomial (with deg 0 := 0 instead of -inf).
    int degree() const {
        if (coefficients.size() == 0)
            return 0;

        int d;
        for (d = coefficients.size() - 1; d > 0; d--) {
            if (coefficients[d] != F::from_int(0))
                break;
        }
        return d;
    }

    // Set the coefficient x_i to c.
    void set(int i, typename F::Element c) {
        if (i + 1 > coefficients.size()) {
            coefficients.resize(i + 1);
        }
        coefficients[i] = c;
    }

    // Construct and return the polynomial x+c.
    static Poly linear(typename F::Element c) {
        Poly result;
        result.set(1, 1);
        result.set(0, c);
        return result;
    }

    // Construct and return the polynomial c.
    static Poly constant(typename F::Element c) {
        Poly result;
        result.set(0, c);
        return result;
    }

    // Construct and return this + rval.
    Poly operator+(const Poly<F>& rval) {
        Poly result = rval;
        if (coefficients.size() > result.coefficients.size())
            result.coefficients.resize(coefficients.size());
        for (int i = 0; i < coefficients.size(); i++) {
            result.coefficients[i] = result.coefficients[i] + coefficients[i];
        }
        return result;
    }

    // Construct and return e * this.
    Poly scalar_multiple(typename F::Element e) {
        Poly result = *this;
        for (int i = 0; i < coefficients.size(); i++) {
            result.coefficients[i] = result.coefficients[i] * e;
        }
        return result;
    }

    // Construct and return this * x^n.
    Poly shift(int n) const {
        Poly result;
        result.coefficients.resize(coefficients.size() + n);
        for (int i = 0; i < coefficients.size(); i++) {
            result.coefficients[i + n] = coefficients[i];
        }
        return result;
    }

    // Construct and return this * rval.
    Poly operator*(const Poly<F>& rval) const {
        Poly result;

        for (int i = 0; i <= degree(); i++) {
            Poly term = rval;
            term = term.scalar_multiple(coefficients[i]);
            term = term.shift(i);
            result = result + term;
        }

        return result;
    }

    // Return true if this is equal to rval.
    bool operator==(const Poly<F>& rval) const {
        if (degree() != rval.degree())
            return false;

        for (int i = 0; i <= degree(); i++) {
            if (coefficients[i] != rval.coefficients[i])
                return false;
        }
        return true;
    }

    // Evaluate this polynomial at posixition x and return the result.
    typename F::Element evaluate(const typename F::Element& x) const {
        typename F::Element result = F::from_int(0);

        for (int i = coefficients.size() - 1; i >= 0; i--) {
            result = result * x + coefficients[i];
        }

        return result;
    }

    // Return a string representation of this polynomial.
    const std::string to_str() const {
        std::ostringstream out;

        if (coefficients.size() == 0) {
            out << "0";
            return out.str();
        }

        int d = degree();
        bool first = true;
        for (int i = d; i >= 0; i--) {
            if (i == 0) {
                if (coefficients[i] == F::from_int(0)) {
                    if (first)
                        out << coefficients[i];
                } else {
                    if (!first)
                        out << " + ";
                    out << coefficients[i];
                }
            } else if (coefficients[i] != F::from_int(0)) {
                if (!first) {
                    out << " + ";
                }
                if (coefficients[i] != F::from_int(1)) {
                    out << "(" << coefficients[i] << ")*";
                }
                out << "X";
                if (i != 1)
                    out << "^" << i;
                first = false;
            }
        }
        return out.str();
    }

    friend std::ostream &operator<<(std::ostream &os, Poly<F> const &p) {
        return os << p.to_str();
    }

    // Construct and return the Lagrange interpolation for the given points.
    static Poly<F> lagrange(std::vector<std::pair<typename F::Element, typename F::Element>> points) {
        Poly<F> result;

        for (int i = 0; i < points.size(); i++) {
            Poly<F> basis = Poly<F>::constant(F::from_int(1));

            for (int j = 0; j < points.size(); j++) {
                if (i != j) {
                    Poly<F> numeratorPoly = Poly<F>::linear(points[j].first);
                    typename F::Element denominator = points[i].first - points[j].first;
                    basis = basis * (numeratorPoly.scalar_multiple(denominator.inv()));
                }
            }

            basis = basis.scalar_multiple(points[i].second);
            spdlog::debug("Lagrange basis {0}: {1}", i, basis.to_str());
            result = result + basis;
        }

        return result;
    }
};

#endif //STOCKFISH_POLY_H
