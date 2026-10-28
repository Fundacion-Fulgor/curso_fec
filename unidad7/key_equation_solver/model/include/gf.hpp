/**
 * @file gf.hpp
 * @brief GF(2^m) arithmetic with exp/log tables, for 1 <= m <= 16.
 *
 * Elements are integers in [0, 2^m - 1] where bit i is the coefficient of x^i, as in galois.py.
 */

#pragma once

#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace fec {

using Symbol = std::uint32_t;

class GaloisField {
public:
    /**
     * @brief Build GF(2^m) with alpha = x.
     * @param m Field degree, from 1 to 16.
     * @param poly Primitive polynomial without the x^m term (x^8 + x^4 + x^3 + x^2 + 1 -> 0x1D).
     * @throws std::invalid_argument If m or poly are out of range, or x is not a generator.
     */
    GaloisField(unsigned m, std::uint32_t poly) : m_(m), poly_(poly) {
        if (m < 1 || m > 16) {
            throw std::invalid_argument("m must be in [1, 16], got " + std::to_string(m));
        }
        order_ = Symbol{1} << m;
        if (poly >= order_) {
            throw std::invalid_argument("poly does not fit in m bits (the x^m term is implicit)");
        }
        if (m > 1 && !(poly & 1u)) {
            throw std::invalid_argument("P(x) is divisible by x");
        }
        build_tables();
    }

    unsigned m() const { return m_; }
    std::uint32_t poly() const { return poly_; }
    Symbol order() const { return order_; }
    Symbol n() const { return order_ - 1; }

    /**
     * @brief Validate a field element.
     * @throws std::invalid_argument If a is outside [0, 2^m - 1].
     */
    Symbol check(Symbol a) const {
        if (a >= order_) {
            throw std::invalid_argument(std::to_string(a) + " is not in GF(2^" + std::to_string(m_) + ")");
        }
        return a;
    }

    /** @brief Product a * b; the arguments are not validated. */
    Symbol mul(Symbol a, Symbol b) const {
        if (a == 0 || b == 0) {
            return 0;
        }
        return exp_[(log_[a] + log_[b]) % n()];
    }

    /**
     * @brief Multiplicative inverse.
     * @throws std::domain_error If a == 0.
     */
    Symbol inverse(Symbol a) const {
        if (check(a) == 0) {
            throw std::domain_error("0 has no multiplicative inverse");
        }
        return exp_[(n() - log_[a]) % n()];
    }

    /** @brief alpha^k for any integer k. */
    Symbol exp(long long k) const {
        long long r = k % static_cast<long long>(n());
        return exp_[static_cast<std::size_t>(r < 0 ? r + n() : r)];
    }

    /**
     * @brief Discrete logarithm in base alpha.
     * @throws std::domain_error If a == 0.
     */
    unsigned log(Symbol a) const {
        if (check(a) == 0) {
            throw std::domain_error("log(0) is undefined");
        }
        return log_[a];
    }

private:
    Symbol xtime(Symbol a) const {
        a <<= 1;
        if (a & order_) {
            a ^= order_ | poly_;
        }
        return a;
    }

    void build_tables() {
        exp_.assign(n(), 0);
        log_.assign(order_, 0);
        std::vector<bool> seen(order_, false);
        Symbol elem = 1;
        const Symbol x = xtime(1);    // 1 when m == 1
        for (Symbol i = 0; i < n(); ++i) {
            if (seen[elem]) {
                throw std::invalid_argument("P(x) is not primitive");
            }
            seen[elem] = true;
            exp_[i] = elem;
            log_[elem] = i;
            Symbol prod = 0;    // elem * x, carry-less, reduced
            for (Symbol a = elem, b = x; b; b >>= 1, a = xtime(a)) {
                if (b & 1u) {
                    prod ^= a;
                }
            }
            elem = prod;
        }
        if (elem != 1) {
            throw std::invalid_argument("P(x) is not primitive");
        }
    }

    unsigned m_;
    std::uint32_t poly_;
    Symbol order_ = 0;
    std::vector<Symbol> exp_;
    std::vector<Symbol> log_;
};

}  // namespace fec
