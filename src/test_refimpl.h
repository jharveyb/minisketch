/**********************************************************************
 * Copyright (c) 2026 The minisketch developers                        *
 * Distributed under the MIT software license, see the accompanying   *
 * file LICENSE or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

/** Independent reference implementations of field and polynomial arithmetic,
 * used as oracles by the unit tests (unit_tests_impl.h) and the fuzz targets
 * (src/fuzz/poly_ops.cpp). Deliberately naive: correctness over speed, and no
 * sharing of code paths with the implementations under test beyond the field
 * element type. Polynomials are coefficient vectors, lowest degree first. */

#ifndef _MINISKETCH_TEST_REFIMPL_H_
#define _MINISKETCH_TEST_REFIMPL_H_

#include <cstdint>
#include <vector>

namespace ut {

/** The moduli minisketch uses, as their low parts (modulus minus the x^bits
 *  term), from doc/moduli.md / tests/pyminisketch.py. Hardcoding them here
 *  (rather than reading the LFSR parameter from the field classes) also
 *  verifies each field was generated for the intended modulus. */
constexpr uint64_t REF_LOWMOD[65] = {
    0, 0, 3, 3, 3, 5, 3, 3, 27, 3, 9, 5, 9, 27, 33, 3, 43,
    9, 9, 39, 9, 5, 3, 33, 27, 9, 27, 39, 3, 5, 3, 9, 141,
    1025, 129, 5, 513, 83, 99, 17, 57, 9, 129, 89, 33, 27, 3, 33, 45,
    513, 29, 75, 9, 71, 513, 129, 149, 17, 524289, 149, 3, 39, 536870913, 3, 27
};

/** Carryless "Russian peasant" multiplication in GF(2^bits) with the given
 *  modulus low part. All intermediate values stay below 2^bits, so this works
 *  for bits up to 64. */
inline uint64_t RefMul(uint64_t a, uint64_t b, int bits, uint64_t lowmod) {
    const uint64_t mask = bits == 64 ? ~uint64_t{0} : (uint64_t{1} << bits) - 1;
    uint64_t ret = 0;
    while (b != 0) {
        if (b & 1) ret ^= a;
        uint64_t top = a >> (bits - 1);
        a = (a << 1) & mask;
        if (top) a ^= lowmod;
        b >>= 1;
    }
    return ret;
}

/** The field trace Tr(a) = a + a^2 + a^4 + ... + a^(2^(bits-1)); always 0 or 1. */
inline uint64_t RefTrace(uint64_t a, int bits, uint64_t lowmod) {
    uint64_t t = a, acc = a;
    for (int i = 1; i < bits; ++i) {
        t = RefMul(t, t, bits, lowmod);
        acc ^= t;
    }
    return acc;
}

/** a^e in GF(2^bits). */
inline uint64_t RefPow(uint64_t a, uint64_t e, int bits, uint64_t lowmod) {
    uint64_t ret = 1;
    while (e != 0) {
        if (e & 1) ret = RefMul(ret, a, bits, lowmod);
        a = RefMul(a, a, bits, lowmod);
        e >>= 1;
    }
    return ret;
}

/** Schoolbook polynomial multiplication. */
template<typename F>
std::vector<typename F::Elem> PolyMulRef(const std::vector<typename F::Elem>& a, const std::vector<typename F::Elem>& b, const F& field) {
    if (a.empty() || b.empty()) return {};
    std::vector<typename F::Elem> out(a.size() + b.size() - 1, 0);
    for (size_t i = 0; i < a.size(); ++i) {
        if (a[i] == 0) continue;
        for (size_t j = 0; j < b.size(); ++j) {
            out[i + j] ^= field.Mul(a[i], b[j]);
        }
    }
    while (!out.empty() && out.back() == 0) out.pop_back();
    return out;
}

/** Naive reduction of val modulo a monic polynomial mod (in place). */
template<typename F>
void PolyReduceRef(std::vector<typename F::Elem>& val, const std::vector<typename F::Elem>& mod, const F& field) {
    while (val.size() >= mod.size()) {
        auto lead = val.back();
        val.pop_back();
        if (lead != 0) {
            for (size_t i = 0; i + 1 < mod.size(); ++i) {
                val[val.size() - (mod.size() - 1) + i] ^= field.Mul(lead, mod[i]);
            }
        }
    }
    while (!val.empty() && val.back() == 0) val.pop_back();
}

/** Strip trailing zero coefficients (for comparing polynomials). */
template<typename E>
std::vector<E> Stripped(std::vector<E> poly) {
    while (!poly.empty() && poly.back() == 0) poly.pop_back();
    return poly;
}

/** Build the monic polynomial prod_i (x + roots[i]). */
template<typename F>
std::vector<typename F::Elem> PolyFromRootsRef(const std::vector<typename F::Elem>& roots, const F& field) {
    std::vector<typename F::Elem> poly{1};
    for (auto root : roots) {
        std::vector<typename F::Elem> factor{root, 1};
        poly = PolyMulRef(poly, factor, field);
    }
    return poly;
}

} // namespace ut

#endif
