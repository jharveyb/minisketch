/**********************************************************************
 * Copyright (c) 2026 The minisketch developers                        *
 * Distributed under the MIT software license, see the accompanying   *
 * file LICENSE or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

#ifndef _MINISKETCH_TEST_UTILS_H_
#define _MINISKETCH_TEST_UTILS_H_

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <vector>

#include "int_utils.h"

/** Deterministic seedable RNG for tests: SipHash-2-4 in counter mode.
 *
 * Unlike std::random_device (not reproducible) or std::uniform_int_distribution
 * (implementation-defined output), every value produced by this generator is a
 * pure function of the seed, on every platform. Any test failure can therefore
 * be replayed exactly by re-running with the seed that is printed at startup.
 */
class TestRand {
    uint64_t m_k0;
    uint64_t m_k1;
    uint64_t m_count{0};

public:
    explicit TestRand(uint64_t seed) : m_k0(seed), m_k1(0x746573747365656dULL /* "testseem" */) {}

    /** Generate a uniformly random 64-bit number. */
    uint64_t Rand64() { return SipHash(m_k0, m_k1, m_count++); }

    /** Generate a uniformly random `bits`-bit number (bits in [0, 64]). */
    uint64_t RandBits(int bits) {
        uint64_t r = Rand64();
        return bits >= 64 ? r : r & ((uint64_t{1} << bits) - 1);
    }

    /** Generate a uniformly random number in [0, range). Requires range > 0.
     *  Uses mask + rejection sampling, so there is no modulo bias. */
    uint64_t RandRange(uint64_t range) {
        uint64_t mask = range - 1;
        mask |= mask >> 1;
        mask |= mask >> 2;
        mask |= mask >> 4;
        mask |= mask >> 8;
        mask |= mask >> 16;
        mask |= mask >> 32;
        while (true) {
            uint64_t ret = Rand64() & mask;
            if (ret < range) return ret;
        }
    }

    /** Generate a uniformly random number in [lo, hi] (inclusive). */
    uint64_t RandIncl(uint64_t lo, uint64_t hi) { return lo + RandRange(hi - lo + 1); }
};

/** Determine the test seed to use.
 *
 * Priority: a "--seed=N" argument in argv (N decimal or 0x-prefixed hex), then
 * the MINISKETCH_TEST_SEED environment variable, then a random seed. The chosen
 * seed is printed to stdout so failing runs can be reproduced.
 */
inline uint64_t GetTestSeed(int argc, char** argv) {
    uint64_t seed = 0;
    bool have_seed = false;
    for (int i = 1; i < argc; ++i) {
        if (std::strncmp(argv[i], "--seed=", 7) == 0) {
            seed = std::strtoull(argv[i] + 7, nullptr, 0);
            have_seed = true;
        }
    }
    if (!have_seed) {
        const char* env = std::getenv("MINISKETCH_TEST_SEED");
        if (env && *env) {
            seed = std::strtoull(env, nullptr, 0);
            have_seed = true;
        }
    }
    if (!have_seed) {
        std::random_device rnd;
        seed = (uint64_t{rnd()} << 32) | rnd();
    }
    printf("Test seed: 0x%016llx (reproduce with --seed=0x%016llx or MINISKETCH_TEST_SEED)\n",
           (unsigned long long)seed, (unsigned long long)seed);
    return seed;
}

/* Generators for field elements and polynomials (usable with any field class F
 * exposing the interface used by sketch_impl.h). Polynomials are coefficient
 * vectors, lowest degree first, like everywhere else in this codebase. */

/** Generate a uniformly random field element (possibly zero). */
template<typename F>
typename F::Elem RandElem(TestRand& rng, const F& field) {
    return static_cast<typename F::Elem>(rng.RandBits(field.Bits()));
}

/** Generate a uniformly random nonzero field element. */
template<typename F>
typename F::Elem RandNonzeroElem(TestRand& rng, const F& field) {
    while (true) {
        auto elem = RandElem(rng, field);
        if (elem != 0) return elem;
    }
}

/** Generate a random polynomial of exactly the given size (degree size-1),
 *  with nonzero leading coefficient. Size 0 gives the zero polynomial. */
template<typename F>
std::vector<typename F::Elem> RandPoly(TestRand& rng, const F& field, size_t size) {
    std::vector<typename F::Elem> poly(size);
    for (size_t i = 0; i < size; ++i) poly[i] = RandElem(rng, field);
    if (size > 0) poly.back() = RandNonzeroElem(rng, field);
    return poly;
}

/** Generate a random monic polynomial of exactly the given size (degree size-1).
 *  Requires size >= 1. */
template<typename F>
std::vector<typename F::Elem> RandMonicPoly(TestRand& rng, const F& field, size_t size) {
    auto poly = RandPoly(rng, field, size);
    poly.back() = 1;
    return poly;
}

/** Generate `count` distinct nonzero field elements. Requires count < 2^bits. */
template<typename F>
std::vector<typename F::Elem> RandDistinctElems(TestRand& rng, const F& field, size_t count) {
    std::vector<typename F::Elem> ret;
    ret.reserve(count);
    while (ret.size() < count) {
        auto elem = RandNonzeroElem(rng, field);
        bool duplicate = false;
        for (auto e : ret) duplicate |= (e == elem);
        if (!duplicate) ret.push_back(elem);
    }
    return ret;
}

#endif
