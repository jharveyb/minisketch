/**********************************************************************
 * Copyright (c) 2026 The minisketch developers                        *
 * Distributed under the MIT software license, see the accompanying   *
 * file LICENSE or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

/** Property tests for the internal algorithm stages of libminisketch.
 *
 * Everything here is templated over a field class F (the interface used by
 * sketch_impl.h), so the same properties run against every generic and clmul
 * field instantiation. Field arithmetic is checked against an independent
 * bit-level reference implementation (RefMul below) with hardcoded moduli, and
 * the polynomial algorithms are checked against naive schoolbook references.
 *
 * This is the layer to extend when swapping out the implementation of a decode
 * stage: add a test that the new implementation matches the corresponding
 * *Ref reference (or the old implementation) on the existing generators.
 *
 * This header must be included from a TU that also includes field
 * implementation files from src/fields/ (which bring in sketch_impl.h and
 * define the field classes in the anonymous namespace).
 */

#ifndef _MINISKETCH_UNIT_TESTS_IMPL_H_
#define _MINISKETCH_UNIT_TESTS_IMPL_H_

#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cstdint>
#include <vector>

#include "sketch_impl.h"
#include "test_refimpl.h"
#include "test_utils.h"

/** Test seed, initialized once by the SeedFixture in unit_tests.cpp. */
extern uint64_t g_test_seed;

/** Fail the current Boost test case with a message; use instead of
 *  BOOST_CHECK in hot loops so millions of passing checks cost nothing. */
#define UT_REQUIRE(cond) do { if (!(cond)) BOOST_FAIL("Check failed: " #cond); } while(0)

namespace ut {

/* Reference field/polynomial arithmetic lives in test_refimpl.h (shared with
 * the fuzz targets). */

/* ---------- Field operation properties ---------- */

template<typename F>
void CheckFieldElem(const F& field, uint64_t lowmod, uint64_t a64, uint64_t b64) {
    typedef typename F::Elem Elem;
    const int bits = field.Bits();
    Elem a = static_cast<Elem>(a64), b = static_cast<Elem>(b64);

    // Multiplication (and the precomputed-table Multiplier) against the reference.
    Elem ab = field.Mul(a, b);
    UT_REQUIRE(uint64_t(ab) == RefMul(a64, b64, bits, lowmod));
    typename F::Multiplier mul_a(field, a);
    UT_REQUIRE(mul_a(b) == ab);

    // Squaring and doubling are just special multiplications.
    UT_REQUIRE(uint64_t(field.Sqr(a)) == RefMul(a64, a64, bits, lowmod));
    UT_REQUIRE(uint64_t(field.Mul2(a)) == RefMul(a64, 2, bits, lowmod));

    // Inverses.
    if (a != 0) {
        Elem inv = field.Inv(a);
        UT_REQUIRE(inv != 0);
        UT_REQUIRE(field.Mul(a, inv) == 1);
        UT_REQUIRE(field.Inv(inv) == a);
    }

    // Qrt solves r^2 + r = a whenever a solution exists (i.e. Tr(a) == 0);
    // for Tr(a) == 1 its result is undefined.
    Elem r = field.Qrt(a);
    if ((field.Sqr(r) ^ r) != a) {
        UT_REQUIRE(RefTrace(a64, bits, lowmod) == 1);
    } else {
        UT_REQUIRE(RefTrace(a64, bits, lowmod) == 0);
    }

    // FromUint64/ToUint64 roundtrip on in-range values.
    UT_REQUIRE(field.ToUint64(field.FromUint64(a64)) == a64);
}

template<typename F>
void TestFieldOps(const F& field, uint64_t lowmod, TestRand& rng, size_t iters) {
    const int bits = field.Bits();

    if (bits <= 8) {
        // Small fields: exhaustive over all element pairs.
        for (uint64_t a = 0; a >> bits == 0; ++a) {
            for (uint64_t b = 0; b >> bits == 0; ++b) {
                CheckFieldElem(field, lowmod, a, b);
            }
        }
    } else {
        for (size_t i = 0; i < iters; ++i) {
            CheckFieldElem(field, lowmod, rng.RandBits(bits), rng.RandBits(bits));
        }
        // Edge cases.
        uint64_t max_elem = bits == 64 ? ~uint64_t{0} : (uint64_t{1} << bits) - 1;
        for (uint64_t a : {uint64_t{0}, uint64_t{1}, uint64_t{2}, max_elem}) {
            for (uint64_t b : {uint64_t{0}, uint64_t{1}, uint64_t{2}, max_elem}) {
                CheckFieldElem(field, lowmod, a, b);
            }
        }
    }

    // FromSeed must produce a nonzero element (it is used as the root-finding basis).
    for (size_t i = 0; i < 16; ++i) {
        UT_REQUIRE(field.FromSeed(rng.Rand64()) != 0);
        UT_REQUIRE(field.FromSeed(i) != 0);
    }
}

/** Serialize a random element sequence via BitWriter and read it back via BitReader. */
template<typename F>
void TestFieldSerialization(const F& field, TestRand& rng, size_t iters) {
    const int bits = field.Bits();
    for (size_t i = 0; i < iters; ++i) {
        size_t count = rng.RandRange(17);
        std::vector<typename F::Elem> elems(count);
        for (auto& e : elems) e = RandElem(rng, field);
        std::vector<unsigned char> buf((count * bits + 7) / 8, 0);
        BitWriter writer(buf.data());
        for (auto e : elems) field.Serialize(writer, e);
        writer.Flush();
        BitReader reader(buf.data());
        for (auto e : elems) UT_REQUIRE(field.Deserialize(reader) == e);
    }
}

/* ---------- Polynomial operation properties ---------- */

template<typename F>
std::vector<typename F::Elem> TraceModRefImpl(const std::vector<typename F::Elem>& tmod, typename F::Elem param, const F& field);

template<typename F>
void TestPolyOps(const F& field, TestRand& rng, size_t iters, size_t maxdeg) {
    typedef typename F::Elem Elem;
    for (size_t i = 0; i < iters; ++i) {
        // DivMod/PolyMod: for random val and monic mod, val == div*mod + rem,
        // deg(rem) < deg(mod), and PolyMod produces the same remainder.
        auto mod = RandMonicPoly(rng, field, 1 + rng.RandRange(maxdeg));
        auto val = RandPoly(rng, field, rng.RandRange(2 * maxdeg + 1));
        std::vector<Elem> quot, rem = val;
        DivMod(mod, rem, quot, field);
        UT_REQUIRE(Stripped(rem).size() < mod.size());
        auto recombined = PolyMulRef(quot, mod, field);
        recombined.resize(std::max(recombined.size(), rem.size()), 0);
        for (size_t j = 0; j < rem.size(); ++j) recombined[j] ^= rem[j];
        UT_REQUIRE(Stripped(recombined) == Stripped(val));
        auto rem2 = val;
        PolyMod(mod, rem2, field);
        UT_REQUIRE(rem2 == Stripped(rem));

        // Sqr against the naive product.
        auto poly = RandPoly(rng, field, rng.RandRange(maxdeg + 1));
        auto sqr = poly;
        Sqr(sqr, field);
        UT_REQUIRE(Stripped(sqr) == PolyMulRef(poly, poly, field));

        // MakeMonic: returns the inverse of the old leading coefficient and
        // scales the polynomial by it.
        auto nonmonic = RandPoly(rng, field, 1 + rng.RandRange(maxdeg));
        auto scaled = nonmonic;
        Elem inv = MakeMonic(scaled, field);
        UT_REQUIRE(scaled.back() == 1);
        Elem lead = nonmonic.back();
        if (lead == 1) {
            UT_REQUIRE(inv == 0); // Special value: no scaling was needed.
            UT_REQUIRE(scaled == nonmonic);
        } else {
            UT_REQUIRE(field.Mul(lead, inv) == 1);
            for (size_t j = 0; j < nonmonic.size(); ++j) {
                UT_REQUIRE(scaled[j] == field.Mul(nonmonic[j], inv));
            }
        }

        // GCD: for a = f*g and b = h*g, GCD(a, b) divides both inputs and is
        // divisible by g.
        auto f = RandPoly(rng, field, 1 + rng.RandRange(maxdeg / 2 + 1));
        auto g = RandPoly(rng, field, 1 + rng.RandRange(maxdeg / 2 + 1));
        auto h = RandPoly(rng, field, 1 + rng.RandRange(maxdeg / 2 + 1));
        auto a = PolyMulRef(f, g, field);
        auto b = PolyMulRef(h, g, field);
        auto a_copy = a, b_copy = b;
        GCD(a, b, field);
        auto d = a;
        UT_REQUIRE(!d.empty());
        MakeMonic(d, field);
        PolyReduceRef(a_copy, d, field);
        UT_REQUIRE(a_copy.empty());
        PolyReduceRef(b_copy, d, field);
        UT_REQUIRE(b_copy.empty());
        auto g_monic = g;
        MakeMonic(g_monic, field);
        auto d_copy = d;
        PolyReduceRef(d_copy, g_monic, field);
        UT_REQUIRE(d_copy.empty());

        // TraceMod against a naive reference: sum of (param*x)^(2^i) mod `mod`,
        // computed with schoolbook squarings and naive reduction. The
        // fast_tracemod_reducers API precomputes per modulus (degree >= 3);
        // this exercises the square-table reducer path.
        auto tmod = RandMonicPoly(rng, field, 4 + rng.RandRange(maxdeg));
        Elem param = RandNonzeroElem(rng, field);
        std::vector<Elem> trace;
        TraceMod<F> trace_mod(tmod, field);
        trace_mod.trace(trace, param);
        UT_REQUIRE(Stripped(trace) == TraceModRefImpl(tmod, param, field));
    }
}

/** Naive TraceMod reference: sum of (param*x)^(2^i) mod `mod`. */
template<typename F>
std::vector<typename F::Elem> TraceModRefImpl(const std::vector<typename F::Elem>& tmod, typename F::Elem param, const F& field) {
    typedef typename F::Elem Elem;
    std::vector<Elem> cur{0, param};
    PolyReduceRef(cur, tmod, field);
    std::vector<Elem> acc = cur;
    for (int bit = 1; bit < field.Bits(); ++bit) {
        cur = PolyMulRef(cur, cur, field);
        PolyReduceRef(cur, tmod, field);
        acc.resize(std::max(acc.size(), cur.size()), 0);
        for (size_t j = 0; j < cur.size(); ++j) acc[j] ^= cur[j];
    }
    return Stripped(acc);
}

/** Properties of the fast_tracemod_reducers building blocks: Karatsuba
 *  multiplication, the Newton inverse power series, and the reciprocal
 *  (large-degree) TraceMod path. */
template<typename F>
void TestTraceModReducers(const F& field, TestRand& rng, size_t iters) {
    typedef typename F::Elem Elem;
    for (size_t i = 0; i < iters; ++i) {
        // Karatsuba full/low products against schoolbook, at sizes crossing
        // TRACEMOD_POLYMUL_CUTOFF.
        auto a = RandPoly(rng, field, 1 + rng.RandRange(80));
        auto b = RandPoly(rng, field, 1 + rng.RandRange(80));
        std::vector<Elem> full;
        TraceModPolyMulFull(a, b, full, field);
        UT_REQUIRE(Stripped(full) == PolyMulRef(a, b, field));

        size_t n = 1 + rng.RandRange(80);
        std::vector<Elem> low;
        TraceModPolyMulLow(a, b, low, n, field);
        auto ref = PolyMulRef(a, b, field);
        ref.resize(n, 0);
        UT_REQUIRE(Stripped(low) == Stripped(ref));

        // f(x)*g(x)^2 mod x^n.
        auto fg2 = TraceModMulBySquareLow(a, b, n, field);
        auto ref2 = PolyMulRef(a, PolyMulRef(b, b, field), field);
        ref2.resize(n, 0);
        UT_REQUIRE(Stripped(fg2) == Stripped(ref2));

        // Newton inverse: f * (1/f mod x^k) == 1 mod x^k, for f with f(0)=1.
        auto f = RandPoly(rng, field, 1 + rng.RandRange(40));
        f[0] = 1;
        size_t k = 1 + rng.RandRange(40);
        auto inv = TraceModInvSeries(f, k, field);
        auto prod = PolyMulRef(f, inv, field);
        prod.resize(std::max(prod.size(), k), 0);
        for (size_t j = 0; j < k; ++j) UT_REQUIRE(prod[j] == (j == 0 ? 1 : 0));

        // Repeated use of one reducer object (results must not depend on
        // earlier calls), and SquareAndReduce against the naive reference.
        auto tmod = RandMonicPoly(rng, field, 4 + rng.RandRange(12));
        TraceMod<F> trace_mod(tmod, field);
        Elem p1 = RandNonzeroElem(rng, field), p2 = RandNonzeroElem(rng, field);
        std::vector<Elem> t1, t2, t1_again;
        trace_mod.trace(t1, p1);
        trace_mod.trace(t2, p2);
        trace_mod.trace(t1_again, p1);
        UT_REQUIRE(Stripped(t1) == TraceModRefImpl(tmod, p1, field));
        UT_REQUIRE(Stripped(t2) == TraceModRefImpl(tmod, p2, field));
        UT_REQUIRE(t1_again == t1);

        auto val = RandPoly(rng, field, 1 + rng.RandRange(tmod.size() - 1));
        auto reduced = val;
        trace_mod.SquareAndReduce(reduced);
        auto ref_sq = PolyMulRef(val, val, field);
        PolyReduceRef(ref_sq, tmod, field);
        UT_REQUIRE(Stripped(reduced) == ref_sq);
    }

    // The reciprocal reducer path (degree >= TRACEMOD_TABLE_CUTOFF) against
    // the naive reference; only for one cheap field, as the reference is
    // quadratic in the degree.
    if (field.Bits() == 11) {
        auto tmod = RandMonicPoly(rng, field, TRACEMOD_TABLE_CUTOFF + 2);
        Elem param = RandNonzeroElem(rng, field);
        TraceMod<F> trace_mod(tmod, field);
        std::vector<Elem> out;
        trace_mod.trace(out, param);
        UT_REQUIRE(Stripped(out) == TraceModRefImpl(tmod, param, field));

        // SquareAndReduce through the reciprocal path.
        auto val = RandPoly(rng, field, 1 + rng.RandRange(tmod.size() - 1));
        auto reduced = val;
        trace_mod.SquareAndReduce(reduced);
        auto ref_sq = PolyMulRef(val, val, field);
        PolyReduceRef(ref_sq, tmod, field);
        UT_REQUIRE(Stripped(reduced) == ref_sq);
    }
}

/* ---------- Syndrome and decode-stage properties ---------- */

template<typename F>
void TestSyndromes(const F& field, uint64_t lowmod, TestRand& rng, size_t iters, size_t maxcap) {
    typedef typename F::Elem Elem;
    const int bits = field.Bits();
    for (size_t i = 0; i < iters; ++i) {
        size_t cap = 1 + rng.RandRange(maxcap);

        // Odd syndromes of a single element m are m^1, m^3, m^5, ...
        Elem m = RandNonzeroElem(rng, field);
        std::vector<Elem> osyndromes(cap, 0);
        AddToOddSyndromes(osyndromes, m, field);
        for (size_t j = 0; j < cap; ++j) {
            UT_REQUIRE(uint64_t(osyndromes[j]) == RefPow(m, 2 * j + 1, bits, lowmod));
        }

        // Adding the same element again cancels out.
        AddToOddSyndromes(osyndromes, m, field);
        for (auto s : osyndromes) UT_REQUIRE(s == 0);

        // ReconstructAllSyndromes: with odd syndromes of a set, all_syndromes[k-1]
        // must equal the k'th power sum of the set, for all k (not just odd ones).
        size_t count = rng.RandRange(std::min<uint64_t>(cap, (uint64_t{1} << std::min(bits, 20)) - 1) + 1);
        auto elems = RandDistinctElems(rng, field, count);
        std::vector<Elem> osyn(cap, 0);
        for (auto e : elems) AddToOddSyndromes(osyn, e, field);
        auto all = ReconstructAllSyndromes(osyn, field);
        UT_REQUIRE(all.size() == 2 * cap);
        for (size_t k = 1; k <= all.size(); ++k) {
            uint64_t expect = 0;
            for (auto e : elems) expect ^= RefPow(e, k, bits, lowmod);
            UT_REQUIRE(uint64_t(all[k - 1]) == expect);
        }
    }
}

template<typename F>
void TestBerlekampMassey(const F& field, TestRand& rng, size_t iters, size_t maxroots) {
    typedef typename F::Elem Elem;
    const int bits = field.Bits();
    const uint64_t field_size = bits >= 20 ? (uint64_t{1} << 20) : (uint64_t{1} << bits);

    // All-zero syndromes: the minimal LFSR is the constant polynomial {1}.
    std::vector<Elem> zeros(8, 0);
    auto res = BerlekampMassey(zeros, 4, field);
    UT_REQUIRE(res == std::vector<Elem>{1});

    for (size_t i = 0; i < iters; ++i) {
        // Syndromes generated by a known root set: BM must return the
        // connection polynomial prod_i (1 + m_i*x).
        size_t count = 1 + rng.RandRange(std::min<uint64_t>(maxroots, field_size - 1));
        size_t cap = count + rng.RandRange(4);
        auto roots = RandDistinctElems(rng, field, count);
        std::vector<Elem> osyn(cap, 0);
        for (auto e : roots) AddToOddSyndromes(osyn, e, field);
        auto syndromes = ReconstructAllSyndromes(osyn, field);

        auto poly = BerlekampMassey(syndromes, cap, field);
        std::vector<Elem> expect{1};
        for (auto m : roots) {
            std::vector<Elem> factor{1, m};
            expect = PolyMulRef(expect, factor, field);
        }
        UT_REQUIRE(poly == expect);

        // Asking for a maximum degree below the true one must fail.
        auto rejected = BerlekampMassey(syndromes, count - 1, field);
        UT_REQUIRE(rejected.empty());

        // For random syndromes not derived from any element set: if BM
        // succeeds, its output must satisfy the defining LFSR recurrence over
        // the whole sequence. Note the input must have the Frobenius structure
        // (s_{2k+2} = s_{k+1}^2) that ReconstructAllSyndromes produces; that is
        // the domain BerlekampMassey is used on (and asserts internally).
        std::vector<Elem> osyn_garbage(1 + rng.RandRange(maxroots));
        for (auto& s : osyn_garbage) s = RandElem(rng, field);
        auto garbage = ReconstructAllSyndromes(osyn_garbage, field);
        auto conn = BerlekampMassey(garbage, garbage.size() / 2, field);
        if (!conn.empty()) {
            UT_REQUIRE(conn[0] == 1);
            for (size_t n = conn.size() - 1; n < garbage.size(); ++n) {
                Elem acc = 0;
                for (size_t j = 0; j < conn.size(); ++j) {
                    acc ^= field.Mul(conn[j], garbage[n - j]);
                }
                UT_REQUIRE(acc == 0);
            }
        }
    }
}

template<typename F>
void TestFindRoots(const F& field, uint64_t lowmod, TestRand& rng, size_t iters, size_t maxroots) {
    typedef typename F::Elem Elem;
    const int bits = field.Bits();
    const uint64_t field_size = bits >= 20 ? (uint64_t{1} << 20) : (uint64_t{1} << bits);

    // A constant polynomial has no roots (and this is not a failure).
    UT_REQUIRE(FindRoots(std::vector<Elem>{1}, Elem(1), field).empty());

    for (size_t i = 0; i < iters; ++i) {
        size_t count = 1 + rng.RandRange(std::min<uint64_t>(maxroots, field_size - 1));
        auto roots = RandDistinctElems(rng, field, count);
        auto poly = PolyFromRootsRef(roots, field);
        std::sort(roots.begin(), roots.end());

        // Roots of prod (x + m_i) are exactly {m_i}, for any nonzero basis.
        for (Elem basis : {Elem(1), RandNonzeroElem(rng, field)}) {
            auto found = FindRoots(poly, basis, field);
            std::sort(found.begin(), found.end());
            UT_REQUIRE(found == roots);
        }

        if (bits > 1) {
            // Multiply in an irreducible quadratic x^2 + x + c (irreducible over
            // GF(2^bits) iff Tr(c) == 1): the polynomial is still square-free but
            // no longer fully factorizable, so FindRoots must return nothing.
            Elem c;
            do {
                c = RandNonzeroElem(rng, field);
            } while (RefTrace(c, bits, lowmod) != 1);
            std::vector<Elem> quad{c, 1, 1};
            auto not_factorizable = PolyMulRef(poly, quad, field);
            UT_REQUIRE(FindRoots(not_factorizable, Elem(1), field).empty());
        }
    }
}

/** Reference-free decode-stage roundtrip, cheap enough for large capacities:
 *  the syndromes of `count` random distinct roots must produce a
 *  BerlekampMassey polynomial of degree exactly `count` whose reversal's
 *  roots are exactly the input set. This is the whole decode pipeline at
 *  full capacity, without the quadratic naive references. */
template<typename F>
void TestDecodeStagesLarge(const F& field, TestRand& rng, size_t count) {
    typedef typename F::Elem Elem;
    const int bits = field.Bits();
    const uint64_t field_size = bits >= 20 ? (uint64_t{1} << 20) : (uint64_t{1} << bits);
    count = std::min<uint64_t>(count, field_size - 1);
    auto roots = RandDistinctElems(rng, field, count);
    std::vector<Elem> osyn(count, 0);
    for (auto e : roots) AddToOddSyndromes(osyn, e, field);
    auto syndromes = ReconstructAllSyndromes(osyn, field);
    auto poly = BerlekampMassey(syndromes, count, field);
    UT_REQUIRE(poly.size() == count + 1);
    std::reverse(poly.begin(), poly.end());
    auto found = FindRoots(poly, field.FromSeed(rng.Rand64()), field);
    std::sort(found.begin(), found.end());
    std::sort(roots.begin(), roots.end());
    UT_REQUIRE(found == roots);
}

/** Run every property group against one field instantiation. */
template<typename F>
void RunAllFieldTests(uint64_t seed_offset) {
    F field;
    const uint64_t lowmod = REF_LOWMOD[field.Bits()];
    TestRand rng(g_test_seed + seed_offset);
    TestFieldOps(field, lowmod, rng, 4000);
    TestFieldSerialization(field, rng, 128);
    TestPolyOps(field, rng, 256, 12);
    TestTraceModReducers(field, rng, 24);
    TestSyndromes(field, lowmod, rng, 128, 12);
    TestBerlekampMassey(field, rng, 128, 10);
    TestFindRoots(field, lowmod, rng, 64, 8);
    TestDecodeStagesLarge(field, rng, 256);
}

} // namespace ut

/** Register a Boost test case running all property groups for one field type.
 *  The case name doubles as a deterministic per-field seed offset. */
#define MINISKETCH_FIELD_TEST(name, FieldType) \
    BOOST_AUTO_TEST_CASE(name) { \
        ut::RunAllFieldTests<FieldType>(__LINE__); \
    }

#endif
