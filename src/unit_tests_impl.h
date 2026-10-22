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
        // All pairs of basis elements x^i, x^j. Every multiplier here (shift-and-add,
        // and carry-less multiply followed by the fold reduction) is GF(2)-bilinear,
        // so agreeing with the reference on a basis means agreeing everywhere: this
        // checks Mul, and the reduction recipes, completely (claims C47, C49).
        for (int i = 0; i < bits; ++i) {
            for (int j = 0; j < bits; ++j) {
                CheckFieldElem(field, lowmod, uint64_t{1} << i, uint64_t{1} << j);
                if constexpr (F::LAZY_REDUCE) {
                    typedef typename F::Elem Elem;
                    UT_REQUIRE(field.Reduce(field.MulWide(Elem(uint64_t{1} << i), Elem(uint64_t{1} << j))) ==
                               field.Mul(Elem(uint64_t{1} << i), Elem(uint64_t{1} << j)));
                }
            }
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

    // One scratch pool shared across mixed-size scratch-threaded calls, as the
    // reducer does internally: results must not depend on what earlier calls
    // left behind in the slots. For fields with the additive-FFT tier, sizes
    // extend beyond the FFT cutoff so the same pool serves both tiers.
    {
        const size_t maxsize = TraceModFFTEligible(field, 400, 400) ? 700 : 80;
        TraceModScratch<F> scratch;
        for (int rep = 0; rep < 4; ++rep) {
            auto a = RandPoly(rng, field, 1 + rng.RandRange(maxsize));
            auto b = RandPoly(rng, field, 1 + rng.RandRange(maxsize));
            std::vector<Elem> full;
            TraceModPolyMulFull(a, b, full, field, scratch, 0);
            UT_REQUIRE(Stripped(full) == PolyMulRef(a, b, field));
            size_t n = 1 + rng.RandRange(maxsize);
            std::vector<Elem> low;
            TraceModPolyMulLow(a, b, low, n, field, scratch, 0);
            auto ref = PolyMulRef(a, b, field);
            ref.resize(n, 0);
            UT_REQUIRE(Stripped(low) == Stripped(ref));
            std::vector<Elem> fg2;
            TraceModMulBySquareLow(a, b, fg2, n, field, scratch, 0);
            auto ref2 = PolyMulRef(a, PolyMulRef(b, b, field), field);
            ref2.resize(n, 0);
            UT_REQUIRE(Stripped(fg2) == Stripped(ref2));
        }
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

        // A repeated trace on the same object must be unaffected by the
        // buffers earlier calls left behind.
        std::vector<Elem> out2;
        trace_mod.trace(out2, param);
        UT_REQUIRE(out2 == out);
    }

    // The same reciprocal-path property once through the additive-FFT tier,
    // just past the (lower) table cutoff of FFT-capable fields; the 11-bit
    // case above keeps covering the pure Karatsuba reciprocal path.
    if (field.Bits() == 32) {
        auto tmod = RandMonicPoly(rng, field, TRACEMOD_FFT_TABLE_CUTOFF + 2);
        Elem param = RandNonzeroElem(rng, field);
        TraceMod<F> trace_mod(tmod, field);
        std::vector<Elem> out;
        trace_mod.trace(out, param);
        UT_REQUIRE(Stripped(out) == TraceModRefImpl(tmod, param, field));

        auto val = RandPoly(rng, field, 1 + rng.RandRange(tmod.size() - 1));
        auto reduced = val;
        trace_mod.SquareAndReduce(reduced);
        auto ref_sq = PolyMulRef(val, val, field);
        PolyReduceRef(ref_sq, tmod, field);
        UT_REQUIRE(Stripped(reduced) == ref_sq);

        std::vector<Elem> out2;
        trace_mod.trace(out2, param);
        UT_REQUIRE(out2 == out);
    }

    // The cached-transform (prepared-operand) reducer path. At the degree
    // just past the table cutoff above, the fixed node transform is twice
    // the dynamically chosen size and the cached path correctly stays off;
    // at 5/4 of the cutoff the full-size reductions land in the upper half
    // of the padded transform and the prepared operands engage.
    if (field.Bits() == 32) {
        auto tmod = RandMonicPoly(rng, field, 5 * TRACEMOD_FFT_TABLE_CUTOFF / 4);
        Elem param = RandNonzeroElem(rng, field);
        TraceMod<F> trace_mod(tmod, field);
        std::vector<Elem> out;
        trace_mod.trace(out, param);
        UT_REQUIRE(Stripped(out) == TraceModRefImpl(tmod, param, field));

        // Full-length SquareAndReduce so both prepared products engage.
        auto val = RandPoly(rng, field, tmod.size() - 1);
        auto reduced = val;
        trace_mod.SquareAndReduce(reduced);
        auto ref_sq = PolyMulRef(val, val, field);
        PolyReduceRef(ref_sq, tmod, field);
        UT_REQUIRE(Stripped(reduced) == ref_sq);

        std::vector<Elem> out2;
        trace_mod.trace(out2, param);
        UT_REQUIRE(out2 == out);
    }
}

/** Properties of the subquadratic division and gcd (FastDivMod, FastGCD):
 *  equivalence with the schoolbook DivMod and the quadratic GCD. A small
 *  explicit cutoff is passed so the reciprocal and half-gcd machinery is
 *  exercised with real recursion depth at reference-affordable degrees, on
 *  every field (its correctness does not depend on the field); the
 *  production cutoffs from FastDivModCutoff/FastGCDCutoff are far larger
 *  and are exercised by the decode-level tests. */
template<typename F>
void TestSubquadraticGCD(const F& field, TestRand& rng, size_t iters) {
    typedef typename F::Elem Elem;
    const size_t maxdeg = 320;
    const size_t cutoff = 24;

    for (size_t i = 0; i < iters; ++i) {
        // FastDivMod == DivMod (quotient and raw remainder buffer alike),
        // for random numerators and for exactly-divisible ones.
        size_t dlen = 2 + rng.RandRange(maxdeg);
        auto mod = RandMonicPoly(rng, field, dlen);
        auto val = RandPoly(rng, field, dlen + rng.RandRange(maxdeg));
        std::vector<Elem> q_ref, rem_ref = val, q_fast, rem_fast = val;
        DivMod(mod, rem_ref, q_ref, field);
        FastDivMod(mod, rem_fast, q_fast, field, cutoff);
        UT_REQUIRE(q_fast == q_ref);
        UT_REQUIRE(rem_fast == rem_ref);

        auto quot_in = RandPoly(rng, field, 1 + rng.RandRange(maxdeg));
        auto exact = PolyMulRef(quot_in, mod, field);
        if (!exact.empty()) {
            std::vector<Elem> q2, rem2 = exact;
            FastDivMod(mod, rem2, q2, field, cutoff);
            UT_REQUIRE(Stripped(rem2).empty());
            UT_REQUIRE(Stripped(q2) == Stripped(quot_in));
        }

        // FastGCD == GCD for operands with a planted common factor.
        auto f = RandPoly(rng, field, 1 + rng.RandRange(maxdeg));
        auto g = RandPoly(rng, field, 1 + rng.RandRange(maxdeg / 2 + 1));
        auto h = RandPoly(rng, field, 1 + rng.RandRange(maxdeg));
        auto a = PolyMulRef(f, g, field);
        auto b = PolyMulRef(h, g, field);
        auto a_ref = a, b_ref = b, a_fast = a, b_fast = b;
        GCD(a_ref, b_ref, field);
        FastGCD(a_fast, b_fast, field, cutoff);
        UT_REQUIRE(!a_ref.empty() && !a_fast.empty());
        MakeMonic(a_ref, field);
        MakeMonic(a_fast, field);
        UT_REQUIRE(a_fast == a_ref);

        // A remainder sequence built backwards with occasional large
        // quotients: exercises the non-normal paths of the half-gcd (the
        // internal division steps and the d2 recursion).
        std::vector<Elem> r1 = RandPoly(rng, field, 1 + rng.RandRange(8));
        std::vector<Elem> r0 = RandPoly(rng, field, r1.size() + 1 + rng.RandRange(8));
        while (r0.size() <= maxdeg && !r0.empty()) {
            size_t qdeg = 1 + (rng.RandRange(8) == 0 ? rng.RandRange(20) : 0);
            auto q = RandPoly(rng, field, qdeg + 1);
            if (q.empty() || q.back() == 0) continue;
            auto next = PolyMulRef(q, r0, field);
            for (size_t j = 0; j < r1.size(); ++j) next[j] ^= r1[j];
            r1.swap(r0);
            r0.swap(next);
        }
        if (!r0.empty() && !r1.empty()) {
            auto na_ref = r0, nb_ref = r1, na_fast = r0, nb_fast = r1;
            GCD(na_ref, nb_ref, field);
            FastGCD(na_fast, nb_fast, field, cutoff);
            UT_REQUIRE(!na_ref.empty() && !na_fast.empty());
            MakeMonic(na_ref, field);
            MakeMonic(na_fast, field);
            UT_REQUIRE(na_fast == na_ref);
        }
    }

    // Edge cases: empty side, identical operands, one dividing the other,
    // and coprime pairs reducing to a constant.
    auto a = RandPoly(rng, field, 150);
    std::vector<Elem> empty_b;
    auto a_copy = a;
    FastGCD(a_copy, empty_b, field, cutoff);
    UT_REQUIRE(Stripped(a_copy) == Stripped(a));

    auto same1 = a, same2 = a;
    FastGCD(same1, same2, field, cutoff);
    auto want = a;
    UT_REQUIRE(!want.empty());
    MakeMonic(want, field);
    MakeMonic(same1, field);
    UT_REQUIRE(same1 == want);

    auto factor = RandMonicPoly(rng, field, 140);
    auto multiple = PolyMulRef(RandPoly(rng, field, 60), factor, field);
    if (!multiple.empty()) {
        auto ma = multiple, mb = factor;
        FastGCD(ma, mb, field, cutoff);
        UT_REQUIRE(!ma.empty());
        MakeMonic(ma, field);
        UT_REQUIRE(ma == factor);
    }
}

/** Properties of the additive (Cantor-basis) FFT engine and its
 *  multiplication tier. Engine checks run only for eligible fields (degree a
 *  power of two, >= 16); other fields verify they are gated out. */
template<typename F>
void TestAdditiveFFT(const F& field, TestRand& rng, size_t iters) {
    typedef typename F::Elem Elem;
    const int bits = field.Bits();
    const bool eligible_field = bits >= 16 && (bits & (bits - 1)) == 0;

    // Gating: never below the product-length cutoff (probe just under it),
    // never for non-power-of-2 degrees, never with an unbalanced (short)
    // operand.
    UT_REQUIRE(!TraceModFFTEligible(field, TRACEMOD_FFT_CUTOFF / 2, TRACEMOD_FFT_CUTOFF / 2));
    UT_REQUIRE(!TraceModFFTEligible(field, 1000, TRACEMOD_POLYMUL_CUTOFF));
    UT_REQUIRE(TraceModFFTEligible(field, 400, 400) == eligible_field);
    if (!eligible_field) return;

    const int max_t = std::min(bits, 12);
    AdditiveFFT<F> fft;
    fft.Extend(field, max_t);

    // Cantor basis: beta[0] == 1 and Sqr(beta[j]) ^ beta[j] == beta[j-1].
    const auto& beta = fft.Basis();
    UT_REQUIRE((int)beta.size() == max_t);
    UT_REQUIRE(beta[0] == Elem(1));
    for (size_t j = 1; j < beta.size(); ++j) {
        UT_REQUIRE((field.Sqr(beta[j]) ^ beta[j]) == beta[j - 1]);
    }

    // Anchor property Z_{W_q}(beta[q]) == 1, with the vanishing polynomial
    // evaluated independently (Lucas' theorem + repeated squaring), not via
    // the engine's tables.
    for (int q = 1; q < max_t; ++q) {
        Elem acc = 0, pw = beta[q];
        for (int j = 0; j <= q; ++j) {
            if ((q & j) == j) acc ^= pw;
            pw = field.Sqr(pw);
        }
        UT_REQUIRE(acc == Elem(1));
    }

    // Roundtrip at every dimension.
    for (int m = 1; m <= max_t; ++m) {
        std::vector<Elem> f(size_t(1) << m);
        for (auto& e : f) e = RandElem(rng, field);
        auto orig = f;
        fft.FFT(f, m, field);
        if (m > 1) UT_REQUIRE(f != orig); // Not the identity map (paranoia).
        fft.IFFT(f, m, field);
        UT_REQUIRE(f == orig);
        fft.IFFT(f, m, field);
        fft.FFT(f, m, field);
        UT_REQUIRE(f == orig);
    }

    // Direct evaluation: FFT output j must equal the Horner evaluation at
    // eta_j = XOR of beta[b] over the set bits b of j. This breaks the
    // fft/ifft symmetry that roundtrip checks cannot see through.
    for (int m = 1; m <= std::min(max_t, 6); ++m) {
        size_t n_points = size_t(1) << m;
        std::vector<Elem> f(n_points);
        for (auto& e : f) e = RandElem(rng, field);
        auto eval = f;
        fft.FFT(eval, m, field);
        for (size_t j = 0; j < n_points; ++j) {
            Elem x = 0;
            for (int b = 0; b < m; ++b) {
                if ((j >> b) & 1) x ^= beta[b];
            }
            Elem horner = 0;
            for (size_t i = f.size(); i > 0; --i) horner = field.Mul(horner, x) ^ f[i - 1];
            UT_REQUIRE(eval[j] == horner);
        }
    }

    // Coset transforms: FFT on block b evaluates at eta_(b*2^m + j) (Horner
    // oracle); the novel-basis path (ToNovel, then NovelFFT) gives the same
    // values; every inverse is exact; NovelDerivative matches the formal
    // derivative of the monomial form.
    for (int m = 1; m + 2 <= std::min(max_t, 8); ++m) {
        size_t n_points = size_t(1) << m;
        for (size_t block = 0; block < 4; ++block) {
            std::vector<Elem> f(n_points);
            for (auto& e : f) e = RandElem(rng, field);
            auto eval = f;
            fft.FFT(eval, m, field, block);
            for (size_t j = 0; j < n_points; ++j) {
                Elem x = 0;
                for (int b = 0; b < max_t; ++b) {
                    if (((block * n_points + j) >> b) & 1) x ^= beta[b];
                }
                Elem horner = 0;
                for (size_t i = f.size(); i > 0; --i) horner = field.Mul(horner, x) ^ f[i - 1];
                UT_REQUIRE(eval[j] == horner);
            }
            auto back = eval;
            fft.IFFT(back, m, field, block);
            UT_REQUIRE(back == f);
            auto nov = f;
            fft.ToNovel(nov, m);
            auto nov_coef = nov;
            fft.NovelFFT(nov, m, field, block);
            UT_REQUIRE(nov == eval);
            fft.NovelIFFT(nov, m, field, block);
            UT_REQUIRE(nov == nov_coef);
            fft.FromNovel(nov, m);
            UT_REQUIRE(nov == f);
            std::vector<Elem> deriv(n_points, 0), nov_deriv;
            for (size_t i = 1; i < n_points; i += 2) deriv[i - 1] = f[i];
            fft.ToNovel(deriv, m);
            AdditiveFFT<F>::NovelDerivative(nov_coef, m, nov_deriv);
            UT_REQUIRE(nov_deriv == deriv);
        }
    }

    // MulFull/MulLow against the schoolbook oracle, one engine across mixed
    // sizes (exercises Extend growth and buffer reuse). Boundary product
    // lengths around powers of two and around the cutoff, then random sizes.
    std::vector<std::pair<size_t, size_t>> shapes;
    for (size_t prod : {size_t{63}, size_t{64}, size_t{65}, size_t{127}, size_t{128}, size_t{129},
                        TRACEMOD_FFT_CUTOFF - 1, TRACEMOD_FFT_CUTOFF, TRACEMOD_FFT_CUTOFF + 1}) {
        size_t la = prod / 2 + 1, lb = prod - la + 1;
        shapes.emplace_back(la, lb);
    }
    for (size_t i = 0; i < iters; ++i) {
        shapes.emplace_back(1 + rng.RandRange(300), 1 + rng.RandRange(300));
    }
    for (const auto& shape : shapes) {
        auto a = RandPoly(rng, field, shape.first);
        auto b = RandPoly(rng, field, shape.second);
        auto ref = PolyMulRef(a, b, field);
        std::vector<Elem> full;
        fft.MulFull(a, b, full, field);
        UT_REQUIRE(full == ref);
        size_t n = 1 + rng.RandRange(shape.first + shape.second + 8);
        std::vector<Elem> low;
        fft.MulLow(a, b, low, n, field);
        UT_REQUIRE(low.size() == n);
        auto ref_low = ref;
        ref_low.resize(n, 0);
        UT_REQUIRE(low == ref_low);
    }

    // Edge cases: empty operands, squaring shape (a == b).
    std::vector<Elem> empty_in, out;
    auto a = RandPoly(rng, field, 40);
    fft.MulFull(a, empty_in, out, field);
    UT_REQUIRE(out.empty());
    fft.MulLow(empty_in, a, out, 7, field);
    UT_REQUIRE(out == std::vector<Elem>(7, 0));
    fft.MulFull(a, a, out, field);
    UT_REQUIRE(out == PolyMulRef(a, a, field));

    // Prepared-operand multiplication: Prepare + MulLowPrepared must match
    // MulLow for every truncation length, both at the minimal transform size
    // and at one above it (the oversized case mirrors TraceMod's fixed
    // per-node size), reusing one prepared table across several multiplicands.
    for (size_t i = 0; i < iters; ++i) {
        auto b = RandPoly(rng, field, 1 + rng.RandRange(200));
        size_t max_la = 1 + rng.RandRange(200);
        size_t max_prod = max_la + b.size() - 1;
        int m_min = 1;
        while ((size_t(1) << m_min) < max_prod) ++m_min;
        for (int m = m_min; m <= std::min(m_min + 1, max_t); ++m) {
            std::vector<Elem> b_fft;
            fft.Prepare(b, m, field, b_fft);
            for (int j = 0; j < 3; ++j) {
                auto a2 = RandPoly(rng, field, 1 + rng.RandRange(max_la));
                size_t n = 1 + rng.RandRange(a2.size() + b.size() + 8);
                std::vector<Elem> low, low_ref;
                fft.MulLowPrepared(a2, b_fft, b.size(), m, low, n, field);
                fft.MulLow(a2, b, low_ref, n, field);
                UT_REQUIRE(low == low_ref);
            }
        }
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

/** BerlekampMassey (even-step skipping, early termination with tail
 *  verification) must return exactly what the textbook iteration returns,
 *  including the max_degree failure, on Frobenius-structured inputs: the
 *  sketches of sets far below capacity (the early-termination regime, both
 *  verification tiers), the same with one syndrome perturbed at a random
 *  position (verification must fail, or the iteration must not stop early),
 *  and unstructured random odd syndromes. */
template<typename F>
void TestBerlekampMasseyVsRef(const F& field, TestRand& rng, size_t iters, size_t maxcap) {
    typedef typename F::Elem Elem;
    const int bits = field.Bits();
    const uint64_t field_size = bits >= 20 ? (uint64_t{1} << 20) : (uint64_t{1} << bits);
    for (size_t i = 0; i < iters; ++i) {
        size_t cap = 1 + rng.RandRange(maxcap);
        std::vector<Elem> osyn(cap, 0);
        int kind = rng.RandRange(3);
        if (kind == 2) {
            for (auto& s : osyn) s = RandElem(rng, field);
        } else {
            // Log-uniform count, so small counts (long verified tails) dominate.
            size_t count = (size_t(1) << rng.RandRange(8)) - 1;
            count = std::min<uint64_t>(rng.RandRange(count + 1), field_size - 1);
            for (auto e : RandDistinctElems(rng, field, count)) AddToOddSyndromes(osyn, e, field);
            if (kind == 1) osyn[rng.RandRange(cap)] ^= RandNonzeroElem(rng, field);
        }
        auto syndromes = ReconstructAllSyndromes(osyn, field);
        size_t max_degree = rng.RandRange(2) ? cap : rng.RandRange(cap + 1);
        auto got = BerlekampMassey(syndromes, max_degree, field);
        auto want = BerlekampMasseyRef(syndromes, max_degree, field);
        UT_REQUIRE(got == want);
    }
    // The tail check itself (whose failure would only cost speed, not
    // change BM's output): the odd syndromes of a root set satisfy the
    // recurrence with the squared connection polynomial past the prefix,
    // and a perturbed tail syndrome is rejected - for lengths on both sides
    // of the FFT threshold.
    // Tails span several blocks of each tier: the Karatsuba tier's blocks
    // are ~L+1 positions (interior blocks read u in place), the FFT tier's
    // up to ~3L.
    const size_t fft_min = F::LAZY_REDUCE ? BM_VERIFY_FFT_MIN_DEGREE : BM_VERIFY_FFT_MIN_DEGREE_TABLE;
    for (size_t count : {size_t(3), size_t(17), size_t(40), size_t(BM_VERIFY_KARA_MIN_DEGREE + 6), size_t(100), fft_min + 5}) {
        if (count + 2 >= field_size) continue;
        size_t cap = count + 40 + rng.RandRange(40) + (count >= BM_VERIFY_KARA_MIN_DEGREE ? 4 * count : 0) + (count >= fft_min ? 4 * count : 0);
        auto roots = RandDistinctElems(rng, field, count);
        std::vector<Elem> osyn(cap, 0);
        for (auto e : roots) AddToOddSyndromes(osyn, e, field);
        std::vector<Elem> conn{1};
        for (auto m : roots) {
            std::vector<Elem> factor{1, m};
            conn = PolyMulRef(conn, factor, field);
        }
        size_t j0 = count + 1 + rng.RandRange(8);
        UT_REQUIRE(BerlekampMasseyVerifyTail(conn, osyn, j0, cap, field));
        osyn[j0 + rng.RandRange(cap - j0)] ^= RandNonzeroElem(rng, field);
        UT_REQUIRE(!BerlekampMasseyVerifyTail(conn, osyn, j0, cap, field));
    }
    // A few LFSR lengths at and above the FFT verification threshold (the
    // random counts above stay below it), with and without a perturbation.
    if (field_size - 1 < BM_VERIFY_FFT_MIN_DEGREE + 8) return;
    for (int i = 0; i < 4; ++i) {
        size_t count = BM_VERIFY_FFT_MIN_DEGREE + rng.RandRange(8);
        size_t cap = count + BM_VERIFY_MIN_TAIL + rng.RandRange(64);
        std::vector<Elem> osyn(cap, 0);
        for (auto e : RandDistinctElems(rng, field, count)) AddToOddSyndromes(osyn, e, field);
        if (i & 1) osyn[count + rng.RandRange(cap - count)] ^= RandNonzeroElem(rng, field);
        auto syndromes = ReconstructAllSyndromes(osyn, field);
        UT_REQUIRE(BerlekampMassey(syndromes, cap, field) == BerlekampMasseyRef(syndromes, cap, field));
    }
}

/** Exhaustive version of the BM oracle comparison for tiny fields: every
 *  odd-syndrome vector up to the largest capacity with at most ~2^12
 *  vectors, every max_degree. (Too short for the tail check, which the
 *  random cases cover; this pins the even-step skip and the degree-limit
 *  logic completely.) */
template<typename F>
void TestBerlekampMasseyExhaustive(const F& field) {
    typedef typename F::Elem Elem;
    const int bits = field.Bits();
    if (bits > 4) return;
    const uint64_t q = uint64_t{1} << bits;
    for (size_t cap = 1; ; ++cap) {
        uint64_t total = 1;
        for (size_t i = 0; i < cap; ++i) total *= q;
        if (total > 4096) break;
        std::vector<Elem> osyn(cap);
        for (uint64_t v = 0; v < total; ++v) {
            uint64_t t = v;
            for (size_t i = 0; i < cap; ++i) {
                osyn[i] = Elem(t % q);
                t /= q;
            }
            auto syndromes = ReconstructAllSyndromes(osyn, field);
            for (size_t max_degree = 0; max_degree <= cap; ++max_degree) {
                UT_REQUIRE(BerlekampMassey(syndromes, max_degree, field) == BerlekampMasseyRef(syndromes, max_degree, field));
            }
        }
    }
}

/** FindRoots on degrees that take the Frobenius-orbit path: square-free
 *  split polynomials must yield their roots; a non-split factor or a
 *  repeated root must make it fail. */
template<typename F>
void TestFindRootsOrbit(const F& field, uint64_t lowmod, TestRand& rng, size_t iters, size_t maxroots) {
    typedef typename F::Elem Elem;
    const int bits = field.Bits();
    const uint64_t field_size = bits >= 20 ? (uint64_t{1} << 20) : (uint64_t{1} << bits);
    if (field_size - 1 < ORBIT_MIN_DEGREE + 2) return;
    for (size_t i = 0; i < iters; ++i) {
        size_t count = ORBIT_MIN_DEGREE + rng.RandRange(std::min<uint64_t>(maxroots, field_size - 2 - ORBIT_MIN_DEGREE));
        auto roots = RandDistinctElems(rng, field, count);
        auto poly = PolyFromRootsRef(roots, field);
        std::sort(roots.begin(), roots.end());
        Elem basis = RandNonzeroElem(rng, field);
        auto found = FindRoots(poly, basis, field);
        std::sort(found.begin(), found.end());
        UT_REQUIRE(found == roots);

        Elem c;
        do {
            c = RandNonzeroElem(rng, field);
        } while (RefTrace(c, bits, lowmod) != 1);
        std::vector<Elem> quad{c, 1, 1};
        UT_REQUIRE(FindRoots(PolyMulRef(poly, quad, field), basis, field).empty());

        std::vector<Elem> dup{roots[rng.RandRange(roots.size())], 1};
        UT_REQUIRE(FindRoots(PolyMulRef(poly, dup, field), basis, field).empty());

        // The split parameters must be GF(2)-linearly independent (rank
        // Bits over the bit representation), and each orbit trace must be
        // Tr(beta_k x) mod poly: GF(2)-valued at every root.
        std::vector<std::vector<Elem>> traces;
        std::vector<Elem> betas;
        UT_REQUIRE(ComputeTraceOrbit(poly, basis, 4, field, traces, betas));
        UT_REQUIRE(betas.size() == size_t(bits));
        uint64_t pivot[64] = {0};
        for (auto b : betas) {
            uint64_t v = uint64_t(b);
            int bit = 63;
            for (; bit >= 0; --bit) {
                if (!((v >> bit) & 1)) continue;
                if (!pivot[bit]) break;
                v ^= pivot[bit];
            }
            UT_REQUIRE(bit >= 0); // reduced to zero: dependent
            pivot[bit] = v;
        }
        for (size_t k = 0; k < traces.size(); ++k) {
            for (auto root : roots) {
                Elem val = 0, pw = 1;
                for (auto c : traces[k]) {
                    val ^= field.Mul(c, pw);
                    pw = field.Mul(pw, root);
                }
                UT_REQUIRE(val == 0 || val == 1);
                UT_REQUIRE(val == Elem(RefTrace(uint64_t(field.Mul(betas[k], root)), bits, lowmod)));
            }
        }
    }
}

/** FindRoots on degrees that take the coset Montgomery path (root degree
 *  >= MONT_ORBIT_MIN_DEGREE on FFT-capable fields with a cheap Multiplier):
 *  degree exactly a power of two (deg f == N), roots planted on the
 *  candidate cosets (forcing the later candidates, then the fallback when
 *  all four are hit), a degree deep enough for the second relative trace
 *  and the normal-form carry, and the failure cases. */
template<typename F>
void TestFindRootsMontgomery(const F& field, uint64_t lowmod, TestRand& rng) {
    typedef typename F::Elem Elem;
    const int bits = field.Bits();
    if (!TraceModFFTCapable(field) || !TraceModCheapScalarMul<F>() || bits % 8 != 0) return;
    AdditiveFFT<F> fft;
    fft.Extend(field, std::min(bits, 13));
    const auto& beta = fft.Basis();
    // A point of coset `block` of W_m: eta_(block * 2^m + j).
    auto coset_point = [&](int m, size_t block) {
        size_t idx = block * (size_t(1) << m) + rng.RandRange(uint64_t(1) << m);
        Elem x = 0;
        for (int k = 0; (idx >> k) != 0; ++k) {
            if ((idx >> k) & 1) x ^= beta[k];
        }
        return x;
    };
    struct Case { size_t degree; int planted; };
    for (Case c : {Case{MONT_ORBIT_MIN_DEGREE, 0}, Case{256, 0}, Case{300, 0}, Case{300, 3}, Case{300, 4}, Case{256, 4}, Case{1100, 0}}) {
        if (c.degree < MONT_ORBIT_MIN_DEGREE) continue;
        if (c.degree + 2 > (size_t(1) << std::min(bits, 20))) continue;
        int m = 1;
        while ((size_t(1) << m) < c.degree) ++m;
        if (m + 2 > bits) continue;
        std::vector<Elem> roots;
        // Planted roots on the candidate cosets (blocks 1, 2, 3, then 0).
        for (int k = 0; k < c.planted; ++k) {
            Elem x;
            do {
                x = coset_point(m, (k + 1) % 4);
            } while (x == 0 || std::find(roots.begin(), roots.end(), x) != roots.end());
            roots.push_back(x);
        }
        for (auto x : RandDistinctElems(rng, field, c.degree)) {
            if (roots.size() == c.degree) break;
            if (std::find(roots.begin(), roots.end(), x) == roots.end()) roots.push_back(x);
        }
        UT_REQUIRE(roots.size() == c.degree);
        auto poly = PolyFromRootsRef(roots, field);
        std::sort(roots.begin(), roots.end());
        Elem basis = RandNonzeroElem(rng, field);
        auto found = FindRoots(poly, basis, field);
        std::sort(found.begin(), found.end());
        UT_REQUIRE(found == roots);

        // With the roots' odd power sums (the decoder's syndromes) the large
        // nodes split by power sums (PowerSumSplit); with as few as 2n - 1
        // of them (a full sketch), and with corrupted ones, which the
        // division check must reject in favour of the gcd.
        for (size_t count : {c.degree + 8, c.degree}) {
            std::vector<Elem> osyn(count, 0);
            for (auto r : roots) {
                Elem r2 = field.Sqr(r), pw = r;
                for (auto& o : osyn) {
                    o ^= pw;
                    pw = field.Mul(pw, r2);
                }
            }
            found = FindRoots(poly, basis, field, &osyn);
            std::sort(found.begin(), found.end());
            UT_REQUIRE(found == roots);
            osyn[rng.RandRange(osyn.size())] ^= RandNonzeroElem(rng, field);
            found = FindRoots(poly, basis, field, &osyn);
            std::sort(found.begin(), found.end());
            UT_REQUIRE(found == roots);
        }

        if (c.degree > 300) continue;
        // ComputeTraceOrbit without a tree: every trace from the orbit's
        // carry must be Tr(beta_k x) mod poly.
        std::vector<std::vector<Elem>> traces;
        std::vector<Elem> betas;
        UT_REQUIRE(ComputeTraceOrbit(poly, basis, 8, field, traces, betas));
        UT_REQUIRE(traces.size() == 8);
        for (size_t k = 0; k < traces.size(); ++k) {
            for (size_t r = 0; r < 16; ++r) {
                Elem root = roots[rng.RandRange(roots.size())];
                Elem val = 0, pw = 1;
                for (auto co : traces[k]) {
                    val ^= field.Mul(co, pw);
                    pw = field.Mul(pw, root);
                }
                UT_REQUIRE(val == Elem(RefTrace(uint64_t(field.Mul(betas[k], root)), bits, lowmod)));
            }
        }

        // Not fully factorizable: an irreducible quadratic factor, or a
        // repeated root.
        Elem q;
        do {
            q = RandNonzeroElem(rng, field);
        } while (RefTrace(q, bits, lowmod) != 1);
        std::vector<Elem> quad{q, 1, 1};
        UT_REQUIRE(FindRoots(PolyMulRef(poly, quad, field), basis, field).empty());
        std::vector<Elem> dup{roots[rng.RandRange(roots.size())], 1};
        UT_REQUIRE(FindRoots(PolyMulRef(poly, dup, field), basis, field).empty());
    }
}

/** FindRoots with every root in the kernel of the first 12 split
 *  functionals alpha -> Tr(beta_k alpha) of the Montgomery tree's tower:
 *  the splits of depths 0..11 all fail (every root on one side), which runs
 *  the carry and the power-sum split through their no-split paths and the
 *  tree past the carried traces into its own-trace fallback. */
template<typename F>
void TestFindRootsDeepTree(const F& field, uint64_t lowmod, TestRand& rng) {
    typedef typename F::Elem Elem;
    const int bits = field.Bits();
    if (!TraceModFFTCapable(field) || !TraceModCheapScalarMul<F>() || bits % MONT_TOWER_BITS != 0 || bits < 32) return;
    const size_t degree = 300, nfun = 12;
    const Elem basis = RandNonzeroElem(rng, field);
    // The tower split parameters (see ComputeTraceOrbit).
    std::vector<Elem> betas;
    const Elem z = SubfieldGenerator(field, MONT_TOWER_BITS);
    for (Elem xl = basis; betas.size() < nfun; xl = field.Mul2(xl)) {
        Elem v = xl;
        for (int k = 0; k < MONT_TOWER_BITS && betas.size() < nfun; ++k, v = field.Mul(v, z)) betas.push_back(v);
    }
    // Kernel of alpha -> (Tr(beta_k alpha))_k over the bit basis, by
    // elimination on [functional values | identity].
    std::vector<std::pair<uint64_t, uint64_t>> rows; // (values, combination)
    for (int b = 0; b < bits; ++b) {
        uint64_t val = 0;
        for (size_t k = 0; k < nfun; ++k) val |= RefTrace(uint64_t(field.Mul(betas[k], Elem(uint64_t{1} << b))), bits, lowmod) << k;
        rows.emplace_back(val, uint64_t{1} << b);
    }
    for (size_t k = 0, r = 0; k < nfun; ++k) {
        size_t p = r;
        while (p < rows.size() && !((rows[p].first >> k) & 1)) ++p;
        if (p == rows.size()) continue;
        std::swap(rows[p], rows[r]);
        for (size_t i = 0; i < rows.size(); ++i) {
            if (i != r && ((rows[i].first >> k) & 1)) {
                rows[i].first ^= rows[r].first;
                rows[i].second ^= rows[r].second;
            }
        }
        ++r;
    }
    std::vector<uint64_t> kernel;
    for (const auto& row : rows) {
        if (row.first == 0) kernel.push_back(row.second);
    }
    UT_REQUIRE(kernel.size() == size_t(bits) - nfun);
    std::vector<Elem> roots;
    while (roots.size() < degree) {
        uint64_t a = 0;
        for (auto kv : kernel) {
            if (rng.RandRange(2)) a ^= kv;
        }
        Elem e = Elem(a);
        if (e == 0 || std::find(roots.begin(), roots.end(), e) != roots.end()) continue;
        roots.push_back(e);
    }
    auto poly = PolyFromRootsRef(roots, field);
    std::sort(roots.begin(), roots.end());
    std::vector<Elem> osyn(degree + 8, 0);
    for (auto r : roots) {
        Elem r2 = field.Sqr(r), pw = r;
        for (auto& o : osyn) {
            o ^= pw;
            pw = field.Mul(pw, r2);
        }
    }
    for (bool with_sums : {false, true}) {
        auto found = FindRoots(poly, basis, field, with_sums ? &osyn : nullptr);
        std::sort(found.begin(), found.end());
        UT_REQUIRE(found == roots);
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
/** The tower split parameters (TowerBetas, claim C27) are GF(2)-linearly
 * independent for every subfield degree s in {4, 8} dividing Bits. Scaling by
 * the seed's nonzero basis preserves independence, so basis = 1 covers every
 * seed. Also checks the hypotheses of the Lean proof: z lies in GF(2^s)
 * (z^(2^s) = z) and in no proper subfield of it (z^(2^(s/2)) != z). */
template<typename F>
void TestTowerBetas() {
    typedef typename F::Elem Elem;
    F field;
    const int bits = field.Bits();
    for (int s : {4, 8}) {
        if (bits % s != 0) continue;
        Elem z = SubfieldGenerator(field, s);
        Elem zp = z;
        for (int i = 0; i < s / 2; ++i) zp = field.Sqr(zp);
        UT_REQUIRE(zp != z);
        for (int i = 0; i < s / 2; ++i) zp = field.Sqr(zp);
        UT_REQUIRE(zp == z);
        std::vector<Elem> betas(bits);
        TowerBetas(field, Elem(1), s, betas.data());
        uint64_t pivot[64] = {0};
        for (auto b : betas) {
            uint64_t v = uint64_t(b);
            int bit = 63;
            for (; bit >= 0; --bit) {
                if (!((v >> bit) & 1)) continue;
                if (!pivot[bit]) break;
                v ^= pivot[bit];
            }
            UT_REQUIRE(bit >= 0); // reduced to zero: dependent
            pivot[bit] = v;
        }
    }
}

template<typename F>
void RunAllFieldTests(uint64_t seed_offset) {
    F field;
    const uint64_t lowmod = REF_LOWMOD[field.Bits()];
    TestRand rng(g_test_seed + seed_offset);
    TestFieldOps(field, lowmod, rng, 4000);
    TestFieldSerialization(field, rng, 128);
    TestPolyOps(field, rng, 256, 12);
    // Engine-level FFT diagnostics run before the reducer tests that use the
    // FFT tier internally, so engine bugs are reported at their source.
    TestAdditiveFFT(field, rng, 24);
    TestTraceModReducers(field, rng, 24);
    TestSubquadraticGCD(field, rng, 10);
    TestSyndromes(field, lowmod, rng, 128, 12);
    TestBerlekampMassey(field, rng, 128, 10);
    TestBerlekampMasseyVsRef(field, rng, 200, 160);
    TestBerlekampMasseyExhaustive(field);
    TestFindRoots(field, lowmod, rng, 64, 8);
    TestFindRootsOrbit(field, lowmod, rng, 16, 64);
    TestFindRootsMontgomery(field, lowmod, rng);
    TestFindRootsDeepTree(field, lowmod, rng);
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
