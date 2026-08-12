/**********************************************************************
 * Copyright (c) 2026 The minisketch developers                        *
 * Distributed under the MIT software license, see the accompanying   *
 * file LICENSE or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

/** Fuzz target: coverage-guided property testing of the internal polynomial
 * algorithms (DivMod/PolyMod/GCD/Sqr/TraceMod/BerlekampMassey/FindRoots)
 * against the naive references in test_refimpl.h, for an odd-size (11-bit)
 * and a byte-aligned (32-bit) generic field.
 */

// This TU compiles field implementation files directly to reach the internal
// templates, while the rest of the fuzz binary links against
// libminisketch_verify (which defines the same constructor symbols) — so
// rename them here.
#define ConstructGeneric2Bytes ConstructGeneric2BytesPolyOpsFuzz
#define ConstructGeneric4Bytes ConstructGeneric4BytesPolyOpsFuzz
#include "../fields/generic_2bytes.cpp"
#include "../fields/generic_4bytes.cpp"
#undef ConstructGeneric2Bytes
#undef ConstructGeneric4Bytes

#include <algorithm>
#include <cstdint>
#include <vector>

#include "../test_refimpl.h"
#include "FuzzedDataProvider.h"
#include "fuzz.h"

namespace {

constexpr size_t MAX_DEG = 24;

/** Consume a polynomial with up to maxsize coefficients (trailing zeros
 *  stripped, so the leading coefficient is nonzero or the poly is empty). */
template<typename F>
std::vector<typename F::Elem> ConsumePoly(FuzzedDataProvider& provider, const F& field, size_t maxsize) {
    size_t count = provider.ConsumeIntegralInRange<size_t>(0, maxsize);
    std::vector<typename F::Elem> poly(count);
    const uint64_t max_elem = field.Bits() == 64 ? ~uint64_t{0} : (uint64_t{1} << field.Bits()) - 1;
    for (auto& coef : poly) {
        coef = static_cast<typename F::Elem>(provider.ConsumeIntegralInRange<uint64_t>(0, max_elem));
    }
    while (!poly.empty() && poly.back() == 0) poly.pop_back();
    return poly;
}

template<typename F>
std::vector<typename F::Elem> ConsumeMonicPoly(FuzzedDataProvider& provider, const F& field, size_t minsize, size_t maxsize) {
    auto poly = ConsumePoly(provider, field, maxsize);
    if (poly.size() < minsize) poly.resize(minsize, 0);
    poly.back() = 1;
    return poly;
}

/** Consume up to maxcount distinct nonzero field elements. */
template<typename F>
std::vector<typename F::Elem> ConsumeDistinctElems(FuzzedDataProvider& provider, const F& field, size_t maxcount) {
    const uint64_t max_elem = field.Bits() == 64 ? ~uint64_t{0} : (uint64_t{1} << field.Bits()) - 1;
    size_t count = provider.ConsumeIntegralInRange<size_t>(0, maxcount);
    std::vector<typename F::Elem> ret;
    while (ret.size() < count && provider.remaining_bytes() > 0) {
        auto elem = static_cast<typename F::Elem>(provider.ConsumeIntegralInRange<uint64_t>(1, max_elem));
        if (std::find(ret.begin(), ret.end(), elem) == ret.end()) ret.push_back(elem);
    }
    return ret;
}

template<typename F>
void PolyOpsForField(FuzzedDataProvider& provider, const F& field) {
    typedef typename F::Elem Elem;
    const int bits = field.Bits();
    const uint64_t lowmod = ut::REF_LOWMOD[bits];
    const uint64_t max_elem = bits == 64 ? ~uint64_t{0} : (uint64_t{1} << bits) - 1;

    // DivMod/PolyMod: val == div*mod + rem, deg(rem) < deg(mod), PolyMod agrees.
    {
        auto mod = ConsumeMonicPoly(provider, field, 1, MAX_DEG);
        auto val = ConsumePoly(provider, field, 2 * MAX_DEG);
        std::vector<Elem> quot, rem = val;
        DivMod(mod, rem, quot, field);
        FUZZ_CHECK(ut::Stripped(rem).size() < mod.size());
        auto recombined = ut::PolyMulRef(quot, mod, field);
        recombined.resize(std::max(recombined.size(), rem.size()), 0);
        for (size_t i = 0; i < rem.size(); ++i) recombined[i] ^= rem[i];
        FUZZ_CHECK(ut::Stripped(recombined) == val);
        auto rem2 = val;
        PolyMod(mod, rem2, field);
        FUZZ_CHECK(rem2 == ut::Stripped(rem));
    }

    // Sqr against the naive product.
    {
        auto poly = ConsumePoly(provider, field, MAX_DEG);
        auto sqr = poly;
        Sqr(sqr, field);
        FUZZ_CHECK(ut::Stripped(sqr) == ut::PolyMulRef(poly, poly, field));
    }

    // GCD(f*g, h*g) divides both products and is divisible by g.
    {
        auto f = ConsumePoly(provider, field, MAX_DEG / 2);
        auto g = ConsumePoly(provider, field, MAX_DEG / 2);
        auto h = ConsumePoly(provider, field, MAX_DEG / 2);
        if (!f.empty() && !g.empty() && !h.empty()) {
            auto a = ut::PolyMulRef(f, g, field);
            auto b = ut::PolyMulRef(h, g, field);
            auto a_copy = a, b_copy = b;
            GCD(a, b, field);
            auto d = a;
            FUZZ_CHECK(!d.empty());
            MakeMonic(d, field);
            ut::PolyReduceRef(a_copy, d, field);
            FUZZ_CHECK(a_copy.empty());
            ut::PolyReduceRef(b_copy, d, field);
            FUZZ_CHECK(b_copy.empty());
            auto g_monic = g;
            MakeMonic(g_monic, field);
            auto d_copy = d;
            ut::PolyReduceRef(d_copy, g_monic, field);
            FUZZ_CHECK(d_copy.empty());
        }
    }

    // TraceMod against the naive sum of Frobenius powers.
    {
        auto mod = ConsumeMonicPoly(provider, field, 4, MAX_DEG);
        Elem param = static_cast<Elem>(provider.ConsumeIntegralInRange<uint64_t>(1, max_elem));
        std::vector<Elem> trace;
        TraceMod<F> trace_mod(mod, field);
        trace_mod.trace(trace, param);
        std::vector<Elem> cur{0, param};
        ut::PolyReduceRef(cur, mod, field);
        std::vector<Elem> acc = cur;
        for (int i = 1; i < bits; ++i) {
            cur = ut::PolyMulRef(cur, cur, field);
            ut::PolyReduceRef(cur, mod, field);
            acc.resize(std::max(acc.size(), cur.size()), 0);
            for (size_t j = 0; j < cur.size(); ++j) acc[j] ^= cur[j];
        }
        FUZZ_CHECK(ut::Stripped(trace) == ut::Stripped(acc));

        // Repeated use of one reducer object, and SquareAndReduce against the
        // naive reference.
        std::vector<Elem> trace2;
        trace_mod.trace(trace2, param);
        FUZZ_CHECK(trace2 == trace);
        auto val = ConsumePoly(provider, field, mod.size() - 1);
        auto reduced = val;
        trace_mod.SquareAndReduce(reduced);
        auto ref_sq = ut::PolyMulRef(val, val, field);
        ut::PolyReduceRef(ref_sq, mod, field);
        FUZZ_CHECK(ut::Stripped(reduced) == ref_sq);
    }

    // Additive-FFT tier (engaged only for the 32-bit field here; the 11-bit
    // field is ineligible and keeps exercising the pure Karatsuba path).
    // Behind a coin flip to keep the average exec cheap.
    if (bits == 32 && provider.ConsumeBool()) {
        // fft/ifft roundtrip at a fuzz-chosen dimension.
        {
            AdditiveFFT<F> fft;
            int m = provider.ConsumeIntegralInRange<int>(1, 10);
            std::vector<Elem> f(size_t(1) << m);
            for (auto& e : f) e = static_cast<Elem>(provider.ConsumeIntegralInRange<uint64_t>(0, max_elem));
            auto orig = f;
            fft.Extend(field, m);
            fft.FFT(f, m, field);
            fft.IFFT(f, m, field);
            FUZZ_CHECK(f == orig);
        }

        // Products crossing the FFT cutoff, against the schoolbook oracle;
        // TraceModPolyMulLow must equal the truncated full product.
        auto fa = ConsumePoly(provider, field, 160);
        auto fb = ConsumePoly(provider, field, 160);
        TraceModScratch<F> scratch;
        std::vector<Elem> full;
        TraceModPolyMulFull(fa, fb, full, field, scratch, 0);
        FUZZ_CHECK(ut::Stripped(full) == ut::PolyMulRef(fa, fb, field));
        size_t n = 1 + provider.ConsumeIntegralInRange<size_t>(0, 340);
        std::vector<Elem> low;
        TraceModPolyMulLow(fa, fb, low, n, field, scratch, 0);
        auto ref_low = full;
        ref_low.resize(n, 0);
        FUZZ_CHECK(ut::Stripped(low) == ut::Stripped(ref_low));
    }

    // BerlekampMassey on syndromes of a known root set returns prod (1 + m*x);
    // FindRoots on prod (x + m) returns the set.
    {
        auto roots = ConsumeDistinctElems(provider, field, 12);
        size_t cap = roots.size() + provider.ConsumeIntegralInRange<size_t>(0, 3);
        if (cap > 0) {
            std::vector<Elem> osyn(cap, 0);
            for (auto m : roots) AddToOddSyndromes(osyn, m, field);
            auto syndromes = ReconstructAllSyndromes(osyn, field);
            auto poly = BerlekampMassey(syndromes, cap, field);
            std::vector<Elem> expect{1};
            for (auto m : roots) {
                std::vector<Elem> factor{1, m};
                expect = ut::PolyMulRef(expect, factor, field);
            }
            FUZZ_CHECK(poly == expect);
            if (!roots.empty()) {
                FUZZ_CHECK(BerlekampMassey(syndromes, roots.size() - 1, field).empty());
            }
        }
        if (!roots.empty()) {
            auto poly = ut::PolyFromRootsRef(roots, field);
            Elem basis = static_cast<Elem>(provider.ConsumeIntegralInRange<uint64_t>(1, max_elem));
            auto found = FindRoots(poly, basis, field);
            std::sort(found.begin(), found.end());
            std::sort(roots.begin(), roots.end());
            FUZZ_CHECK(found == roots);

            // Multiplying in an irreducible quadratic (x^2 + x + c, Tr(c) = 1)
            // makes the polynomial not fully factorizable.
            Elem c = static_cast<Elem>(provider.ConsumeIntegralInRange<uint64_t>(1, max_elem));
            if (ut::RefTrace(c, bits, lowmod) == 1) {
                std::vector<Elem> quad{c, 1, 1};
                auto not_fact = ut::PolyMulRef(poly, quad, field);
                FUZZ_CHECK(FindRoots(not_fact, basis, field).empty());
            }
        }
    }

    // Large-degree decode stages, reference-free (the naive references above
    // are quadratic, so they stay small): syndromes of a root set must
    // BerlekampMassey to a polynomial of exactly that degree whose reversal's
    // roots are exactly the set. Root count is log-uniform up to 256.
    {
        size_t want = size_t{1} << provider.ConsumeIntegralInRange<int>(0, 8);
        auto roots = ConsumeDistinctElems(provider, field, want);
        if (!roots.empty()) {
            std::vector<Elem> osyn(roots.size(), 0);
            for (auto m : roots) AddToOddSyndromes(osyn, m, field);
            auto syndromes = ReconstructAllSyndromes(osyn, field);
            auto poly = BerlekampMassey(syndromes, roots.size(), field);
            FUZZ_CHECK(poly.size() == roots.size() + 1);
            std::reverse(poly.begin(), poly.end());
            Elem basis = static_cast<Elem>(provider.ConsumeIntegralInRange<uint64_t>(1, max_elem));
            auto found = FindRoots(poly, basis, field);
            std::sort(found.begin(), found.end());
            std::sort(roots.begin(), roots.end());
            FUZZ_CHECK(found == roots);
        }
    }
}

} // namespace

FUZZ_TARGET(poly_ops) {
    FuzzedDataProvider provider(data, size);
    if (provider.ConsumeBool()) {
        PolyOpsForField(provider, Field11());
    } else {
        PolyOpsForField(provider, Field32());
    }
}
