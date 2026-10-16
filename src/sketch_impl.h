/**********************************************************************
 * Copyright (c) 2018 Pieter Wuille, Greg Maxwell, Gleb Naumenko      *
 * Distributed under the MIT software license, see the accompanying   *
 * file LICENSE or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

#ifndef _MINISKETCH_SKETCH_IMPL_H_
#define _MINISKETCH_SKETCH_IMPL_H_

#include <deque>
#include <memory>
#include <random>

#include "util.h"
#include "sketch.h"
#include "int_utils.h"
#include "additive_fft.h"

/* Algorithm crossover cutoffs. Each can be overridden at build time with
 * -DMINISKETCH_<NAME>=<value> (e.g. for cross-machine tuning via
 * tools/bench_ab.sh CXXFLAGS_B); the defaults were tuned on Zen 5. */
#ifndef MINISKETCH_TRACEMOD_POLYMUL_CUTOFF
#define MINISKETCH_TRACEMOD_POLYMUL_CUTOFF 24
#endif
#ifndef MINISKETCH_TRACEMOD_TABLE_CUTOFF
#define MINISKETCH_TRACEMOD_TABLE_CUTOFF 512
#endif
#ifndef MINISKETCH_TRACEMOD_FFT_TABLE_CUTOFF
#define MINISKETCH_TRACEMOD_FFT_TABLE_CUTOFF 704
#endif
#ifndef MINISKETCH_TRACEMOD_FFT_CUTOFF
#define MINISKETCH_TRACEMOD_FFT_CUTOFF 128
#endif
#ifndef MINISKETCH_FAST_DIVMOD_CUTOFF
#define MINISKETCH_FAST_DIVMOD_CUTOFF 1024
#endif
#ifndef MINISKETCH_HGCD_CUTOFF
#define MINISKETCH_HGCD_CUTOFF 4096
#endif
#ifndef MINISKETCH_HGCD_ITER_CROSSOVER
#define MINISKETCH_HGCD_ITER_CROSSOVER 24
#endif
/* Divisor degree below which normal Euclid steps (on lazily-reducing
 * fields) use the inversion-free pseudo-remainder. */
#ifndef MINISKETCH_GCD_INVFREE_MAX
#define MINISKETCH_GCD_INVFREE_MAX 128
#endif
/* Divisor length above which DivMod and GCD use lazy reduction (on fields
 * that support it). */
#ifndef MINISKETCH_DIVMOD_LAZY_MIN
#define MINISKETCH_DIVMOD_LAZY_MIN 4
#endif
/* Node degree up to which RecFindRoots computes traces with SmallTrace. */
#ifndef MINISKETCH_SMALL_TRACE_MAX
#define MINISKETCH_SMALL_TRACE_MAX 8
#endif
/* Berlekamp-Massey early termination: a tail verification is attempted after
 * a zero discrepancy once at least MIN_PREFIX syndromes are processed, the
 * processed count exceeds twice the LFSR length by GAP, and at least MIN_TAIL
 * syndromes remain. */
#ifndef MINISKETCH_BM_VERIFY_GAP
#define MINISKETCH_BM_VERIFY_GAP 2
#endif
#ifndef MINISKETCH_BM_VERIFY_MIN_PREFIX
#define MINISKETCH_BM_VERIFY_MIN_PREFIX 0
#endif
/* Minimum LFSR length for the additive-FFT tail verification, for fields
 * with a cheap, lazily-reducing scalar multiplication (CLMUL; below it, the
 * lazy direct / Karatsuba middle-product checks are faster, measured on
 * Haswell and Cortex-A76) and for table-based ones (generic, whose
 * pointwise FFT products are slow one-off multiplications; below it, a
 * direct Multiplier-table check). */
#ifndef MINISKETCH_BM_VERIFY_FFT_MIN_DEGREE
#define MINISKETCH_BM_VERIFY_FFT_MIN_DEGREE 512
#endif
#ifndef MINISKETCH_BM_VERIFY_FFT_MIN_DEGREE_TABLE
#define MINISKETCH_BM_VERIFY_FFT_MIN_DEGREE_TABLE 40
#endif
/* Root finding via the Frobenius orbit of the root polynomial (TraceOrbit):
 * engaged for root degrees >= MIN_DEGREE; precomputes the traces of the
 * first log2(Bits) + EXTRA_DEPTHS split depths; a node of degree d uses them
 * when NODE_FACTOR * deg(root) <= d * Bits. */
#ifndef MINISKETCH_ORBIT_MIN_DEGREE
#define MINISKETCH_ORBIT_MIN_DEGREE 8
#endif
#ifndef MINISKETCH_ORBIT_EXTRA_DEPTHS
#define MINISKETCH_ORBIT_EXTRA_DEPTHS 2
#endif
#ifndef MINISKETCH_ORBIT_NODE_FACTOR
#define MINISKETCH_ORBIT_NODE_FACTOR 2
#endif
/* Root degree from which the Frobenius orbit uses coset Montgomery
 * reduction (MontgomeryOrbit) on fields with the FFT tier and a cheap
 * scalar Multiplier. */
#ifndef MINISKETCH_MONT_ORBIT_MIN_DEGREE
#define MINISKETCH_MONT_ORBIT_MIN_DEGREE 224
#endif
/* Node degree up to which the root-finding tree carries the relative
 * traces in normal form instead of Montgomery form (MontCarry). */
/* Subfield degree of the Montgomery tree's tower split parameters: each
 * carried relative trace serves this many depths, at this many minus one
 * modular squarings per node (1, 2, 4 or 8). */
#ifndef MINISKETCH_MONT_TOWER_BITS
#define MINISKETCH_MONT_TOWER_BITS 4
#endif
#ifndef MINISKETCH_SMALL_CARRY_MAX
#define MINISKETCH_SMALL_CARRY_MAX 16
#endif
/* Node degree from which Montgomery-form split traces are normalized
 * before the gcd (MontTowerTrace). */
/* Node degree from which Montgomery-form nodes divide out a split factor
 * pointwise on their coset instead of by long division. */
#ifndef MINISKETCH_COSET_DIV_MIN
#define MINISKETCH_COSET_DIV_MIN 64
#endif
/* Node degree from which Montgomery-form nodes split by power sums
 * (PowerSumSplit) instead of a gcd when the power sums are known. */
#ifndef MINISKETCH_PSUM_SPLIT_MIN
#define MINISKETCH_PSUM_SPLIT_MIN 256
#endif
#ifndef MINISKETCH_MONT_NORMALIZE_MIN
#define MINISKETCH_MONT_NORMALIZE_MIN 512
#endif
/* Lazy Karatsuba middle-product tail check (lazily-reducing fields): used
 * for LFSR lengths in [KARA_MIN_DEGREE, FFT_MIN_DEGREE), below which the
 * direct check wins (from ~32 on Haswell, ~96 on Cortex-A76); direct
 * kernels at middle-product sizes <= MP_BASE. */
#ifndef MINISKETCH_BM_VERIFY_KARA_MIN_DEGREE
#define MINISKETCH_BM_VERIFY_KARA_MIN_DEGREE 64
#endif
#ifndef MINISKETCH_MP_BASE
#define MINISKETCH_MP_BASE 16
#endif
#ifndef MINISKETCH_BM_VERIFY_MIN_TAIL
#define MINISKETCH_BM_VERIFY_MIN_TAIL 32
#endif

static const size_t TRACEMOD_POLYMUL_CUTOFF = MINISKETCH_TRACEMOD_POLYMUL_CUTOFF;
static const size_t TRACEMOD_TABLE_CUTOFF = MINISKETCH_TRACEMOD_TABLE_CUTOFF;
static const size_t TRACEMOD_FFT_TABLE_CUTOFF = MINISKETCH_TRACEMOD_FFT_TABLE_CUTOFF;
static const size_t TRACEMOD_FFT_CUTOFF = MINISKETCH_TRACEMOD_FFT_CUTOFF;
static const size_t FAST_DIVMOD_CUTOFF = MINISKETCH_FAST_DIVMOD_CUTOFF;
static const size_t HGCD_CUTOFF = MINISKETCH_HGCD_CUTOFF;
static const size_t HGCD_ITER_CROSSOVER = MINISKETCH_HGCD_ITER_CROSSOVER;
static const size_t BM_VERIFY_GAP = MINISKETCH_BM_VERIFY_GAP;
static const size_t SMALL_TRACE_MAX = MINISKETCH_SMALL_TRACE_MAX;
static const size_t DIVMOD_LAZY_MIN = MINISKETCH_DIVMOD_LAZY_MIN;
static const size_t GCD_INVFREE_MAX = MINISKETCH_GCD_INVFREE_MAX;
static const size_t BM_VERIFY_MIN_PREFIX = MINISKETCH_BM_VERIFY_MIN_PREFIX;
static const size_t BM_VERIFY_MIN_TAIL = MINISKETCH_BM_VERIFY_MIN_TAIL;
static const size_t BM_VERIFY_FFT_MIN_DEGREE = MINISKETCH_BM_VERIFY_FFT_MIN_DEGREE;
static const size_t BM_VERIFY_FFT_MIN_DEGREE_TABLE = MINISKETCH_BM_VERIFY_FFT_MIN_DEGREE_TABLE;
static const size_t BM_VERIFY_KARA_MIN_DEGREE = MINISKETCH_BM_VERIFY_KARA_MIN_DEGREE;
static const size_t MP_BASE = MINISKETCH_MP_BASE;
static const size_t ORBIT_MIN_DEGREE = MINISKETCH_ORBIT_MIN_DEGREE;
static const size_t ORBIT_EXTRA_DEPTHS = MINISKETCH_ORBIT_EXTRA_DEPTHS;
static const size_t ORBIT_NODE_FACTOR = MINISKETCH_ORBIT_NODE_FACTOR;
static const size_t MONT_ORBIT_MIN_DEGREE = MINISKETCH_MONT_ORBIT_MIN_DEGREE;
static const size_t SMALL_CARRY_MAX = MINISKETCH_SMALL_CARRY_MAX;
static const size_t MONT_NORMALIZE_MIN = MINISKETCH_MONT_NORMALIZE_MIN;
static const size_t COSET_DIV_MIN = MINISKETCH_COSET_DIV_MIN;
static const size_t PSUM_SPLIT_MIN = MINISKETCH_PSUM_SPLIT_MIN;

/** Compute the remainder of a polynomial division of val by mod, putting the result in mod. */
template<typename F>
void PolyMod(const std::vector<typename F::Elem>& mod, std::vector<typename F::Elem>& val, const F& field) {
    size_t modsize = mod.size();
    CHECK_SAFE(modsize > 0 && mod.back() == 1);
    if (val.size() < modsize) return;
    CHECK_SAFE(val.back() != 0);
    while (val.size() >= modsize) {
        auto term = val.back();
        val.pop_back();
        if (term != 0) {
            typename F::Multiplier mul(field, term);
            for (size_t x = 0; x < mod.size() - 1; ++x) {
                val[val.size() - modsize + 1 + x] ^= mul(mod[x]);
            }
        }
    }
    while (val.size() > 0 && val.back() == 0) val.pop_back();
}

/** Division in dot-product form for lazily-reducing fields: a (length na)
 * divided by b (degree db, inv = b[db]^-1, 1 for monic b) gives the
 * quotient q (length na - db, requires na > db) and remainder r (length
 * db). Every coefficient is one register-accumulated unreduced dot product
 * reduced once:
 *     q_i = inv * (a_{db+i} + sum_{t>i} q_t b_{db+i-t}),
 *     r_j = a_j + sum_t q_t b_{j-t}.
 * No accumulator array is kept (cf. the row-update form of DivMod), so
 * the only memory traffic is reading the operands. q and r may not alias
 * a or b. */
template<typename F>
void LazyDivRem(const typename F::Elem* a, size_t na, const typename F::Elem* b, size_t db, typename F::Elem inv, typename F::Elem* q, typename F::Elem* r, const F& field) {
    typedef typename F::Wide Wide;
    const size_t nq = na - db;
    for (size_t i = nq; i-- > 0; ) {
        Wide acc = field.WideFrom(a[db + i]);
        const size_t tmax = std::min(nq - 1, i + db);
        for (size_t t = i + 1; t <= tmax; ++t) acc = field.WideXor(acc, field.MulWide(q[t], b[db + i - t]));
        auto v = field.Reduce(acc);
        q[i] = inv == 1 ? v : field.Mul(v, inv);
    }
    for (size_t j = 0; j < db; ++j) {
        Wide acc = field.WideFrom(a[j]);
        const size_t tmax = std::min(nq - 1, j);
        for (size_t t = 0; t <= tmax; ++t) acc = field.WideXor(acc, field.MulWide(q[t], b[j - t]));
        r[j] = field.Reduce(acc);
    }
}

/** Compute the quotient of a polynomial division of val by mod, putting the quotient in div and the remainder in val. */
template<typename F>
void DivMod(const std::vector<typename F::Elem>& mod, std::vector<typename F::Elem>& val, std::vector<typename F::Elem>& div, const F& field) {
    size_t modsize = mod.size();
    CHECK_SAFE(mod.size() > 0 && mod.back() == 1);
    if (val.size() < mod.size()) {
        div.clear();
        return;
    }
    CHECK_SAFE(val.back() != 0);
    div.resize(val.size() - mod.size() + 1);
    if (F::LAZY_REDUCE && modsize > DIVMOD_LAZY_MIN) {
        const size_t dm = modsize - 1;
        std::vector<typename F::Elem> rem(dm);
        LazyDivRem(val.data(), val.size(), mod.data(), dm, typename F::Elem(1), div.data(), rem.data(), field);
        val.swap(rem);
        return;
    }
    while (val.size() >= modsize) {
        auto term = val.back();
        div[val.size() - modsize] = term;
        val.pop_back();
        if (term != 0) {
            typename F::Multiplier mul(field, term);
            for (size_t x = 0; x < mod.size() - 1; ++x) {
                val[val.size() - modsize + 1 + x] ^= mul(mod[x]);
            }
        }
    }
}

/** Make a polynomial monic. */
template<typename F>
typename F::Elem MakeMonic(std::vector<typename F::Elem>& a, const F& field) {
    CHECK_SAFE(a.back() != 0);
    if (a.back() == 1) return 0;
    auto inv = field.Inv(a.back());
    typename F::Multiplier mul(field, inv);
    a.back() = 1;
    for (size_t i = 0; i < a.size() - 1; ++i) {
        a[i] = mul(a[i]);
    }
    return inv;
}

/** Compute the GCD of two polynomials, putting the result in a. b will be cleared. */
template<typename F>
void GCD(std::vector<typename F::Elem>& a, std::vector<typename F::Elem>& b, const F& field) {
    if (a.size() < b.size()) std::swap(a, b);
    if (F::LAZY_REDUCE && b.size() > DIVMOD_LAZY_MIN) {
        // Euclid without per-step normalization: the divisor's inverse
        // leading coefficient is folded into the quotient terms (one
        // multiplication each) instead of rescaling the divisor to monic,
        // and each step is one lazy dot-product division (LazyDivRem).
        // The result is made monic once at the end.
        std::vector<typename F::Elem> q, r;
        bool divided = false;
        while (b.size() > 1) {
            divided = true;
            const size_t db = b.size() - 1;
            if (a.size() == db + 2 && db < GCD_INVFREE_MAX) {
                // Normal step (quotient degree 1), inversion-free: with
                // beta = lc(b), the pseudo-remainder
                //   beta^2*a - (beta*a_n*x + beta*a_{n-1} - a_n*b_{m-1})*b
                // is beta^2 times the remainder; the gcd is unaffected and
                // the result is normalized once at the end. One more
                // product per coefficient instead of a field inversion.
                typedef typename F::Wide Wide;
                const auto beta = b[db], an = a[db + 1], an1 = a[db];
                const auto p1 = field.Mul(beta, an);
                const auto p0 = field.Mul(beta, an1) ^ field.Mul(an, b[db - 1]);
                const auto beta2 = field.Sqr(beta);
                r.resize(db);
                r[0] = field.Reduce(field.WideXor(field.MulWide(beta2, a[0]), field.MulWide(p0, b[0])));
                for (size_t j = 1; j < db; ++j) {
                    Wide acc = field.WideXor(field.MulWide(beta2, a[j]), field.MulWide(p0, b[j]));
                    r[j] = field.Reduce(field.WideXor(acc, field.MulWide(p1, b[j - 1])));
                }
                a.swap(r);
            } else if (a.size() > db) {
                q.resize(a.size() - db);
                r.resize(db);
                LazyDivRem(a.data(), a.size(), b.data(), db, field.Inv(b.back()), q.data(), r.data(), field);
                a.swap(r);
            }
            while (!a.empty() && a.back() == 0) a.pop_back();
            std::swap(a, b);
        }
        if (b.size() == 1) {
            // A nonzero constant remainder: coprime.
            a.assign(1, 1);
            return;
        }
        // a is the last divisor (the unchanged input if there was none, as
        // in the plain loop below, which only normalizes divisors).
        if (divided) MakeMonic(a, field);
        return;
    }
    while (b.size() > 0) {
        if (b.size() == 1) {
            a.resize(1);
            a[0] = 1;
            return;
        }
        MakeMonic(b, field);
        PolyMod(b, a, field);
        std::swap(a, b);
    }
}

/** Square a polynomial. */
template<typename F>
void Sqr(std::vector<typename F::Elem>& poly, const F& field) {
    if (poly.size() == 0) return;
    poly.resize(poly.size() * 2 - 1);
    for (size_t i = 0; i < poly.size(); ++i) {
        auto x = poly.size() - i - 1;
        poly[x] = (x & 1) ? 0 : field.Sqr(poly[x / 2]);
    }
}

/** Compute a*b mod x^n: the low n coefficients of the product in F[x]. */
template<typename F>
void PolyMulLowNaive(const std::vector<typename F::Elem>& a, const std::vector<typename F::Elem>& b, std::vector<typename F::Elem>& out, size_t n, const F& field) {
    typedef typename F::Elem Elem;
    out.assign(n, 0);
    size_t na = std::min(a.size(), n);
    if (F::LAZY_REDUCE) {
        // Convolution order: each output coefficient is one unreduced
        // accumulation, reduced once.
        size_t nb = std::min(b.size(), n);
        if (na == 0 || nb == 0) return;
        size_t n_out = std::min(n, na + nb - 1);
        for (size_t k = 0; k < n_out; ++k) {
            size_t lo = k + 1 > nb ? k + 1 - nb : 0, hi = std::min(k + 1, na);
            typename F::Wide acc = field.WideZero();
            for (size_t i = lo; i < hi; ++i) acc = field.WideXor(acc, field.MulWide(a[i], b[k - i]));
            out[k] = field.Reduce(acc);
        }
        return;
    }
    for (size_t i = 0; i < na; ++i) {
        Elem coeff = a[i];
        if (coeff == 0) continue;
        size_t max_j = std::min(b.size(), n - i);
        if (coeff == 1) {
            for (size_t j = 0; j < max_j; ++j) out[i + j] ^= b[j];
        } else {
            typename F::Multiplier mul(field, coeff);
            for (size_t j = 0; j < max_j; ++j) out[i + j] ^= mul(b[j]);
        }
    }
}

/** Whether the field supports the additive-FFT tier at all: its degree must
 * be a power of two of at least 16 (the Cantor special basis requirement;
 * smaller power-of-two fields never reach the product-length cutoff). */
template<typename F>
bool TraceModFFTCapable(const F& field) {
    int bits = field.Bits();
    return bits >= 16 && (bits & (bits - 1)) == 0;
}

/** Whether the field implementation's scalar Multiplier is a plain wrapper
 * (the CLMUL fields, which store just the element) rather than a
 * precomputation table (the generic fields, which pay a table build per
 * use). Several cutoffs below depend on this distinction. */
template<typename F>
constexpr bool TraceModCheapScalarMul() {
    return sizeof(typename F::Multiplier) <= 2 * sizeof(typename F::Elem);
}

/** The degree at which TraceMod switches from the square-table reducer to
 * the reciprocal reducer. For field implementations with a cheap scalar
 * Multiplier (CLMUL) the additive-FFT tier speeds up the reciprocal path,
 * but the lazily-reducing, streaming square-table reducer is faster still,
 * so the crossover sits above the original value (measured: ~700 on
 * Haswell, above 768 on Cortex-A76; the table far exceeds L2 there but is
 * streamed sequentially). The generic fields keep the original value:
 * their reciprocal path pays slow one-off products in the FFT. */
template<typename F>
size_t TraceModTableCutoff(const F& field) {
    return (TraceModFFTCapable(field) && TraceModCheapScalarMul<F>()) ? TRACEMOD_FFT_TABLE_CUTOFF : TRACEMOD_TABLE_CUTOFF;
}

/** Whether the additive-FFT multiplication tier applies to a product of
 * (truncated) operand lengths la and lb: the field must support it, the
 * product must be balanced enough and long enough, and 2^t >= la + lb - 1
 * evaluation points must exist with t <= field.Bits() (a real constraint
 * only for the 16-bit field; the fallback is Karatsuba, whose recursive
 * halves may re-enter this tier). */
template<typename F>
bool TraceModFFTEligible(const F& field, size_t la, size_t lb) {
    int bits = field.Bits();
    if (!TraceModFFTCapable(field)) return false;
    if (la == 0 || lb == 0) return false;
    if (std::min(la, lb) <= TRACEMOD_POLYMUL_CUTOFF) return false;
    size_t prod = la + lb - 1;
    if (prod < TRACEMOD_FFT_CUTOFF) return false;
    int t = 1;
    while ((size_t(1) << t) < prod && t < 64) ++t;
    return t <= bits;
}

/** Reusable scratch buffers for the TraceMod polynomial multiplications below.
 * Slot (level, idx) is dedicated to one temporary of one recursion level, so a
 * callee never touches its caller's buffers, and buffer capacity persists
 * across calls. A deque keeps references to existing slots valid while deeper
 * recursion levels append new ones. Also owns the additive-FFT engine, whose
 * basis/twiddle tables amortize over all multiplications sharing the scratch
 * (one TraceMod node performs 2*(Bits-1) reductions with it). */
template<typename F>
struct TraceModScratch {
    static const size_t SLOTS = 9;
    std::deque<std::vector<typename F::Elem>> slots;
    AdditiveFFT<F> fft;
    std::vector<typename F::Elem>& Get(size_t level, size_t idx) {
        CHECK_SAFE(idx < SLOTS);
        size_t i = level * SLOTS + idx;
        while (i >= slots.size()) slots.emplace_back();
        return slots[i];
    }
};

template<typename F>
void TraceModPolyMulFull(const std::vector<typename F::Elem>& a, const std::vector<typename F::Elem>& b, std::vector<typename F::Elem>& out, const F& field, TraceModScratch<F>& scratch, size_t level);

/** Compute a*b mod x^n: the low n coefficients of the product in F[x].
 * Uses Karatsuba with the characteristic-2 identity
 * (a0+a1)*(b0+b1) + a0*b0 + a1*b1 for the middle term.
 * Temporaries live in scratch slots of the given recursion level, so repeated
 * calls with the same scratch are allocation-free in steady state; out must
 * not alias a, b, or a scratch slot. Peak live polynomial storage is bounded
 * by about 11*n field elements, including recursive calls. */
template<typename F>
void TraceModPolyMulLow(const std::vector<typename F::Elem>& a, const std::vector<typename F::Elem>& b, std::vector<typename F::Elem>& out, size_t n, const F& field, TraceModScratch<F>& scratch, size_t level) {
    typedef typename F::Elem Elem;
    if (n == 0 || a.empty() || b.empty()) {
        out.assign(n, 0);
        return;
    }
    if (TraceModFFTEligible(field, std::min(a.size(), n), std::min(b.size(), n))) {
        scratch.fft.MulLow(a, b, out, n, field);
        return;
    }
    if (n <= TRACEMOD_POLYMUL_CUTOFF || std::min(a.size(), b.size()) <= TRACEMOD_POLYMUL_CUTOFF) {
        PolyMulLowNaive(a, b, out, n, field);
        return;
    }

    size_t split = (n + 1) / 2;
    size_t a0_len = std::min(a.size(), split);
    size_t b0_len = std::min(b.size(), split);
    auto& a0 = scratch.Get(level, 0);
    auto& b0 = scratch.Get(level, 1);
    auto& z0 = scratch.Get(level, 2);
    a0.assign(a.begin(), a.begin() + a0_len);
    b0.assign(b.begin(), b.begin() + b0_len);
    TraceModPolyMulFull(a0, b0, z0, field, scratch, level + 1);
    if (z0.size() > n) z0.resize(n);

    out.assign(n, 0);
    for (size_t i = 0; i < z0.size(); ++i) out[i] = z0[i];
    if (split >= n) return;

    size_t mid_n = n - split;
    auto& a1 = scratch.Get(level, 3);
    auto& b1 = scratch.Get(level, 4);
    a1.clear();
    b1.clear();
    if (a.size() > split) a1.assign(a.begin() + split, a.begin() + std::min(a.size(), n));
    if (b.size() > split) b1.assign(b.begin() + split, b.begin() + std::min(b.size(), n));
    if (a1.empty() && b1.empty()) return;

    if (!a1.empty() && !b1.empty()) {
        size_t asum_len = std::max(a0.size(), a1.size());
        size_t bsum_len = std::max(b0.size(), b1.size());
        auto& asum = scratch.Get(level, 5);
        auto& bsum = scratch.Get(level, 6);
        asum.assign(asum_len, 0);
        bsum.assign(bsum_len, 0);
        for (size_t i = 0; i < asum_len; ++i) {
            if (i < a0.size()) asum[i] ^= a0[i];
            if (i < a1.size()) asum[i] ^= a1[i];
        }
        for (size_t i = 0; i < bsum_len; ++i) {
            if (i < b0.size()) bsum[i] ^= b0[i];
            if (i < b1.size()) bsum[i] ^= b1[i];
        }
        auto& zsum = scratch.Get(level, 7);
        auto& z2 = scratch.Get(level, 8);
        TraceModPolyMulLow(asum, bsum, zsum, mid_n, field, scratch, level + 1);
        TraceModPolyMulLow(a1, b1, z2, mid_n, field, scratch, level + 1);
        for (size_t i = 0; i < mid_n; ++i) {
            Elem v = (i < zsum.size() ? zsum[i] : Elem(0)) ^ (i < z2.size() ? z2[i] : Elem(0)) ^ (i < z0.size() ? z0[i] : Elem(0));
            out[split + i] ^= v;
        }
    } else {
        auto& tmp = scratch.Get(level, 7);
        if (!a1.empty()) {
            TraceModPolyMulLow(a1, b0, tmp, mid_n, field, scratch, level + 1);
            for (size_t i = 0; i < tmp.size(); ++i) out[split + i] ^= tmp[i];
        }
        if (!b1.empty()) {
            TraceModPolyMulLow(a0, b1, tmp, mid_n, field, scratch, level + 1);
            for (size_t i = 0; i < tmp.size(); ++i) out[split + i] ^= tmp[i];
        }
    }
}

/** Compute the full product a*b in F[x] using the same characteristic-2 Karatsuba identity.
 * Same scratch-slot discipline as TraceModPolyMulLow. Peak live polynomial
 * storage is bounded by about 10*n field elements, including recursive calls. */
template<typename F>
void TraceModPolyMulFull(const std::vector<typename F::Elem>& a, const std::vector<typename F::Elem>& b, std::vector<typename F::Elem>& out, const F& field, TraceModScratch<F>& scratch, size_t level) {
    typedef typename F::Elem Elem;
    if (a.empty() || b.empty()) {
        out.clear();
        return;
    }
    if (TraceModFFTEligible(field, a.size(), b.size())) {
        scratch.fft.MulFull(a, b, out, field);
        return;
    }
    if (std::min(a.size(), b.size()) <= TRACEMOD_POLYMUL_CUTOFF) {
        PolyMulLowNaive(a, b, out, a.size() + b.size() - 1, field);
        return;
    }
    size_t split = (std::max(a.size(), b.size()) + 1) / 2;
    auto& a0 = scratch.Get(level, 0);
    auto& b0 = scratch.Get(level, 1);
    auto& a1 = scratch.Get(level, 2);
    auto& b1 = scratch.Get(level, 3);
    a0.assign(a.begin(), a.begin() + std::min(a.size(), split));
    b0.assign(b.begin(), b.begin() + std::min(b.size(), split));
    a1.clear();
    b1.clear();
    if (a.size() > split) a1.assign(a.begin() + split, a.end());
    if (b.size() > split) b1.assign(b.begin() + split, b.end());

    auto& z0 = scratch.Get(level, 4);
    auto& z1 = scratch.Get(level, 5);
    auto& z2 = scratch.Get(level, 6);
    z1.clear();
    z2.clear();
    TraceModPolyMulFull(a0, b0, z0, field, scratch, level + 1);
    if (!a1.empty() && !b1.empty()) TraceModPolyMulFull(a1, b1, z2, field, scratch, level + 1);

    size_t asum_len = std::max(a0.size(), a1.size());
    size_t bsum_len = std::max(b0.size(), b1.size());
    auto& asum = scratch.Get(level, 7);
    auto& bsum = scratch.Get(level, 8);
    asum.assign(asum_len, 0);
    bsum.assign(bsum_len, 0);
    for (size_t i = 0; i < asum_len; ++i) {
        if (i < a0.size()) asum[i] ^= a0[i];
        if (i < a1.size()) asum[i] ^= a1[i];
    }
    for (size_t i = 0; i < bsum_len; ++i) {
        if (i < b0.size()) bsum[i] ^= b0[i];
        if (i < b1.size()) bsum[i] ^= b1[i];
    }
    while (!asum.empty() && asum.back() == 0) asum.pop_back();
    while (!bsum.empty() && bsum.back() == 0) bsum.pop_back();
    if (!asum.empty() && !bsum.empty()) TraceModPolyMulFull(asum, bsum, z1, field, scratch, level + 1);

    out.assign(a.size() + b.size() - 1, 0);
    for (size_t i = 0; i < z0.size(); ++i) out[i] ^= z0[i];
    for (size_t i = 0; i < z2.size(); ++i) out[2 * split + i] ^= z2[i];
    size_t mid_len = std::max(z1.size(), std::max(z0.size(), z2.size()));
    for (size_t i = 0; i < mid_len && split + i < out.size(); ++i) {
        Elem v = (i < z1.size() ? z1[i] : Elem(0)) ^ (i < z0.size() ? z0[i] : Elem(0)) ^ (i < z2.size() ? z2[i] : Elem(0));
        out[split + i] ^= v;
    }
    while (!out.empty() && out.back() == 0) out.pop_back();
}

/** Compute f(x)*g(x)^2 mod x^n into out. Writing y=x^2, g(x)^2=h(y), and
 * f(x)=f_even(y)+x*f_odd(y), so the product is f_even(y)*h(y) + x*f_odd(y)*h(y).
 * out must not alias f, g, or a scratch slot. */
template<typename F>
void TraceModMulBySquareLow(const std::vector<typename F::Elem>& f, const std::vector<typename F::Elem>& g, std::vector<typename F::Elem>& out, size_t n, const F& field, TraceModScratch<F>& scratch, size_t level) {
    out.clear();
    if (n == 0 || f.empty() || g.empty()) return;
    size_t even_n = (n + 1) / 2;
    size_t odd_n = n / 2;
    auto& h = scratch.Get(level, 0);
    h.resize(std::min(g.size(), even_n));
    for (size_t i = 0; i < h.size(); ++i) h[i] = field.Sqr(g[i]);
    while (!h.empty() && h.back() == 0) h.pop_back();
    if (h.empty()) return;

    auto& f_even = scratch.Get(level, 1);
    auto& f_odd = scratch.Get(level, 2);
    f_even.clear();
    f_odd.clear();
    for (size_t i = 0; i < f.size() && i < n; ++i) {
        if (i & 1) f_odd.push_back(f[i]);
        else f_even.push_back(f[i]);
    }
    while (!f_even.empty() && f_even.back() == 0) f_even.pop_back();
    while (!f_odd.empty() && f_odd.back() == 0) f_odd.pop_back();

    auto& even_prod = scratch.Get(level, 3);
    auto& odd_prod = scratch.Get(level, 4);
    TraceModPolyMulLow(f_even, h, even_prod, even_n, field, scratch, level + 1);
    TraceModPolyMulLow(f_odd, h, odd_prod, odd_n, field, scratch, level + 1);
    out.assign(n, 0);
    for (size_t i = 0; i < even_prod.size() && 2 * i < n; ++i) out[2 * i] = even_prod[i];
    for (size_t i = 0; i < odd_prod.size() && 2 * i + 1 < n; ++i) out[2 * i + 1] = odd_prod[i];
    while (!out.empty() && out.back() == 0) out.pop_back();
}

/** Compute 1/f mod x^n into out, assuming f(0)=1. In characteristic 2,
 * Newton iteration is g' = f*g^2 mod x^(2m): if f*g = 1+e, then f*g' = (f*g)^2 = 1+e^2.
 * Uses scratch levels 0 and up; out must not alias f or a scratch slot. */
template<typename F>
void TraceModInvSeries(const std::vector<typename F::Elem>& f, size_t n, std::vector<typename F::Elem>& out, const F& field, TraceModScratch<F>& scratch) {
    CHECK_SAFE(n > 0 && !f.empty() && f[0] == 1);
    out.assign(1, 1);
    auto& next_g = scratch.Get(0, 5);
    size_t m = 1;
    while (m < n) {
        size_t next = std::min(2 * m, n);
        TraceModMulBySquareLow(f, out, next_g, next, field, scratch, 1);
        if (next_g.empty()) next_g.push_back(0);
        out.swap(next_g);
        m = next;
    }
    out.resize(n, 0);
}

/* Scratch-free convenience wrappers (used by the unit and fuzz test layers). */
template<typename F>
void TraceModPolyMulLow(const std::vector<typename F::Elem>& a, const std::vector<typename F::Elem>& b, std::vector<typename F::Elem>& out, size_t n, const F& field) {
    TraceModScratch<F> scratch;
    TraceModPolyMulLow(a, b, out, n, field, scratch, 0);
}

template<typename F>
void TraceModPolyMulFull(const std::vector<typename F::Elem>& a, const std::vector<typename F::Elem>& b, std::vector<typename F::Elem>& out, const F& field) {
    TraceModScratch<F> scratch;
    TraceModPolyMulFull(a, b, out, field, scratch, 0);
}

template<typename F>
std::vector<typename F::Elem> TraceModMulBySquareLow(const std::vector<typename F::Elem>& f, const std::vector<typename F::Elem>& g, size_t n, const F& field) {
    TraceModScratch<F> scratch;
    std::vector<typename F::Elem> out;
    TraceModMulBySquareLow(f, g, out, n, field, scratch, 0);
    return out;
}

template<typename F>
std::vector<typename F::Elem> TraceModInvSeries(const std::vector<typename F::Elem>& f, size_t n, const F& field) {
    TraceModScratch<F> scratch;
    std::vector<typename F::Elem> out;
    TraceModInvSeries(f, n, out, field, scratch);
    return out;
}

/** The FastDivMod engage threshold for a field: the measured win requires
 * the FFT tier and a cheap scalar Multiplier (the strong schoolbook
 * baseline is hard to beat otherwise); ineligible fields always take the
 * schoolbook path. */
template<typename F>
size_t FastDivModCutoff(const F& field) {
    if (!TraceModFFTCapable(field) || !TraceModCheapScalarMul<F>()) return std::numeric_limits<size_t>::max();
    return FAST_DIVMOD_CUTOFF;
}

/** Subquadratic DivMod, same contract as DivMod: divide val by monic mod,
 * quotient into div (leading coefficient nonzero), remainder into val (size
 * mod.size() - 1, not stripped). Uses the reciprocal method of
 * ReciprocalReduce as a one-shot division: the reversed quotient is
 * rev(val) * InvSeries(rev(mod)) mod x^k, and the remainder follows from one
 * more low product. Falls back to the schoolbook DivMod when the smaller of
 * quotient length and divisor degree is at most `cutoff` (pass
 * FastDivModCutoff(field); tests pass smaller values to exercise the
 * reciprocal path at small sizes). */
template<typename F>
NOINLINE void FastDivMod(const std::vector<typename F::Elem>& mod, std::vector<typename F::Elem>& val, std::vector<typename F::Elem>& div, const F& field, size_t cutoff) {
    typedef typename F::Elem Elem;
    size_t m = mod.size();
    CHECK_SAFE(m > 0 && mod.back() == 1);
    if (val.size() < m) {
        div.clear();
        return;
    }
    size_t k = val.size() - m + 1;
    if (std::min(k, m - 1) <= cutoff) {
        DivMod(mod, val, div, field);
        return;
    }
    TraceModScratch<F> scratch;
    std::vector<Elem> rev_mod(m);
    for (size_t i = 0; i < m; ++i) rev_mod[i] = mod[m - 1 - i];
    std::vector<Elem> inv;
    TraceModInvSeries(rev_mod, k, inv, field, scratch);
    std::vector<Elem> rev_val(k);
    for (size_t i = 0; i < k; ++i) rev_val[i] = val[val.size() - 1 - i];
    std::vector<Elem> q_rev;
    TraceModPolyMulLow(rev_val, inv, q_rev, k, field, scratch, 0);
    div.resize(k);
    for (size_t i = 0; i < k; ++i) div[k - 1 - i] = q_rev[i];
    std::vector<Elem> prod;
    TraceModPolyMulLow(div, mod, prod, m - 1, field, scratch, 0);
    val.resize(m - 1);
    for (size_t i = 0; i < m - 1; ++i) val[i] ^= prod[i];
}

/* ---------- Subquadratic GCD (half-gcd) ----------
 *
 * FastGCD computes gcd like GCD, but reduces large inputs with the classical
 * recursive half-gcd: the transition matrix of a stretch of the Euclidean
 * remainder sequence depends only on the top coefficients of the operands
 * (see von zur Gathen & Gerhard ch. 11; the recursion shape follows NTL's
 * GF2EX HalfGCD), so it can be computed at half size recursively and applied
 * to the full-length pair with four multiplications, giving O(M(d) log d)
 * total. Only the gcd itself is needed (no Bezout coefficients), so no
 * cofactor back-substitution is performed. */

/** Remove trailing zero coefficients. */
template<typename E>
void PolyStrip(std::vector<E>& a) {
    while (!a.empty() && a.back() == 0) a.pop_back();
}

/** a ^= b, growing a as needed; result stripped. */
template<typename E>
void PolyAddInto(std::vector<E>& a, const std::vector<E>& b) {
    if (b.size() > a.size()) a.resize(b.size(), 0);
    for (size_t i = 0; i < b.size(); ++i) a[i] ^= b[i];
    PolyStrip(a);
}

/** out = a*b (stripped), tolerating empty operands. out must not alias. */
template<typename F>
void PolyMulStripped(const std::vector<typename F::Elem>& a, const std::vector<typename F::Elem>& b, std::vector<typename F::Elem>& out, const F& field, TraceModScratch<F>& scratch) {
    if (a.empty() || b.empty()) {
        out.clear();
        return;
    }
    TraceModPolyMulFull(a, b, out, field, scratch, 0);
    PolyStrip(out);
}

/** Transforms of the operands of a 2x2 matrix-vector or matrix-matrix
 * product, one per distinct operand, at a size fitting products of length
 * max_prod (PolyDot2FromFFT then forms each sum of two products with one
 * interpolation): 6 + 2 transforms instead of 4 * 3 for the matrix-vector
 * product, 8 + 4 instead of 8 * 3 for the matrix product. Returns false
 * (nothing computed) unless the FFT tier applies to every product. */
template<typename F>
bool PolyMat22FFT(const std::vector<const std::vector<typename F::Elem>*>& ops, size_t max_prod, const F& field, TraceModScratch<F>& scratch, std::vector<std::vector<typename F::Elem>>& tr) {
    if (!TraceModFFTCapable(field) || !TraceModCheapScalarMul<F>()) return false;
    if (max_prod < TRACEMOD_FFT_CUTOFF) return false;
    for (auto o : ops) {
        if (o->size() <= TRACEMOD_POLYMUL_CUTOFF) return false;
    }
    int m = 1;
    while ((size_t(1) << m) < max_prod) ++m;
    if (m > field.Bits()) return false;
    tr.resize(ops.size());
    for (size_t i = 0; i < ops.size(); ++i) scratch.fft.Prepare(*ops[i], m, field, tr[i]);
    return true;
}

/** out = a0*b0 + a1*b1 from transforms (size 2^m), interpolated and
 * stripped. */
template<typename F>
void PolyDot2FromFFT(const std::vector<typename F::Elem>& a0, const std::vector<typename F::Elem>& b0, const std::vector<typename F::Elem>& a1, const std::vector<typename F::Elem>& b1, std::vector<typename F::Elem>& out, const F& field, TraceModScratch<F>& scratch) {
    const size_t n = a0.size();
    int m = 0;
    while ((size_t(1) << m) < n) ++m;
    out.resize(n);
    for (size_t i = 0; i < n; ++i) out[i] = field.Reduce(field.WideXor(field.MulWide(a0[i], b0[i]), field.MulWide(a1[i], b1[i])));
    scratch.fft.IFFT(out, m, field);
    PolyStrip(out);
}

/** Quotient and remainder for a not-necessarily-monic divisor: den must be
 * nonempty with den.back() != 0; num (stripped) becomes the stripped
 * remainder and q the quotient. Schoolbook; used for the single division
 * steps inside the half-gcd, whose quotients are almost always tiny. */
template<typename F>
void DivModNonMonic(const std::vector<typename F::Elem>& den, std::vector<typename F::Elem>& num, std::vector<typename F::Elem>& q, const F& field) {
    typedef typename F::Elem Elem;
    CHECK_SAFE(!den.empty() && den.back() != 0);
    if (num.size() < den.size()) {
        q.clear();
        return;
    }
    Elem ilead = den.back() == 1 ? Elem(1) : field.Inv(den.back());
    q.assign(num.size() - den.size() + 1, 0);
    while (num.size() >= den.size()) {
        Elem term = ilead == 1 ? num.back() : field.Mul(num.back(), ilead);
        q[num.size() - den.size()] = term;
        num.pop_back();
        if (term != 0) {
            typename F::Multiplier mul(field, term);
            for (size_t x = 0; x + 1 < den.size(); ++x) {
                num[num.size() - den.size() + 1 + x] ^= mul(den[x]);
            }
        }
    }
    PolyStrip(num);
}

/** 2x2 matrix of polynomials: the accumulated transition of a stretch of
 * the remainder sequence, with (u', v')^T = M * (u, v)^T. */
template<typename F>
struct PolyMat22 {
    std::vector<typename F::Elem> e[2][2];
    void SetIdentity() {
        e[0][0].assign(1, 1);
        e[0][1].clear();
        e[1][0].clear();
        e[1][1].assign(1, 1);
    }
};

/** Fold one division step with quotient q into M: the pair update
 * (u, v) -> (v, u - q*v) corresponds to rows (m0, m1) -> (m1, m0 + q*m1)
 * in characteristic 2. */
template<typename F>
void PolyMat22Step(PolyMat22<F>& mat, const std::vector<typename F::Elem>& q, const F& field, TraceModScratch<F>& scratch) {
    std::vector<typename F::Elem> tmp;
    for (int j = 0; j < 2; ++j) {
        PolyMulStripped(q, mat.e[1][j], tmp, field, scratch);
        PolyAddInto(mat.e[0][j], tmp);
        mat.e[0][j].swap(mat.e[1][j]);
    }
}

/** (u, v) <- M * (u, v), results stripped. */
template<typename F>
void PolyMat22Apply(const PolyMat22<F>& mat, std::vector<typename F::Elem>& u, std::vector<typename F::Elem>& v, const F& field, TraceModScratch<F>& scratch) {
    std::vector<typename F::Elem> t00, t01, t10, t11;
    if (!u.empty() && !v.empty() && !mat.e[0][0].empty() && !mat.e[0][1].empty() && !mat.e[1][0].empty() && !mat.e[1][1].empty()) {
        size_t lm = std::max(std::max(mat.e[0][0].size(), mat.e[0][1].size()), std::max(mat.e[1][0].size(), mat.e[1][1].size()));
        size_t max_prod = lm + std::max(u.size(), v.size()) - 1;
        std::vector<std::vector<typename F::Elem>> tr;
        if (PolyMat22FFT<F>({&mat.e[0][0], &mat.e[0][1], &mat.e[1][0], &mat.e[1][1], &u, &v}, max_prod, field, scratch, tr)) {
            PolyDot2FromFFT(tr[0], tr[4], tr[1], tr[5], u, field, scratch);
            PolyDot2FromFFT(tr[2], tr[4], tr[3], tr[5], v, field, scratch);
            return;
        }
    }
    PolyMulStripped(mat.e[0][0], u, t00, field, scratch);
    PolyMulStripped(mat.e[0][1], v, t01, field, scratch);
    PolyMulStripped(mat.e[1][0], u, t10, field, scratch);
    PolyMulStripped(mat.e[1][1], v, t11, field, scratch);
    PolyAddInto(t00, t01);
    PolyAddInto(t10, t11);
    u.swap(t00);
    v.swap(t10);
}

/** out = a * b (matrix product). out must not alias a or b. */
template<typename F>
void PolyMat22Mul(PolyMat22<F>& out, const PolyMat22<F>& a, const PolyMat22<F>& b, const F& field, TraceModScratch<F>& scratch) {
    std::vector<typename F::Elem> tmp;
    bool all = true;
    size_t la = 0, lb = 0;
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            all = all && !a.e[i][j].empty() && !b.e[i][j].empty();
            la = std::max(la, a.e[i][j].size());
            lb = std::max(lb, b.e[i][j].size());
        }
    }
    std::vector<std::vector<typename F::Elem>> tr;
    if (all && PolyMat22FFT<F>({&a.e[0][0], &a.e[0][1], &a.e[1][0], &a.e[1][1], &b.e[0][0], &b.e[0][1], &b.e[1][0], &b.e[1][1]}, la + lb - 1, field, scratch, tr)) {
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) PolyDot2FromFFT(tr[2 * i], tr[4 + j], tr[2 * i + 1], tr[6 + j], out.e[i][j], field, scratch);
        }
        return;
    }
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            PolyMulStripped(a.e[i][0], b.e[0][j], out.e[i][j], field, scratch);
            PolyMulStripped(a.e[i][1], b.e[1][j], tmp, field, scratch);
            PolyAddInto(out.e[i][j], tmp);
        }
    }
}

/** Iterative base case: Euclidean steps on (u, v) with matrix accumulation,
 * until v is empty or deg(v) <= stop_deg. u and v are consumed. */
template<typename F>
void IterHalfGCD(PolyMat22<F>& mat, std::vector<typename F::Elem>& u, std::vector<typename F::Elem>& v, size_t stop_deg, const F& field, TraceModScratch<F>& scratch) {
    mat.SetIdentity();
    std::vector<typename F::Elem> q;
    while (!v.empty() && v.size() - 1 > stop_deg) {
        DivModNonMonic(v, u, q, field);
        u.swap(v);
        PolyMat22Step(mat, q, field, scratch);
    }
}

/** Compute the transition matrix M of the remainder sequence of (u, v) —
 * requires deg(u) > deg(v) — up to the first pair (u', v') with
 * deg(v') <= deg(u) - d_red. Only the top 2*d_red - 1 coefficients of u
 * and v influence M (the classical half-gcd truncation lemma), which is
 * what makes the half-size recursion valid. u and v are unchanged; the
 * caller applies M to them (or to the full-length pair they were
 * truncated from). */
template<typename F>
void HalfGCDMatrix(PolyMat22<F>& mat, const std::vector<typename F::Elem>& u_in, const std::vector<typename F::Elem>& v_in, size_t d_red, const F& field, TraceModScratch<F>& scratch) {
    typedef typename F::Elem Elem;
    CHECK_SAFE(!u_in.empty() && (v_in.empty() || v_in.size() < u_in.size()));
    if (v_in.empty() || v_in.size() - 1 + d_red <= u_in.size() - 1) {
        mat.SetIdentity();
        return;
    }
    // Work on the top 2*d_red - 1 coefficients only.
    size_t deg_u = u_in.size() - 1;
    size_t shift = deg_u >= 2 * d_red - 2 ? deg_u - (2 * d_red - 2) : 0;
    std::vector<Elem> u(u_in.begin() + shift, u_in.end());
    std::vector<Elem> v;
    if (v_in.size() > shift) v.assign(v_in.begin() + shift, v_in.end());
    size_t du = u.size() - 1;
    if (d_red <= HGCD_ITER_CROSSOVER) {
        IterHalfGCD(mat, u, v, du - d_red, field, scratch);
        return;
    }
    // First half of the reduction, recursively at half size.
    size_t d1 = (d_red + 1) / 2;
    PolyMat22<F> m1;
    HalfGCDMatrix(m1, u, v, d1, field, scratch);
    PolyMat22Apply(m1, u, v, field, scratch);
    if (v.empty() || v.size() - 1 + d_red <= du) {
        mat = m1;
        return;
    }
    // One explicit division step to cross the halfway point.
    std::vector<Elem> q;
    DivModNonMonic(v, u, q, field);
    u.swap(v);
    PolyMat22Step(m1, q, field, scratch);
    if (v.empty() || v.size() - 1 + d_red <= du) {
        mat = m1;
        return;
    }
    // Second half: the remaining degrees to shave, again at half size. The
    // sub-call's target is relative to its own leading degree, so the
    // remaining reduction is measured against the current u, giving a final
    // deg(v) <= (u.size() - 1) - d2 == du - d_red as required.
    size_t d2 = u.size() - 1 + d_red - du;
    PolyMat22<F> m2;
    HalfGCDMatrix(m2, u, v, d2, field, scratch);
    PolyMat22Mul(mat, m2, m1, field, scratch);
}

/** The FastGCD engage threshold for a field, on the same eligibility terms
 * as FastDivModCutoff. The measured crossover against the quadratic loop is
 * high: the lazily-reducing, normalization-free Euclid loop (GCD) beats the
 * half-gcd's matrix machinery at every degree up to 4096 on an i7-4770 and
 * a Pi 5 (by ~10% at 2048 and 4096). */
template<typename F>
size_t FastGCDCutoff(const F& field) {
    if (!TraceModFFTCapable(field) || !TraceModCheapScalarMul<F>()) return std::numeric_limits<size_t>::max();
    return HGCD_CUTOFF;
}

/** Compute the GCD of a and b like GCD (result in a, b destroyed, inputs
 * stripped), using half-gcd reductions while operands are longer than
 * `cutoff` (pass FastGCDCutoff(field); tests pass smaller values to
 * exercise the half-gcd at small sizes), and the quadratic loop below. */
template<typename F>
NOINLINE void FastGCD(std::vector<typename F::Elem>& a, std::vector<typename F::Elem>& b, const F& field, size_t cutoff) {
    typedef typename F::Elem Elem;
    PolyStrip(a);
    PolyStrip(b);
    if (a.size() < b.size()) std::swap(a, b);
    if (b.size() <= cutoff) {
        GCD(a, b, field);
        return;
    }
    TraceModScratch<F> scratch;
    PolyMat22<F> mat;
    std::vector<Elem> q;
    while (!b.empty() && b.size() > cutoff) {
        size_t deg_a = a.size() - 1;
        size_t d_red = (deg_a + 1) / 2;
        if (b.size() == a.size() || b.size() - 1 + d_red <= deg_a) {
            // Equal degrees, or the gap already exceeds the half-gcd target:
            // take one division step (subquadratic for large quotients).
            MakeMonic(b, field);
            FastDivMod(b, a, q, field, cutoff);
            PolyStrip(a);
            std::swap(a, b);
            continue;
        }
        HalfGCDMatrix(mat, a, b, d_red, field, scratch);
        PolyMat22Apply(mat, a, b, field, scratch);
        // The half-gcd target guarantees progress: deg(b) is now at most
        // deg_a - d_red < deg_a / 2 + 1 (or b is empty).
        CHECK_SAFE(b.empty() || b.size() - 1 + d_red <= deg_a);
    }
    GCD(a, b, field);
}

/** Compute repeated TraceMod operations with a fixed modulus.
 * For smaller degrees it precomputes rows x^e mod mod(x) for even e >= deg(mod), so
 * squaring sum a_i*x^i can be reduced by adding a_i^2*(x^(2i) mod mod). For larger
 * degrees, division by monic mod(x) uses reciprocal series: reverse the high part of val,
 * multiply by 1/reverse(mod) mod x^k to get the reversed quotient, then subtract q*mod;
 * subtraction is XOR in characteristic 2. */
template<typename F>
class TraceMod {
    typedef typename F::Elem Elem;
    const std::vector<Elem>& mod;
    const F& field;
    bool use_square_table;

    size_t d, first_even;
    std::vector<Elem> rows;
    // Lazy-reduction layout of rows (see the constructor).
    size_t num_rows = 0;
    std::vector<Elem> blk, tail;

    std::vector<Elem> rev_mod, inv;

    // Cached 2^fft_m-point transforms of the two per-node-constant reducer
    // operands (inv and mod), or fft_m == 0 when the FFT tier does not apply
    // to this node. With them, a full-size ReciprocalReduce costs four
    // transforms instead of six.
    int fft_m = 0;
    std::vector<Elem> inv_fft, mod_fft;

    // Reused buffers for SquareReduce/ReciprocalReduce, and the scratch pool
    // for the Karatsuba multiplications, hoisted here so the trace loop is
    // allocation-free in steady state.
    std::vector<Elem> sqr_tmp;
    std::vector<Elem> sqr_coeffs;
    std::vector<Elem> rev_val, q_rev, quot, prod;
    TraceModScratch<F> scratch;

public:
    TraceMod(const std::vector<Elem>& mod_in, const F& field_in) : mod(mod_in), field(field_in), use_square_table(mod_in.size() - 1 < TraceModTableCutoff(field_in)), d(mod_in.size() - 1), first_even(0) {
        CHECK_SAFE(!mod.empty() && mod.back() == 1);
        // RecFindRoots handles degree 1 and 2 before constructing this object.
        CHECK_SAFE(d >= 3);
        if (use_square_table) {
            size_t max_e = 2 * d - 2;
            first_even = (d & 1) ? d + 1 : d;
            const size_t nrows = first_even <= max_e ? (max_e - first_even) / 2 + 1 : 0;
            // Lazily-reducing fields store the rows directly in the
            // column-block-major layout their reducer walks (all rows per
            // block of 8 output columns: contiguous instead of one cache line
            // per row 8*d bytes apart), with the d % 8 tail columns in a
            // small row-major array; other fields row-major.
            const size_t nb = d / 8, tl = d % 8;
            if (F::LAZY_REDUCE) {
                num_rows = nrows;
                blk.assign(nb * nrows * 8, 0);
                tail.assign(nrows * tl, 0);
            } else {
                rows.assign(nrows * d, 0);
            }
            std::vector<Elem> rem(d, 0);
            for (size_t j = 0; j < d; ++j) rem[j] = mod[j];
            for (size_t e = d; e <= max_e; ++e) {
                if ((e & 1) == 0 && e >= first_even) {
                    const size_t r = (e - first_even) / 2;
                    if (F::LAZY_REDUCE) {
                        for (size_t b = 0; b < nb; ++b) std::copy(rem.begin() + 8 * b, rem.begin() + 8 * b + 8, blk.begin() + (b * nrows + r) * 8);
                        std::copy(rem.begin() + 8 * nb, rem.end(), tail.begin() + r * tl);
                    } else {
                        std::copy(rem.begin(), rem.end(), rows.begin() + r * d);
                    }
                }
                if (e == max_e) break;
                Elem fold = rem[d - 1];
                for (size_t j = d - 1; j > 0; --j) rem[j] = rem[j - 1];
                rem[0] = 0;
                if (fold != 0) {
                    if (fold == 1) {
                        for (size_t j = 0; j < d; ++j) rem[j] ^= mod[j];
                    } else {
                        typename F::Multiplier mul(field, fold);
                        for (size_t j = 0; j < d; ++j) rem[j] ^= mul(mod[j]);
                    }
                }
            }
        } else {
            size_t m = mod.size();
            rev_mod.resize(m);
            for (size_t i = 0; i < m; ++i) rev_mod[i] = mod[m - 1 - i];
            TraceModInvSeries(rev_mod, m - 1, inv, field, scratch);
            // Every product in ReciprocalReduce has length at most 2*d, so a
            // fixed 2^fft_m >= 2*d transform size fits them all; the two
            // one-time Prepare transforms amortize over the (Bits-1)-round
            // trace loop and the fully-factorizable test.
            if (TraceModFFTCapable(field) && 2 * d - 1 >= TRACEMOD_FFT_CUTOFF) {
                int t = 1;
                while ((size_t(1) << t) < 2 * d && t < 64) ++t;
                if (t <= field.Bits()) {
                    fft_m = t;
                    scratch.fft.Prepare(inv, fft_m, field, inv_fft);
                    scratch.fft.Prepare(mod, fft_m, field, mod_fft);
                }
            }
        }
    }

private:

    const Elem* Row(size_t e) const {
        CHECK_SAFE(e >= first_even && ((e - first_even) & 1) == 0);
        return rows.data() + ((e - first_even) / 2) * d;
    }

    void SquareReduce(std::vector<Elem>& val) {
        sqr_tmp.assign(d, 0);
        if (F::LAZY_REDUCE) {
            // Accumulate the row contributions unreduced (one reduction per
            // output coefficient instead of one per row entry), in column
            // blocks of 8 whose accumulators stay in registers; the table
            // (column-block-major) is streamed once per call.
            typedef typename F::Wide Wide;
            sqr_coeffs.assign(num_rows, 0);
            size_t used = 0;
            for (size_t i = 0; i < val.size(); ++i) {
                if (val[i] == 0) continue;
                Elem coeff = field.Sqr(val[i]);
                size_t e = 2 * i;
                if (e < d) {
                    sqr_tmp[e] ^= coeff;
                } else {
                    size_t r = (e - first_even) / 2;
                    sqr_coeffs[r] = coeff;
                    used = r + 1;
                }
            }
            if (used) {
                const Elem* c = sqr_coeffs.data();
                const size_t nb = d / 8, tl = d % 8;
                for (size_t b = 0; b < nb; ++b) {
                    const Elem* p = blk.data() + b * num_rows * 8;
                    Wide a0 = field.WideZero(), a1 = a0, a2 = a0, a3 = a0, a4 = a0, a5 = a0, a6 = a0, a7 = a0;
                    for (size_t r = 0; r < used; ++r, p += 8) {
                        Elem cr = c[r];
                        a0 = field.WideXor(a0, field.MulWide(cr, p[0]));
                        a1 = field.WideXor(a1, field.MulWide(cr, p[1]));
                        a2 = field.WideXor(a2, field.MulWide(cr, p[2]));
                        a3 = field.WideXor(a3, field.MulWide(cr, p[3]));
                        a4 = field.WideXor(a4, field.MulWide(cr, p[4]));
                        a5 = field.WideXor(a5, field.MulWide(cr, p[5]));
                        a6 = field.WideXor(a6, field.MulWide(cr, p[6]));
                        a7 = field.WideXor(a7, field.MulWide(cr, p[7]));
                    }
                    Elem* out = sqr_tmp.data() + 8 * b;
                    out[0] ^= field.Reduce(a0);
                    out[1] ^= field.Reduce(a1);
                    out[2] ^= field.Reduce(a2);
                    out[3] ^= field.Reduce(a3);
                    out[4] ^= field.Reduce(a4);
                    out[5] ^= field.Reduce(a5);
                    out[6] ^= field.Reduce(a6);
                    out[7] ^= field.Reduce(a7);
                }
                for (size_t t = 0; t < tl; ++t) {
                    Wide a = field.WideZero();
                    for (size_t r = 0; r < used; ++r) a = field.WideXor(a, field.MulWide(c[r], tail[r * tl + t]));
                    sqr_tmp[8 * nb + t] ^= field.Reduce(a);
                }
            }
            val.swap(sqr_tmp);
            while (!val.empty() && val.back() == 0) val.pop_back();
            return;
        }
        for (size_t i = 0; i < val.size(); ++i) {
            if (val[i] == 0) continue;
            Elem coeff = field.Sqr(val[i]);
            size_t e = 2 * i;
            if (e < d) {
                sqr_tmp[e] ^= coeff;
            } else {
                const Elem* row = Row(e);
                if (coeff == 1) {
                    for (size_t j = 0; j < d; ++j) sqr_tmp[j] ^= row[j];
                } else {
                    typename F::Multiplier mul(field, coeff);
                    for (size_t j = 0; j < d; ++j) sqr_tmp[j] ^= mul(row[j]);
                }
            }
        }
        // The swap donates val's old buffer back as the next call's scratch.
        val.swap(sqr_tmp);
        while (!val.empty() && val.back() == 0) val.pop_back();
    }

    /** Whether MulLowPrepared beats the dynamic multiplication path for a
     * product of (truncated) operand lengths la and lb: exactly when the
     * dynamic path would run the FFT at the node transform size anyway, so
     * the cached operand transform saves one of its three transforms. For
     * anything smaller the dynamic path (half-size FFT, Karatsuba or naive)
     * is cheaper than a fixed full-size transform. */
    bool CachedMulEligible(size_t la, size_t lb) const {
        if (fft_m == 0) return false;
        if (!TraceModFFTEligible(field, la, lb)) return false;
        return la + lb - 1 > (size_t(1) << (fft_m - 1));
    }

    void ReciprocalReduce(std::vector<Elem>& val) {
        size_t m = mod.size();
        if (val.size() < m) return;
        size_t k = val.size() - m + 1;
        rev_val.resize(k);
        for (size_t i = 0; i < k; ++i) rev_val[i] = val[val.size() - 1 - i];
        if (CachedMulEligible(k, std::min(inv.size(), k))) {
            scratch.fft.MulLowPrepared(rev_val, inv_fft, inv.size(), fft_m, q_rev, k, field);
        } else {
            TraceModPolyMulLow(rev_val, inv, q_rev, k, field, scratch, 0);
        }
        quot.resize(k);
        for (size_t i = 0; i < k; ++i) quot[k - 1 - i] = q_rev[i];
        while (!quot.empty() && quot.back() == 0) quot.pop_back();
        size_t rem_len = m - 1;
        if (CachedMulEligible(quot.size(), std::min(mod.size(), rem_len))) {
            scratch.fft.MulLowPrepared(quot, mod_fft, mod.size(), fft_m, prod, rem_len, field);
        } else {
            TraceModPolyMulLow(quot, mod, prod, rem_len, field, scratch, 0);
        }
        for (size_t i = 0; i < prod.size(); ++i) val[i] ^= prod[i];
        val.resize(rem_len);
        while (!val.empty() && val.back() == 0) val.pop_back();
    }

public:
    /** val := val^2 mod mod. Requires val to be reduced already (val.size() <= deg(mod)). */
    void SquareAndReduce(std::vector<Elem>& val) {
        CHECK_SAFE(val.size() <= d);
        if (use_square_table) {
            SquareReduce(val);
        } else {
            Sqr(val, field);
            ReciprocalReduce(val);
        }
    }

    void trace(std::vector<Elem>& out, const Elem& param) {
        if (use_square_table) {
            out.resize(2);
            out[0] = 0;
            out[1] = param;
            for (int i = 0; i < field.Bits() - 1; ++i) {
                SquareReduce(out);
                if (out.size() < 2) out.resize(2);
                // SquareReduce already reduced modulo mod, so add rather than replace the x coefficient.
                out[1] ^= param;
                while (!out.empty() && out.back() == 0) out.pop_back();
            }
        } else {
            out.reserve(mod.size() * 2);
            out.resize(2);
            out[0] = 0;
            out[1] = param;
            for (int i = 0; i < field.Bits() - 1; ++i) {
                Sqr(out, field);
                if (out.size() < 2) out.resize(2);
                out[1] = param;
                ReciprocalReduce(out);
            }
        }
    }
};

/** An element generating the subfield GF(2^sub) of GF(2^Bits) (sub a
 * power of two dividing Bits): the relative trace sum_l y^(2^(sub*l)) of
 * y = x, x^2, ... maps into GF(2^sub); the first image outside the largest
 * proper subfield GF(2^(sub/2)) has degree sub over GF(2), so its powers
 * 1..w^(sub-1) are a GF(2)-basis of GF(2^sub). */
template<typename F>
typename F::Elem SubfieldGenerator(const F& field, int sub) {
    typedef typename F::Elem Elem;
    if (sub == 1) return 1;
    Elem y = 1;
    while (true) {
        y = field.Mul2(y);
        Elem w = 0, t = y;
        for (int l = 0; l < field.Bits() / sub; ++l) {
            w ^= t;
            for (int s = 0; s < sub; ++s) t = field.Sqr(t);
        }
        Elem wh = w;
        for (int s = 0; s < sub / 2; ++s) wh = field.Sqr(wh);
        if (wh != w) return w;
    }
}

template<typename F>
typename F::Elem SubfieldGenerator8(const F& field) { return SubfieldGenerator(field, 8); }

/* ---------- Montgomery arithmetic on additive-FFT cosets ----------
 *
 * For a monic g of degree <= N = 2^m and a coset C = W_m + eta_(b*2^m) of
 * the additive-FFT subspace (AdditiveFFT::FFT with block b), let Z_C be the
 * vanishing polynomial of C: in the Cantor basis
 *     Z_C(x) = sum_(j subset of m) x^(2^j) + eta_b
 * (the subspace polynomial S^m(x), S(x) = x^2 + x, which maps
 * eta_(b*2^m + j) to eta_b), monic of degree N with formal derivative
 * Z_C' = 1. If g has no root on C, Z_C is invertible mod g, and a residue r
 * is held in Montgomery form r*Z_C mod g as a polynomial of degree < N,
 * represented by its N values on C. The Montgomery reduction of a
 * polynomial a,
 *     q = a / g mod Z_C          (pointwise on C),
 *     r = (a + q*g) / Z_C        (exact; r = a / Z_C mod g),
 * has deg r < N if deg a < 2N, and deg r < deg g if deg a < N. On C both
 * sides of h = a + q*g = r * Z_C vanish; differentiating,
 * h' = r' * Z_C + r * Z_C' gives r(w) = h'(w) = a'(w) + q'(w)*g(w) +
 * q(w)*g'(w) for w in C. Only q' (and a', unless a is a square, whose
 * derivative is 0 in characteristic 2) is not available pointwise:
 * interpolating q and evaluating its formal derivative are the two size-N
 * transforms of a modular squaring, against roughly eight for the
 * reciprocal (Barrett) reducer, whose products have length up to 2*deg g.
 *
 * Halving the coset: C is the union of its halves C0 (block 2b) and C1
 * (block 2b+1) of W_(m-1), with Z_C = Z_C0 * Z_C1 and Z_C0 = 1 on C1. So for
 * a factor g of the modulus with deg g <= N/2, the Montgomery form v of a
 * residue T on C (v = T*Z_C mod g as well) gives the one on C1 as
 *     v / Z_C0 mod g = (v + q*g) / Z_C0,  q = v / g mod Z_C0,
 * whose values on C1 are v + q*g pointwise: two size-N/2 transforms carry
 * a residue from a node of the root-finding tree to its child. */

/** out = 1 / in elementwise (Montgomery's batch inversion); false if some
 * element is 0. */
template<typename F>
bool BatchInvert(const std::vector<typename F::Elem>& in, std::vector<typename F::Elem>& out, const F& field) {
    typedef typename F::Elem Elem;
    const size_t n = in.size();
    out.resize(n);
    Elem prod = 1;
    for (size_t j = 0; j < n; ++j) {
        if (in[j] == 0) return false;
        out[j] = prod;
        prod = field.Mul(prod, in[j]);
    }
    prod = field.Inv(prod);
    for (size_t j = n; j-- > 0; ) {
        out[j] = field.Mul(out[j], prod);
        prod = field.Mul(prod, in[j]);
    }
    return true;
}

/** Montgomery arithmetic modulo a monic g on one coset (see above). */
template<typename F>
class MontCoset {
    typedef typename F::Elem Elem;
    const AdditiveFFT<F>& fft;
    const F& field;
    std::vector<Elem> ga, inv_ga, dga, q, da, tmp;

public:
    const int m;
    const size_t block, n;

    MontCoset(const AdditiveFFT<F>& fft_in, const F& field_in, int m_in, size_t block_in) : fft(fft_in), field(field_in), m(m_in), block(block_in), n(size_t(1) << m_in) {}

    /** Z_C(0) = eta_b, the coset's offset under S^m. */
    static Elem Offset(const AdditiveFFT<F>& fft, size_t block) {
        Elem e = 0;
        for (size_t k = 0; block >> k; ++k) {
            if ((block >> k) & 1) e ^= fft.Basis()[k];
        }
        return e;
    }

    /** out = the values of g (monic, degree <= 2^m) on coset (m, block). */
    static void Eval(const AdditiveFFT<F>& fft, const F& field, int m, size_t block, const std::vector<Elem>& g, std::vector<Elem>& out) {
        const size_t n = size_t(1) << m, d = g.size() - 1;
        CHECK_SAFE(d <= n);
        out.assign(n, 0);
        std::copy(g.begin(), g.begin() + std::min(d + 1, n), out.begin());
        if (d == n) {
            // Reduce mod Z_C first: subtract the monic Z_C.
            for (int j = 0; j < m; ++j) {
                if ((m & j) == j) out[size_t(1) << j] ^= 1;
            }
            out[0] ^= Offset(fft, block);
        }
        fft.ToNovel(out, m);
        fft.NovelFFT(out, m, field, block);
    }

    /** Set the modulus g (monic, degree <= 2^m), optionally with its values
     *  on the coset precomputed. False if g has a root on the coset. */
    bool SetModulus(const std::vector<Elem>& g, const std::vector<Elem>* g_values = nullptr) {
        if (g_values) {
            ga = *g_values;
        } else {
            Eval(fft, field, m, block, g, ga);
        }
        if (!BatchInvert(ga, inv_ga, field)) return false;
        dga.assign(n, 0);
        for (size_t i = 1; i < g.size(); i += 2) dga[i - 1] = g[i];
        fft.ToNovel(dga, m);
        fft.NovelFFT(dga, m, field, block);
        return true;
    }

    const std::vector<Elem>& Values() const { return ga; }

    /** out = the values of the formal derivative of the polynomial of degree
     *  < 2^m with values v. */
    void Derivative(const std::vector<Elem>& v, std::vector<Elem>& out) {
        tmp = v;
        fft.NovelIFFT(tmp, m, field, block);
        AdditiveFFT<F>::NovelDerivative(tmp, m, out);
        fft.NovelFFT(out, m, field, block);
    }

    /** v := v / Z_C mod g, values in and out. With square set, v must hold
     *  a square of degree < 2^(m+1) (the result has degree < 2^m), else any
     *  polynomial of degree < 2^m (the result has degree < deg g). */
    void Reduce(std::vector<Elem>& v, bool square) {
        if (!square) Derivative(v, da);
        for (size_t j = 0; j < n; ++j) v[j] = field.Mul(v[j], inv_ga[j]);
        Derivative(v, q);
        for (size_t j = 0; j < n; ++j) v[j] = field.Reduce(field.WideXor(field.MulWide(q[j], ga[j]), field.MulWide(v[j], dga[j])));
        if (!square) {
            for (size_t j = 0; j < n; ++j) v[j] ^= da[j];
        }
    }

    /** v := v^2 / Z_C mod g: Montgomery squaring. */
    void Square(std::vector<Elem>& v) {
        for (size_t j = 0; j < n; ++j) v[j] = field.Sqr(v[j]);
        Reduce(v, true);
    }

    /** Values v of the Montgomery form of a residue to its coefficients
     *  (normal form, degree < deg g = d): v is destroyed. */
    void ToNormal(std::vector<Elem>& v, size_t d, std::vector<Elem>& out) {
        Reduce(v, false);
        fft.NovelIFFT(v, m, field, block);
        fft.FromNovel(v, m);
#ifdef MINISKETCH_VERIFY
        for (size_t j = d; j < n; ++j) CHECK_SAFE(v[j] == 0);
#endif
        out.assign(v.begin(), v.begin() + d);
        PolyStrip(out);
    }

    /** Values v to the coefficients (degree < 2^m, stripped) of the same
     *  polynomial, without reduction: v is destroyed. */
    void ToCoefficients(std::vector<Elem>& v, std::vector<Elem>& out) {
        fft.NovelIFFT(v, m, field, block);
        fft.FromNovel(v, m);
        out.swap(v);
        PolyStrip(out);
    }
};

/** Shared state of the Montgomery root-finding tree: the transform tables
 * and the generator z of the tower subfield GF(2^MONT_TOWER_BITS). */
template<typename F>
struct MontTree {
    AdditiveFFT<F> fft;
    typename F::Elem z;
};

/** Tower subfield GF(2^MONT_TOWER_BITS) of the Montgomery root-finding
 * tree's split parameters, the number of split depths it carries relative
 * traces for, and the number of those traces. */
static const int MONT_TOWER_BITS = MINISKETCH_MONT_TOWER_BITS;
static const int CARRY_DEPTHS = 16;
static const int CARRY_TRACES = CARRY_DEPTHS / MONT_TOWER_BITS;

/** The relative traces T_l = Tr_(GF(2^Bits)/GF(2^s))(basis * xi^l * x)
 * (s = MONT_TOWER_BITS, xi = Mul2's generator, l < CARRY_TRACES) modulo a
 * node polynomial g of the root-finding tree. With the tower split
 * parameters beta_(sl+k) = basis * xi^l * z^k (z generating GF(2^s)), the
 * split trace of depth sl + k at the node is
 * Tr(beta_(sl+k) * x) = Tr_(GF(2^s)/GF(2))(z^k * T_l) =
 * sum_(r<s) (z^k T_l)^(2^r) mod g: s - 1 modular squarings. A smaller s
 * trades squarings for more carried traces. Large nodes
 * hold the T_l in Montgomery form on a coset (see above); nodes of degree
 * at most SMALL_CARRY_MAX hold them in normal form (coefficients), carried
 * to factors by plain division. T_l is empty once no descendant needs
 * it. */
template<typename F>
struct MontCarry {
    const MontTree<F>* tree = nullptr;
    bool normal = false;
    int m = 0;
    size_t block = 0;
    std::vector<typename F::Elem> t[CARRY_TRACES]; //!< Montgomery values on the coset, or coefficients if normal
    std::vector<typename F::Elem> ga;              //!< the node polynomial on the coset, if computed
};

/** Carry the Montgomery-form traces of a parent node to its factor g
 * (monic, 3 <= deg g < parent's degree) for the traces l >= l_min: halve
 * the coset while deg g fits in the half. gv, if given, holds g's values
 * on the parent's coset (both halves of the first halving). */
template<typename F>
void MontCarryToFactor(const MontCarry<F>& parent, const std::vector<typename F::Elem>& g, int l_min, MontCarry<F>& child, const F& field, const std::vector<typename F::Elem>* gv = nullptr) {
    typedef typename F::Elem Elem;
    const AdditiveFFT<F>& fft = parent.tree->fft;
    const size_t d = g.size() - 1;
    child.tree = parent.tree;
    child.normal = false;
    child.m = parent.m;
    child.block = parent.block;
    for (int l = 0; l < CARRY_TRACES; ++l) {
        if (l >= l_min) {
            child.t[l] = parent.t[l];
        } else {
            child.t[l].clear();
        }
    }
    child.ga.clear();
    std::vector<Elem> g0, inv_g0, q;
    while (child.m > 1 && d <= (size_t(1) << (child.m - 1))) {
        const int mh = child.m - 1;
        const size_t h = size_t(1) << mh, b0 = 2 * child.block, b1 = b0 + 1;
        if (gv && child.m == parent.m) {
            g0.assign(gv->begin(), gv->begin() + h);
            child.ga.assign(gv->begin() + h, gv->end());
        } else {
            MontCoset<F>::Eval(fft, field, mh, b0, g, g0);
            MontCoset<F>::Eval(fft, field, mh, b1, g, child.ga);
        }
        // g divides the root polynomial, which has no root on the root
        // coset, of which every descendant coset is a subset.
        bool ok = BatchInvert(g0, inv_g0, field);
        CHECK_SAFE(ok);
        (void)ok;
        for (int l = l_min; l < CARRY_TRACES; ++l) {
            auto& v = child.t[l];
            if (v.empty()) continue;
            q.resize(h);
            for (size_t j = 0; j < h; ++j) q[j] = field.Mul(v[j], inv_g0[j]);
            fft.NovelIFFT(q, mh, field, b0);
            fft.NovelFFT(q, mh, field, b1);
            for (size_t j = 0; j < h; ++j) q[j] = v[h + j] ^ field.Mul(q[j], child.ga[j]);
            v.swap(q);
        }
        child.m = mh;
        child.block = b1;
    }
    if (gv && child.m == parent.m) child.ga = *gv;
}

/** out = sum_(r<MONT_TOWER_BITS) (zk * t)^(2^r) mod g for t of degree
 * < deg g, all in coefficient form (for small g). */
template<typename F>
void TowerTraceSmall(const std::vector<typename F::Elem>& g, const std::vector<typename F::Elem>& t, typename F::Elem zk, std::vector<typename F::Elem>& out, const F& field) {
    typedef typename F::Elem Elem;
    const size_t n = g.size() - 1;
    std::vector<Elem> y(n, 0), sq(2 * n - 1);
    for (size_t j = 0; j < t.size(); ++j) y[j] = field.Mul(zk, t[j]);
    out = y;
    for (int r = 1; r < MONT_TOWER_BITS; ++r) {
        for (size_t j = 0; j < n; ++j) {
            sq[2 * j] = field.Sqr(y[j]);
            if (j + 1 < n) sq[2 * j + 1] = 0;
        }
        for (size_t e = 2 * n - 2; e >= n; --e) {
            const Elem c = sq[e];
            if (c == 0) continue;
            typename F::Multiplier mul(field, c);
            for (size_t k = 0; k < n; ++k) sq[e - n + k] ^= mul(g[k]);
        }
        for (size_t j = 0; j < n; ++j) {
            y[j] = sq[j];
            out[j] ^= y[j];
        }
    }
    PolyStrip(out);
}

/** out = sum_(r<MONT_TOWER_BITS) (zk * T)^(2^r) mod g from T's Montgomery values v on
 * mc's coset (mc's modulus g). Below degree MONT_NORMALIZE_MIN, out is
 * instead the Montgomery form itself (degree < 2^m, congruent to the trace
 * times Z_C, a unit mod g): its gcd with g, the only use of a split
 * trace, is the same, and the extra length costs the gcd's first division
 * less than the reduction would at those sizes. force_normal always
 * normalizes. */
template<typename F>
void MontTowerTrace(MontCoset<F>& mc, const std::vector<typename F::Elem>& v, size_t deg_g, typename F::Elem zk, std::vector<typename F::Elem>& out, const F& field, bool force_normal = false) {
    typedef typename F::Elem Elem;
    std::vector<Elem> y(mc.n), acc(mc.n);
    for (size_t j = 0; j < mc.n; ++j) acc[j] = y[j] = field.Mul(zk, v[j]);
    for (int r = 1; r < MONT_TOWER_BITS; ++r) {
        mc.Square(y);
        for (size_t j = 0; j < mc.n; ++j) acc[j] ^= y[j];
    }
    if (force_normal || deg_g >= MONT_NORMALIZE_MIN) {
        mc.ToNormal(acc, deg_g, out);
    } else {
        mc.ToCoefficients(acc, out);
    }
}

/** Frobenius orbit sums by Montgomery reduction on a coset (see above).
 *
 * With X_i = x^(2^i) mod f, accumulates (on the coset, in Montgomery form)
 * while x is squared Bits times the root's split trace
 * trace0 = Tr(basis * x) = sum_i basis^(2^i) * X_i, returned in normal
 * form, and the relative traces T_l = sum_(i = 0 mod s) (basis * xi^l)^(2^i)
 * * X_i (see MontCarry), returned in carry (whose tree is left to the caller)
 * with fft as its transform engine, which this function extends. Deeper
 * split traces come from the T_l (MontTowerTrace). X_Bits == x (full
 * factorizability) is tested by reducing the difference of the Montgomery
 * forms.
 *
 * Requires f to have no root on the coset; the four cosets of W_m in
 * W_(m+2) are tried in turn. Returns 1 (factorizable, acc filled), 0 (not
 * factorizable), or -1 when not applicable (caller falls back). */
template<typename F>
int MontgomeryOrbit(const std::vector<typename F::Elem>& poly, typename F::Elem basis, const F& field, std::vector<typename F::Elem>& trace0, AdditiveFFT<F>& fft, MontCarry<F>& carry) {
    typedef typename F::Elem Elem;
    typedef typename F::Wide Wide;
    const size_t d = poly.size() - 1;
    const int bits = field.Bits();
    int m = 1;
    while ((size_t(1) << m) < d) ++m;
    if (m + 2 > bits) return -1;
    const size_t n = size_t(1) << m;
    fft.Extend(field, m + 2);
    std::unique_ptr<MontCoset<F>> mc;
    for (size_t b : {1, 2, 3, 0}) {
        mc.reset(new MontCoset<F>(fft, field, m, b));
        if (mc->SetModulus(poly)) break;
        mc.reset();
    }
    if (!mc) return -1;

    // Montgomery form of x: x * Z_C mod f.
    std::vector<Elem> ra;
    {
        std::vector<Elem> c(n + 2, 0), quot;
        for (int j = 0; j <= m; ++j) {
            if ((m & j) == j) c[(size_t(1) << j) + 1] = 1;
        }
        c[1] ^= MontCoset<F>::Offset(fft, mc->block);
        FastDivMod(poly, c, quot, field, 32);
        c.resize(n, 0);
        ra.swap(c);
        fft.ToNovel(ra, m);
        fft.NovelFFT(ra, m, field, mc->block);
    }
    const std::vector<Elem> ra0 = ra;

    // Accumulators of trace0 = sum_i basis^(2^i) X_i and of the relative
    // traces T_l = sum_(i = 0 mod s) (basis * xi^l)^(2^i) X_i.
    // Only the traces for the depths a degree-d tree reaches with nodes of
    // degree >= 3 (about log2(d) + 3; deeper nodes, if any, compute their
    // own traces).
    const int n_traces = std::min<int>(CARRY_TRACES, (m + 3 + MONT_TOWER_BITS - 1) / MONT_TOWER_BITS);
    std::vector<Wide> wide_all(n, field.WideZero()), wide_t(n_traces * n, field.WideZero());
    Elem bcoef = basis, tcoef[CARRY_TRACES];
    tcoef[0] = basis;
    for (int l = 1; l < CARRY_TRACES; ++l) tcoef[l] = field.Mul2(tcoef[l - 1]);
    for (int i = 0; i < bits; ++i) {
        Wide* w = wide_all.data();
        for (size_t j = 0; j < n; ++j) w[j] = field.WideXor(w[j], field.MulWide(bcoef, ra[j]));
        if (i % MONT_TOWER_BITS == 0) {
            for (int l = 0; l < n_traces; ++l) {
                w = wide_t.data() + l * n;
                for (size_t j = 0; j < n; ++j) w[j] = field.WideXor(w[j], field.MulWide(tcoef[l], ra[j]));
            }
        }
        bcoef = field.Sqr(bcoef);
        for (int l = 0; l < CARRY_TRACES; ++l) tcoef[l] = field.Sqr(tcoef[l]);
        mc->Square(ra);
    }
    // X_Bits == x iff the Montgomery forms agree mod f, i.e. iff the
    // reduction of their difference vanishes (on the coset, hence
    // identically, as its degree is < d).
    for (size_t j = 0; j < n; ++j) ra[j] ^= ra0[j];
    mc->Reduce(ra, false);
    for (size_t j = 0; j < n; ++j) {
        if (ra[j] != 0) return 0;
    }

    carry.normal = false;
    carry.m = m;
    carry.block = mc->block;
    for (int l = 0; l < CARRY_TRACES; ++l) {
        carry.t[l].clear();
        if (l >= n_traces) continue;
        carry.t[l].resize(n);
        for (size_t j = 0; j < n; ++j) carry.t[l][j] = field.Reduce(wide_t[l * n + j]);
    }
    std::vector<Wide>().swap(wide_t);
    carry.ga = mc->Values();
    for (size_t j = 0; j < n; ++j) ra[j] = field.Reduce(wide_all[j]);
    mc->ToNormal(ra, d, trace0);
    return 1;
}

/** The tower scheme's traces from the eight bucket sums acc[r] (length
 * d each): traces[k] = sum_r (z^k)^(2^r) * acc[r], stripped. */
template<typename F>
void CombineTowerTraces(const std::vector<std::vector<typename F::Elem>>& acc, size_t depths, typename F::Elem z, const F& field, std::vector<std::vector<typename F::Elem>>& traces) {
    typedef typename F::Elem Elem;
    typedef typename F::Wide Wide;
    Elem zk = 1;
    for (size_t k = 0; k < depths; ++k) {
        Elem w[8];
        w[0] = zk;
        for (int r = 1; r < 8; ++r) w[r] = field.Sqr(w[r - 1]);
        auto& t = traces[k];
        const size_t d = t.size();
        if (F::LAZY_REDUCE) {
            for (size_t j = 0; j < d; ++j) {
                Wide s = field.WideZero();
                for (int r = 0; r < 8; ++r) s = field.WideXor(s, field.MulWide(w[r], acc[r][j]));
                t[j] = field.Reduce(s);
            }
        } else {
            for (int r = 0; r < 8; ++r) {
                typename F::Multiplier mul(field, w[r]);
                for (size_t j = 0; j < d; ++j) t[j] ^= mul(acc[r][j]);
            }
        }
        zk = field.Mul(zk, z);
    }
    for (auto& t : traces) PolyStrip(t);
}

/** Compute, for the first depths split levels of the root finding
 * recursion on poly (monic, degree >= 3), the split traces
 * traces[k] = Tr(beta_k * x) mod poly, the split parameters betas[k] for
 * all Bits depths, and test full factorizability.
 *
 * All traces come from one pass over the Frobenius orbit of x mod poly:
 * with X_i = x^(2^i) mod poly, Tr(beta*x) = sum_i beta^(2^i) * X_i, so each
 * trace is a linear combination accumulated while the orbit is generated by
 * Bits-1 modular squarings - the cost of a single trace computation. The
 * orbit is consumed as it is produced, keeping the working set at
 * O(depths) * deg(poly) field elements.
 *
 * Split parameters: the betas only need to be GF(2)-linearly independent
 * (so that the splits along any path eventually separate all roots). When
 * 8 divides Bits they are chosen from a tower basis,
 * beta_(8l+j) = basis * x^l * z^j with z generating GF(2^8): for the first
 * eight depths beta_k = basis * z^k with z^k in GF(2^8), whose Frobenius
 * powers repeat with period 8, so
 *     Tr(beta_k x) = sum_(r<8) (z^k)^(2^r) * Z_r,
 *     Z_r = sum_(i = r mod 8) basis^(2^i) * X_i,
 * and the accumulation costs (Bits + 8*depths) * deg(poly) products
 * instead of Bits * depths * deg(poly). Otherwise beta_k = basis * x^k.
 *
 * One more squaring gives X_Bits = x^(2^Bits) mod poly, which equals x iff
 * poly divides x^(2^Bits) - x, i.e. iff poly is square-free with all roots in
 * the field: the return value. This decides factorizability exactly at the
 * root, replacing the deferred trace^2 + trace test.
 *
 * From degree MONT_ORBIT_MIN_DEGREE on fields with the FFT tier and a cheap
 * Multiplier, the orbit runs in coset Montgomery form instead
 * (MontgomeryOrbit), with the tower over GF(2^MONT_TOWER_BITS) in place of
 * GF(2^8) (so betas differ). Given tree and carry, only traces[0] is
 * returned and carry receives the relative traces from which the
 * root-finding tree computes the deeper ones (MontCarry); without them all
 * depths traces are returned. */
template<typename F>
bool ComputeTraceOrbit(const std::vector<typename F::Elem>& poly, typename F::Elem basis, size_t depths, const F& field, std::vector<std::vector<typename F::Elem>>& traces, std::vector<typename F::Elem>& betas, MontTree<F>* tree = nullptr, MontCarry<F>* carry = nullptr) {
    typedef typename F::Elem Elem;
    typedef typename F::Wide Wide;
    const size_t d = poly.size() - 1;
    const int bits = field.Bits();
    CHECK_SAFE(d >= 3 && poly.back() == 1);
    if (carry) carry->tree = nullptr;
    traces.assign(depths, std::vector<Elem>(d, 0));
    betas.resize(bits);
    const bool tower = (bits % 8) == 0 && depths <= 8;
    // Number of accumulators and their per-step coefficients coef[a] =
    // (accumulator's beta)^(2^i): one per depth, or the basis^(2^i) of the
    // eight residue buckets of the tower scheme.
    const size_t nacc = tower ? 8 : depths;
    Elem z = 0;
    if (tower) {
        z = SubfieldGenerator8(field);
        Elem xl = basis;
        for (int l = 0; l < bits / 8; ++l) {
            Elem v = xl;
            for (int j = 0; j < 8; ++j) {
                betas[8 * l + j] = v;
                v = field.Mul(v, z);
            }
            xl = field.Mul2(xl);
        }
    } else {
        Elem beta = basis;
        for (int k = 0; k < bits; ++k) {
            betas[k] = beta;
            beta = field.Mul2(beta);
        }
    }
    std::vector<std::vector<Elem>> acc;
    if (tower && bits % MONT_TOWER_BITS == 0 && d >= MONT_ORBIT_MIN_DEGREE && TraceModFFTCapable(field) && TraceModCheapScalarMul<F>()) {
        // The Montgomery tree's tower: beta_(sl+k) = basis * xi^l * z^k
        // with z generating GF(2^s), s = MONT_TOWER_BITS.
        std::vector<Elem> mbetas(bits);
        const Elem mz = SubfieldGenerator(field, MONT_TOWER_BITS);
        Elem xl = basis;
        for (int l = 0; l < bits / MONT_TOWER_BITS; ++l) {
            Elem v = xl;
            for (int j = 0; j < MONT_TOWER_BITS; ++j) {
                mbetas[MONT_TOWER_BITS * l + j] = v;
                v = field.Mul(v, mz);
            }
            xl = field.Mul2(xl);
        }
        MontTree<F> local_tree;
        MontCarry<F> local_carry;
        MontTree<F>& mtree = (tree && carry) ? *tree : local_tree;
        MontCarry<F>& mcarry = (tree && carry) ? *carry : local_carry;
        int res = MontgomeryOrbit(poly, basis, field, traces[0], mtree.fft, mcarry);
        if (res == 0) return false;
        if (res == 1) {
            betas.swap(mbetas);
            mtree.z = mz;
            mcarry.tree = &mtree;
            if (tree && carry) {
                // The root-finding tree computes deeper traces from the carry.
                traces.resize(1);
            } else {
                MontCoset<F> mc(mtree.fft, field, mcarry.m, mcarry.block);
                bool ok = mc.SetModulus(poly, &mcarry.ga);
                CHECK_SAFE(ok);
                (void)ok;
                for (size_t k = 1; k < depths; ++k) {
                    CHECK_SAFE(k < size_t(CARRY_DEPTHS));
                    Elem zk = 1;
                    for (size_t r = 0; r < k % MONT_TOWER_BITS; ++r) zk = field.Mul(zk, mz);
                    MontTowerTrace(mc, mcarry.t[k / MONT_TOWER_BITS], d, zk, traces[k], field, true);
                }
            }
            return true;
        }
    }
    std::vector<Elem> coef(depths);
    for (size_t k = 0; k < depths; ++k) coef[k] = betas[k];
    Elem bcoef = basis; // basis^(2^i), the tower scheme's single coefficient
    TraceMod<F> trace_mod(poly, field);
    std::vector<Elem> orbit{0, 1};
    std::vector<Wide> wide;
    if (F::LAZY_REDUCE) {
        wide.assign(nacc * d, field.WideZero());
    } else {
        acc.assign(nacc, std::vector<Elem>(d, 0));
    }
    for (int i = 0; i < bits; ++i) {
        // Accumulators receiving X_i, with their coefficients.
        size_t a_lo = tower ? size_t(i % 8) : 0, a_hi = tower ? a_lo + 1 : depths;
        for (size_t a = a_lo; a < a_hi; ++a) {
            Elem c = tower ? bcoef : coef[a];
            if (F::LAZY_REDUCE) {
                Wide* w = wide.data() + a * d;
                for (size_t j = 0; j < orbit.size(); ++j) w[j] = field.WideXor(w[j], field.MulWide(c, orbit[j]));
            } else {
                typename F::Multiplier mul(field, c);
                auto& t = acc[a];
                for (size_t j = 0; j < orbit.size(); ++j) t[j] ^= mul(orbit[j]);
            }
        }
        if (tower) {
            bcoef = field.Sqr(bcoef);
        } else {
            for (size_t k = 0; k < depths; ++k) coef[k] = field.Sqr(coef[k]);
        }
        trace_mod.SquareAndReduce(orbit);
    }
    if (F::LAZY_REDUCE) {
        acc.assign(nacc, std::vector<Elem>(d, 0));
        for (size_t a = 0; a < nacc; ++a) {
            for (size_t j = 0; j < d; ++j) acc[a][j] = field.Reduce(wide[a * d + j]);
        }
    }
    if (tower) {
        CombineTowerTraces(acc, depths, z, field, traces);
    } else {
        for (size_t k = 0; k < depths; ++k) traces[k].swap(acc[k]);
        for (auto& t : traces) PolyStrip(t);
    }
    return orbit.size() == 2 && orbit[0] == 0 && orbit[1] == 1;
}

/** Tr(param * x) mod poly for small monic poly (3 <= degree <= 8), on
 * fixed-size local arrays: the per-node TraceMod machinery (tables,
 * scratch buffers, vector bookkeeping per squaring) costs far more than
 * the ~d^2/2 products per squaring at these degrees. Rows x^e mod poly for
 * e in [d, 2d-2] are built once; each of the Bits-1 steps squares the
 * coefficients and folds the high ones in through the rows (lazily
 * reduced where the field supports it). */
template<typename F>
void SmallTrace(const std::vector<typename F::Elem>& poly, typename F::Elem param, std::vector<typename F::Elem>& out, const F& field) {
    typedef typename F::Elem Elem;
    typedef typename F::Wide Wide;
    const size_t d = poly.size() - 1;
    CHECK_SAFE(d >= 3 && d <= SMALL_TRACE_MAX && poly.back() == 1);
    Elem rows[SMALL_TRACE_MAX - 1][SMALL_TRACE_MAX];
    Elem rem[SMALL_TRACE_MAX];
    for (size_t j = 0; j < d; ++j) rem[j] = poly[j];
    // rows[e - d] = x^e mod poly.
    for (size_t e = d; e + 1 < 2 * d; ++e) {
        for (size_t j = 0; j < d; ++j) rows[e - d][j] = rem[j];
        Elem fold = rem[d - 1];
        for (size_t j = d - 1; j > 0; --j) rem[j] = rem[j - 1] ^ field.Mul(fold, poly[j]);
        rem[0] = field.Mul(fold, poly[0]);
    }
    Elem x[SMALL_TRACE_MAX] = {0};
    x[1] = param;
    for (int it = 0; it < field.Bits() - 1; ++it) {
        Elem nx[SMALL_TRACE_MAX] = {0};
        if (F::LAZY_REDUCE) {
            Wide acc[SMALL_TRACE_MAX];
            for (size_t j = 0; j < d; ++j) acc[j] = field.WideZero();
            for (size_t i = 0; i < d; ++i) {
                Elem c = field.Sqr(x[i]);
                if (2 * i < d) {
                    nx[2 * i] ^= c;
                } else {
                    const Elem* row = rows[2 * i - d];
                    for (size_t j = 0; j < d; ++j) acc[j] = field.WideXor(acc[j], field.MulWide(c, row[j]));
                }
            }
            for (size_t j = 0; j < d; ++j) nx[j] ^= field.Reduce(acc[j]);
        } else {
            for (size_t i = 0; i < d; ++i) {
                Elem c = field.Sqr(x[i]);
                if (2 * i < d) {
                    nx[2 * i] ^= c;
                } else if (c != 0) {
                    typename F::Multiplier mul(field, c);
                    const Elem* row = rows[2 * i - d];
                    for (size_t j = 0; j < d; ++j) nx[j] ^= mul(row[j]);
                }
            }
        }
        nx[1] ^= param;
        for (size_t j = 0; j < d; ++j) x[j] = nx[j];
    }
    out.assign(x, x + d);
    PolyStrip(out);
}

template<typename F>
std::vector<typename F::Elem> BerlekampMasseyOdd(const std::vector<typename F::Elem>& osyn, size_t n_syn, size_t max_degree, const F& field);
template<typename F>
void LazyMiddleProduct(const typename F::Elem* c, const typename F::Elem* w, size_t n, typename F::Wide* out, typename F::Elem* es, typename F::Wide* wsc, const F& field);

/** q = g / h on mc's coset (gv: g's values there, deg_g = deg g), for monic
 * h: hv receives h's values, qv q's. Returns false unless h divides g: the
 * interpolated quotient must be monic of degree deg_g - deg h, which with
 * deg_g <= 2^m forces q * h = g. */
template<typename F>
bool CosetDivide(const AdditiveFFT<F>& fft, MontCoset<F>& mc, const std::vector<typename F::Elem>& gv, size_t deg_g, const std::vector<typename F::Elem>& h, std::vector<typename F::Elem>& hv, std::vector<typename F::Elem>& q, std::vector<typename F::Elem>& qv, const F& field) {
    typedef typename F::Elem Elem;
    const size_t deg_h = h.size() - 1;
    if (deg_h > deg_g) return false;
    MontCoset<F>::Eval(fft, field, mc.m, mc.block, h, hv);
    std::vector<Elem> inv;
    if (!BatchInvert(hv, inv, field)) return false;
    qv.resize(mc.n);
    for (size_t j = 0; j < mc.n; ++j) qv[j] = field.Mul(gv[j], inv[j]);
    std::vector<Elem> v = qv;
    mc.ToCoefficients(v, q);
    return q.size() == deg_g - deg_h + 1 && q.back() == 1;
}

/** Split a Montgomery-form node g (monic, degree n, on mc's coset with
 * values gv) by its split trace t (normal form: t = 0 or 1 at every root)
 * using the power sums q[k-1] = Q_k = sum over roots of alpha^k (k <=
 * q.size()) instead of a gcd.
 *
 * The roots A with t = 1 have the power sums
 *     P_k = sum_(alpha in A) alpha^k = sum_alpha t(alpha) alpha^k
 *         = sum_(j<n) t_j Q_(j+k),
 * a Hankel product, and t(alpha)^2 = t(alpha) gives P_(2k) = P_k^2: the
 * syndrome structure, so BerlekampMasseyOdd recovers prod_(A) (1 - alpha y)
 * from the odd P_k alone, at the cost of an |A|-root decode - against a
 * Euclidean gcd running from degree n down to n - |A|. With 2M >= n power
 * sums one of A and its complement B (power sums Q_k - P_k) has at most n/2
 * roots and is recovered; the other side follows by division on the coset,
 * whose exactness check also rejects a BM result from too few terms.
 *
 * Returns 1 with g0/g1 (the factors with t = 0 and t = 1), their coset
 * values v0/v1 and power sums q0/q1 (2M terms each); 0 if all roots lie on
 * one side; -1 if neither side was recovered (the caller uses the gcd). */
template<typename F>
int PowerSumSplit(const AdditiveFFT<F>& fft, MontCoset<F>& mc, const std::vector<typename F::Elem>& g, const std::vector<typename F::Elem>& gv, const std::vector<typename F::Elem>& t, const std::vector<typename F::Elem>& q, std::vector<typename F::Elem>& g0, std::vector<typename F::Elem>& g1, std::vector<typename F::Elem>& v0, std::vector<typename F::Elem>& v1, std::vector<typename F::Elem>& q0, std::vector<typename F::Elem>& q1, const F& field) {
    typedef typename F::Elem Elem;
    typedef typename F::Wide Wide;
    const size_t n = g.size() - 1;
    if (q.size() + 1 < 2 * n) return -1;
    // 2M power sums of each side: enough for the smaller side (2M >= n),
    // and for the children's own splits (which need twice their degree)
    // unless the split is unbalanced by more than ~4 standard deviations.
    size_t sq = 1;
    while (sq * sq < n) ++sq;
    size_t two_m = std::min(q.size() + 1 - n, n + 4 * sq + 32) & ~size_t(1);
    if (two_m < n) return -1;
    const size_t m_odd = two_m / 2;
    // Odd P_k, k = 2i + 1: sum_j t_j Q_(j+2i+1). Split j by parity:
    //     P_(2i+1) = sum_a t_(2a) Q_(2(a+i)+1) + sum_a t_(2a+1) Q_(2(a+i)+2),
    // two middle products of half length against the odd and the even
    // power sums, by transposed Karatsuba in blocks of h outputs; outputs
    // past the last full block as direct dot products.
    size_t h = (n + 1) / 2;
    h += h & 1;
    const size_t nb = m_odd / h, wlen = (nb + 1) * h + h;
    std::vector<Elem> ce(h, 0), co(h, 0), we(wlen, 0), wo(wlen, 0), es(4 * h);
    for (size_t j = 0; j < t.size(); ++j) ((j & 1) ? co : ce)[j / 2] = t[j];
    for (size_t s = 0; s < wlen && 2 * s < q.size(); ++s) we[s] = q[2 * s];
    for (size_t s = 0; s < wlen && 2 * s + 1 < q.size(); ++s) wo[s] = q[2 * s + 1];
    std::vector<Wide> oute(h), outo(h), wsc(2 * h);
    std::vector<Elem> p_odd(m_odd);
    for (size_t b = 0; b < nb; ++b) {
        LazyMiddleProduct(ce.data(), we.data() + b * h, h, oute.data(), es.data(), wsc.data(), field);
        LazyMiddleProduct(co.data(), wo.data() + b * h, h, outo.data(), es.data(), wsc.data(), field);
        for (size_t s = 0; s < h; ++s) p_odd[b * h + s] = field.Reduce(field.WideXor(oute[s], outo[s]));
    }
    for (size_t i = nb * h; i < m_odd; ++i) {
        Wide acc = field.WideZero();
        for (size_t a = 0; a < h; ++a) acc = field.WideXor(acc, field.WideXor(field.MulWide(ce[a], we[a + i]), field.MulWide(co[a], wo[a + i])));
        p_odd[i] = field.Reduce(acc);
    }
    std::vector<Elem> hv;
    for (int side = 1; side >= 0; --side) {
        std::vector<Elem> osyn = p_odd;
        if (side == 0) {
            for (size_t i = 0; i < m_odd; ++i) osyn[i] ^= q[2 * i];
        }
        std::vector<Elem> conn = BerlekampMasseyOdd(osyn, two_m, m_odd, field);
        if (conn.empty()) continue;
        if (conn.size() == 1) return 0;     // this side has no roots
        if (conn.size() == n + 1) return 0; // all roots on this side
        std::vector<Elem>& gs = side ? g1 : g0;
        std::vector<Elem>& go = side ? g0 : g1;
        std::vector<Elem>& vs = side ? v1 : v0;
        std::vector<Elem>& vo = side ? v0 : v1;
        gs.assign(conn.rbegin(), conn.rend());
        if (!CosetDivide(fft, mc, gv, n, gs, vs, go, vo, field)) continue;
        // Power sums of both sides.
        q1.resize(two_m);
        for (size_t i = 0; i < m_odd; ++i) q1[2 * i] = p_odd[i];
        for (size_t k = 2; k <= two_m; k += 2) q1[k - 1] = field.Sqr(q1[k / 2 - 1]);
        q0.resize(two_m);
        for (size_t k = 0; k < two_m; ++k) q0[k] = q[k] ^ q1[k];
        return 1;
    }
    return -1;
}

/** One step of the root finding algorithm; finds roots of stack[pos] and adds them to roots. Stack elements >= pos are destroyed.
 *
 * It operates on a stack of polynomials. The polynomial operated on is `stack[pos]`, where elements of `stack` with index higher
 * than `pos` are used as scratch space.
 *
 * `stack[pos]` is assumed to be square-free polynomial. If `fully_factorizable` is true, it is also assumed to have no irreducible
 * factors of degree higher than 1.

 * This implements the Berlekamp trace algorithm, plus an efficient test to fail fast in
 * case the polynomial cannot be fully factored.
 */
template<typename F>
bool RecFindRoots(std::vector<std::vector<typename F::Elem>>& stack, size_t pos, std::vector<typename F::Elem>& roots, bool fully_factorizable, int depth, typename F::Elem randv, const F& field, const std::vector<std::vector<typename F::Elem>>* orbit = nullptr, size_t orbit_degree = 0, const std::vector<typename F::Elem>* betas = nullptr, const MontCarry<F>* carry = nullptr, const std::vector<typename F::Elem>* psum = nullptr) {
    auto& ppoly = stack[pos];
    // We assert ppoly.size() > 1 (instead of just ppoly.size() > 0) to additionally exclude
    // constants polynomials because
    //  - ppoly is not constant initially (this is ensured by FindRoots()), and
    //  - we never recurse on a constant polynomial.
    CHECK_SAFE(ppoly.size() > 1 && ppoly.back() == 1);
    /* 1st degree input: constant term is the root. */
    if (ppoly.size() == 2) {
        roots.push_back(ppoly[0]);
        return true;
    }
    /* 2nd degree input: use direct quadratic solver. */
    if (ppoly.size() == 3) {
        CHECK_RETURN(ppoly[1] != 0, false); // Equations of the form (x^2 + a) have two identical solutions; contradicts square-free assumption. */
        auto input = field.Mul(ppoly[0], field.Sqr(field.Inv(ppoly[1])));
        auto root = field.Qrt(input);
        if ((field.Sqr(root) ^ root) != input) {
            CHECK_SAFE(!fully_factorizable);
            return false; // No root found.
        }
        auto sol = field.Mul(root, ppoly[1]);
        roots.push_back(sol);
        roots.push_back(sol ^ ppoly[1]);
        return true;
    }
    /* 3rd degree input and more: recurse further. */
    if (pos + 3 > stack.size()) {
        // Allocate memory if necessary.
        stack.resize((pos + 3) * 2);
    }
    auto& poly = stack[pos];
    auto& tmp = stack[pos + 1];
    auto& trace = stack[pos + 2];
    trace.clear();
    tmp.clear();
    // The node modulus on the carry's coset, if the carry is in Montgomery
    // form and the node needs it.
    std::unique_ptr<MontCoset<F>> node_mc;
    auto get_node_mc = [&]() -> MontCoset<F>& {
        if (!node_mc) {
            node_mc.reset(new MontCoset<F>(carry->tree->fft, field, carry->m, carry->block));
            bool ok = node_mc->SetModulus(poly, carry->ga.empty() ? nullptr : &carry->ga);
            CHECK_SAFE(ok);
            (void)ok;
        }
        return *node_mc;
    };
    const size_t node_degree = poly.size() - 1;
    // Split by power sums instead of a gcd (PowerSumSplit) when the node's
    // power sums are known and it is in Montgomery form; its traces must
    // then be normalized.
    const bool use_psum = psum && carry && carry->tree && !carry->normal && node_degree >= PSUM_SPLIT_MIN && psum->size() + 1 >= 2 * node_degree;
    bool psum_split = false;
    std::vector<typename F::Elem> trace_v, tmp_v, psum0, psum1, ps_g0, ps_g1;
    // Limit TraceMod lifetime to free precomputed table after trace loop.
    {
        // Built only when a trace is not available from the orbit.
        std::unique_ptr<TraceMod<F>> trace_mod;
        for (int iter = 0;; ++iter) {
            if (betas) {
                // Split parameters chosen by ComputeTraceOrbit, one per depth.
                CHECK_RETURN(size_t(depth) < betas->size(), false);
                randv = (*betas)[depth];
            }
            // Compute the polynomial (trace(x*randv) mod poly(x)) symbolically,
            // and put the result in `trace`.
            MINISKETCH_COZ_BEGIN("findroots-trace");
            const bool at_root = orbit && size_t(depth) < orbit->size() && poly.size() - 1 == orbit_degree;
            if (carry && carry->tree && !at_root && depth < CARRY_DEPTHS && !carry->t[depth / MONT_TOWER_BITS].empty()) {
                // The relative trace T_l (l = depth / s) mod poly is carried
                // down the tree; the depth's split parameter is
                // basis * xi^l * z^(depth % s) (tower scheme).
                const int l = depth / MONT_TOWER_BITS;
                typename F::Elem zk = 1;
                for (int k = 0; k < depth % MONT_TOWER_BITS; ++k) zk = field.Mul(zk, carry->tree->z);
                if (carry->normal) {
                    TowerTraceSmall(poly, carry->t[l], zk, trace, field);
                } else {
                    MontTowerTrace(get_node_mc(), carry->t[l], poly.size() - 1, zk, trace, field, use_psum);
                }
            } else if (orbit && size_t(depth) < orbit->size() && ORBIT_NODE_FACTOR * orbit_degree <= (poly.size() - 1) * field.Bits()) {
                // randv is the orbit's split parameter for this depth, so the
                // trace mod the root polynomial reduces to the one mod poly
                // (a factor of it): one division instead of Bits-1 modular
                // squarings.
                trace = (*orbit)[depth];
                FastDivMod(poly, trace, tmp, field, FastDivModCutoff(field));
                PolyStrip(trace);
            } else {
                if (fully_factorizable && poly.size() - 1 <= SMALL_TRACE_MAX) {
                    SmallTrace(poly, randv, trace, field);
                } else {
                    if (!trace_mod) trace_mod.reset(new TraceMod<F>(poly, field));
                    trace_mod->trace(trace, randv);
                }
            }
            MINISKETCH_COZ_END("findroots-trace");

            if (iter >= 1 && !fully_factorizable) {
                // If the polynomial cannot be factorized completely (it has an
                // irreducible factor of degree higher than 1), we want to avoid
                // the case where this is only detected after trying all BITS
                // independent split attempts fail (see the assert below).
                //
                // Observe that if we call y = randv*x, it is true that:
                //
                //   trace = y + y^2 + y^4 + y^8 + ... y^(FIELDSIZE/2) mod poly
                //
                // Due to the Frobenius endomorphism, this means:
                //
                //   trace^2 = y^2 + y^4 + y^8 + ... + y^FIELDSIZE mod poly
                //
                // Or, adding them up:
                //
                //   trace + trace^2 = y + y^FIELDSIZE mod poly.
                //                   = randv*x + randv^FIELDSIZE*x^FIELDSIZE
                //                   = randv*x + randv*x^FIELDSIZE
                //                   = randv*(x + x^FIELDSIZE).
                //     (all mod poly)
                //
                // x + x^FIELDSIZE is the polynomial which has every field element
                // as root once. Whenever x + x^FIELDSIZE is multiple of poly,
                // this means it only has unique first degree factors. The same
                // holds for its constant multiple randv*(x + x^FIELDSIZE) =
                // trace + trace^2.
                //
                // We use this test to quickly verify whether the polynomial is
                // fully factorizable after already having computed a trace.
                // We don't invoke it immediately; only when splitting has failed
                // at least once, which avoids it for most polynomials that are
                // fully factorizable (or at least pushes the test down the
                // recursion to factors which are smaller and thus faster).
                CHECK_SAFE(trace_mod); // Orbit traces are only passed down fully factorizable trees.
                tmp = trace;
                trace_mod->SquareAndReduce(tmp);
                // (trace^2 + trace) mod poly == (trace^2 mod poly) + trace,
                // as deg(trace) < deg(poly).
                if (tmp.size() < trace.size()) tmp.resize(trace.size(), 0);
                for (size_t i = 0; i < trace.size(); ++i) {
                    tmp[i] ^= trace[i];
                }
                while (tmp.size() && tmp.back() == 0) tmp.pop_back();

                // Whenever the test fails, we can immediately abort the root
                // finding. Whenever it succeeds, we can remember and pass down
                // the information that it is in fact fully factorizable, avoiding
                // the need to run the test again.
                if (tmp.size() != 0) return false;
                fully_factorizable = true;
            }

            if (fully_factorizable) {
                // Every successful iteration of this algorithm splits the input
                // polynomial further into buckets, each corresponding to a subset
                // of 2^(BITS-depth) roots. If after depth splits the degree of
                // the polynomial is >= 2^(BITS-depth), something is wrong.
                CHECK_RETURN(field.Bits() - depth >= std::numeric_limits<decltype(poly.size())>::digits ||
                    (poly.size() - 2) >> (field.Bits() - depth) == 0, false);
            }

            depth++;
            // In every iteration we multiply randv by 2. As a result, the set
            // of randv values forms a GF(2)-linearly independent basis of splits.
            randv = field.Mul2(randv);
            if (use_psum) {
                MINISKETCH_COZ_BEGIN("findroots-psum");
                MontCoset<F>& mc = get_node_mc();
                int r = PowerSumSplit(carry->tree->fft, mc, poly, mc.Values(), trace, *psum, ps_g0, ps_g1, trace_v, tmp_v, psum0, psum1, field);
                MINISKETCH_COZ_END("findroots-psum");
                if (r == 1) {
                    psum_split = true;
                    break;
                }
                if (r == 0) continue;
            }
            tmp = poly;
            MINISKETCH_COZ_BEGIN("findroots-gcd");
            FastGCD(trace, tmp, field, FastGCDCutoff(field));
            MINISKETCH_COZ_END("findroots-gcd");
            if (trace.size() != poly.size() && trace.size() > 1) break;
        }
    }
    MINISKETCH_COZ_BEGIN("findroots-divmod");
    if (psum_split) {
        trace.swap(ps_g0);
        tmp.swap(ps_g1);
    } else if (carry && carry->tree && !carry->normal && node_degree >= COSET_DIV_MIN) {
        // Values of the two factors on the node's coset, for Montgomery-form
        // nodes: poly = trace * tmp exactly, so tmp is poly / trace pointwise
        // there (trace has no root on the coset, as a factor of the root
        // polynomial), one interpolation instead of a long division; the
        // values also seed both factors' carries.
        MakeMonic(trace, field);
        MontCoset<F>& mc = get_node_mc();
        bool ok = CosetDivide(carry->tree->fft, mc, mc.Values(), node_degree, trace, trace_v, tmp, tmp_v, field);
        CHECK_SAFE(ok);
        (void)ok;
    } else {
        MakeMonic(trace, field);
        FastDivMod(trace, poly, tmp, field, FastDivModCutoff(field));
    }
    MINISKETCH_COZ_END("findroots-divmod");
    // At this point, the stack looks like [... (poly) tmp trace], and we want to recursively
    // find roots of trace and tmp (= poly/trace). As we don't care about poly anymore, move
    // trace into its position first.
    std::swap(poly, trace);
    // Carry the relative traces still needed to both factors.
    MontCarry<F> carry_poly, carry_tmp;
    if (carry && carry->tree && depth < CARRY_DEPTHS) {
        const int l_min = depth / MONT_TOWER_BITS;
        std::vector<typename F::Elem> node_t[CARRY_TRACES];
        bool have_node_t = carry->normal;
        auto to_factor = [&](const std::vector<typename F::Elem>& g, const std::vector<typename F::Elem>& gv, MontCarry<F>& out) {
            if (g.size() - 1 < 3) return;
            if (!carry->normal && g.size() - 1 > SMALL_CARRY_MAX) {
                MontCarryToFactor(*carry, g, l_min, out, field, gv.empty() ? nullptr : &gv);
                return;
            }
            if (!have_node_t) {
                // The node's traces in normal form, once.
                for (int l = l_min; l < CARRY_TRACES; ++l) {
                    if (carry->t[l].empty()) continue;
                    std::vector<typename F::Elem> v = carry->t[l];
                    get_node_mc().ToNormal(v, node_degree, node_t[l]);
                }
                have_node_t = true;
            }
            out.tree = carry->tree;
            out.normal = true;
            for (int l = 0; l < CARRY_TRACES; ++l) {
                out.t[l].clear();
                if (l < l_min || carry->t[l].empty()) continue;
                out.t[l] = carry->normal ? carry->t[l] : node_t[l];
                PolyMod(g, out.t[l], field);
                // An all-zero residue is a valid trace input; keep a marker
                // coefficient so that "empty" keeps meaning "not carried".
                if (out.t[l].empty()) out.t[l].assign(1, 0);
            }
        };
        to_factor(poly, trace_v, carry_poly);
        to_factor(tmp, tmp_v, carry_tmp);
    }
    // Release this node's coset state before recursing: summed along the
    // recursion path it would otherwise dominate the peak memory.
    node_mc.reset();
    std::vector<typename F::Elem>().swap(trace_v);
    std::vector<typename F::Elem>().swap(tmp_v);
    // Now the stack is [... (trace) tmp ...]. First we factor tmp (at pos = pos+1), and then
    // we factor trace (at pos = pos).
    if (!RecFindRoots(stack, pos + 1, roots, fully_factorizable, depth, randv, field, orbit, orbit_degree, betas, &carry_tmp, psum_split ? &psum1 : nullptr)) return false;
    carry_tmp = MontCarry<F>();
    std::vector<typename F::Elem>().swap(psum1);
    // The stack position pos contains trace, the polynomial with all of poly's roots which (after
    // multiplication with randv) have trace 0. This is never the case for irreducible factors
    // (which always end up in tmp), so we can set fully_factorizable to true when recursing.
    bool ret = RecFindRoots(stack, pos, roots, true, depth, randv, field, orbit, orbit_degree, betas, &carry_poly, psum_split ? &psum0 : nullptr);
    // Because of the above, recursion can never fail here.
    CHECK_SAFE(ret);
    return ret;
}

/** Returns the roots of a fully factorizable polynomial
 *
 * This function assumes that the input polynomial is square-free
 * and not the zero polynomial (represented by an empty vector).
 *
 * In case the square-free polynomial is not fully factorizable, i.e., it
 * has fewer roots than its degree, the empty vector is returned.
 */
template<typename F>
std::vector<typename F::Elem> FindRoots(const std::vector<typename F::Elem>& poly, typename F::Elem basis, const F& field, const std::vector<typename F::Elem>* odd_power_sums = nullptr) {
    std::vector<typename F::Elem> roots;
    CHECK_RETURN(poly.size() != 0, {});
    CHECK_RETURN(basis != 0, {});
    if (poly.size() == 1) return roots; // No roots when the polynomial is a constant.
    roots.reserve(poly.size() - 1);
    std::vector<std::vector<typename F::Elem>> stack = {poly};

    if (poly.size() - 1 >= std::max<size_t>(ORBIT_MIN_DEGREE, 3)) {
        // Enough depths for the nodes whose degree passes RecFindRoots' orbit
        // test while nodes roughly halve per level; deeper (smaller) nodes
        // compute their own traces.
        size_t depths = ORBIT_EXTRA_DEPTHS;
        for (int b = field.Bits(); b > 1; b >>= 1) ++depths;
        depths = std::min<size_t>(depths, field.Bits());
        std::vector<std::vector<typename F::Elem>> traces;
        std::vector<typename F::Elem> betas;
        MINISKETCH_COZ_BEGIN("findroots-orbit");
        MontTree<F> tree;
        MontCarry<F> carry;
        bool split = ComputeTraceOrbit(poly, basis, depths, field, traces, betas, &tree, &carry);
        MINISKETCH_COZ_END("findroots-orbit");
        if (!split) return {};
        // The roots' power sums Q_k (k = 1..2 * odd_power_sums->size()) from
        // the odd ones (Q_2k = Q_k^2), if given.
        // Only built when the root is large enough to split by them, and
        // only as far as the root's split can use (PowerSumSplit).
        std::vector<typename F::Elem> psum;
        const size_t d = poly.size() - 1;
        const bool use_psum = odd_power_sums && carry.tree && d >= PSUM_SPLIT_MIN;
        if (use_psum) {
            size_t sq = 1;
            while (sq * sq < d) ++sq;
            psum.resize(std::min(2 * odd_power_sums->size(), 2 * d + 4 * sq + 64));
            for (size_t k = 1; k <= psum.size(); ++k) psum[k - 1] = (k & 1) ? (*odd_power_sums)[k / 2] : field.Sqr(psum[k / 2 - 1]);
        }
        if (!RecFindRoots(stack, 0, roots, true, 0, basis, field, &traces, d, &betas, &carry, use_psum ? &psum : nullptr)) return {};
        CHECK_RETURN(poly.size() - 1 == roots.size(), {});
        return roots;
    }

    // Invoke the recursive factorization algorithm.
    if (!RecFindRoots(stack, 0, roots, false, 0, basis, field)) {
        // Not fully factorizable.
        return {};
    }
    CHECK_RETURN(poly.size() - 1 == roots.size(), {});
    return roots;
}

/** Unreduced middle product for lazily-reducing fields:
 * out[s] = sum_(k<n) c[k] * w[s+k] for s < n (w has 2n-1 entries), with
 * transposed Karatsuba while n is even and above MP_BASE:
 *     low  = alpha + MP(c0, w0 + w1),  high = alpha + MP(c1, w1 + w2),
 *     alpha = MP(c0 + c1, w1),
 * (c0, c1 the halves of c; w0, w1, w2 the length-(n-1) windows of w at
 * offsets 0, n/2, n), i.e. three half-size middle products instead of
 * four; below that, direct dot products. es must hold 4n elements. */
template<typename F>
void LazyMiddleProduct(const typename F::Elem* c, const typename F::Elem* w, size_t n, typename F::Wide* out, typename F::Elem* es, typename F::Wide* wsc, const F& field) {
    typedef typename F::Wide Wide;
    if (n <= MP_BASE || (n & 1)) {
        for (size_t s = 0; s < n; ++s) {
            Wide acc = field.WideZero();
            for (size_t k = 0; k < n; ++k) acc = field.WideXor(acc, field.MulWide(c[k], w[s + k]));
            out[s] = acc;
        }
        return;
    }
    const size_t h = n / 2, wl = 2 * h - 1;
    typename F::Elem* cs = es;
    typename F::Elem* ws = es + h;
    typename F::Elem* rest = es + h + wl;
    Wide* gamma = wsc;
    Wide* wrest = wsc + h;
    for (size_t k = 0; k < h; ++k) cs[k] = c[k] ^ c[h + k];
    // alpha into out[0..h), beta = MP(c0, w0 + w1) into out[h..n), then
    // low = alpha + beta and high = alpha (+ gamma below).
    LazyMiddleProduct(cs, w + h, h, out, rest, wrest, field);
    for (size_t t = 0; t < wl; ++t) ws[t] = w[t] ^ w[h + t];
    LazyMiddleProduct(c, ws, h, out + h, rest, wrest, field);
    for (size_t s = 0; s < h; ++s) {
        Wide beta = out[h + s];
        out[h + s] = out[s];
        out[s] = field.WideXor(out[s], beta);
    }
    // gamma = MP(c1, w1 + w2), added to the high half.
    for (size_t t = 0; t < wl; ++t) ws[t] = w[h + t] ^ w[2 * h + t];
    LazyMiddleProduct(c + h, ws, h, gamma, rest, wrest, field);
    for (size_t s = 0; s < h; ++s) out[h + s] = field.WideXor(out[h + s], gamma[s]);
}

/** Tail check of the early-terminating Berlekamp-Massey, on the odd
 * syndromes only: whether u[j] + sum_(a=1..L) conn[a]^2 * u[j-a] == 0 for all
 * j in [j0, j_end), where u[j] = S_(2j+1) and conn has degree L <= j0.
 *
 * The odd power sums u[j] = sum zeta * (zeta^2)^j are an LFSR sequence with
 * connection polynomial conn^2 (coefficients squared: roots squared), and
 * past the converged prefix this is equivalent to conn generating every
 * remaining syndrome (see BerlekampMassey). Each block of positions is one
 * product conn^2 * window: with the additive FFT (conn^2 transformed once,
 * one forward and one inverse transform per block), a direct
 * Multiplier-table check (table-based fields), or - on lazily-reducing
 * fields - a direct check at one carry-less product per term (small L) or
 * lazy Karatsuba middle products (LazyMiddleProduct, moderate L). Stops at
 * the first violated block. */
template<typename F>
bool BerlekampMasseyVerifyTail(const std::vector<typename F::Elem>& conn, const std::vector<typename F::Elem>& u, size_t j0, size_t j_end, const F& field) {
    typedef typename F::Elem Elem;
    const size_t L = conn.size() - 1;
    CHECK_SAFE(j0 >= L && j_end <= u.size());
    if (j0 >= j_end) return true;
    std::vector<Elem> conn2(L + 1);
    for (size_t a = 0; a <= L; ++a) conn2[a] = field.Sqr(conn[a]);
    if (F::LAZY_REDUCE && L >= BM_VERIFY_KARA_MIN_DEGREE && L < BM_VERIFY_FFT_MIN_DEGREE) {
        // Blocks of n positions as lazy Karatsuba middle products of the
        // reversed, zero-padded conn^2 with windows of u, where
        // n = b * 2^k >= L+1 with b <= MP_BASE (so the recursion halves
        // exactly down to direct kernels).
        typedef typename F::Wide Wide;
        size_t b = L + 1, k = 0;
        while (b > MP_BASE) {
            b = (b + 1) / 2;
            ++k;
        }
        const size_t n = b << k;
        std::vector<Elem> crev(n, 0), win(2 * n - 1), es(4 * n);
        for (size_t a = 0; a <= L; ++a) crev[n - 1 - a] = conn2[a];
        std::vector<Wide> out(n), wsc(2 * n);
        for (size_t p = j0; p < j_end; p += n) {
            // The window is u[p - n + 1 .. p + n - 1]; interior blocks read
            // it in place, edge blocks from a zero-padded copy.
            const Elem* wp;
            if (p >= n - 1 && p + n - 1 < j_end) {
                wp = u.data() + (p - (n - 1));
            } else {
                for (size_t t = 0; t < 2 * n - 1; ++t) {
                    size_t idx = p + t;
                    win[t] = (idx >= n - 1 && idx - (n - 1) < j_end) ? u[idx - (n - 1)] : Elem(0);
                }
                wp = win.data();
            }
            LazyMiddleProduct(crev.data(), wp, n, out.data(), es.data(), wsc.data(), field);
            const size_t cnt = std::min(n, j_end - p);
            for (size_t s2 = 0; s2 < cnt; ++s2) {
                if (field.Reduce(out[s2]) != 0) return false;
            }
        }
        return true;
    }
    if (F::LAZY_REDUCE && L < BM_VERIFY_FFT_MIN_DEGREE) {
        // Four positions at a time, so every coefficient load feeds four
        // accumulators.
        typedef typename F::Wide Wide;
        size_t j = j0;
        for (; j + 4 <= j_end; j += 4) {
            Wide a0 = field.WideFrom(u[j]), a1 = field.WideFrom(u[j + 1]);
            Wide a2 = field.WideFrom(u[j + 2]), a3 = field.WideFrom(u[j + 3]);
            for (size_t a = 1; a <= L; ++a) {
                Elem c = conn2[a];
                a0 = field.WideXor(a0, field.MulWide(c, u[j - a]));
                a1 = field.WideXor(a1, field.MulWide(c, u[j + 1 - a]));
                a2 = field.WideXor(a2, field.MulWide(c, u[j + 2 - a]));
                a3 = field.WideXor(a3, field.MulWide(c, u[j + 3 - a]));
            }
            if (field.Reduce(a0) | field.Reduce(a1) | field.Reduce(a2) | field.Reduce(a3)) return false;
        }
        for (; j < j_end; ++j) {
            Wide acc = field.WideFrom(u[j]);
            for (size_t a = 1; a <= L; ++a) acc = field.WideXor(acc, field.MulWide(conn2[a], u[j - a]));
            if (field.Reduce(acc)) return false;
        }
        return true;
    }
    const bool use_fft = TraceModFFTCapable(field) && L >= (TraceModCheapScalarMul<F>() ? BM_VERIFY_FFT_MIN_DEGREE : BM_VERIFY_FFT_MIN_DEGREE_TABLE);
    if (!F::LAZY_REDUCE && !use_fft) {
        // Table-based fields: one Multiplier per coefficient for the whole
        // tail (instead of per block), applied coefficient by coefficient
        // to L1-sized blocks of positions so each table stays hot.
        std::vector<typename F::Multiplier> mul;
        mul.reserve(L);
        for (size_t a = 1; a <= L; ++a) mul.emplace_back(field, conn2[a]);
        std::vector<Elem> acc;
        const size_t B = 512;
        for (size_t p = j0; p < j_end; p += B) {
            const size_t b = std::min(B, j_end - p);
            acc.assign(u.begin() + p, u.begin() + p + b);
            for (size_t a = 1; a <= L; ++a) {
                if (conn2[a] == 0) continue;
                const Elem* src = u.data() + p - a;
                for (size_t t = 0; t < b; ++t) acc[t] ^= mul[a - 1](src[t]);
            }
            for (size_t t = 0; t < b; ++t) {
                if (acc[t] != 0) return false;
            }
        }
        return true;
    }
    // A block covers j in [p, p+B); its window u[p-L .. p+B) times conn^2
    // puts position j at product index j - p + L.
    int fft_m = 0;
    size_t block;
    if (use_fft) {
        // The product (length B + 2L) must fit the transform without
        // wraparound; 2^fft_m >= 4(L+1) keeps at least half of it useful.
        fft_m = 1;
        while ((size_t(1) << fft_m) < 4 * (L + 1)) ++fft_m;
        if (fft_m > field.Bits()) fft_m = 0;
    }
    block = fft_m ? (size_t(1) << fft_m) - 2 * L : std::max<size_t>(2 * L, 64);

    TraceModScratch<F> scratch;
    std::vector<Elem> win, prod, conn2_fft;
    if (fft_m) scratch.fft.Prepare(conn2, fft_m, field, conn2_fft);
    const size_t n_points = size_t(1) << fft_m;
    for (size_t p = j0; p < j_end; p += block) {
        size_t b = std::min(block, j_end - p);
        win.assign(u.begin() + (p - L), u.begin() + (p + b));
        if (fft_m) {
            prod.assign(n_points, 0);
            std::copy(win.begin(), win.end(), prod.begin());
            scratch.fft.FFT(prod, fft_m, field);
            for (size_t i = 0; i < n_points; ++i) prod[i] = field.Mul(prod[i], conn2_fft[i]);
            scratch.fft.IFFT(prod, fft_m, field);
        } else {
            TraceModPolyMulFull(conn2, win, prod, field, scratch, 0);
        }
        for (size_t i = L; i < L + b; ++i) {
            if (i < prod.size() && prod[i] != 0) return false;
        }
    }
    return true;
}

/** Find the minimal LFSR (connection polynomial, constant term 1) generating
 * the n_syn syndromes S_1..S_(n_syn) of a sketch, given only its odd
 * syndromes osyn[j] = S_(2j+1) (the even ones are S_(2k) = S_k^2), or {} if
 * the degree would exceed max_degree.
 *
 * The Frobenius structure S_(2k) = S_k^2 - S(x)^2 = S_even(x) with
 * S(x) = sum S_k x^k - makes this cheaper than the textbook iteration while
 * returning exactly the same polynomial:
 *
 *  - The discrepancy at every even-indexed syndrome S_(2k) is identically
 *    zero, so those steps are skipped. This is the classical binary-BCH
 *    simplification of Berlekamp's algorithm (Berlekamp, "Algebraic Coding
 *    Theory", 1968), and holds for every Frobenius-structured sequence, not
 *    just genuine power sums (the unit tests compare against the
 *    unsimplified iteration on random inputs).
 *
 *  - Early termination: once the LFSR length L has stopped growing (after
 *    a zero discrepancy with at least 2L+2 syndromes processed), the
 *    current connection polynomial conn is checked against the whole tail
 *    at once (BerlekampMasseyVerifyTail). With C(x) = conn*S, the structure
 *    gives conn*C = C^2 + conn^2 * S_odd, so at odd positions k
 *        [x^k] (conn^2 * S_odd) = sum_i conn_i C_(k-i),
 *    and at even positions
 *        C_(2j) = C_j^2 + [x^2j](conn * (C mod x^(L+1))).
 *    By induction over the positions past the prefix (all beyond 2L): C
 *    vanishes on the tail iff the odd syndromes satisfy the recurrence with
 *    conn^2 there. So the full iteration would see zero discrepancies
 *    throughout the tail and return the same polynomial (including the
 *    max_degree failure, as a failed check just resumes the iteration).
 *    This replaces the serial O(c*L) tail by O(c log L) (FFT tier) work, or
 *    for moderate L on lazily-reducing fields by an independent-position
 *    check at one carry-less product per term. The even syndromes are only
 *    computed for the prefix actually walked. */
template<typename F>
std::vector<typename F::Elem> BerlekampMasseyOdd(const std::vector<typename F::Elem>& osyn, size_t n_syn, size_t max_degree, const F& field) {
    typedef typename F::Elem Elem;
    CHECK_SAFE(n_syn <= 2 * osyn.size());
    // Generic fields: one Multiplier table per processed syndrome. A deque
    // grows without reallocation copies (the prefix length is not known in
    // advance, and a vector's doubling would transiently hold both buffers).
    std::deque<typename F::Multiplier> table;
    std::vector<Elem> syn, current, prev, tmp;
    current.reserve(n_syn / 2 + 1);
    prev.reserve(n_syn / 2 + 1);
    tmp.reserve(n_syn / 2 + 1);
    current.resize(1);
    current[0] = 1;
    prev.resize(1);
    prev[0] = 1;
    Elem b = 1, b_inv = 1;
    bool b_have_inv = true;
    // Earliest number of processed syndromes at which a tail verification
    // is attempted; pushed back geometrically after a failed attempt.
    size_t next_verify = BM_VERIFY_MIN_PREFIX;

    for (size_t n = 0; n != n_syn; ++n) {
        // syn[n] = S_(n+1): odd ones from the sketch, S_(2k) = S_k^2.
        syn.push_back((n & 1) ? field.Sqr(syn[n / 2]) : osyn[n / 2]);
        // Lazily-reducing fields need no per-syndrome Multiplier.
        if (!F::LAZY_REDUCE) table.emplace_back(field, syn[n]);
        // Even-indexed S have zero discrepancy.
        if (n & 1) continue;
        auto discrepancy = syn[n];
        if (F::LAZY_REDUCE) {
            typename F::Wide acc = field.WideZero();
            for (size_t i = 1; i < current.size(); ++i) acc = field.WideXor(acc, field.MulWide(syn[n - i], current[i]));
            discrepancy ^= field.Reduce(acc);
        } else {
            for (size_t i = 1; i < current.size(); ++i) discrepancy ^= table[n - i](current[i]);
        }
        if (discrepancy != 0) {
            int x = static_cast<int>(n + 1 - (current.size() - 1) - (prev.size() - 1));
            if (!b_have_inv) {
                b_inv = field.Inv(b);
                b_have_inv = true;
            }
            bool swap = 2 * (current.size() - 1) <= n;
            if (swap) {
                if (prev.size() + x - 1 > max_degree) return {}; // We'd exceed maximum degree
                tmp = current;
                current.resize(prev.size() + x);
            }
            typename F::Multiplier mul(field, field.Mul(discrepancy, b_inv));
            for (size_t i = 0; i < prev.size(); ++i) current[i + x] ^= mul(prev[i]);
            if (swap) {
                std::swap(prev, tmp);
                b = discrepancy;
                b_have_inv = false;
            }
        } else {
            // A zero discrepancy with 2L + BM_VERIFY_GAP <= n + 1 processed
            // syndromes: the LFSR has likely converged. The next unprocessed
            // odd syndrome is S_(n+3) = osyn[n/2 + 1].
            size_t done = n + 1;
            size_t L = current.size() - 1;
            if (done >= next_verify && done >= 2 * L + BM_VERIFY_GAP && n_syn - done >= BM_VERIFY_MIN_TAIL) {
                if (BerlekampMasseyVerifyTail(current, osyn, n / 2 + 1, (n_syn + 1) / 2, field)) break;
                next_verify = done + done / 2;
            }
        }
    }
    CHECK_RETURN(current.size() && current.back() != 0, {});
    return current;
}

/** BerlekampMasseyOdd on a full Frobenius-structured syndrome sequence
 *  (syndromes[k-1] = S_k, as produced by ReconstructAllSyndromes). */
template<typename F>
std::vector<typename F::Elem> BerlekampMassey(const std::vector<typename F::Elem>& syndromes, size_t max_degree, const F& field) {
    std::vector<typename F::Elem> osyn;
    osyn.reserve(syndromes.size() / 2 + 1);
    for (size_t k = 0; k < syndromes.size(); k += 2) osyn.push_back(syndromes[k]);
    return BerlekampMasseyOdd(osyn, syndromes.size(), max_degree, field);
}

template<typename F>
std::vector<typename F::Elem> ReconstructAllSyndromes(const std::vector<typename F::Elem>& odd_syndromes, const F& field) {
    std::vector<typename F::Elem> all_syndromes;
    all_syndromes.resize(odd_syndromes.size() * 2);
    for (size_t i = 0; i < odd_syndromes.size(); ++i) {
        all_syndromes[i * 2] = odd_syndromes[i];
        all_syndromes[i * 2 + 1] = field.Sqr(all_syndromes[i]);
    }
    return all_syndromes;
}

template<typename F>
void AddToOddSyndromes(std::vector<typename F::Elem>& osyndromes, typename F::Elem data, const F& field) {
    auto sqr = field.Sqr(data);
    typename F::Multiplier mul(field, sqr);
    for (auto& osyndrome : osyndromes) {
        osyndrome ^= data;
        data = mul(data);
    }
}

// Batched variants of AddToOddSyndromes that walk the odd-syndrome vector once
// for K elements at a time (K in {2, 4, 8}). Each builds K Multipliers up front
// and then issues K independent self-multiplications per syndrome iteration.
// The K mul chains are independent across the K lanes, so the CPU pipelines
// them; the per-syndrome read+write is amortized Kx compared to K serial calls.
//
// We keep three explicit K-instantiations rather than a single C++11 template
// over K because Multiplier has no default constructor — a generic version
// requires placement-new on aligned_storage, which adds boilerplate without
// adding flexibility (only K = 1, 2, 4, 8 are wired into the dispatch).
template<typename F>
void AddToOddSyndromesBatch2(std::vector<typename F::Elem>& osyndromes,
                              typename F::Elem e0, typename F::Elem e1,
                              const F& field) {
    typename F::Multiplier mul0(field, field.Sqr(e0));
    typename F::Multiplier mul1(field, field.Sqr(e1));
    auto d0 = e0;
    auto d1 = e1;
    for (auto& osyndrome : osyndromes) {
        osyndrome ^= d0 ^ d1;
        d0 = mul0(d0);
        d1 = mul1(d1);
    }
}

template<typename F>
void AddToOddSyndromesBatch4(std::vector<typename F::Elem>& osyndromes,
                              typename F::Elem e0, typename F::Elem e1,
                              typename F::Elem e2, typename F::Elem e3,
                              const F& field) {
    typename F::Multiplier mul0(field, field.Sqr(e0));
    typename F::Multiplier mul1(field, field.Sqr(e1));
    typename F::Multiplier mul2(field, field.Sqr(e2));
    typename F::Multiplier mul3(field, field.Sqr(e3));
    auto d0 = e0;
    auto d1 = e1;
    auto d2 = e2;
    auto d3 = e3;
    for (auto& osyndrome : osyndromes) {
        osyndrome ^= d0 ^ d1 ^ d2 ^ d3;
        d0 = mul0(d0);
        d1 = mul1(d1);
        d2 = mul2(d2);
        d3 = mul3(d3);
    }
}

template<typename F>
void AddToOddSyndromesBatch8(std::vector<typename F::Elem>& osyndromes,
                              typename F::Elem e0, typename F::Elem e1,
                              typename F::Elem e2, typename F::Elem e3,
                              typename F::Elem e4, typename F::Elem e5,
                              typename F::Elem e6, typename F::Elem e7,
                              const F& field) {
    typename F::Multiplier mul0(field, field.Sqr(e0));
    typename F::Multiplier mul1(field, field.Sqr(e1));
    typename F::Multiplier mul2(field, field.Sqr(e2));
    typename F::Multiplier mul3(field, field.Sqr(e3));
    typename F::Multiplier mul4(field, field.Sqr(e4));
    typename F::Multiplier mul5(field, field.Sqr(e5));
    typename F::Multiplier mul6(field, field.Sqr(e6));
    typename F::Multiplier mul7(field, field.Sqr(e7));
    auto d0 = e0;
    auto d1 = e1;
    auto d2 = e2;
    auto d3 = e3;
    auto d4 = e4;
    auto d5 = e5;
    auto d6 = e6;
    auto d7 = e7;
    for (auto& osyndrome : osyndromes) {
        // Tree-reduce the eight XORs so the reduction takes 3 dependent
        // cycles instead of 7, leaving more room for the 8 muls to retire.
        auto x01 = d0 ^ d1;
        auto x23 = d2 ^ d3;
        auto x45 = d4 ^ d5;
        auto x67 = d6 ^ d7;
        osyndrome ^= (x01 ^ x23) ^ (x45 ^ x67);
        d0 = mul0(d0);
        d1 = mul1(d1);
        d2 = mul2(d2);
        d3 = mul3(d3);
        d4 = mul4(d4);
        d5 = mul5(d5);
        d6 = mul6(d6);
        d7 = mul7(d7);
    }
}

template<typename F>
std::vector<typename F::Elem> FullDecode(const std::vector<typename F::Elem>& osyndromes, const F& field) {
    auto asyndromes = ReconstructAllSyndromes<typename F::Elem>(osyndromes, field);
    auto poly = BerlekampMassey(asyndromes, field);
    std::reverse(poly.begin(), poly.end());
    return FindRoots(poly, field);
}

template<typename F>
class SketchImpl final : public Sketch
{
    const F m_field;
    std::vector<typename F::Elem> m_syndromes;
    typename F::Elem m_basis;
    // Batch size used by AddBatch(). One of {1, 2, 4, 8}; default 4 matches
    // the original single-helper implementation. Settable at runtime via
    // minisketch_set_batch_size() so callers (and the bench tool) can sweep.
    uint32_t m_batch_size = 4;

public:
    template<typename... Args>
    SketchImpl(int implementation, int bits, const Args&... args) : Sketch(implementation, bits), m_field(args...) {
#ifdef MINISKETCH_FUZZ_DETERMINISTIC
        // Fuzz builds only: never draw entropy, so that fuzz input -> behavior
        // is a pure function. Fuzz targets override the basis per input with
        // SetSeed(); this fixed value only covers construction itself.
        m_basis = m_field.FromSeed(0x6d696e69736b6574 /* "minisket" */);
#else
        std::random_device rng;
        std::uniform_int_distribution<uint64_t> dist;
        m_basis = m_field.FromSeed(dist(rng));
#endif
    }

    size_t Syndromes() const override { return m_syndromes.size(); }
    void Init(size_t count) override { m_syndromes.assign(count, 0); }

    void Add(uint64_t val) override
    {
        auto elem = m_field.FromUint64(val);
        AddToOddSyndromes(m_syndromes, elem, m_field);
    }

    int SetBatchSize(uint32_t batch_size) override
    {
        if (batch_size != 1 && batch_size != 2 && batch_size != 4 && batch_size != 8) return 0;
        m_batch_size = batch_size;
        return 1;
    }

    void AddBatch(const uint64_t* elements, size_t count) override
    {
        size_t i = 0;
        if (m_batch_size == 8) {
            while (i + 8 <= count) {
                auto e0 = m_field.FromUint64(elements[i]);
                auto e1 = m_field.FromUint64(elements[i + 1]);
                auto e2 = m_field.FromUint64(elements[i + 2]);
                auto e3 = m_field.FromUint64(elements[i + 3]);
                auto e4 = m_field.FromUint64(elements[i + 4]);
                auto e5 = m_field.FromUint64(elements[i + 5]);
                auto e6 = m_field.FromUint64(elements[i + 6]);
                auto e7 = m_field.FromUint64(elements[i + 7]);
                AddToOddSyndromesBatch8(m_syndromes, e0, e1, e2, e3, e4, e5, e6, e7, m_field);
                i += 8;
            }
        } else if (m_batch_size == 4) {
            while (i + 4 <= count) {
                auto e0 = m_field.FromUint64(elements[i]);
                auto e1 = m_field.FromUint64(elements[i + 1]);
                auto e2 = m_field.FromUint64(elements[i + 2]);
                auto e3 = m_field.FromUint64(elements[i + 3]);
                AddToOddSyndromesBatch4(m_syndromes, e0, e1, e2, e3, m_field);
                i += 4;
            }
        } else if (m_batch_size == 2) {
            while (i + 2 <= count) {
                auto e0 = m_field.FromUint64(elements[i]);
                auto e1 = m_field.FromUint64(elements[i + 1]);
                AddToOddSyndromesBatch2(m_syndromes, e0, e1, m_field);
                i += 2;
            }
        }
        for (; i < count; ++i) {
            auto elem = m_field.FromUint64(elements[i]);
            AddToOddSyndromes(m_syndromes, elem, m_field);
        }
    }

    void Serialize(unsigned char* ptr) const override
    {
        BitWriter writer(ptr);
        for (const auto& val : m_syndromes) {
            m_field.Serialize(writer, val);
        }
        writer.Flush();
    }

    void Deserialize(const unsigned char* ptr) override
    {
        BitReader reader(ptr);
        for (auto& val : m_syndromes) {
            val = m_field.Deserialize(reader);
        }
    }

    int Decode(int max_count, uint64_t* out) const override
    {
        MINISKETCH_COZ_BEGIN("decode-berlekamp-massey");
        // The sketch stores the odd syndromes; BM derives the even ones
        // (S_2k = S_k^2) only for the prefix it walks.
        auto poly = BerlekampMasseyOdd(m_syndromes, 2 * m_syndromes.size(), max_count, m_field);
        MINISKETCH_COZ_END("decode-berlekamp-massey");
        if (poly.size() == 0) return -1;
        if (poly.size() == 1) return 0;
        if ((int)poly.size() > 1 + max_count) return -1;
        std::reverse(poly.begin(), poly.end());
        MINISKETCH_COZ_BEGIN("decode-findroots");
        // The sketch's syndromes are the power sums of the roots poly is
        // expected to have; the root finder uses them for splitting.
        auto roots = FindRoots(poly, m_basis, m_field, &m_syndromes);
        MINISKETCH_COZ_END("decode-findroots");
        if (roots.size() == 0) return -1;

        for (const auto& root : roots) {
            *(out++) = m_field.ToUint64(root);
        }
        return static_cast<int>(roots.size());
    }

    size_t Merge(const Sketch* other_sketch) override
    {
        // Sad cast. This is safe only because the caller code in minisketch.cpp checks
        // that implementation and field size match.
        const SketchImpl* other = static_cast<const SketchImpl*>(other_sketch);
        m_syndromes.resize(std::min(m_syndromes.size(), other->m_syndromes.size()));
        for (size_t i = 0; i < m_syndromes.size(); ++i) {
            m_syndromes[i] ^= other->m_syndromes[i];
        }
        return m_syndromes.size();
    }

    void SetSeed(uint64_t seed) override
    {
        if (seed == (uint64_t)-1) {
            m_basis = 1;
        } else {
            m_basis = m_field.FromSeed(seed);
        }
    }
};

#endif
