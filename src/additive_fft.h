/**********************************************************************
 * Copyright (c) 2026 The minisketch developers                        *
 * Distributed under the MIT software license, see the accompanying   *
 * file LICENSE or http://www.opensource.org/licenses/mit-license.php.*
 **********************************************************************/

#ifndef _MINISKETCH_ADDITIVE_FFT_H_
#define _MINISKETCH_ADDITIVE_FFT_H_

#include <vector>

#include "util.h"

/** Additive (Cantor-basis) FFT engine over GF(2^Bits), usable when Bits is a
 * power of two (the caller gates on this; see TraceModFFTEligible).
 *
 * The transform evaluates a polynomial of degree < 2^t at the 2^t points of
 * the subspace spanned by a Cantor special basis: beta_0 = 1 and
 * beta_j^2 + beta_j = beta_{j-1}, built with the field's Qrt operation. The
 * chain cannot fail in fields of power-of-two degree: the image of x^2 + x
 * is the trace-zero hyperplane, and every chain element stays inside it.
 *
 * Algorithm: the "Cantor algorithm" as presented by Badakhshan, Samanta and
 * Gong, "Accelerating Post-quantum Secure zkSNARKs by Optimizing Additive
 * FFT" (SAC 2025; doc/sac2025-2-paper17_optimizing_additive_fft.pdf),
 * following Cantor (1989). Each of the t rounds divides the current blocks
 * by a vanishing polynomial Z_{W_i}(x) = S^i(x) with S(x) = x^2 + x. In the
 * Cantor basis these polynomials have all coefficients in GF(2), nonzero
 * exactly at degrees 2^j with binomial(i, j) odd (Lucas' theorem), so the
 * divisions are pure XOR; the only field multiplications are one
 * multiply-accumulate per high-half position, (1/2)*2^t*t per transform.
 * As the evaluation shift is 0 (the subspace itself), the twiddle factor of
 * block i is independent of the round: the XOR of beta_{j+1} over the set
 * bits j of i.
 *
 * Output order: after FFT, index j holds the evaluation at
 * eta_j = XOR of beta_b over the set bits b of j. IFFT is the exact inverse.
 * For multiplication (evaluate, multiply pointwise, interpolate) the point
 * order is irrelevant.
 *
 * The Frobenius truncation of Li et al., "Frobenius Additive Fast Fourier
 * Transform" (ISSAC 2018; doc/3208976.3208998_frobenius_additive_fft.pdf)
 * does not apply here: it requires polynomial coefficients in GF(2), while
 * the decoder's coefficients are full field elements. That paper's review
 * of the Cantor construction was used as a cross-check of this one.
 */
template<typename F>
class AdditiveFFT {
public:
    typedef typename F::Elem Elem;

    /** Grow the tables to dimension t (2^t evaluation points). Requires
     *  t <= field.Bits() with field.Bits() a power of two. Idempotent. */
    void Extend(const F& field, int t) {
        if (t <= max_t) return;
        CHECK_SAFE(t <= field.Bits() && (field.Bits() & (field.Bits() - 1)) == 0);
        if (beta.empty()) beta.push_back(1);
        while (beta.size() < (size_t)t) {
            Elem next = field.Qrt(beta.back());
            CHECK_SAFE((field.Sqr(next) ^ next) == beta.back());
            beta.push_back(next);
        }
        while (zeta.size() < (size_t)t) {
            // zeta[q]: distances from the leading term x^(2^q) of Z_{W_q} to
            // its lower terms x^(2^j), present iff binomial(q, j) is odd,
            // i.e. iff the bits of j are a subset of the bits of q.
            size_t q = zeta.size();
            zeta.emplace_back();
            for (size_t j = 0; j < q; ++j) {
                if ((q & j) == j) zeta.back().push_back((size_t(1) << q) - (size_t(1) << j));
            }
        }
        int s = 0;
        if (twiddle.empty()) twiddle.push_back(0);
        while ((size_t(1) << s) < twiddle.size()) ++s;
        while (twiddle.size() < (size_t(1) << (t - 1))) {
            // twiddle[i] = XOR of beta[j + 1] over the set bits j of i.
            size_t half = twiddle.size();
            for (size_t i = 0; i < half; ++i) twiddle.push_back(twiddle[i] ^ beta[s + 1]);
            ++s;
        }
        max_t = t;
    }

    int MaxDim() const { return max_t; }

    /** The Cantor basis built so far (element j spans dimension j); exposed
     *  for tests and diagnostics. */
    const std::vector<Elem>& Basis() const { return beta; }

    /** In-place transform of f (size 2^m, m <= MaxDim()): coefficients to
     *  evaluations at eta_0..eta_{2^m - 1}. */
    void FFT(std::vector<Elem>& f, int m, const F& field) const {
        CHECK_SAFE(m >= 1 && m <= max_t && f.size() == (size_t(1) << m));
        for (int r = 0; r < m; ++r) {
            int p = m - r;
            size_t h = size_t(1) << (p - 1);
            const std::vector<size_t>& zs = zeta[p - 1];
            for (size_t i = 0; i < (size_t(1) << r); ++i) {
                size_t off = i << p;
                Elem tw = twiddle[i];
                // Divide the block by Z_{W_{p-1}} (XOR cascades) while
                // applying the twiddle multiply-accumulate into the low half.
                // Decreasing k is essential: cascades from higher k can land
                // on lower positions of the high half, and must do so before
                // those positions are read.
                if (tw == 0) {
                    for (size_t k = off + 2 * h; k-- > off + h; ) {
                        Elem c = f[k];
                        if (c == 0) continue;
                        for (size_t z : zs) f[k - z] ^= c;
                    }
                } else if (h >= 4) {
                    typename F::Multiplier mul(field, tw);
                    for (size_t k = off + 2 * h; k-- > off + h; ) {
                        Elem c = f[k];
                        if (c == 0) continue;
                        for (size_t z : zs) f[k - z] ^= c;
                        f[k - h] ^= mul(c);
                    }
                } else {
                    for (size_t k = off + 2 * h; k-- > off + h; ) {
                        Elem c = f[k];
                        if (c == 0) continue;
                        for (size_t z : zs) f[k - z] ^= c;
                        f[k - h] ^= field.Mul(tw, c);
                    }
                }
                // Now the low half holds f0 (= remainder + twiddle * quotient)
                // and the high half holds the quotient q; f1 = f0 + q.
                for (size_t k = off; k < off + h; ++k) f[k + h] ^= f[k];
            }
        }
    }

    /** Exact inverse of FFT. */
    void IFFT(std::vector<Elem>& f, int m, const F& field) const {
        CHECK_SAFE(m >= 1 && m <= max_t && f.size() == (size_t(1) << m));
        for (int r = m - 1; r >= 0; --r) {
            int p = m - r;
            size_t h = size_t(1) << (p - 1);
            const std::vector<size_t>& zs = zeta[p - 1];
            for (size_t i = 0; i < (size_t(1) << r); ++i) {
                size_t off = i << p;
                Elem tw = twiddle[i];
                // Undo f1 = f0 + q: the high half becomes q again.
                for (size_t k = off; k < off + h; ++k) f[k + h] ^= f[k];
                // Undo the fused division and twiddle multiply-accumulate.
                // Increasing k mirrors the decreasing forward pass: each k
                // reads f[k] before any later operation writes into it.
                if (tw == 0) {
                    for (size_t k = off + h; k < off + 2 * h; ++k) {
                        Elem c = f[k];
                        if (c == 0) continue;
                        for (size_t z : zs) f[k - z] ^= c;
                    }
                } else if (h >= 4) {
                    typename F::Multiplier mul(field, tw);
                    for (size_t k = off + h; k < off + 2 * h; ++k) {
                        Elem c = f[k];
                        if (c == 0) continue;
                        for (size_t z : zs) f[k - z] ^= c;
                        f[k - h] ^= mul(c);
                    }
                } else {
                    for (size_t k = off + h; k < off + 2 * h; ++k) {
                        Elem c = f[k];
                        if (c == 0) continue;
                        for (size_t z : zs) f[k - z] ^= c;
                        f[k - h] ^= field.Mul(tw, c);
                    }
                }
            }
        }
    }

    /** out = a*b with trailing zeros stripped. out must not alias a or b. */
    void MulFull(const std::vector<Elem>& a, const std::vector<Elem>& b, std::vector<Elem>& out, const F& field) {
        if (a.empty() || b.empty()) {
            out.clear();
            return;
        }
        size_t prod = a.size() + b.size() - 1;
        Transform(a, a.size(), b, b.size(), prod, field);
#ifdef MINISKETCH_VERIFY
        // Wraparound tripwire: coefficients at or above the product length
        // must have interpolated to zero.
        for (size_t i = prod; i < buf_a.size(); ++i) CHECK_SAFE(buf_a[i] == 0);
#endif
        out.assign(buf_a.begin(), buf_a.begin() + prod);
        while (!out.empty() && out.back() == 0) out.pop_back();
    }

    /** Transform v into its 2^m-point evaluation table for reuse as the
     *  fixed operand of MulLowPrepared calls. Requires v.size() <= 2^m; out
     *  is caller-owned and must not alias a transform buffer. */
    void Prepare(const std::vector<Elem>& v, int m, const F& field, std::vector<Elem>& out) {
        Extend(field, m);
        CHECK_SAFE(v.size() <= (size_t(1) << m));
        out.assign(size_t(1) << m, 0);
        std::copy(v.begin(), v.end(), out.begin());
        FFT(out, m, field);
    }

    /** out = a*b mod x^n with out.size() == n exactly (the MulLow contract),
     *  where b_fft is Prepare(b, m)'s output and b_size is b's length.
     *  Requires a.size() + b_size - 1 <= 2^m, so the product cannot wrap
     *  around the evaluation set. Costs two transforms instead of MulLow's
     *  three when b's transform is reused across calls. out must not alias
     *  a or b_fft. */
    void MulLowPrepared(const std::vector<Elem>& a, const std::vector<Elem>& b_fft, size_t b_size, int m, std::vector<Elem>& out, size_t n, const F& field) {
        CHECK_SAFE(m >= 1 && m <= max_t && b_fft.size() == (size_t(1) << m));
        CHECK_SAFE(!a.empty() && b_size > 0 && a.size() + b_size - 1 <= (size_t(1) << m));
        size_t n_points = size_t(1) << m;
        buf_a.assign(n_points, 0);
        std::copy(a.begin(), a.end(), buf_a.begin());
        FFT(buf_a, m, field);
        for (size_t i = 0; i < n_points; ++i) buf_a[i] = field.Mul(buf_a[i], b_fft[i]);
        IFFT(buf_a, m, field);
        size_t prod = a.size() + b_size - 1;
#ifdef MINISKETCH_VERIFY
        // Wraparound tripwire, as in MulFull: coefficients at or above the
        // product length must have interpolated to zero.
        for (size_t i = prod; i < buf_a.size(); ++i) CHECK_SAFE(buf_a[i] == 0);
#endif
        out.assign(n, 0);
        size_t copy = std::min(prod, n);
        for (size_t i = 0; i < copy; ++i) out[i] = buf_a[i];
    }

    /** out = a*b mod x^n, with out.size() == n exactly (not stripped, like
     *  PolyMulLowNaive). out must not alias a or b. */
    void MulLow(const std::vector<Elem>& a, const std::vector<Elem>& b, std::vector<Elem>& out, size_t n, const F& field) {
        size_t la = std::min(a.size(), n), lb = std::min(b.size(), n);
        if (n == 0 || la == 0 || lb == 0) {
            out.assign(n, 0);
            return;
        }
        size_t prod = la + lb - 1;
        Transform(a, la, b, lb, prod, field);
        out.assign(n, 0);
        size_t copy = std::min(prod, n);
        for (size_t i = 0; i < copy; ++i) out[i] = buf_a[i];
    }

private:
    /** Compute a[0..la) * b[0..lb) into buf_a (as coefficients), using
     *  transforms of size 2^m with 2^m >= prod = la + lb - 1. */
    void Transform(const std::vector<Elem>& a, size_t la, const std::vector<Elem>& b, size_t lb, size_t prod, const F& field) {
        int m = 1;
        while ((size_t(1) << m) < prod) ++m;
        Extend(field, m);
        size_t n_points = size_t(1) << m;
        buf_a.assign(n_points, 0);
        buf_b.assign(n_points, 0);
        for (size_t i = 0; i < la; ++i) buf_a[i] = a[i];
        for (size_t i = 0; i < lb; ++i) buf_b[i] = b[i];
        FFT(buf_a, m, field);
        FFT(buf_b, m, field);
        for (size_t i = 0; i < n_points; ++i) buf_a[i] = field.Mul(buf_a[i], buf_b[i]);
        IFFT(buf_a, m, field);
    }

    std::vector<Elem> beta;                 //!< Cantor basis: beta[0] = 1, Sqr(beta[j]) ^ beta[j] == beta[j-1].
    std::vector<Elem> twiddle;              //!< twiddle[i] = XOR of beta[j+1] over set bits j of i; size 2^(max_t-1).
    std::vector<std::vector<size_t>> zeta;  //!< Per-dimension sparse offsets of the vanishing polynomials.
    std::vector<Elem> buf_a, buf_b;         //!< Transform buffers, reused across calls.
    int max_t = 0;
};

#endif
