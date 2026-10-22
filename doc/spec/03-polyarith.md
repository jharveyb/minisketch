# 3. Polynomial arithmetic

Everything here is informative: the decoder needs exact polynomial
products, quotients, remainders and gcds over K, and any correct algorithm
gives the same decode result. This chapter describes the algorithms and
tier choices of `src/sketch_impl.h` and `src/additive_fft.h` so that a port
can match their performance. Cutoff values are collected in chapter 7.

## Representation

A polynomial is a vector of field elements, lowest degree first. Most
functions keep vectors **stripped** (no trailing zeros; the zero polynomial
is the empty vector). A few return fixed-length, unstripped results (noted
below), and remainders of length deg(mod) are sometimes left unstripped.

Two predicates select tiers throughout:

- **FFT-capable**: b ∈ {16, 32, 64} (`TraceModFFTCapable`: b ≥ 16 and a
  power of two; the Cantor basis needs b = 2^ℓ, and 2, 4, 8 never reach
  the product-length cutoff).
- **Cheap multiplier**: CLMUL fields (chapter 1).

## Schoolbook products

`PolyMulLowNaive(a, b, n)` (`src/sketch_impl.h:342`) returns the low n
coefficients of a·b (length exactly n, unstripped). CLMUL fields compute it
in convolution order, out[k] = Reduce(⊕_i MulWide(a_i, b_(k−i))), one
reduction per output. Generic fields accumulate rows: for each nonzero a_i,
out[i + j] ^= Multiplier(a_i)(b_j), skipping the table build when a_i = 1.

`Sqr(p)` squares in place: coefficient i of p^2 is Sqr(p_(i/2)) for even i
and 0 for odd i (Frobenius).

## Karatsuba (characteristic 2)

`TraceModPolyMulFull(a, b)` (`src/sketch_impl.h:533`), full product,
stripped:

```
mul_full(a, b):
    if a or b empty: return []
    if fft_eligible(|a|, |b|): return fft.MulFull(a, b)
    if min(|a|, |b|) <= 24: return naive(a, b)          # TRACEMOD_POLYMUL_CUTOFF
    s = ceil(max(|a|, |b|) / 2)
    a0, a1 = a[0:s], a[s:];  b0, b1 = b[0:s], b[s:]
    z0 = mul_full(a0, b0)
    z2 = mul_full(a1, b1)                    # if both nonempty, else 0
    z1 = mul_full(strip(a0+a1), strip(b0+b1))
    return strip(z0 + x^s (z1 + z0 + z2) + x^(2s) z2)
```

The middle term uses (a0 + a1)(b0 + b1) + a0b0 + a1b1 = a0b1 + a1b0, where
"+" is XOR (C39c).

`TraceModPolyMulLow(a, b, n)` (`src/sketch_impl.h:454`) returns a·b mod x^n
(length exactly n):

```
mul_low(a, b, n):
    if n == 0 or a, b empty: return zeros(n)
    if fft_eligible(min(|a|,n), min(|b|,n)): return fft.MulLow(a, b, n)
    if n <= 24 or min(|a|, |b|) <= 24: return naive_low(a, b, n)
    s = ceil(n / 2)
    a0, b0 = a[0:s], b[0:s]
    out = (mul_full(a0, b0) truncated to n), zero-padded to n
    if s >= n: return out
    a1, b1 = a[s:n], b[s:n]                  # truncated to n
    if a1, b1 both nonempty:
        out[s:] += mul_low(a0+a1, b0+b1, n-s) + mul_low(a1, b1, n-s) + z0[0:n-s]
    else:
        out[s:] += mul_low(a1, b0, n-s)  (if a1 nonempty)
        out[s:] += mul_low(a0, b1, n-s)  (if b1 nonempty)
    return out
```

(The a1·b1·x^(2s) term lies entirely above x^n.)

`fft_eligible(la, lb)` (`TraceModFFTEligible`, `src/sketch_impl.h:411`):
FFT-capable, la, lb > 0, min(la, lb) > 24, la + lb − 1 ≥ 128
(`TRACEMOD_FFT_CUTOFF`), and the transform dimension
t = max(1, ⌈log2(la + lb − 1)⌉) satisfies t ≤ b (binding only for b = 16).
Note this applies to generic fields too (their FFT is slower but still
asymptotically better).

Temporaries come from `TraceModScratch` (9 slots per recursion level in a
deque, reused across calls, and one `AdditiveFFT` engine whose tables
amortize across calls).

## Multiplication by a square, Newton inverse series

`TraceModMulBySquareLow(f, g, n)` (`src/sketch_impl.h:600`) computes
f·g^2 mod x^n. With y = x^2, g(x)^2 = h(y), h_i = g_i^2, and
f = f_even(y) + x·f_odd(y), so the result interleaves
f_even·h mod y^⌈n/2⌉ (even coefficients) and f_odd·h mod y^⌊n/2⌋ (odd
coefficients): two half-size products.

`TraceModInvSeries(f, n)` (`src/sketch_impl.h:636`) computes 1/f mod x^n for
f(0) = 1 by Newton iteration in characteristic 2 (C39d): if f·g = 1 + e then
f·(f g^2) = (1 + e)^2 = 1 + e^2.

```
inv_series(f, n):
    g = [1]; m = 1
    while m < n:
        m2 = min(2m, n)
        g = f * g^2 mod x^m2          # mul_by_square_low
        m = m2
    return g zero-padded to n
```

## Transposed Karatsuba middle product (C17)

The middle product of c (length n) and w (length 2n − 1) is
MP(c, w)_s = Σ_{k<n} c_k w_(s+k), s < n. `LazyMiddleProduct`
(`src/sketch_impl.h:2414`) returns it unreduced (as `Wide` values):

```
mp(c, w, n):
    if n <= 16 or n odd:                       # MP_BASE
        return [⊕_k MulWide(c[k], w[s+k]) for s < n]
    h = n/2
    c0, c1 = c[0:h], c[h:n]
    W0, W1, W2 = w[0:2h-1], w[h:3h-1], w[2h:4h-1]   # length-(n-1) windows
    alpha = mp(c0 + c1, W1, h)
    beta  = mp(c0, W0 + W1, h)
    gamma = mp(c1, W1 + W2, h)
    return (alpha + beta) ++ (alpha + gamma)
```

Three half-size middle products instead of four, with only XORs of the
inputs; every output is reduced once by the caller (~33 instead of 65
carry-less products per output at n = 64). Users: the BM tail check
(chapter 4) and the power-sum split (chapter 5). `middle_product_kara` in
`tests/pyminisketch.py` is the same recursion.

## The additive FFT

`AdditiveFFT` (`src/additive_fft.h`) evaluates polynomials over affine
GF(2)-subspaces of K. It requires b to be a power of two; the decoder uses
it only for FFT-capable fields.

### Cantor basis (C38)

β_0 = 1 and β_j = Qrt(β_(j−1)), i.e. β_j^2 + β_j = β_(j−1), for
1 ≤ j < t (`Extend(t)`, t ≤ b). For b = 2^ℓ the equation is always
solvable (Tr(β_(j−1)) = 0 for j < b) and the β_j are GF(2)-independent
(C29, C38). The basis depends on which Qrt solution is taken (bit 0 clear,
chapter 1); a different choice gives a different point order only.

The point with index i is η_i = ⊕_{k ∈ bits(i)} β_k. W_m = {η_i : i < 2^m}
is the subspace spanned by β_0..β_(m−1), and the coset of W_m with block
number B is C(m, B) = {η_(B·2^m + j) : j < 2^m} = η_(B·2^m) + W_m.

### Subspace (vanishing) polynomials (C29, C39a)

S(x) = x^2 + x, S^0 = x, S^(i+1) = S(S^i). Then

```
S^i(x) = Σ_{j : (i AND j) == j} x^(2^j)
```

(coefficients in GF(2), by Lucas' theorem), S^i is GF(2)-linear, monic of
degree 2^i, (S^i)' = 1, S^i vanishes exactly on W_i, and
S^i(β_(j+i)) = β_j. Hence the vanishing polynomial of C(m, B) is

```
Z_C(x) = S^m(x) + S^m(η_(B·2^m)) = S^m(x) + η_B,   Z_C' = 1.
```

`zeta[q]` lists the offsets 2^q − 2^j for the lower terms x^(2^j) of S^q
(j a proper bit-subset of q), and `twiddle[i] = ⊕_{k ∈ bits(i)} β_(k+1)`
(size 2^(t−1)), which equals S^(p−1)(η_(i·2^p)) for every p: the value of
S^(p−1) on the first half of block i at any level.

### Novel basis (Lin–Chung–Han)

X_k(x) = Π_{j ∈ bits(k)} S^j(x), k < 2^m, is a basis of polynomials of
degree < 2^m. Conversions, in place on a length-2^m vector:

```
ToNovel(f, m):                           # monomial -> novel
    for r in 0 .. m-1:
        p = m - r; h = 2^(p-1)
        for each block [off, off + 2h):
            for k from off+2h-1 down to off+h:     # decreasing k is essential
                c = f[k]
                if c: for z in zeta[p-1]: f[k - z] ^= c
FromNovel: the same rounds in reverse order (r = m-1 .. 0), k increasing.
```

Each block step divides the block by S^(p−1) (leading term x^h, the lower
terms XORed in by the offsets): the high half becomes the quotient, the low
half the remainder. No multiplications.

### Butterflies and coset transforms

```
NovelFFT(f, m, B):                       # novel coeffs -> values on C(m, B)
    for r in 0 .. m-1:
        h = 2^(m-r-1)
        for i in 0 .. 2^r - 1:
            lo = f[2h*i : 2h*i+h];  hi = f[2h*i+h : 2h*(i+1)]
            tw = twiddle[i + (B << r)]
            lo[k] ^= tw * hi[k];  hi[k] ^= lo[k]     for k < h   (tw = 0: only hi ^= lo)
NovelIFFT: rounds r = m-1 .. 0, per butterfly hi[k] ^= lo[k]; lo[k] ^= tw * hi[k].
```

Output index j holds the value at η_(B·2^m + j). With f = g0 + S^(p−1)·g1
on a block, S^(p−1) takes the value tw on the first half of the block's
points and tw + 1 on the second, giving g0 + tw·g1 and g0 + (tw + 1)·g1
(`novel_fft` in `tests/pyminisketch.py` is the recursive form). The coset
B is a block of a larger transform, so the same twiddle table serves:
requires (B + 1)·2^m ≤ 2^t.

`FFT(f, m, B)` = NovelFFT(ToNovel(f), m, B) and `IFFT` its inverse. The
code fuses the two: in round r it divides each block by S^(p−1) and applies
the twiddle multiply-accumulate in the same descending loop, then
hi ^= lo. The division and the butterfly commute (both are linear and act
blockwise), so this equals the composition. Cost: ½·N·log2 N
multiplications per size-N transform. Generic fields cache Multiplier
tables for twiddles 0..63 (`TW_MUL_CACHE`), the large blocks of the early
rounds; a performance detail only.

### Formal derivative in the novel basis (C34)

Since (S^j)' = 1, X_k' = Σ_{j ∈ bits(k)} X_(k − 2^j). `NovelDerivative`
(`src/additive_fft.h:312`):

```
out[i] = ⊕_{j < m, bit j of i clear} f[i + 2^j]
```

one XOR pass, no multiplications.

### Products

- `MulFull(a, b)`: t = max(1, ⌈log2(|a| + |b| − 1)⌉); FFT both
  zero-padded to 2^t, multiply pointwise, IFFT, take |a| + |b| − 1
  coefficients, strip. The product length must not exceed 2^t (no
  wraparound).
- `MulLow(a, b, n)`: the same on a, b truncated to n, result truncated or
  zero-padded to exactly n.
- `Prepare(v, m)`: the 2^m-point FFT of v, kept by the caller.
- `MulLowPrepared(a, v_fft, |v|, m, n)`: FFT(a), pointwise product with
  v_fft, IFFT; requires |a| + |v| − 1 ≤ 2^m. Two transforms instead of
  three when v is reused.

## Division

- `PolyMod(mod, val)`: schoolbook remainder by a monic divisor, stripped.
- `DivMod(mod, val, div)` (`src/sketch_impl.h:221`): quotient and remainder
  by a monic divisor. CLMUL fields with deg(mod) ≥ 4 (`modsize >
  DIVMOD_LAZY_MIN = 4`) use `LazyDivRem`; otherwise schoolbook rows with a
  Multiplier per quotient coefficient. The remainder has length deg(mod).
- `MakeMonic(p)`: multiply by the inverse leading coefficient; returns that
  inverse (0 if already monic).

`LazyDivRem(a, b, inv)` (`src/sketch_impl.h:201`), division in dot-product
form: for a of length n_a, b of degree d_b with inv = b_(d_b)^(−1), the
quotient q (length n_q = n_a − d_b) and remainder r (length d_b) are

```
q_i = inv * Reduce( WideFrom(a_(d_b+i)) ⊕ ⊕_{t=i+1}^{min(n_q-1, i+d_b)} MulWide(q_t, b_(d_b+i-t)) ),  i = n_q-1 .. 0
r_j = Reduce( WideFrom(a_j) ⊕ ⊕_{t=0}^{min(n_q-1, j)} MulWide(q_t, b_(j-t)) ),             j < d_b
```

(skip the multiplication by inv when inv = 1). Each coefficient is one
register accumulation; no accumulator array.

`FastDivMod(mod, val, cutoff)` (`src/sketch_impl.h:700`), reciprocal
division by a monic divisor of length m: with k = |val| − m + 1 (quotient
length), if min(k, m − 1) ≤ cutoff use `DivMod`; else

```
inv   = inv_series(rev(mod), k)
q_rev = rev(val)[0:k] * inv mod x^k          # top k coefficients of val, reversed
div   = rev(q_rev)
val   = (val + div * mod mod x^(m-1))[0:m-1] # remainder, length m-1, unstripped
```

The cutoff is `FastDivModCutoff` = 1024 for FFT-capable cheap fields, ∞
otherwise (chapter 7).

## Greatest common divisor

All gcd routines return the gcd up to a scalar factor; callers make it
monic. If the second operand is 0 the result is the first operand,
unchanged.

### Plain Euclid

`GCD(a, b)` (`src/sketch_impl.h:266`), generic fields or when
|b| ≤ DIVMOD_LAZY_MIN = 4 at entry:

```
if |a| < |b|: swap
while b nonzero:
    if deg b == 0: return 1
    make b monic; a = a mod b; swap(a, b)
return a
```

### Normalization-free lazy Euclid (CLMUL fields, |b| > 4 at entry)

The divisor is never made monic; its inverse leading coefficient is folded
into the quotient by `LazyDivRem(a, b, Inv(lc(b)))`. Normal steps
(deg a = deg b + 1) with deg b < 128 (`GCD_INVFREE_MAX`) use the
inversion-free pseudo-remainder (C39e): with β = lc(b), a_n = a_(d_b+1),
a_(n−1) = a_(d_b),

```
p1 = β·a_n
p0 = β·a_(n-1) + a_n·b_(d_b-1)
r  = β²·a − (p1·x + p0)·b        # = β² (a mod b); degree < d_b
r_0 = Reduce(MulWide(β², a_0) ⊕ MulWide(p0, b_0))
r_j = Reduce(MulWide(β², a_j) ⊕ MulWide(p0, b_j) ⊕ MulWide(p1, b_(j-1))),  1 <= j < d_b
```

The loop runs while deg b ≥ 1. If it ends with b a nonzero constant the
inputs are coprime (result 1); otherwise the last divisor a is made monic
once (only if at least one division happened).

### No half-gcd

libminisketch uses only the quadratic `GCD` above. A recursive half-gcd
was implemented and removed: at full fill it only engaged for decodes of
more than about 8192 elements and gained 4–6% there
(`doc/optimizations.md`, "Half-gcd removed"). A port need not implement
one.

## Batch inversion

`BatchInvert(in)` (`src/sketch_impl.h:1350`), Montgomery's trick: prefix
products, one `Inv`, then a backward pass; three multiplications per
element. Returns false if any input is 0 (the root finder uses this to
detect a root on a coset).
