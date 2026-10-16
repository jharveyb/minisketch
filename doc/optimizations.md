# Decode optimization notes

This document records the cost model of `minisketch_decode`, the algorithmic
optimizations applied to it, and — mainly — the investigated-but-deferred
follow-ups, so that future work can pick them up without redoing the analysis.
See [math.md](math.md) for the underlying mathematics of the decoder and
[testing.md](testing.md) for the workflow used when swapping out a decode
stage.

## Where decode time goes

Decoding a sketch of capacity `c` containing `n` differences over
GF(2<sup>m</sup>) runs three stages (`src/sketch_impl.h`):

1. **Syndrome reconstruction** — the even syndromes s<sub>2k</sub> =
   s<sub>k</sub>², one squaring each. Since 2026-10 only for the
   Berlekamp-Massey prefix (~2n), not all 2c.
2. **Berlekamp-Massey** — O(n²) for the LFSR itself (only the odd steps;
   see "Frobenius-structured Berlekamp-Massey" below) plus an
   independent-position check of the remaining O(c) recurrence constraints
   once the LFSR has converged: O(c·n) single carry-less products for
   small n, O(c·n<sup>0.585</sup>) with lazy Karatsuba middle products,
   O(c log n) with the additive FFT for large n. Before 2026-10 this stage
   was a serial O(c·n) loop over all 2c syndromes, which dominated decode
   whenever n ≪ c.
3. **Root finding** (`FindRoots`/`RecFindRoots`, the Berlekamp Trace
   Algorithm). Every split needs the trace Tr(βx) mod the node polynomial,
   i.e. the Frobenius orbit x<sup>2<sup>i</sup></sup> mod f for i < m: `m−1`
   squarings + reductions. Since 2026-10 that orbit is computed **once, at
   the root** (degree n), and the traces of the top ~log₂m+2 split levels
   are linear combinations of it; deeper nodes need one division each.
   Cost: ~`m·n²/2` carry-less products for the root chain (square-table
   reducer, lazy reduction) plus ~`n²` per orbit level for the node
   divisions, plus the per-node GCD and division (~`n²` per level).
   Previously every node ran its own chain, ~`m·n²` in total (level `k`
   costing ~`m·n²/2^k`).

   For roots of degree ≥ 224 on CLMUL/PMULL fields (the "coset Montgomery"
   series below), the chain is m modular squarings at two size-N
   transforms each (N = 2<sup>⌈log₂n⌉</sup>), the tree carries relative
   traces down in Montgomery form (O(N log N) per level), and nodes of
   degree ≥ 256 split by power sums and a half-size Berlekamp–Massey
   instead of a gcd. The quadratic stages left are the decode
   Berlekamp–Massey and the splits' BMs (~35% at 2048/1536), and the
   gcds of the nodes below degree 256.

For n ≳ 64 below that crossover the root chain is the dominant cost of
decode (~50–60%), and on Haswell it runs at the PCLMULQDQ throughput
limit; see "Remaining" below for why it is hard to shorten.

A note on split retries: a failed trace split at degree `d` has probability
~2<sup>1−d</sup>, so retries are essentially free above degree ~8 and only
matter near the leaves. (This kills some otherwise-attractive optimizations;
see "Demoted ideas" below.)

## Implemented in this series

* **Reducer-powered non-factorizability test.** The fail-fast test
  (`trace² + trace mod poly`) reuses the per-node `TraceMod` reducer
  (`SquareAndReduce`) instead of a schoolbook `PolyMod`, making it ~2× cheaper
  below `TRACEMOD_TABLE_CUTOFF` and subquadratic above. Impact is small by
  design — the test runs at most once per decode (a pass propagates down the
  recursion as `fully_factorizable`, a fail aborts) — so this trims one
  schoolbook reduction from failed decodes; slow-unit replays measure within
  noise. Kept because it is strictly less work per call.
* **Allocation hoisting in the reducer hot path.** `SquareReduce`,
  `ReciprocalReduce`, and the Karatsuba multiplication kernels reuse
  per-object scratch buffers (`TraceModScratch`) instead of allocating
  temporaries per call (previously `Bits-1` allocations per trace and ~6 per
  Karatsuba level). The trace loop is allocation-free in steady state.

### Measured impact

Benchmarks: `bench <syndromes> <errors> <iters>` (per-run best over iters,
then min over 3 invocations), CLMUL implementation on x86-64. The b0 column
is the pre-series baseline (the TraceMod-reducer branch rebased onto the test
stack); b2 is after both changes above.

| syndromes/errors | bits | b0 (ms) | b2 (ms) | delta |
|---|---|---|---|---|
| 64 / 8      | 64 | 0.00953 | 0.00951 | -0.2% |
| 128 / 16    | 64 | 0.03463 | 0.03470 | +0.2% |
| 150 / 150   | 64 | 1.994 | 1.910 | -4.2% |
| 500 / 500   | 64 | 20.175 | 19.884 | -1.4% |
| 600 / 600   | 64 | 29.294 | 28.172 | -3.8% |
| 1024 / 1024 | 64 | 73.247 | 67.978 | -7.2% |

The improvement grows with capacity (the reciprocal-path Karatsuba
temporaries dominate the saved allocations); the 32-bit field and the generic
implementation show the same shape (−1% to −5.6%, largest at 1024/1024). No
configuration regressed beyond run-to-run noise.

## Measured negative result: explicit cubic/quartic root solvers

The README long listed *"explicit formulas for the roots of polynomials of
higher degree than 2"* as a TODO. This was implemented and benchmarked — and
turned out to be a **net regression at every field size**, so it was not
merged. The working implementation, its property tests (including an
exhaustive tiny-field enumeration against a distinct-root-count oracle), and
the analysis live on the branch `experiment-affine-solvers`.

The method (Berlekamp-Rumsey-Solomon): `L(x) = x⁴ + a·x² + b·x` is
GF(2)-linear in the element representation, so the roots of an affine quartic
`L(x) + c` are the solutions of an m×m GF(2) linear system (bit-matrix
Gauss-Jordan). Cubics lift to affine quartics via multiplication by
`(x + a₂)`; general quartics reduce via the shift `x → x + sqrt(a₁/a₃)` plus
a reciprocal substitution. `RecFindRoots` then stops recursing at degree 4
instead of degree 2.

Decode time, 64 syndromes / 8 errors, min-of-3-runs (`bench 64 8 50`):

| bits | generic before | generic after | clmul before | clmul after |
|---|---|---|---|---|
| 8  | 0.00318 ms | 0.00322 ms | 0.00317 ms | 0.00320 ms |
| 11 | 0.00444 | 0.00483 | 0.00342 | 0.00379 |
| 16 | 0.00788 | 0.00852 | 0.00434 | 0.00507 |
| 24 | 0.01239 | 0.01363 | 0.00528 | 0.00701 |
| 32 | 0.01578 | 0.01861 | 0.00568 | 0.00937 |
| 48 | 0.03669 | 0.04257 | 0.00789 | 0.01566 |
| 64 | 0.04784 | 0.05962 | 0.00953 | 0.02345 |

(150/150 at 64 bits: 1.99 ms → 2.34 ms. The regression grows with field
width and never flips sign; there is no field size worth gating on.)

Why it loses: the square-table `TraceMod` reducer already reduces a
degree-3/4 node with `(Bits−1)·d²/2` *table-row* multiplications — about 500
cheap operations at m=64 — while the solver pays `m·(2 Sqr + 2 Mul)` to build
the linear system plus O(m²) transpose/elimination word-ops per node. The
system build alone (~4.5k cycles at m=64) exceeds the whole trace-path cost,
before any elimination work. The literature that reports wins for this method
(e.g. Biswas-Herbert 2009, "BTZ") compares against a *plain* Berlekamp trace
algorithm with schoolbook per-squaring reduction; against a precomputed
reducer the trade flips. The conclusion generalizes: per-node linear algebra
in the *field representation* cannot amortize at degrees this small — the
reducer's per-*polynomial* precomputation is the stronger primitive.

What would change the calculus: a per-*decode* (not per-node) precomputation
reusable across all small-degree nodes, or field sizes small enough that the
m×m solve is trivial — neither applies to the 32/64-bit fields that matter
here.

## Implemented: the additive-FFT multiplication tier

The bottleneck quantity for everything below is `M(d)`, the cost of
multiplying degree-`d` polynomials with coefficients in GF(2<sup>m</sup>).
The classical `O(d log d)` route — roots of unity — does not exist in
characteristic 2 (the multiplicative group has odd order); the
characteristic-2 answer is the **additive FFT**, which evaluates over a
GF(2)-linear subspace instead of a multiplicative subgroup.

`src/additive_fft.h` implements the **Cantor-basis additive FFT** ("Cantor
algorithm", after Cantor 1989, in the presentation of Badakhshan-Samanta-
Gong, SAC 2025): the transform is `t` rounds of division by vanishing
polynomials Z<sub>W<sub>i</sub></sub>(x) = S<sup>i</sup>(x), S(x) = x²+x,
whose coefficients in the Cantor special basis are in GF(2) and *sparse*
(2<sup>wt(i)</sup> terms, by Lucas' theorem) — so divisions are pure XOR and
the only field multiplications are one twiddle multiply-accumulate per
high-half position: ½·N·log₂N per transform. The basis is built with the
existing `Qrt` field operation (β₀ = 1, β<sub>j</sub> = Qrt(β<sub>j−1</sub>);
the chain provably cannot fail in fields of power-of-two degree), and with
evaluation shift 0 the per-block twiddles collapse to one round-independent
array of basis-element XORs. `TraceModPolyMulLow/Full` dispatch to the
engine (owned by `TraceModScratch`, so tables amortize over a node's
2·(Bits−1) reductions) when `TraceModFFTEligible`; `ReciprocalReduce`, the
Newton inverse, and `SquareAndReduce` accelerate transparently.

Scope and tuned constants:

* The Cantor special basis requires field degree 2<sup>ℓ</sup> — the tier
  covers the 16/32/64-bit fields (both implementations). Other fields keep
  Karatsuba; extending to them needs the general-basis Gao-Mateer
  construction (future work — notably the committed worst-case slow-units
  use 54/56/62-bit fields and are untouched by this tier).
* `TRACEMOD_FFT_CUTOFF = 128` (product length; response is nearly flat in
  64..512, so this is not a sensitive knob).
* The square-table/reciprocal crossover halves to
  `TRACEMOD_FFT_TABLE_CUTOFF = 256` — but **only for field implementations
  whose scalar `Multiplier` is a plain wrapper** (the CLMUL fields). The
  generic fields pay a Multiplier table build per FFT block, which keeps
  their crossover at 512; giving them the lower cutoff measurably regressed
  500-1024-syndrome decodes by ~10-25%. `TraceModTableCutoff` distinguishes
  the two by the Multiplier size.

Measured decode time (interleaved A/B vs the pre-FFT branch point,
best-of-iters, min of 3 invocations, otherwise-idle machine):

| syndromes/errors | 64-bit CLMUL | 64-bit generic | 32-bit CLMUL |
|---|---|---|---|
| 150/150   | 1.92 ms (parity) | parity | parity |
| 500/500   | 19.57 → 13.96 ms (1.40×) | parity | 1.36× |
| 1024/1024 | 68.66 → 44.01 ms (1.56×) | 1.16× | 1.50× |
| 2048/2048 | 228.90 → 119.94 ms (1.91×) | 1.24× | 1.68× |
| 4096/4096 | 763.16 → 337.93 ms (2.26×) | 1.45× | 1.84× |

The speedup grows with capacity, as the FFT's `M(d) ≈ 3d·log d` pulls away
from Karatsuba's `d^1.585`; the remaining quadratic terms (GCD, division,
Berlekamp-Massey — see below) now bound the curve. Benchmark hygiene note:
these runs must be taken on an otherwise-idle machine — concurrent IDE
indexing was observed to inflate short rows by ~2.5×, and sub-millisecond
rows show bimodal scheduling jitter that only warm min-of-many runs see
through.

### Cached node transforms ("Phase 2")

Each `ReciprocalReduce` multiplies by the same two per-node-constant
operands: the reversed-modulus inverse series (`inv`) and the modulus
itself (`mod`). Their transforms are prepared once per `TraceMod` node at
a fixed size 2<sup>fft_m</sup> ≥ 2·deg(mod) — which fits every product in
a reduction — and reused through a prepared-operand engine entry point
(`AdditiveFFT::Prepare` / `MulLowPrepared`): a full-size reduction then
costs **4 transforms instead of 6**. The prepared path engages exactly
when the dynamic path would run the FFT at the node size anyway
(`CachedMulEligible`), so the short-value early rounds of a trace keep
their cheaper half-size-FFT/Karatsuba handling, and the two one-time
`Prepare` transforms amortize over the ~2·(Bits−1) reductions of a node.

Measured on top of the FFT tier (3-way interleaved base/FFT/cached, same
protocol; 100% error rows):

| syndromes/errors | 64-bit CLMUL | 64-bit generic | 32-bit CLMUL |
|---|---|---|---|
| 500/500   | 13.81 → 12.34 ms (1.12×) | parity (square-table mode) | 1.09× |
| 1024/1024 | 42.35 → 36.91 ms (1.15×) | 1.25× | 1.13× |
| 2048/2048 | 122.89 → 100.58 ms (1.22×) | 1.23× | 1.11× |
| 4096/4096 | 342.48 → 286.80 ms (1.19×) | 1.19× | 1.14× |

Cumulatively vs the pre-FFT baseline (64-bit CLMUL): 1.58× at 500/500,
1.87× at 1024/1024, 2.26× at 2048/2048, **2.67× at 4096/4096**.

A deeper restructuring that keeps the trace value in the *evaluation
domain* (squaring is pointwise there: `f²(p) = f(p)²` in char 2) was
analyzed and found to be a wash: the quotient needed for reduction must be
built from the *reversed* high coefficients of the squared value, and
reversal is not representable in the evaluation domain, so each round
still needs the same 4 transforms (IFFT of the squared evaluations, FFT of
the reversed value, IFFT of the quotient series, FFT of the quotient) —
the only saving is the already-cheap coefficient-wise squaring pass.
Recorded here so it is not re-derived.

### Decode cost at partial error counts

Decode cost tracks the *actual* number of differences `e`, not the
capacity: Berlekamp-Massey is adaptive O(e·n) (as of this measurement;
see the 2026-10 section for the current BM), and the locator degree —
which drives the whole root-finding stage and hence the FFT tier — equals
`e`. Measured with errors at a fraction of the syndrome count (64-bit
CLMUL, same 3-way protocol; base → FFT → cached, speedup = base/cached):

| syndromes | errors | base (ms) | +FFT | +cached | speedup |
|---|---|---|---|---|---|
| 4096 | 4096 (100%) | 765.16 | 342.48 | 286.80 | 2.67× |
| 4096 | 3072 (75%)  | 483.77 | 228.09 | 190.84 | 2.53× |
| 4096 | 2048 (50%)  | 239.19 | 133.87 | 112.28 | 2.13× |
| 4096 | 1024 (25%)  | 77.23  | 52.34  | 45.42  | 1.70× |
| 4096 | 410 (10%)   | 16.96  | 14.73  | 13.69  | 1.24× |
| 2048 | 2048 (100%) | 226.81 | 122.89 | 100.58 | 2.26× |
| 2048 | 1536 (75%)  | 143.07 | 82.59  | 68.92  | 2.08× |
| 2048 | 1024 (50%)  | 71.82  | 47.71  | 38.41  | 1.87× |
| 2048 | 512 (25%)   | 21.81  | 15.81  | 13.80  | 1.58× |
| 2048 | 204 (10%)   | 4.46   | 4.45   | 4.41   | parity |
| 1024 | 1024 (100%) | 69.12  | 42.35  | 36.91  | 1.87× |
| 1024 | 768 (75%)   | 42.29  | 29.69  | 25.39  | 1.67× |
| 1024 | 512 (50%)   | 20.52  | 15.05  | 12.38  | 1.66× |
| 1024 | 256 (25%)   | 5.87   | 5.20   | 4.69   | 1.25× |
| 1024 | 102 (10%)   | 1.14   | 1.15   | 1.15   | parity |
| 500  | 500 (100%)  | 19.47  | 13.81  | 12.34  | 1.58× |
| 500  | 375 (75%)   | 11.29  | 9.87   | 8.77   | 1.29× |
| 500  | ≤250 (≤50%) | —      | —      | —      | parity |

The speedup is a function of `e` with mild dilution from the O(e·n) BM
term as `n` grows (e.g. e = 1024: 1.87× at n = 1024/2048, 1.70× at
n = 4096); below the reciprocal-reducer crossover (deg < 256 on CLMUL
fields) the FFT tier never engages and every configuration is at parity —
no regression anywhere, in either field implementation. The e = 256 row
engaging (1.25×) while e = 250 does not is the table cutoff, working as
intended.

Remaining follow-ups on this tier:

* **General-basis Gao-Mateer** for non-power-of-two fields.
* **Coset-shifted truncated multiplication** to soften the power-of-two
  padding cliff: for product length 2<sup>t</sup> + c (c small) the engine
  currently pads to 2<sup>t+1</sup>, doubling transform cost. The CRT
  trick of Chen et al. (SFAFFT, §4.2) computes `fg mod q` with a
  size-2<sup>t</sup> transform — wraparound *is* reduction modulo the
  vanishing polynomial of the evaluation set — plus `fg mod x^c` with a
  small Karatsuba product, and recombines; the recombination divides by
  the sparse vanishing polynomial, nearly free XOR. Caveat: it needs
  gcd(Z, x<sup>c</sup>) = 1, and the θ=0 subspace contains 0 (Z has zero
  constant term), so this requires a coset-shifted transform (evaluate on
  β<sub>t</sub> + W<sub>t</sub>; twiddles gain a per-round offset). Payoff
  is up to ~2× on transform cost right above each power of two, shrinking
  to nothing at the next one.
* The **Frobenius additive FFT** (Li et al., ISSAC 2018) does *not* apply:
  its factor-d truncation requires coefficients in GF(2), while the
  decoder's are full field elements.

### Strided Frobenius additive FFT (2026): assessed, core inapplicable

Chen et al., *Strided Frobenius Additive FFT and its Application to HQC*
(2026; `doc/2026-1588_strided_frobenius.pdf`) reframes additive FFTs
ring-theoretically and adds two generalizations: *strided* transforms
(treat x<sup>k</sup> as the variable; stop log₂k butterfly rounds early
and do the pointwise step as length-k polynomial products) and
*incomplete* transforms (its Theorem 1 formalizes that multiplication only
needs the quotient-algebra degree to exceed the product degree — the same
fact `MulFull`'s wraparound tripwire relies on). Assessment:

* The **Frobenius core does not transfer**, for the same reason as the
  ISSAC 2018 paper: evaluation on orbit representatives needs GF(2)
  coefficients so that f(σ²) = f(σ)².
* **Striding is field-agnostic but near-neutral here**: every round of our
  transform costs N/2 twiddle-macs regardless of depth, so truncating
  log₂k rounds saves (3/2)·N·log₂k multiplications across a product's
  three transforms while the pointwise step grows from N to
  (N/k)·(M(k)+k−1) — net ~0.5N saved at k=2, ~2N at k=4, i.e. 3–6% of FFT
  multiplications. The paper's headline wins (butterflies dropping into a
  smaller field, byte-aligned basis conversion, sparse CRT moduli) are
  specific to the F₂[x]/HQC setting.
* The **CRT-with-x^c truncation** (§4.2) is the one transferable idea —
  recorded as the coset-shifted follow-up above.

References:

* Badakhshan, Samanta, Gong, *Accelerating Post-quantum Secure zkSNARKs by
  Optimizing Additive FFT* (SAC 2025;
  `doc/sac2025-2-paper17_optimizing_additive_fft.pdf`) — primary source for
  the implemented algorithm.
* Li, Chen, Kuo, Cheng, Yang, *Frobenius Additive Fast Fourier Transform*
  (ISSAC 2018; `doc/3208976.3208998_frobenius_additive_fft.pdf`,
  <https://arxiv.org/abs/1802.03932>) — Cantor-construction cross-check.
* Chen, Chien, Chiu, Huang, Lin, Peng, Yang, *Strided Frobenius Additive
  FFT and its Application to HQC* (2026;
  `doc/2026-1588_strided_frobenius.pdf`) — assessed above; source of the
  CRT-truncation follow-up.
* Gao, Mateer, *Additive Fast Fourier Transforms over Finite Fields*
  (<https://www.math.clemson.edu/~sgao/papers/GM10.pdf>) — the general-basis
  construction for the non-power-of-two follow-up.
* Bernstein, *Multiplication of polynomials over F₂* survey
  (<https://cr.yp.to/f2mult.html>)

## Implemented: subquadratic division and GCD

With the FFT tier and cached transforms in place, the quadratic
`GCD(trace, poly)` (once per split attempt at every `RecFindRoots` node)
and the quadratic `DivMod` (once per successful split) were the next
structural bottleneck: ~2d² multiplications summed over the recursion, at
capacity 4096 comparable to the whole FFT-accelerated trace stage. Both
are now subquadratic (`src/sketch_impl.h`), following the algorithm choice
of the mid-2026 literature review: **classic recursive half-gcd, gcd-only
— not jumpdivstep** (van der Hoeven arXiv:2212.12389 for the constants;
NTL's `GF2EX` for the recursion shape; polynomial divstep is ~2× behind in
the same multiplication model and stays relevant only for a hypothetical
constant-time decoder).

* **FastDivMod** divides by a monic divisor with the reciprocal method the
  reducer already uses — Newton inverse series of the reversed divisor to
  quotient-length precision, one low product for the reversed quotient,
  one more for the remainder.
* **FastGCD** reduces large pairs with the half-gcd: the transition matrix
  of a stretch of the remainder sequence depends only on the top
  2·d_red−1 coefficients of the pair (the classical truncation lemma), so
  `HalfGCDMatrix` computes it recursively at half size — first-half
  matrix, one explicit division step, second-half matrix, composed with
  2×2 polynomial matrix products — and the driver applies it to the
  full-length pair with four multiplications, halving deg(b) per outer
  round: O(M(d) log d), all products through the additive-FFT tier. No
  Bezout coefficients are computed (the split step needs only the gcd).
  The quadratic loop remains as the base case and as the recursion's
  iterative bottom (`HGCD_ITER_CROSSOVER = 24`).

### What the measurements actually said

The textbook/NTL crossover estimates (degree 40–340) were off by an order
of magnitude here, for the same reason the affine root solvers lost
(above): **the baseline is strong**. The schoolbook remainder steps are
Multiplier-row operations, allocation-free and cache-friendly — a whole
degree-1024 quadratic gcd costs ~1–2 ms — while the half-gcd pays FFT
setup, matrix products, and allocation churn. Interleaved sweeps with
disabled-path control rows (which calibrated ±5% measurement noise on
this hardware, and exposed one systematic effect worth recording):

* Engage cutoffs at or below degree 512 regressed 1024-syndrome decode by
  ~20%; `HGCD_CUTOFF = 1024` is neutral at 1024/2048 syndromes and wins
  at 4096.
* The **generic implementation regressed at every tested cutoff** (its
  FFT butterflies pay Multiplier table builds), so both fast paths gate
  on the cheap-scalar-Multiplier heuristic (`TraceModCheapScalarMul`,
  shared with the table cutoff): generic fields keep the quadratic loop.
* `FAST_DIVMOD_CUTOFF = 1024`: division alone measured neutral below
  that, positive in combination with the half-gcd at 4096.
* **Inlining find**: merely instantiating the fast paths at their
  `RecFindRoots` call sites cost a consistent ~5% at 1024 syndromes even
  though they never executed there — a control build with reverted call
  sites (templates never instantiated) matched the old branch exactly,
  attributing the loss to code bloat in the hot recursion. `NOINLINE` on
  the two entry points (new macro in `util.h`) restored parity. Moral:
  when an A/B shows a regression the algorithm cannot explain, check the
  code layout before the code.

Final measured decode (interleaved, alternating order, min over runs, vs
the cached-transform tip):

| syndromes/errors | 64-bit CLMUL | 32-bit CLMUL | 64-bit generic |
|---|---|---|---|
| 1024/1024 | parity | parity | parity (gated) |
| 2048/2048 | parity | — | parity (gated) |
| 4096/4096 | 284.8 → 264.4 ms (1.08×) | 172.6 → 151.6 ms (1.14×) | parity (gated) |

Cumulative vs the pre-FFT baseline at 4096/4096: **765 → 264 ms (2.89×)
on 64-bit CLMUL**, 350 → 152 ms (2.31×) on 32-bit CLMUL.

### Remaining in this area

* (Superseded: after the 2026-10 lazy-reduction series the plain Euclid
  loop beats the half-gcd up to degree 4096 and `HGCD_CUTOFF` is 4096;
  see below.) **Constant-factor half-gcd work** is the lever to push the
  crossover down again: middle products for the matrix-vector
  applications, transform caching across the 2×2 matrix recursion (the
  vdH / eprint 2024/644 trick; `Prepare`/`MulLowPrepared` is the hook),
  and a scratch pool for the per-call matrix/temporary allocations.
* **Hybrid fast Berlekamp-Massey** via the Dornstetter/key-equation
  equivalence (partial XGCD at degree 2t; modern form: PM-Basis
  approximants, practical in Neiger's PML at degrees 10³–10⁵): only wins
  when the difference count approaches capacity, since the current BM is
  adaptive O(e·n); gate on measured e ≈ n profiles.
* What is still worth taking from safegcd: libsecp256k1's batched
  transition-matrix trick (`doc/safegcd_implementation.md`) as a
  memory-traffic optimization, and divstep's data-independent shape if a
  constant-time decoder is ever required (no published
  divstep-over-GF(2<sup>m</sup>)[x] implementation appears to exist).

References: <https://arxiv.org/abs/2212.12389>,
<https://eprint.iacr.org/2019/266>, <https://eprint.iacr.org/2024/644>,
NTL `GF2EX` HalfGCD (<https://github.com/libntl/ntl>),
<https://github.com/vneiger/pml>,
<https://github.com/bitcoin-core/secp256k1/blob/master/doc/safegcd_implementation.md>,
<https://github.com/sipa/safegcd-bounds>.

## Implemented: Frobenius-structured BM, shared trace orbit, lazy reduction (2026-10)

Target: high capacity (1024+) with a low-to-medium number of differences
(32–256+). In that regime the old decoder spent most of its time in a
Berlekamp-Massey pass over all 2c syndromes and in per-node trace chains.
Everything below is exact: decode returns the same result, including the
same failure behaviour, as before.

Gentler explanations of the same ideas: [decode_explained.md](decode_explained.md)
(undergraduate discrete math), [decode_coding_theory.md](decode_coding_theory.md)
(coding theory) and [math_decode.md](math_decode.md) (in the style of
[math.md](math.md)).

### Frobenius-structured Berlekamp-Massey

The syndromes are power sums in characteristic 2, so s<sub>2k</sub> =
s<sub>k</sub>² (equivalently S(x)² = S<sub>even</sub>(x)). Two consequences:

* **Even steps have zero discrepancy.** This is the classical binary-BCH
  simplification (Berlekamp 1968). It holds for every Frobenius-structured
  sequence, not only genuine sketches, so BM skips those steps. Verified
  exhaustively against a textbook oracle over GF(4)…GF(256), plus the
  property tests.
* **Early termination with an exact tail check.** Once the LFSR length L
  stops growing, the connection polynomial Λ is checked against all
  remaining *odd*-indexed constraints at once. The even ones then follow
  from S² = S<sub>even</sub>:
  [x<sup>2j</sup>](ΛS) = ([x<sup>j</sup>](ΛS))² + [x<sup>2j</sup>](Λ·(ΛS mod x<sup>L+1</sup>)),
  and both terms vanish for j > L once Λ generates the prefix (induction
  over j, as checking starts past 2L). So the full iteration would see zero
  discrepancies through the tail and return the same Λ. That includes
  rejecting a sketch with more than `max_count` differences, since a
  failed check just resumes the iteration.

* **The tail only needs the sketch itself.** With C = ΛS the same
  structure gives Λ·C = C² + Λ²·S<sub>odd</sub>. Past the prefix, C
  vanishes on the tail iff the raw sketch values (the odd syndromes
  u<sub>j</sub> = s<sub>2j+1</sub> = Σζ·(ζ²)<sup>j</sup>) satisfy the
  recurrence whose connection polynomial is Λ with each coefficient
  squared (roots squared). So decode no longer reconstructs all 2c
  syndromes (`BerlekampMasseyOdd` derives s<sub>2k</sub> = s<sub>k</sub>²
  only for the prefix it walks), and the tail check reads half the memory.

The tail check (`BerlekampMasseyVerifyTail`) has these tiers on
lazily-reducing fields:
* **Direct** (L < 64). One carry-less product per term, no serial
  dependency, and each coefficient load shared across 4 positions.
* **Lazy Karatsuba middle products** (64 ≤ L < 512). Each block of n
  positions is a middle product of the reversed Λ² with a window of u.
  The transposed-Karatsuba identity low = α + MP(c₀, w₀+w₁),
  high = α + MP(c₁, w₁+w₂), α = MP(c₀+c₁, w₁) needs three half-size
  middle products instead of four, with only XORs of the input windows,
  and every output is reduced once. That is ~33 products per position at
  L = 64 instead of 65. It beats the direct check from L ≈ 32 on the i7
  and ≈ 96 on the Pi, and the FFT up to L ≈ 384 on both.
* **FFT** (L ≥ 512). One product per block of Λ² with a window of u,
  with Λ²'s transform prepared once: two transforms per block.

Generic fields use the FFT from L = 40, and below that a direct check
that builds each coefficient's Multiplier table once for the whole tail.

BM goes from a serial O(c·n) to O(n²) + O(c·n) independent products, and
it builds Multipliers only for the processed prefix. On the generic
fields that was a 2 KB table per syndrome: 16 MB at capacity 4096, a real
cost on small caches.

For table-based (generic) fields, the additive FFT now builds the
Multiplier tables of its 64 lowest twiddles once per engine instead of once
per block per transform. That halved generic BM time on a Pi 5.

### One Frobenius orbit for the whole root-finding tree

All nodes at depth k split with the same parameter β<sub>k</sub>
(basis·2<sup>k</sup>, or the tower basis below), and every node g divides
the root polynomial f, so
Tr(β<sub>k</sub>x) mod g = (Tr(β<sub>k</sub>x) mod f) mod g.
`ComputeTraceOrbit` generates X<sub>i</sub> = x<sup>2<sup>i</sup></sup> mod f
once (m−1 modular squarings, the cost of one trace) and accumulates
Tr(β<sub>k</sub>x) mod f = Σ<sub>i</sub> β<sub>k</sub><sup>2<sup>i</sup></sup> X<sub>i</sub>
for the first log₂m+2 depths on the fly. Only those accumulators stay live,
not the m orbit elements. A node of degree d<sub>g</sub> with
2·deg f ≤ d<sub>g</sub>·m then needs a single division instead of a
reducer build plus m−1 squarings; smaller nodes keep the old path.

**Tower-basis split parameters.** The split parameters only need to be
GF(2)-linearly independent. When 8 divides m they are
β<sub>8l+j</sub> = basis·x<sup>l</sup>·z<sup>j</sup>, with z generating the
subfield GF(2⁸), found via the relative trace. For the first eight
depths β<sub>k</sub> = basis·z<sup>k</sup>, and (z<sup>k</sup>)<sup>2<sup>i</sup></sup>
depends only on i mod 8. So the chain feeds eight residue buckets
Z<sub>r</sub> = Σ<sub>i≡r</sub> basis<sup>2<sup>i</sup></sup>X<sub>i</sub>,
and each depth's trace is Σ<sub>r</sub>(z<sup>k</sup>)<sup>2<sup>r</sup></sup>Z<sub>r</sub>.
That is (m + 8·depths)·d products instead of m·depths·d, i.e. 128·d
instead of 512·d at m = 64, and gives +11–14% at n = 16…64.

One more squaring gives x<sup>2<sup>m</sup></sup> mod f, which equals x
iff f is square-free and splits completely. That is an **exact
factorizability test at the root**: undecodable sketches fail before any
splitting. (This differs from the "Frobenius-power precomputation for
trace retries" demoted below: that cached per node for retries, while this
shares the root's orbit across all nodes of a depth.)

Why the root chain stays: any split needs a nonconstant GF(2)-valued
function of the roots. Such a function is a trace and involves every
Frobenius power. The routes around the chain all cost more than m cheap
squarings at these degrees, because squaring is linear in characteristic
2: subfield/relative traces still need m powers in total, and
modular-composition doubling (Brent–Kung) or Hankel solves from the
sketch's extra syndromes cost more. The power-sum (dual) representation,
where squaring is free at even indices, is circular by the transposition
principle.

### Lazy modular reduction (CLMUL/PMULL fields)

Reduction modulo the field polynomial is GF(2)-linear, so a sum of
carry-less products can be XOR-accumulated in 128 bits and reduced once.
That is one CLMUL per term instead of three. The CLMUL `GenField` now
exposes `Wide`/`MulWide`/`WideXor`/`WideFrom`/`Reduce` (`Mul` is
unchanged); the generic fields expose trivial versions with
`LAZY_REDUCE = false` and keep their Multiplier loops. Uses:

* **Square-table reducer** (the root chain and small-node traces):
  register-blocked over 8 output columns, with the table stored
  column-block-major so each block streams contiguous memory. The
  row-major table had one cache line per row 8·d bytes apart, more
  concurrent streams than the prefetchers track. On Haswell this loop now
  runs at ~2.1 cycles per product, the PCLMULQDQ throughput limit.
* **BM discrepancy and tail check**, the orbit accumulation, and the
  schoolbook product (convolution order, register accumulator).
* **Division in dot-product form** (`LazyDivRem`): each quotient and
  remainder coefficient is one register-accumulated dot product, reduced
  once, with no accumulator array. DivMod and Euclid both use it.
* **Normalization-free Euclid.** The divisor's inverse leading coefficient
  is folded into the quotient instead of rescaling the divisor every step.
  Below degree 128, normal steps (quotient degree 1) use the
  inversion-free pseudo-remainder β²a − (βa<sub>n</sub>x + βa<sub>n−1</sub> −
  a<sub>n</sub>b<sub>m−1</sub>)b, trading one product per coefficient for
  a field inversion.

The cheaper schoolbook primitives moved several crossovers:
* Table/FFT reducer: 256 → 704 for CLMUL fields (~700 on Haswell, above
  768 on Cortex-A76).
* `HGCD_CUTOFF`: 1024 → 4096; the lazy Euclid beats the half-gcd at every
  degree up to 4096.
* BM tail check: direct < 64 ≤ Karatsuba < 512 ≤ FFT on CLMUL fields;
  FFT from L ≥ 40 on generic fields.
* Shared orbit: from root degree 8 (was 16).

### Measured

Interleaved A/B, min over 3 rounds (`bench --decode-only`), 64-bit field
unless noted. "Before" is the tip before this series and "after" is the
end of the series, all commits included.

| capacity/differences | i7-4770 CLMUL | Pi 5 (A76) PMULL |
|---|---|---|
| 1024/16 | 0.123 → 0.034 ms (3.6×) | 0.224 → 0.046 ms (4.9×) |
| 1024/32 | 0.318 → 0.077 ms (4.1×) | 0.589 → 0.104 ms (5.7×) |
| 1024/64 | 0.900 → 0.201 ms (4.5×) | 1.697 → 0.286 ms (5.9×) |
| 1024/128 | 2.80 → 0.66 ms (4.3×) | 5.29 → 0.90 ms (5.9×) |
| 1024/256 | 8.18 → 2.41 ms (3.4×) | 14.90 → 3.03 ms (4.9×) |
| 1024/512 | 22.3 → 9.0 ms (2.5×) | 39.9 → 10.8 ms (3.7×) |
| 2048/64 | 1.154 → 0.228 ms (5.1×) | 2.174 → 0.329 ms (6.6×) |
| 2048/256 | 9.19 → 2.50 ms (3.7×) | 16.92 → 3.16 ms (5.4×) |
| 4096/32 | 0.697 → 0.133 ms (5.2×) | 1.344 → 0.168 ms (8.0×) |
| 4096/64 | 1.669 → 0.285 ms (5.9×) | 3.172 → 0.414 ms (7.7×) |
| 4096/128 | 4.34 → 0.78 ms (5.5×) | 8.29 → 1.09 ms (7.6×) |
| 4096/256 | 11.10 → 2.64 ms (4.2×) | 20.84 → 3.29 ms (6.3×) |
| 8192/64 | 2.744 → 0.400 ms (6.9×) | 5.159 → 0.582 ms (8.9×) |
| 8192/256 | 15.13 → 2.90 ms (5.2×) | 28.66 → 3.66 ms (7.8×) |
| 1024/1024 | 66.9 → 21.2 ms (3.2×) | 120.5 → 35.0 ms (3.4×) |
| 32-bit 1024/64 | 0.778 → 0.165 ms (4.7×) | 0.694 → 0.204 ms (3.4×) |
| 32-bit 4096/256 | 9.89 → 1.86 ms (5.3×) | 9.08 → 2.16 ms (4.2×) |
| generic 1024/64 | 3.93 → 1.96 ms (2.0×) | 5.85 → 2.70 ms (2.2×) |
| generic 4096/64 | 9.56 → 2.92 ms (3.3×) | 12.76 → 4.18 ms (3.1×) |

Fuller sketches (50–100% of capacity) gain less: the BM tail is short
and the root degree reaches the FFT reducer (≥ 704). Same protocol:

| capacity/differences | i7-4770 CLMUL | Pi 5 (A76) PMULL |
|---|---|---|
| 1024/512 (50%) | 22.1 → 7.9 ms (2.8×) | 40.0 → 10.0 ms (4.0×) |
| 1024/768 (75%) | 44.7 → 17.0 ms (2.6×) | 82.1 → 28.3 ms (2.9×) |
| 2048/1024 (50%) | 71.1 → 21.6 ms (3.3×) | 129.4 → 35.7 ms (3.6×) |
| 2048/1536 (75%) | 124.5 → 44.3 ms (2.8×) | 225.2 → 72.6 ms (3.1×) |
| 4096/2048 (50%) | 201.7 → 60.4 ms (3.3×) | 366.5 → 97.2 ms (3.8×) |
| 4096/3072 (75%) | 335.8 → 126.0 ms (2.7×) | 600.5 → 202.4 ms (3.0×) |
| 4096/4096 (100%) | 504.5 → 177.1 ms (2.8×) | 881.4 → 279.0 ms (3.2×) |
| 32-bit 1024/768 (75%) | 30.1 → 11.0 ms (2.7×) | 28.8 → 13.1 ms (2.2×) |
| generic 1024/768 (75%) | 208.7 → 117.5 ms (1.8×) | 339.5 → 195.7 ms (1.7×) |

**Memory.** Peak heap during one decode (global `operator new`
instrumented, 64-bit field), before → after:

| capacity/differences | CLMUL | generic |
|---|---|---|
| 1024/32 | 57 → 27 KB | 4.2 MB → 236 KB |
| 1024/256 | 107 → 349 KB | 4.2 → 1.3 MB |
| 1024/512 | 299 KB → 1.21 MB | 4.2 → 2.4 MB |
| 1024/1024 | 333 → 458 KB | 4.2 → 4.3 MB |
| 4096/64 | 229 → 109 KB | 16.9 MB → 574 KB |
| 4096/4096 | 1.37 → 1.81 MB | — |

The old BM built a 2 KB Multiplier table per syndrome on the generic
fields; it now builds them only for the walked prefix, in a deque. CLMUL
decodes with root degree in [256, 704) keep the square table (d²/2
elements, ~2 MB at the top of that range) that the raised table/FFT
crossover uses there instead of the FFT reducer. That is the one real
memory cost of the series, traded for speed on the small-cache machines.
Allocation counts per decode dropped 1.2–4× throughout.

Failed decodes (more differences than capacity) got faster too, from the
exact root-level factorizability test. On the i7, capacity 1024 with 1100
differences went 28.4 → 12.6 ms, and capacity 4096 with 4200 went
192 → 76 ms.

Stage split at 1024/256 on the i7 afterwards: the root orbit (chain,
accumulation, exact test) ~57%, node divisions and small-node traces
~24%, Euclid ~9%, BM ~7%. At 8192/128 BM (the tail check) is ~39%.

Per-commit gains for the later commits, against the preceding commit on
i7-4770 / Pi 5:
* Tower basis: 1024/32 +12% / +12%, 1024/64 +11% / +7%.
* Odd-syndrome BM: 4096/16 +27% / +45%, 4096/32 +14% / +28%,
  8192/64 +9% / +24%.
* Orbit from root degree 8 (was 16): 1024/12 +43% / +51%.
* Karatsuba middle-product tail check: 4096/128 +12% / +2%,
  8192/128 +22% / +5%.

### Measured negative results in this series

* **Lazy-Karatsuba reciprocal reducer for the chain.** The Frobenius-sparse
  first Barrett product splits into two half-size products, with unreduced
  Karatsuba down to a schoolbook base of 8–32. It was slower than the
  streaming lazy table at every degree on the i7 (n=128 1.11 vs 0.80 ms,
  n=256 3.29 vs 2.66 ms, n=1024 29.5 vs 22.4 ms). With base b,
  M(n) = (n/b)<sup>1.585</sup>b², so at d = 256 the three products cost
  ~34.5k CLMULs (b = 16), against d²/2 = 32k for the table plus no
  recursion overhead. The old non-lazy Karatsuba reciprocal path (FFT tier
  disabled) was worse still.
* **Reciprocal (FFT) reducer on generic fields**, even with the cached
  twiddle tables: the pointwise products are slow one-off GFMuls. At
  1024/256 on the Pi, 22.0 vs 15.1 ms for the table.
* **Orbit node-factor / depth retune** after the lazy division: flat
  within ±3%.
* **Operand order for SSE PCLMULQDQ** (put the loop-invariant operand
  second, to avoid a register copy per product without VEX encoding): no
  measurable effect. The copies are eliminated at register rename; their
  apparent cost in `perf annotate` was skid.

### Probabilistic tail verification (not implemented: changes guarantees)

At high capacity and small n the exact O(c·n) (or O(c log n)) tail check
is a large share of decode: about half at 4096/32 before the odd-syndrome
change. A **secret random linear combination** of the tail constraints,
Σ<sub>k</sub> r<sup>k</sup>c<sub>k</sub> = Σ<sub>i</sub> Λ<sub>i</sub>r<sup>i</sup>(P(N−i) − P(t−i))
with prefix sums P(j) = Σ<sub>j'<j</sub> r<sup>j'</sup>s<sub>j'</sub>,
costs O(c + n), roughly 2× decode at 4096/32. But it accepts a sketch
whose tail does not match with probability up to N/2<sup>m</sup> per
random r (Schwartz–Zippel): ~2<sup>−51</sup> at m = 64, but 2<sup>−19</sup>
at m = 32. That weakens the decoder's guarantee against adversarially
built overfull sketches (low-weight BCH codewords make prefix-only
matches constructible), so it would need an explicit API-level opt-in
and several rounds for small fields. Exact verification is Ω(c log n)-ish
(it is a polynomial product), so this is the only route to O(c).

### Remaining in this area

* **Lehmer-style blocked Euclid**: compute k quotient steps from the top
  2k coefficients and apply the composed 2×2 matrix in one lazy pass. That
  cuts the remainder reductions (half of the remaining GCD cost) by ~k/2.
* **Lazy square-table construction**: the shift-and-fold builds 3 CLMULs
  per entry today; a two-step fold with lazy accumulation needs ~2.
* **Orbit-node divisions at n ≥ 2048** would benefit from fast (FFT)
  remainders at the top levels (`FastDivMod` currently never engages
  there).
* The **chain itself** is at the hardware throughput bound on Haswell.
  Shortening it needs a cheaper Frobenius than modular squaring; see
  above for why the obvious routes do not pay.

## Implemented: coset Montgomery root finding and power-sum splits (2026-10, high fill)

This series targets large, full sketches: capacity 1024–4096 with 50–100%
of it in use, so root degrees d ≥ 256. Before it, root finding at those
degrees cost (per decode) the root chain on the reciprocal (Barrett) FFT
reducer, whose products need transforms of size ≥ 2d; node traces as
divisions of the root's traces (~d·n per node, ~d² per tree level); own
Bits−1-squaring traces for the nodes past the orbit depths; and a
Euclidean gcd plus a long division per split (~1.5·d² in total). All
items below are in `src/sketch_impl.h` unless noted.

### Montgomery squaring on an FFT coset

For the root f (monic, degree d), pick N = 2<sup>m</sup> ≥ d and a coset
C = W<sub>m</sub> + η of the additive-FFT subspace, with vanishing
polynomial Z<sub>C</sub>: in the Cantor basis
Z<sub>C</sub>(x) = Σ<sub>j⊆m</sub> x<sup>2<sup>j</sup></sup> + const,
degree N, formal derivative **Z<sub>C</sub>′ = 1**. Residues are held in
Montgomery form r·Z<sub>C</sub> mod f, as a polynomial of degree < N given
by its N values on C (`MontCoset`). The reduction of a product a is
q = a/f mod Z<sub>C</sub> (pointwise on C) and r = (a + q·f)/Z<sub>C</sub>.
On C both sides of h = a + q·f = r·Z<sub>C</sub> vanish. Differentiating,
h′ = r′Z<sub>C</sub> + r gives **r(w) = a′(w) + q′(w)f(w) + q(w)f′(w)** on C,
and for a square a = r² the term a′ = 2rr′ vanishes in characteristic 2.
So a modular squaring is pointwise work plus the values of q′: one
interpolation of q and one evaluation of its formal derivative, **two
size-N transforms**. The reciprocal reducer needs about eight of size 2N
for the same squaring. f must have no root on C (1/f there); the four
cosets of W<sub>m</sub> in W<sub>m+2</sub> are tried in turn, with the old
path as the fallback (exercised in the unit tests by planting roots on
all four).

The orbit (`MontgomeryOrbit`) runs the Bits squarings this way. It
accumulates, on the coset, only the root trace Σ basis<sup>2<sup>i</sup></sup>X<sub>i</sub>
(normalized once at the end) and the relative traces described below. The
exact factorizability test X<sub>Bits</sub> = x becomes the vanishing of
the reduced difference of two Montgomery forms. `AdditiveFFT` gained coset
transforms: the block index is a twiddle offset, since a coset of
W<sub>m</sub> is one block of a larger transform.

### Novel-basis transforms

The Cantor FFT is the Lin–Chung–Han butterfly transform with the
monomial→novel basis conversion (repeated division by
S<sup>j</sup>(x), S = x² + x) fused in as XOR cascades. Between the
interpolation and the evaluation of a Montgomery squaring nothing needs
monomial coefficients, and the formal derivative is cheap in the novel
basis X<sub>k</sub> = Π<sub>j∈bits(k)</sub> S<sup>j</sup>: since
(S<sup>j</sup>)′ = 1, X<sub>k</sub>′ = Σ<sub>j∈bits(k)</sub> X<sub>k−2<sup>j</sup></sub>,
which is one XOR pass. `NovelFFT`/`NovelIFFT` are plain butterflies, 1.55× faster
than the fused transform (2.2 vs 3.6 ns per butterfly on the devbox);
conversions (`ToNovel`/`FromNovel`) happen only at the edges.

### The Montgomery remainder tree

The tree's nodes are factors of f, and a Montgomery form mod f is also one
mod any factor with the same Z<sub>C</sub>. Shrinking to the child's size
uses the halves C = C<sub>0</sub> ∪ C<sub>1</sub>:
Z<sub>C</sub> = Z<sub>C0</sub>·Z<sub>C1</sub> and Z<sub>C0</sub> = 1 on
C<sub>1</sub>. So the child's form on C<sub>1</sub> is v + q·g with
q = v/g mod Z<sub>C0</sub>, which is **two half-size transforms per carried
residue** (`MontCarryToFactor`). The tree carries relative traces this way
instead of dividing the root's traces by every node.

**Tower split parameters, GF(16).** With z generating GF(2<sup>s</sup>)
and split parameters β<sub>sl+k</sub> = basis·ξ<sup>l</sup>·z<sup>k</sup>
(ξ the field generator), the split trace of depth sl + k at any node is
Tr<sub>GF(2<sup>s</sup>)/GF(2)</sub>(z<sup>k</sup>T<sub>l</sub>), with
T<sub>l</sub> = Tr<sub>GF(2<sup>Bits</sup>)/GF(2<sup>s</sup>)</sub>(basis·ξ<sup>l</sup>x)
mod g: **s − 1 squarings** of a carried trace. The orbit accumulates
T<sub>l</sub> for free (every s-th step). Carrying one trace to both
children costs two full-size transforms while one squaring costs two per
node, so a smaller subfield, which needs more carried traces but fewer squarings, wins.
s = 4 (`MONT_TOWER_BITS`, up to four traces for 16 depths) beats s = 8 by
5–8% on the i7 and the Pi, and s = 2 is in between. Carrying traces down to
the deep small nodes removes their 63-squaring own traces (`SmallTrace`),
which were 4–10% of decode at high fill. Only the traces for the
log₂N + 2 levels a degree-d tree reaches with nodes of degree ≥ 3 are
carried (`CARRY_DEPTH_MARGIN`): three up to N = 2048. That is 0–6% faster
than always carrying four, and the rare deeper node computes its own
trace.

**Small nodes.** Below degree 16 (`SMALL_CARRY_MAX`) the traces are
carried in normal form (division by the factor, coefficient-domain
squarings: `TowerTraceSmall`), avoiding the per-node coset setup. Below
degree 512 (`MONT_NORMALIZE_MIN`) a node's split trace is left in
Montgomery form, a unit multiple mod g: the gcd with g, its only use, is
the same, and the gcd's first division is cheaper than the reduction at
those sizes.

**Division on the coset.** After a split g = g<sub>0</sub>·g<sub>1</sub>,
g<sub>1</sub>'s values are g/g<sub>0</sub> pointwise on the node's coset
(`CosetDivide`). One evaluation, a batch inversion and one interpolation
replace the long division. The quotient must come out monic of the right
degree, which with deg g ≤ N proves divisibility. Both factors' values are
also the halves their carries start from.

### Power-sum splits: BM instead of the gcd

The decoder passes the sketch's syndromes, the power sums
Q<sub>k</sub> = Σ α<sup>k</sup> of the roots, to `FindRoots`. For a node g
with split trace t (0 or 1 at every root), the roots A with t = 1 have
power sums
P<sub>k</sub> = Σ<sub>α∈A</sub> α<sup>k</sup> = Σ<sub>j<n</sub> t<sub>j</sub>Q<sub>j+k</sub>,
a Hankel product (`PowerSumSplit`: for odd k, two half-size
transposed-Karatsuba middle products split by the parity of j). Since
t² = t on the roots, P<sub>2k</sub> = P<sub>k</sub>². These are syndromes,
and `BerlekampMasseyOdd` recovers Π<sub>A</sub>(1 − αy) from the odd ones.
An |A| ≈ n/2-root decode costs roughly half of the Euclidean gcd from
degree n down to n/2 that it replaces: ≈0.45·n² against ≈1.09·n² ns
measured at n = 1536 on the devbox. The other factor follows by coset division, whose
exactness check rejects a wrong BM result (too few power sums, or power
sums that do not belong to g). With 2M ≥ n power sums one of A and its
complement has at most n/2 roots, so the complement (power sums
Q<sub>k</sub> − P<sub>k</sub>) is tried next, then the gcd. Children
inherit P and Q − P, so no node recomputes power sums. Nodes of degree
≥ 256 in Montgomery form use it.

### Field arithmetic (`src/fields/clmul_common_impl.h`)

* CLMUL fields above 32 bits square by one carry-less product instead of
  the linear-map tables (64-bit: 2.0 vs 3.4 ns x86-64, 3.7 vs 6.1 ns
  Cortex-A76). At 32 bits they are even and the tables stay.
* On ARM the 64-bit reduction folds the high word with shifts in
  general registers, avoiding two dependent PMULLs with their NEON
  transfers: Mul 4.0 → 3.0 ns, latency 10.0 → 6.4 ns on the A76, 4–10%
  on decode. x86-64 keeps PCLMULQDQ (the shift form is slower there,
  2× on Haswell).

### Measured

Pinned (`tools/bench_host.sh`: performance governor, turbo off, min = max:
i7-4770 at 3.4 GHz, Pi 5 at 2.4 GHz), interleaved, min of 7 rounds,
CLMUL/PMULL implementation unless noted; starting point: the 2026-10
"Frobenius-structured BM, shared trace orbit, lazy reduction" series.

| capacity/differences | i7-4770 | Pi 5 (A76) |
|---|---|---|
| 1024/512 (50%) | 9.08 → 4.14 ms (2.19×) | 10.05 → 5.84 ms (1.72×) |
| 1024/768 (75%) | 19.34 → 8.01 ms (2.42×) | 28.34 → 11.19 ms (2.53×) |
| 1024/1024 (100%) | 23.94 → 10.12 ms (2.37×) | 34.43 → 14.02 ms (2.46×) |
| 2048/1024 (50%) | 24.66 → 11.00 ms (2.24×) | 35.45 → 15.09 ms (2.35×) |
| 2048/1536 (75%) | 50.85 → 21.44 ms (2.37×) | 72.49 → 30.19 ms (2.40×) |
| 2048/2048 (100%) | 67.65 → 28.15 ms (2.40×) | 94.35 → 40.50 ms (2.33×) |
| 4096/2048 (50%) | 69.35 → 30.46 ms (2.28×) | 96.75 → 42.58 ms (2.27×) |
| 4096/3072 (75%) | 144.1 → 59.6 ms (2.42×) | 202.2 → 83.3 ms (2.43×) |
| 4096/4096 (100%) | 203.5 → 87.1 ms (2.34×) | 279.3 → 123.2 ms (2.27×) |
| 32-bit 1024/768 | 12.62 → 8.23 ms (1.53×) | 13.06 → 7.52 ms (1.74×) |
| 32-bit 2048/1536 | 36.03 → 23.48 ms (1.53×) | 36.50 → 20.85 ms (1.75×) |
| 32-bit 2048/2048 | 52.25 → 33.10 ms (1.58×) | 51.95 → 28.79 ms (1.80×) |
| 1024/256 (25%) | 2.50 → 1.61 ms (1.55×) | 2.79 → 2.19 ms (1.27×) |
| 1024/102, 2048/204 (10%) | 1.01× | 1.10–1.12× |
| 1024/32, 1024/64, 4096/64 | 1.03× | 1.07–1.13× |
| 32-bit 1024/64 | 0.99× | 0.98× |
| generic 1024/768 (64-, 32-bit) | 1.01–1.02× | 1.05× |
| 1024/1100, 4096/4200 (over capacity: fail) | 3.22×, 2.25× | 3.39×, 2.35× |

The 32-bit field gains less: its chain is half as long, so the quadratic
stages (BM and the splits' BMs, the small nodes' gcds), which this series
changes least, are a larger share. Raw data:
`doc/benchmark_results/high_fill/`.

Peak heap of one decode (instrumented `operator new`, 64-bit / 32-bit
field) is lower everywhere: the d²/2 square table of the old orbit is
gone, and the coset vectors are O(N) per tree level.

| capacity/differences | 64-bit | 32-bit |
|---|---|---|
| 1024/256 | 340 → 88 KB | 185 → 51 KB |
| 1024/512 | 1183 → 180 KB | 622 → 97 KB |
| 1024/768 | 358 → 263 KB | 225 → 148 KB |
| 1024/1024 | 446 → 325 KB | 284 → 189 KB |
| 2048/1536 | 710 → 573 KB | 446 → 335 KB |
| 4096/3072 | 1414 → 1130 KB | 888 → 669 KB |
| 4096/4096 | 1766 → 1458 KB | 1124 → 729 KB |

Allocation counts per decode roughly doubled (per-node coset vectors);
reusing the carry's buffers measured no difference on x86-64.

### Measured negative results in this series

* **Two-coset Montgomery with 4 transforms.** Representing residues on two
  cosets A and B, reducing pointwise on B and moving q (A→B) and r (B→A),
  costs four size-N transforms per squaring. The derivative trick needs two.
* **Normalizing all eight orbit buckets.** Each costs ~5 transforms; with
  the carry only the root trace needs normal form (−20% orbit).
* **Phased XOR cascades.** Every offset z of the division by
  S<sup>q</sup> satisfies h/2 ≤ z < h. So the cascades out of a block's
  top quarter never touch it, and two dependence-free passes replace the
  sequential top-down division. Neutral end to end: the conversion is
  bound by its ~13 XORs per element at N = 2048, not by the dependences.
* **Half-gcd with shared transforms.** The 2×2 matrix-vector and matrix
  products now transform each distinct operand once (8 instead of 12, 12
  instead of 24 transforms), which made HGCD ~10% faster than the
  quadratic loop at 3072 but still a tie at 1024–1536. The cutoff stays
  at 4096.
* **32-bit reduction in general registers.** On x86-64, +1% at high fill
  and +5–7% at low fill on Haswell, but −8% on Zen 3. On the A76, −17%.
  Not adopted.
* **Blocked Berlekamp–Massey / Lehmer Euclid (analysis).** Simulating k
  steps on short discrepancy windows and applying the 2×2 matrix with
  Karatsuba products costs about the classical ~4 carry-less products per
  coefficient-step at block sizes where the block simulation stays cheap.
  It pays only in a full divide-and-conquer form with fast products, where
  the half-gcd's constants apply.
* **Fused Berlekamp–Massey update.** Accumulating the next odd step's
  discrepancy in the same pass that writes the updated polynomial (one
  pass over it per step instead of two) measured no gain on x86-64
  (−2% for the 32-bit field). The loop is bound by its carry-less
  products, not its memory traffic.
* **Peeling splits (analysis).** The gcd work telescopes:
  Σ(n<sub>parent</sub>² − n<sub>gcd</sub>²) over the gcds computed. So
  extracting large gcds (small complements by coset division) could cut it
  by up to a third. But the needed indicator polynomials cost more
  Montgomery products than that saves.

### Remaining in this area

* **Truncated (union-of-cosets) Montgomery transforms.** N is the next
  power of two above the node degree: at 75% fill (768, 1536, ...) a
  quarter of every transform is padding. A coset union
  C<sub>1</sub> ∪ C<sub>2</sub> (1024 + 512) works with
  Z = Z<sub>1</sub>Z<sub>2</sub>, Z′ = Z<sub>1</sub> + Z<sub>2</sub> known
  pointwise, at the cost of truncated transforms.
* **Divide-and-conquer Berlekamp–Massey.** BM is now the largest
  quadratic stage at high fill: the decode BM plus the power-sum splits'
  BMs, ~35% at 2048/1536.
* **Allocation pooling** for the per-node coset vectors (2–4% estimated,
  more on the Pi).

## Demoted ideas (investigated, not worth it)

* **Frobenius-power precomputation for trace retries** — cache
  `T_i = x^(2^i) mod poly` per node so retried splits cost `m·d` instead of
  `m·d²`. Sounds attractive, but split retries have probability
  ~2<sup>1−d</sup> at degree d, so above degree ~8 there is nothing to
  amortize, and at the retry-prone degrees (3, 4) the square-table reducer
  already makes each retry nearly free. Estimated ≤10% in the small-n
  regime, and only before the reducer existed.
* **Toom-3 for the Karatsuba tier** — evaluation points and interpolation
  divisions exist in GF(2<sup>m</sup>) for m ≥ 3, but the extra
  interpolation cost makes it marginal over Karatsuba in practice; an
  additive FFT is the better use of complexity budget.
* **Explicit quintic+ solvers** — no closed form exists for degree ≥ 5
  (the affine trick stops at linearized-polynomial degrees 2<sup>k</sup>;
  degree-8 affine multiples exist but the general reduction to them is no
  longer cheap). The Fedorenko-Trifonov affine-decomposition family targets
  full-domain evaluation (Chien-style, feasible only for small fields), not
  the large-m fields minisketch cares about.
