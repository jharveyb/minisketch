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

1. **Syndrome reconstruction** — O(c) squarings.
2. **Berlekamp-Massey** — adaptive O(c·n) multiplications: the inner loops
   only touch the current LFSR length, so the cost tracks the *actual* number
   of differences, not the capacity.
3. **Root finding** (`RecFindRoots`, the Berlekamp Trace Algorithm) — the
   dominant stage. Each node of the recursion at degree `d` computes a
   symbolic trace: `Bits-1` iterations of squaring + reduction modulo the
   node's polynomial. With schoolbook reduction that is ~`m·d²`
   multiplications per node; the per-node GCD and division are a further
   ~`d²`, i.e. `1/m` of the trace cost. Summing over the (roughly binary)
   recursion tree gives ~`2·m·n²` total multiplications, **dominated by the
   top 2–3 recursion levels**: level `k` costs ~`m·n²/2^k`.

Two consequences shape everything below:

* Speedups to the *top-of-tree* reduction cost translate almost 1:1 into
  total decode time for large `n`.
* Speedups to *small-degree nodes* only matter when `n` itself is small
  (shallow tree) — but that is exactly the low-latency regime that set
  reconciliation protocols like Erlay live in.

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
capacity: Berlekamp-Massey is adaptive O(e·n), and the locator degree —
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

## Next: subquadratic GCD (roadmap)

With the FFT tier and cached transforms in place, the quadratic GCD is the
next structural bottleneck: `GCD(trace, poly)` runs once per split attempt
at every `RecFindRoots` node (~d² multiplications at degree d, ~2d² summed
over the recursion — at capacity 4096 already comparable to the whole
FFT-accelerated trace stage), plus one quadratic `DivMod` per successful
split. Berlekamp-Massey stays adaptive O(e·n) and only matters when e ≈ n.

Algorithm choice (literature review, mid-2026):

* **Classic recursive HGCD, gcd-only — not jumpdivstep.** The
  root-splitting step needs only the gcd itself, not Bezout coefficients,
  which drops the cofactor reconstruction. van der Hoeven, *Optimizing the
  half-gcd algorithm* (arXiv:2212.12389) gives ~5.5·M(d)·log₂d in the FFT
  model for *normal* remainder sequences — the expected case over
  GF(2<sup>m</sup>), where sequences are normal with probability
  1 − O(d/2<sup>m</sup>) — via middle products and transform caching
  across the 2×2-matrix recursion. Both tricks fit the additive FFT; the
  `Prepare`/`MulLowPrepared` entry point is the needed hook. (An earlier
  note here claimed the paper's FFT-model constants were "not reachable in
  characteristic 2" — true before the additive-FFT tier existed, obsolete
  now.)
* Polynomial **divstep** (safegcd §3–7, eprint 2019/266) processes exactly
  2d−1 steps, cannot exploit degree drops, and lands ~2× behind optimized
  HGCD in the same multiplication model; its constant-time regularity is
  irrelevant to decode (the decoder is variable-time by design, and BM is
  a Euclid special case per the paper — see below). The one measured
  jumpdivstep implementation (Jian-Wang-Yang-Chen, eprint 2024/644, NTRU
  Prime inversion) crosses over against *quadratic divstep* at degree
  ~653 and confirms transform caching as the main practical lever. Keep
  divstep in reserve only for a hypothetical constant-time decoder (no
  published divstep-over-GF(2<sup>m</sup>)[x] implementation appears to
  exist), along with libsecp256k1's batched transition-matrix trick
  (`doc/safegcd_implementation.md`) as a memory-traffic optimization.
* **Crossover evidence**: NTL's `GF2EX` (the same coefficient-ring shape)
  switches Euclid→HGCD around degree 40 with recursion base 40; FLINT's
  `nmod_poly` at ~340 over cheaper word-size prime fields. Expect ours in
  the 100–500 range; sweep it.
* **Honest ceiling**: at d = 4096, ~5.5·M(d)·log₂d ≈ 11M multiplications
  (M(4096) ≈ 168k with the FFT tier) vs ~16.7M quadratic — ~1.5× on the
  gcd itself by the naive constant, growing with capacity and with the
  constant-factor tricks; the structural point is that decode stops being
  quadratic overall.
* **Fast BM** via the Dornstetter/key-equation equivalence (one partial
  XGCD at degree 2t; modern formulation: PM-Basis minimal approximants,
  practical in Neiger's PML library at degrees 10³–10⁵) only wins when
  e ≈ n, since the current BM is adaptive O(e·n): the pragmatic hybrid
  keeps adaptive BM and restarts with the HGCD key-equation solver when
  the LFSR length crosses a tuned threshold. Strictly after the gcd work.

Implementation order (branch `subquadratic-gcd`):

1. Subquadratic `DivMod` after each successful split, using the existing
   reciprocal machinery (`TraceModInvSeries` + two `MulLow`s) — the same
   algorithm `ReciprocalReduce` already uses, as a one-shot division with
   quotient output. Small, low-risk, immediately measurable.
2. Recursive gcd-only HGCD over `TraceModPolyMulFull` with a swept
   degree-cutoff fallback to the quadratic loop; van der Hoeven
   middle-product and transform-caching constants. Tests per `testing.md`:
   property tests against the quadratic GCD reference, a poly_ops fuzz
   mirror, planted mutations. Acceptance: no regression at any capacity,
   measured win at 2048+.
3. Hybrid fast-BM key-equation solver, gated on measured e ≈ n profiles.

References: <https://arxiv.org/abs/2212.12389>,
<https://eprint.iacr.org/2019/266>, <https://eprint.iacr.org/2024/644>,
NTL `GF2EX` HalfGCD (<https://github.com/libntl/ntl>),
<https://github.com/vneiger/pml>,
<https://github.com/bitcoin-core/secp256k1/blob/master/doc/safegcd_implementation.md>,
<https://github.com/sipa/safegcd-bounds>.

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
