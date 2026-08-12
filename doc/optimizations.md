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
indexing was observed to inflate short rows by ~2.5×.

Remaining follow-ups on this tier:

* **Cached node transforms ("Phase 2")**: `inv` and `mod` are fixed per
  TraceMod node, so their transforms can be computed once and reused across
  all `Bits−1` reductions, cutting each reduction from 6 transforms to 4
  (bounded ~33% of the FFT time; needs a pretransformed-operand entry point
  and a fixed per-node transform size). Additionally, in the evaluation
  domain squaring is pointwise (`f²(p) = f(p)²` in char 2), which a deeper
  restructuring of the trace loop could exploit.
* **General-basis Gao-Mateer** for non-power-of-two fields.
* The **Frobenius additive FFT** (Li et al., ISSAC 2018) does *not* apply:
  its factor-d truncation requires coefficients in GF(2), while the
  decoder's are full field elements.

References:

* Badakhshan, Samanta, Gong, *Accelerating Post-quantum Secure zkSNARKs by
  Optimizing Additive FFT* (SAC 2025;
  `doc/sac2025-2-paper17_optimizing_additive_fft.pdf`) — primary source for
  the implemented algorithm.
* Li, Chen, Kuo, Cheng, Yang, *Frobenius Additive Fast Fourier Transform*
  (ISSAC 2018; `doc/3208976.3208998_frobenius_additive_fft.pdf`,
  <https://arxiv.org/abs/1802.03932>) — Cantor-construction cross-check.
* Gao, Mateer, *Additive Fast Fourier Transforms over Finite Fields*
  (<https://www.math.clemson.edu/~sgao/papers/GM10.pdf>) — the general-basis
  construction for the non-power-of-two follow-up.
* Bernstein, *Multiplication of polynomials over F₂* survey
  (<https://cr.yp.to/f2mult.html>)

## Deferred: safegcd (Bernstein-Yang divstep) for Berlekamp-Massey

The safegcd paper (Bernstein-Yang 2019, eprint 2019/266) is half about
polynomials: its §3–7 define the `divstep` iteration for `k[[x]]` over an
arbitrary field k — GF(2<sup>m</sup>) applies directly — and the paper itself
notes that **Berlekamp-Massey is a special case of a Euclid/half-gcd
computation** with reversed coefficients (with one input a power of x). So a
divstep-based BM is well-defined and would slot into `BerlekampMassey`'s
place, with the error locator falling out of the accumulated 2×2 transition
matrix.

Why it is deferred — the arithmetic at minisketch sizes goes the wrong way:

* Plain `divstepsx` is Θ(steps × length) = Θ(c²) for capacity c: it always
  processes full-length inputs. The current BM is adaptive Θ(c·n) (n =
  actual LFSR length), which is *better* whenever n ≪ c — the common case.
  A straight port would be a 2–8× regression.
* The subquadratic variant `jumpdivstepsx` (§5.3) is elegant — unlike
  classical HGCD it needs **no polynomial division subroutine at all**, just
  2×2 matrix products of half precision — but its constants are HGCD-class.
  With Karatsuba-only `M(d)` its crossover against schoolbook sat around
  degree 10⁴; with the implemented additive-FFT `M(d)` the break-even
  estimate drops to a few thousand — at the edge of real sketch sizes.
  Relevant because the FFT changed the cost balance: at capacity 4096 the
  per-node quadratic GCD is now comparable to the whole FFT-accelerated
  trace computation of the top node, so subquadratic GCD/BM is the next
  structural bottleneck for large-capacity decode.

What *is* worth taking from safegcd today:

* The **batched transition matrix** trick from libsecp256k1's implementation
  (`doc/safegcd_implementation.md`): N divsteps depend only on the bottom N
  coefficients, so one can compute an N-step 2×2 matrix from a prefix and
  apply it to the full-length polynomials in one pass — a constant-factor
  memory-traffic optimization that does not need subquadratic multiplication.
* Divstep's data-independent control flow is the natural shape for a
  **constant-time decoder**, should side-channel-uniform decoding ever become
  a requirement (the current decoder is variable-time and instead randomizes
  the trace basis per sketch). No published divstep-over-GF(2<sup>m</sup>)[x]
  implementation appears to exist, so this would be new work.

References: <https://eprint.iacr.org/2019/266>,
<https://github.com/bitcoin-core/secp256k1/blob/master/doc/safegcd_implementation.md>,
iteration bounds <https://github.com/sipa/safegcd-bounds>.

## Deferred: Half-GCD

The README's remaining TODO. Classical Knuth-Schönhage-Moenck HGCD computes
polynomial GCDs in `O(M(d) log d)`; the best published constants (van der
Hoeven, *Optimizing the half-gcd algorithm*, 2022) are ~10–22·`M(d)`·log₂d
depending on normality of the remainder sequence, and the favorable "FFT
model" constants in that paper are not reachable in characteristic 2. Two
independent reasons to defer:

* With Karatsuba `M(d)`, break-even against the current quadratic GCD is
  around degree 3×10⁴ — an order of magnitude above the largest practical
  sketch capacities.
* Before the additive FFT, GCD was not even the bottleneck: ~1/m of a node's
  trace computation. The FFT changed this — at capacity 4096 the quadratic
  GCD is now comparable to the top node's whole trace stage — so a
  subquadratic GCD is the next structural lever for large capacities.

With the additive FFT in place, `jumpdivstepsx` above is likely the
better-shaped candidate for this role (no division subroutine; uniform
structure); the break-even remains to be measured, estimated at a few
thousand syndromes.

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
