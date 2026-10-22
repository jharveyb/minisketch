# 7. Tuning parameters

All algorithm crossovers are compile-time macros at the top of
`src/sketch_impl.h` (lines 19–142), each overridable with
`-DMINISKETCH_<NAME>=<value>`. None of them changes a decode result (C45:
they are measurements, and every tier is exact); they only move work
between equivalent algorithms. A port can start from these defaults.

## Field kinds

Two properties of the field select tiers:

- **FFT-capable**: b ∈ {16, 32, 64} (power of two, ≥ 16), any
  implementation.
- **Cheap multiplier** (`TraceModCheapScalarMul`): the CLMUL
  implementations 1 and 2, whose `Multiplier` is the bare element; the
  generic implementation 0 builds a lookup table per scalar. CLMUL fields
  are also the lazily-reducing ones (`LAZY_REDUCE`).

| tier | generic, any b | generic, b ∈ {16,32,64} | CLMUL, other b | CLMUL, b ∈ {16,32,64} |
|---|---|---|---|---|
| additive-FFT products (`fft_eligible`) | - | yes | - | yes |
| square-table / reciprocal reducer cutoff | 512 | 512 | 512 | 704 |
| cached reducer transforms | - | yes | - | yes |
| FastDivMod | - | - | - | yes |
| lazy dot products, inversion-free Euclid | - | - | yes | yes |
| BM tail: direct / Karatsuba / FFT | table direct | table < 40 ≤ FFT | direct < 64 ≤ Karatsuba < 512 ≤ Karatsuba blocks | direct < 64 ≤ Karatsuba < 512 ≤ FFT |
| shared orbit, tower betas (8 ∣ b) | yes | yes | yes | yes |
| Montgomery orbit, carries, coset division, power-sum splits | - | - | - | implementation 1 only (impl 2 has no b ∈ {16,32,64}) |

## Cutoff macros

"Hardware" reports what the code comment (or, where marked, the section of
`doc/optimizations.md` that introduced the value) says it was measured on.

| macro | default | selects | trade-off | hardware |
|---|---|---|---|---|
| `TRACEMOD_POLYMUL_CUTOFF` | 24 | Karatsuba recursion bottoms out in schoolbook when n ≤ 24 or min operand length ≤ 24; FFT products also require operands longer than 24 | recursion overhead vs d^2 | not stated |
| `TRACEMOD_TABLE_CUTOFF` | 512 | TraceMod: square table for node degree < 512, reciprocal (Barrett) reducer above; for generic and non-FFT-capable fields | table d^2/2 elements, streamed, vs Barrett products | optimizations.md: generic fields regressed 10–25% with a lower value |
| `TRACEMOD_FFT_TABLE_CUTOFF` | 704 | the same cutoff for FFT-capable cheap fields | the lazy streaming table beats the FFT reducer up to here | "~700 on Haswell, above 768 on Cortex-A76" |
| `TRACEMOD_FFT_CUTOFF` | 128 | minimum product length for the additive-FFT tier; TraceMod caches node transforms when 2d − 1 ≥ 128 | transform setup vs Karatsuba | optimizations.md: "nearly flat in 64..512" |
| `FAST_DIVMOD_CUTOFF` | 1024 | FastDivMod uses the reciprocal method when min(quotient length, divisor degree) > 1024 (FFT-capable cheap fields only) | Newton inverse + 2 FFT products vs lazy schoolbook | optimizations.md (pre-lazy measurements) |
| `GCD_INVFREE_MAX` | 128 | lazy Euclid uses the inversion-free pseudo-remainder for normal steps with divisor degree < 128 | one extra product per coefficient vs one inversion | not stated |
| `DIVMOD_LAZY_MIN` | 4 | DivMod and GCD use lazy dot-product division when the divisor has more than 4 coefficients (CLMUL) | | not stated |
| `SMALL_TRACE_MAX` | 8 | nodes of degree ≤ 8 on fully factorizable paths use SmallTrace | TraceMod setup vs d^2/2 products | not stated |
| `BM_VERIFY_GAP` | 2 | a tail check is attempted once done ≥ 2L + 2 | earlier attempts vs wasted checks | not stated |
| `BM_VERIFY_MIN_PREFIX` | 0 | earliest done count for a tail check | | not stated |
| `BM_VERIFY_MIN_TAIL` | 32 | a tail check needs at least 32 unprocessed syndromes | check overhead vs plain BM steps | not stated |
| `BM_VERIFY_FFT_MIN_DEGREE` | 512 | FFT tail check from L = 512 on cheap multipliers | FFT vs lazy Karatsuba | "measured on Haswell and Cortex-A76" |
| `BM_VERIFY_FFT_MIN_DEGREE_TABLE` | 40 | FFT tail check from L = 40 on generic FFT-capable fields | slow pointwise products vs table direct | same comment |
| `BM_VERIFY_KARA_MIN_DEGREE` | 64 | lazy Karatsuba tail check for 64 ≤ L < 512 (CLMUL) | ~n^0.585 vs n products per position | "from ~32 on Haswell, ~96 on Cortex-A76" |
| `MP_BASE` | 16 | middle products of size ≤ 16 (or odd) are direct | | not stated |
| `ORBIT_MIN_DEGREE` | 8 | shared Frobenius orbit for root degree ≥ 8 (else the small path) | orbit setup vs per-node traces | optimizations.md: from 16 to 8 gave +43%/+51% at 1024/12 (i7/Pi 5) |
| `ORBIT_EXTRA_DEPTHS` | 2 | orbit traces for min(b, 2 + floor(log2 b)) depths | accumulators vs nodes' own traces | optimizations.md: retune flat ±3% |
| `ORBIT_NODE_FACTOR` | 2 | a node uses the orbit trace (one division) when 2·deg f ≤ deg g · b | division by a long trace vs own b − 1 squarings | same |
| `MONT_ORBIT_MIN_DEGREE` | 224 | Montgomery orbit from root degree 224 (FFT-capable cheap fields, 8 ∣ b, 4 ∣ b) | 2 size-N transforms per squaring vs table/Barrett | optimizations.md (i7-4770, Pi 5) |
| `MONT_TOWER_BITS` | 4 | subfield GF(2^s) of the Montgomery betas; s − 1 squarings per carried trace; 16/s carried traces (s ∈ {1, 2, 4, 8}, s ∣ b) | squarings per node vs carried traces | optimizations.md: s = 4 beats s = 8 by 5–8% on the i7 and the Pi |
| `CARRY_DEPTH_MARGIN` | 2 | carried traces: min(4, ⌈(m + 2)/4⌉) for root coset size 2^m | carrying a trace vs deep nodes' own traces | optimizations.md: 0–6% vs always four |
| `SMALL_CARRY_MAX` | 16 | factors of degree 3..16 get normal-form carries (division + schoolbook squarings) | coset setup vs schoolbook | not stated |
| `COSET_DIV_MIN` | 64 | Montgomery nodes of degree ≥ 64 divide out the split factor on the coset | interpolation vs long division | not stated |
| `PSUM_SPLIT_MIN` | 256 | power-sum splits for Montgomery nodes of degree ≥ 256 (and root psum only from degree 256) | half-size BM vs Euclid | optimizations.md: ≈0.45 n^2 vs ≈1.09 n^2 ns at n = 1536 ("devbox") |
| `MONT_NORMALIZE_MIN` | 512 | Montgomery split traces are normalized (ToNormal) from node degree 512; below, the unit multiple goes to the gcd | one reduction vs a longer first gcd division | not stated |

**Tuning hardware.** The defaults were tuned on an i7-4770 (Haswell), a
Raspberry Pi 5 (Cortex-A76) and a Zen 5 machine; the measurements cited in
the macro comments and `doc/optimizations.md` are from the first two.

## Internal constants (not macros)

| constant | value | where |
|---|---|---|
| `CARRY_DEPTHS` | 16 | depths for which relative traces can be carried |
| `CARRY_TRACES` | 16 / MONT_TOWER_BITS = 4 | number of carried relative traces |
| Montgomery coset order | blocks 1, 2, 3, 0 | MontgomeryOrbit |
| x·Z_C mod f | `FastDivMod` with cutoff 32 | MontgomeryOrbit |
| power-sum budget | 2M = min(\|q\| + 1 − n, n + 4⌈√n⌉ + 32), even | PowerSumSplit |
| root power sums | min(2c, 2d + 4⌈√d⌉ + 64) | FindRoots |
| PowerSumSplit block | h = ⌈n/2⌉ rounded up to even | PowerSumSplit |
| tail check, table tier block | 512 positions | BerlekampMasseyVerifyTail |
| tail check, FFT size | least 2^m ≥ 4(L + 1), block 2^m − 2L | BerlekampMasseyVerifyTail |
| tail check, Karatsuba-block size | max(2L, 64) | BerlekampMasseyVerifyTail |
| tail check retry | next attempt at done ≥ 1.5 × done | BerlekampMasseyOdd |
| `TW_MUL_CACHE` | 64 twiddle Multiplier tables (generic fields) | AdditiveFFT |
| `TraceModScratch::SLOTS` | 9 scratch vectors per recursion level | Karatsuba |
| add batch size | 1, 2, 4 (default) or 8 | AddBatch |

## Retuning

`CLAUDE.md` describes the procedure: pin the CPU (`tools/bench_host.sh`),
then `CXXFLAGS_B=-DMINISKETCH_<NAME>=<v> tools/bench_ab.sh . .` for an
interleaved A/B; unpinned numbers are not comparable.
`-DMINISKETCH_STAGE_TIMING` prints per-stage times, with timer overhead
that inflates tiny per-node calls. `doc/optimizations.md` records the
measured crossovers and the negative results; check it before retuning.
