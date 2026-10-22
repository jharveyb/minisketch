# minisketch implementation specification

This directory specifies libminisketch (PinSketch set-reconciliation sketches
over GF(2^b), 2 ≤ b ≤ 64) in enough detail to write an independent
implementation in another language, such as Rust or Go, that:

1. is **byte-compatible**: identical serializations, identical element
   encodings, identical capacity arithmetic;
2. **decodes exactly the same sets** and fails on exactly the same inputs as
   libminisketch, for every input, including adversarial and overfull
   sketches;
3. **performs comparably**, because it uses the same algorithms, tiers and
   cutoffs.

The spec describes the code at the tip of the `verification` branch (line
numbers as of commit `e42e84fc`). Where it cites code, it uses
`src/file.h:line` once per algorithm. Claim IDs such as
(C30) refer to [claims.md](claims.md), the registry of every mathematical
statement the decoder relies on, with its proof status.

## Normative and informative content

- **Normative** (a conforming implementation must match it): the field
  element encoding (chapter 1, "Element representation"), sketch contents,
  add, merge, serialization and deserialization (chapter 2), the capacity
  functions (chapter 2), and the decode *result*: the set returned, or
  failure (chapter 6, "Conformance").
- **Informative** (matching it gives the same performance profile, but any
  algorithm with the same result conforms): every algorithm in chapters 3–5
  and 7, the internal representation used for arithmetic (including the
  internal modulus of the 58-bit CLMUL field), the split parameters, the
  evaluation points of the FFT, and the output order of decoded elements.

## Conformance requirement

For every field size b, capacity c, input sketch (any c·b bits, not only
sketches of real sets) and `max_elements`, decode must return:

- the same **set** as the textbook decoder described in chapter 6
  (Berlekamp–Massey over all 2c syndromes, followed by Berlekamp's trace
  root finding), or
- **failure** exactly when the textbook decoder fails.

The **order** of the returned elements is unspecified. In libminisketch it
depends on the per-sketch seed, the implementation number (generic or
CLMUL), the root degree (which picks the split-parameter scheme) and the
field's internal representation; see chapter 6, "Output order".

Serialization is byte-exact: chapter 2 gives the bit layout.

The decoder's claim to equal the textbook decoder is (C46). It is the
conjunction of the claims cited in chapters 4–6, most of them proved in Lean
(see [claims.md](claims.md) for which).

## Chapters

| Chapter | Contents |
|---|---|
| [01-field.md](01-field.md) | GF(2^b): encoding, moduli, multiplication, CLMUL reduction, squaring, Qrt, inversion, lazy reduction |
| [02-sketch.md](02-sketch.md) | Sketch contents, add, merge, serialization, seeds, implementations, capacity functions |
| [03-polyarith.md](03-polyarith.md) | Polynomial arithmetic: Karatsuba, middle products, additive FFT, Newton inverse, division, gcd, half-gcd |
| [04-bm.md](04-bm.md) | `BerlekampMasseyOdd`: odd steps only, early termination, exact tail check |
| [05-roots.md](05-roots.md) | `FindRoots` / `RecFindRoots`: orbit, tower bases, coset Montgomery tree, power-sum splits, every fallback |
| [06-decode.md](06-decode.md) | The decode pipeline, call graph, failure conditions, output order |
| [07-tuning.md](07-tuning.md) | Every `MINISKETCH_*` cutoff and internal constant |
| [08-porting.md](08-porting.md) | Notes for Rust and Go ports, what matters at which sizes, testing a port |

Read chapters 1, 2 and 6 for a minimal conforming implementation (textbook
decoder). Chapters 3–5 and 7 describe the fast decoder.

## Related material

- [claims.md](claims.md): every mathematical claim, its source, how it is
  checked, and its status (`proved`, `stated`, `tested`, `open`).
- `verification/`: the Lean 4 proofs referenced by claims.md
  (`verification/README.md`; `python3 verification/check_claims.py` checks
  the registry).
- `tests/pyminisketch.py`: a pure-Python implementation. Its first half is
  the textbook sketch and decoder, the normative reference for decode
  results. Its section "Optimized decode (oracle)" implements each fast
  technique of chapters 4–5 in a few lines (`bm_odd`, `middle_product_kara`,
  `frobenius_orbit`, `tower_betas`, `cantor_basis`, `novel_fft`,
  `mont_square`, `mont_carry_to_factor`, `coset_divide`,
  `montgomery_orbit`, `power_sum_split`, `find_roots_fast`, `decode_fast`)
  and is tested against the textbook functions. It is the best pseudocode
  companion to this spec.
- `tests/golden_vectors.json` and `tests/differential.py`: fixed
  serialization vectors and a differential test against the C library
  (chapter 8).
- Background: [../math.md](../math.md), [../math_decode.md](../math_decode.md),
  [../decode_explained.md](../decode_explained.md),
  [../decode_coding_theory.md](../decode_coding_theory.md),
  [../optimizations.md](../optimizations.md) (measurements and negative
  results), [../moduli.md](../moduli.md).

## Notation

K = GF(2^b). Polynomials are coefficient vectors, lowest degree first, with
no trailing zeros unless stated ("stripped"); the zero polynomial is the
empty vector. Addition is XOR. S_k is the k-th power sum (syndrome) of the
set; the sketch stores u_j = S_{2j+1} for j < c. Λ(x) = Π(1 − m x) is the
locator, f = rev(Λ) = Π(x − m) the root polynomial. Tr is the absolute
trace K → GF(2). "Lazy" means "accumulate unreduced carry-less products and
reduce once" (chapter 1).

## Discrepancies found

Where the code and a document disagree, this spec follows the code.

1. **Tuning hardware.** The header comment of the cutoff block in
   `src/sketch_impl.h:19` says "the defaults were tuned on Zen 5".
   `CLAUDE.md` says the cutoffs were tuned on an i7-4770 and a Pi 5
   (Cortex-A76), and the individual macro comments and
   `doc/optimizations.md` cite Haswell (i7-4770), Cortex-A76 (Pi 5) and a
   "devbox". No measurement in the repository is attributed to Zen 5 (the
   only AMD figure is one negative result on Zen 3, and `tools/bench_host.sh`
   mentions a mobile Zen 5 machine; the "devbox" may be that machine, which
   could not be confirmed). Chapter 7 lists what each comment says.
2. **Stale cutoffs in `doc/optimizations.md`.** The section "Implemented:
   the additive-FFT multiplication tier" gives
   `TRACEMOD_FFT_TABLE_CUTOFF = 256`; the code has 704 (the later "Lazy
   modular reduction" section records the 256 → 704 change). The section
   "Implemented: subquadratic division and GCD" gives `HGCD_CUTOFF = 1024`;
   the code has 4096 (marked as superseded in that section's "Remaining").
3. **Carried traces at N = 2048.** `doc/optimizations.md` ("The Montgomery
   remainder tree") says three relative traces are carried "up to
   N = 2048". The code (`src/sketch_impl.h:1673`) carries
   min(4, ⌈(m + 2)/4⌉) traces for N = 2^m: three for m ≤ 10 (N ≤ 1024),
   four for m = 11, 12, 13, 14 (N = 2048 … 16384).
4. **Tail-check tiers.** `doc/optimizations.md` and
   `doc/decode_coding_theory.md` describe the BM tail check as "FFT
   (L ≥ 512)" on CLMUL fields and "FFT from L = 40" on generic fields. In
   the code the FFT tier exists only for b ∈ {16, 32, 64}. On CLMUL fields
   of other sizes, L ≥ 512 uses blockwise Karatsuba full products; on
   generic fields of other sizes every L uses the direct Multiplier-table
   check (chapter 4).
5. **Misplaced comment.** In `src/sketch_impl.h:115` the comment "Node
   degree from which Montgomery-form split traces are normalized before the
   gcd (MontTowerTrace)" stands above `MINISKETCH_COSET_DIV_MIN`; it
   describes `MINISKETCH_MONT_NORMALIZE_MIN`, defined further down.
6. **`minisketch_add_uint64s` documentation.** `include/minisketch.h` says
   the function "processes elements in groups of four"; the group size is
   configurable (1, 2, 4 or 8, default 4) with `minisketch_set_batch_size`,
   as the next comment in the same header says. The result is the same for
   every batch size.
7. **Oracle simplifications** (not errors; noted so that a reader comparing
   them is not surprised). `tests/pyminisketch.py`'s `bm_odd` uses a
   different early-termination schedule from the C++ (it verifies at
   `processed = n + 2 ≥ 2L` and retries at `processed + processed//2 + 2`,
   with no `BM_VERIFY_GAP` or `BM_VERIFY_MIN_TAIL`), and `find_roots_fast`
   keeps the GF(2^8) tower split parameters on the Montgomery path where the
   C++ switches to the GF(16) tower. Both return the same results (C48, C25).

## Not determined

- Whether the degree-2 case `x^2 + a` (a repeated root) can be reached in a
  decode. The code returns failure there (`CHECK_RETURN`, which aborts in
  `MINISKETCH_VERIFY` builds, so the authors believe it unreachable); a port
  should return failure.
- Whether the `CHECK_RETURN(current.back() != 0)` at the end of
  `BerlekampMasseyOdd` (a connection polynomial with Λ_L = 0) can trigger.
  C37 (proved in Lean, for every BM length L) says a decode that would
  otherwise succeed never has Λ_L = 0; a port should still return failure
  in that case, as the C++ does. Whether the check can trigger at all
  (BM ending with Λ_L = 0 on some input) was not determined; it is then
  just a decode failure.
- `minisketch_decode` converts `max_elements` to a C `int`; behaviour for
  `max_elements ≥ 2^31` is not specified by the API and is not specified
  here.
