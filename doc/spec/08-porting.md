# 8. Porting notes (Rust, Go)

Informative. What a port needs, in what order, and how to test it.

## A minimal conforming port

Chapters 1, 2 and 6 suffice: any correct GF(2^b) arithmetic with the
canonical moduli, the sketch operations and serialization of chapter 2, the
exact capacity functions, and the textbook decoder (BM over all 2c
syndromes, the X^(2^b) = x test, Berlekamp trace splitting). The first half
of `tests/pyminisketch.py` is such an implementation in about 300 lines. It
is byte-compatible and decodes the same sets; it is only slower: BM is
O(c·n) serial work, and root finding costs about b·n^2 per tree level.

Then add, roughly in order of payoff (measurements from
`doc/optimizations.md`, "Measured" sections; i7-4770 / Pi 5, 64-bit field,
capacity/differences):

| step | chapter | regime where it matters | measured effect |
|---|---|---|---|
| carry-less multiply and lazy reduction | 1 | everything | ~3× fewer multiplies in dot products; root chain at the PCLMULQDQ throughput limit |
| odd-step BM + early termination + tail check | 4 | n ≪ c | part of the 2026-10 series: 3.6–8.9× at 1024–8192 capacity, 16–256 differences; odd-syndrome tail alone +27–45% at 4096/16 |
| shared orbit from degree 8, tower betas | 5 | n ≥ 8 | orbit from 8 (was 16): +43–51% at 1024/12; GF(2^8) tower +7–12% at 1024/32–64 |
| square-table reducer (column-blocked) | 5 | 8 ≤ n < 704 | the dominant root-finding stage at moderate n |
| Karatsuba / FFT tail tiers | 4 | n ≥ 64 with c ≫ n | Karatsuba tail +12–22% at 4096–8192/128 (i7) |
| additive FFT products, Barrett reducer | 3, 5 | n ≥ 704 (b ∈ {16,32,64}) | pre-lazy: 1.4–2.3× at 500–4096 full |
| Montgomery orbit, carries, CosetDivide, power-sum splits | 5 | n ≥ 224, b ∈ {16,32,64} | 2.2–2.5× at 1024–4096, 50–100% full (32-bit: 1.5–1.8×) |
| FastDivMod | 3 | degrees > 1024 | small; safe to omit |

## Carry-less multiplication

**Rust.** `core::arch::x86_64::_mm_clmulepi64_si128` (enable with
`#[target_feature(enable = "pclmulqdq")]`, detect with
`is_x86_feature_detected!("pclmulqdq")`), and on AArch64
`core::arch::aarch64::vmull_p64` (target feature `aes`, detect with
`is_aarch64_feature_detected!("pmull")`). Use `u128` (or `__m128i` /
`uint64x2_t`) for the unreduced `Wide` type. Monomorphize the field per b
(const generics or macros), as the C++ does with templates: the reduction
shifts and constants must be compile-time constants to be cheap. Dispatch
once at sketch construction, like `EnableClmul` (chapter 2), and keep the
hot loops inside `#[target_feature]` functions so the intrinsics inline.

**Go.** There is no carry-less multiply intrinsic. Options:

- Go assembly (`PCLMULQDQ` on amd64, `VPMULL` on arm64; the standard
  library's GCM code is a model), with `golang.org/x/sys/cpu`
  (`cpu.X86.HasPCLMULQDQ`, `cpu.ARM64.HasPMULL`) for detection. Assembly
  functions are not inlined, so a per-multiplication call is too slow:
  write whole kernels in assembly (the lazy dot product, the square-table
  reduction, the pointwise loops of the Montgomery squaring, the BM
  discrepancy and tail check).
- A table-based field (implementation 0's design): `Multiplier` tables of
  the map v ↦ a·v by 4-bit chunks (16 entries per chunk, ⌈b/4⌉ chunks) for
  loops that reuse a scalar, shift-and-add for one-off products. This is
  also the fallback for CPUs without CLMUL.
- `math/bits` has no carry-less operations; a software 64×64 carry-less
  multiply (bit-sliced with integer multiplies and masks, as in BearSSL's
  `ghash_ctmul64`) is a portable middle ground.

Represent `Wide` as two `uint64` (hi, lo); the reduction recipes of
chapter 1 are a few shifts and XORs on them.

## Field implementation choices

- Use the canonical modulus for every b, including 58 (the C++ CLMUL
  field's internal modulus is an artifact of a lane-shift constraint;
  chapter 1). `ReduceTrinomial`'s two-fold form works for every trinomial
  modulus; `ReduceWithClMul` for the rest (2·deg MOD − 2 < b holds for all
  canonical moduli).
- Squaring: one carry-less product for b > 32; a linear-map table is as
  fast at b ≤ 32.
- Qrt: a precomputed linear map (b images) or the half-trace for odd b;
  return the solution with bit 0 clear if you want libminisketch's Cantor
  basis and root order.
- Inversion: any; it is rare (batch inversion and inversion-free Euclid
  avoid it in loops).

## Tiers for table-based fields

If a port keeps a table-based field, mirror the C++ distinctions (chapter
7): the generic fields use Multiplier rows instead of lazy dot products,
the square/Barrett reducer cutoff stays at 512, FFT tail checks start at
L = 40, and FastDivMod and the whole Montgomery path are off
(each measured slower there: the FFT's pointwise products are one-off
multiplications). The BM builds one Multiplier per processed syndrome,
only for the walked prefix.

## Memory and scratch layout

- **Square table.** d^2/2 elements for a degree-d node (rows x^e mod g for
  even e in [d, 2d − 2]). For lazy fields store it column-block-major
  (blocks of 8 output columns, all rows of a block contiguous) and keep 8
  `Wide` accumulators in registers; the row-major layout had more
  concurrent streams than the prefetchers track.
- **Scratch reuse.** The Karatsuba recursion uses fixed per-level scratch
  slots and one FFT engine per reducer, so the trace loop is
  allocation-free in steady state; Rust can do the same with a
  `Vec<Vec<Elem>>` arena indexed by level, Go with per-reducer slices.
- **Recursion.** The root-finding recursion depth is at most b (C24), so
  plain recursion is fine. Release a node's coset state (values of g and
  g', inverse values, carried traces) before recursing; summed along the
  recursion path it would otherwise dominate peak memory.
- **Peak heap** (64-bit field, measured): about 27 KB at 1024/32 and
  0.3–1.5 MB at 1024–4096 full capacity on CLMUL. Generic fields: up to a
  few MB from Multiplier tables.

## What can be omitted

Decode is at least quadratic in n (the serial BM part is O(n^2)), so many
techniques are constant factors:

- **Safe to omit** (constant factors, or never engaged at practical
  sizes): FastDivMod (engages only above degree 1024), cached
  reducer transforms, the inversion-free Euclid step, SmallTrace and
  TowerTraceSmall (use the general paths), batched add, the ARM
  shift-fold, PowerSumSplit (the gcd gives the same split at about twice
  the cost), CosetDivide (use long division).
- **Changes the asymptotic profile** if omitted: the early-terminating BM
  with a subquadratic tail check (otherwise O(c·n) serial for n ≪ c; the
  direct tail check alone is O(c·n) but parallel and lazy), and, for large
  n, the FFT-based root chain (Montgomery or Barrett with FFT products),
  which replaces ~b·n^2 by ~b·n log n for the Frobenius orbit.
- **Required for conformance:** nothing beyond chapters 1, 2 and 6. Every
  optional technique must, if implemented, keep the exact result; the
  claims cited in chapters 4 and 5 are the conditions (for example the
  CosetDivide monic-degree check, C33, and PowerSumSplit's fallback, C35).

## Testing a port

1. **Golden vectors.** `tests/golden_vectors.json`: a list of objects
   `{bits, capacity, elements, serialized_hex}` (28 vectors,
   b ∈ {2, 3, 8, 11, 16, 32, 64}, c ∈ {1, 2, 4, 8}, elements sorted).
   Adding the elements must give exactly `serialized_hex`, and decoding
   the deserialized bytes must give `elements` (after sorting).
   `pyminisketch.TestGoldenVectors` and `src/test.cpp TestGoldenVectors`
   check the same file.
2. **Differential testing.** `tests/differential.py` loads libminisketch
   via ctypes and, for seeded random cases (b ∈ 2..16 ∪ {32, 64},
   c ≤ 16, up to 2c random elements, so overfull sketches included),
   checks against pyminisketch: byte-identical serialization,
   cross-deserialization, identical decode results (sorted, including
   failure), merge equal to XOR of serializations, for every
   implementation. A port can be wired in the same way (FFI, or a small
   CLI that reads cases on stdin) and compared to the C library or to
   pyminisketch.
3. **Capacity functions.** Compare `compute_capacity` for b = 2..64,
   n = 0..128 and fpbits ∈ {0, 1, 255, 256, random ≤ 256}, and
   `compute_max_elements` on random (c, fpbits), against the big-integer
   definitions `compute_capacity` / `compute_max_elements` in
   pyminisketch (`check_capacity_functions`).
4. **Exhaustive tiny fields.** For GF(4)…GF(32) at small capacity,
   decode every possible syndrome vector and compare with the textbook
   decoder; check that exactly C(2^b − 1, i) vectors decode to i elements
   (C11) and that every decoded set re-encodes to its input and has no 0
   (C37). pyminisketch's `TestOracle.test_decode_fast_exhaustive` and
   `test_bm_odd_exhaustive` do this for the oracle.
5. **Large-degree tiers.** Small random cases never reach the fast tiers.
   Include cases that do, and check each tier actually ran: BM tail with
   L ≥ 64 and L ≥ 512 (and ≥ 40 on table fields), capacity ≫ n; root
   degree ≥ 224 (Montgomery), ≥ 256 (power-sum splits), ≥ 704 (Barrett
   reducer), on b ∈ {16, 32, 64}; roots planted on the first Montgomery
   coset candidates to force the coset fallback (`TestFindRootsMontgomery`
   plants roots on all four to force the standard orbit); roots chosen in
   the kernel of the first split functionals to force long chains of
   trivial splits (`TestFindRootsDeepTree`). `tests/conformance.py` (being
   added at the time of writing) compares libminisketch's internal stages
   against the oracle on such cases and counts which tiers ran.
6. **Adversarial and overfull inputs.** Random byte strings, sketches of
   more than c elements, `max_elements` below and above the true size,
   and sketches whose first syndromes match a small set but whose tail
   does not (the case the exact tail check exists for). The expected
   result is always the textbook decoder's.
7. **The oracle as pseudocode.** For each fast technique, the
   corresponding function in pyminisketch's "Optimized decode (oracle)"
   section (`bm_odd`, `middle_product_kara`, `montgomery_orbit`,
   `mont_square`, `mont_carry_to_factor`, `coset_divide`,
   `power_sum_split`, `find_roots_fast`, `decode_fast`) is a short,
   tested reference to unit-test a port's component against.
