# 1. The field GF(2^b)

Every sketch works over K = GF(2^b) for one b in 2..64. This chapter defines
the element encoding (normative) and the arithmetic the decoder uses
(informative: any correct field arithmetic gives the same decode results).

Sources: `src/fields/clmul_common_impl.h` (CLMUL/PMULL fields),
`src/fields/generic_common_impl.h` (table fields), `src/int_utils.h`
(shift-and-add multiplication, inversion, SipHash, bit I/O),
`src/lintrans.h` (table-driven GF(2)-linear maps), the per-size instance
files `src/fields/{clmul,generic}_{1..8}byte(s).cpp`, and
`doc/gen_params.sage`, which generated them.

## Element representation (normative)

An element is the polynomial a(x) = Σ_{k<b} a_k x^k over GF(2), reduced
modulo the field's canonical modulus x^b + MOD(x). Its integer encoding is
Σ a_k 2^k, an integer in [0, 2^b). Addition is XOR. A `uint64_t` element
given to the API is reduced to its low b bits (`FromUint64` masks); 0 is
not a valid set element (adding it is a no-op, chapter 2).

The canonical moduli are the minimal-weight, lexicographically first
irreducible polynomials of each degree (C7; `doc/moduli.md`,
`GF2_MODULI` in `tests/pyminisketch.py`). The table gives MOD as an
integer (x^b + MOD with MOD = Σ m_k 2^k), and which CLMUL implementations
exist for each b (chapter 2, "Implementations"); implementation 0 (generic)
exists for every b.

| b | modulus | MOD | impl 1 (CLMUL) | impl 2 (CLMUL_TRI) |
|---|---|---|---|---|
| 2 | x^2 + x + 1 | 3 | - | yes |
| 3 | x^3 + x + 1 | 3 | - | yes |
| 4 | x^4 + x + 1 | 3 | - | yes |
| 5 | x^5 + x^2 + 1 | 5 | yes | yes |
| 6 | x^6 + x + 1 | 3 | - | yes |
| 7 | x^7 + x + 1 | 3 | - | yes |
| 8 | x^8 + x^4 + x^3 + x + 1 | 27 | yes | - |
| 9 | x^9 + x + 1 | 3 | - | yes |
| 10 | x^10 + x^3 + 1 | 9 | yes | yes |
| 11 | x^11 + x^2 + 1 | 5 | yes | yes |
| 12 | x^12 + x^3 + 1 | 9 | yes | yes |
| 13 | x^13 + x^4 + x^3 + x + 1 | 27 | yes | - |
| 14 | x^14 + x^5 + 1 | 33 | yes | yes |
| 15 | x^15 + x + 1 | 3 | - | yes |
| 16 | x^16 + x^5 + x^3 + x + 1 | 43 | yes | - |
| 17 | x^17 + x^3 + 1 | 9 | yes | yes |
| 18 | x^18 + x^3 + 1 | 9 | yes | yes |
| 19 | x^19 + x^5 + x^2 + x + 1 | 39 | yes | - |
| 20 | x^20 + x^3 + 1 | 9 | yes | yes |
| 21 | x^21 + x^2 + 1 | 5 | yes | yes |
| 22 | x^22 + x + 1 | 3 | - | yes |
| 23 | x^23 + x^5 + 1 | 33 | yes | yes |
| 24 | x^24 + x^4 + x^3 + x + 1 | 27 | yes | - |
| 25 | x^25 + x^3 + 1 | 9 | yes | yes |
| 26 | x^26 + x^4 + x^3 + x + 1 | 27 | yes | - |
| 27 | x^27 + x^5 + x^2 + x + 1 | 39 | yes | - |
| 28 | x^28 + x + 1 | 3 | - | yes |
| 29 | x^29 + x^2 + 1 | 5 | yes | yes |
| 30 | x^30 + x + 1 | 3 | - | yes |
| 31 | x^31 + x^3 + 1 | 9 | yes | yes |
| 32 | x^32 + x^7 + x^3 + x^2 + 1 | 141 | yes | - |
| 33 | x^33 + x^10 + 1 | 1025 | yes | yes |
| 34 | x^34 + x^7 + 1 | 129 | yes | yes |
| 35 | x^35 + x^2 + 1 | 5 | yes | yes |
| 36 | x^36 + x^9 + 1 | 513 | yes | yes |
| 37 | x^37 + x^6 + x^4 + x + 1 | 83 | yes | - |
| 38 | x^38 + x^6 + x^5 + x + 1 | 99 | yes | - |
| 39 | x^39 + x^4 + 1 | 17 | yes | yes |
| 40 | x^40 + x^5 + x^4 + x^3 + 1 | 57 | yes | - |
| 41 | x^41 + x^3 + 1 | 9 | yes | yes |
| 42 | x^42 + x^7 + 1 | 129 | yes | yes |
| 43 | x^43 + x^6 + x^4 + x^3 + 1 | 89 | yes | - |
| 44 | x^44 + x^5 + 1 | 33 | yes | yes |
| 45 | x^45 + x^4 + x^3 + x + 1 | 27 | yes | - |
| 46 | x^46 + x + 1 | 3 | - | yes |
| 47 | x^47 + x^5 + 1 | 33 | yes | yes |
| 48 | x^48 + x^5 + x^3 + x^2 + 1 | 45 | yes | - |
| 49 | x^49 + x^9 + 1 | 513 | yes | yes |
| 50 | x^50 + x^4 + x^3 + x^2 + 1 | 29 | yes | - |
| 51 | x^51 + x^6 + x^3 + x + 1 | 75 | yes | - |
| 52 | x^52 + x^3 + 1 | 9 | yes | yes |
| 53 | x^53 + x^6 + x^2 + x + 1 | 71 | yes | - |
| 54 | x^54 + x^9 + 1 | 513 | yes | yes |
| 55 | x^55 + x^7 + 1 | 129 | yes | yes |
| 56 | x^56 + x^7 + x^4 + x^2 + 1 | 149 | yes | - |
| 57 | x^57 + x^4 + 1 | 17 | yes | yes |
| 58 | x^58 + x^19 + 1 | 524289 | yes (internal MOD 99, see below) | yes |
| 59 | x^59 + x^7 + x^4 + x^2 + 1 | 149 | yes | - |
| 60 | x^60 + x + 1 | 3 | - | yes |
| 61 | x^61 + x^5 + x^2 + x + 1 | 39 | yes | - |
| 62 | x^62 + x^29 + 1 | 536870913 | yes | yes |
| 63 | x^63 + x + 1 | 3 | - | yes |
| 64 | x^64 + x^4 + x^3 + x + 1 | 27 | yes | - |

Implementation 1 is missing exactly where MOD = 3 (x^b + x + 1, where the
trinomial reduction is strictly better); implementation 2 exists exactly
where the canonical modulus is a trinomial. Every implementation uses the
canonical modulus for its arithmetic except implementation 1 at b = 58.

## The 58-bit CLMUL field: a different internal modulus

`Field58` in `src/fields/clmul_8bytes.cpp:47` computes internally modulo
y^58 + y^6 + y^5 + y + 1 (MOD = 99), not the canonical x^58 + x^19 + 1.
`gen_params.sage` (`pick_modulus`, style `CLMUL`) prefers moduli with
b + deg(MOD) ≤ 66, so that `ReduceWithClMul` can use 64-bit lane shifts;
58 + 19 = 77 fails that, so it took the first such alternative.

The two representations are related by a field isomorphism, applied at
every boundary:

- `LOAD` (canonical → internal) is the GF(2)-linear map with
  LOAD(2^i) = ν^i, where ν = `0x77dd188d5d600` (internal encoding) is a
  root of x^58 + x^19 + 1 in GF(2)[y]/(y^58 + y^6 + y^5 + y + 1). So
  LOAD(a(x)) = a(ν).
- `SAVE` (internal → canonical) is its inverse.
- The tables are `LOAD_TABLE_58` and `SAVE_TABLE_58` (58 entries each,
  image of each bit).

`FromUint64`, `Deserialize` and `FromSeed` apply LOAD; `ToUint64` and
`Serialize` apply SAVE. All other operations stay internal. The
consequences:

- Serializations and decoded sets are identical to the canonical
  representation (the map is an isomorphism, so power sums, Λ and roots
  correspond).
- `Mul2` multiplies by y (the internal generator), which is a different
  field element from canonical 2, and `Qrt` picks its solution by the
  internal bit 0. So split parameters, the Cantor basis and the
  seed-derived basis differ from implementations 0 and 2, which changes
  only the output order.

A port should use the canonical modulus everywhere. The trinomial recipe
below handles x^58 + x^19 + 1 (2·19 − 2 < 58).

No other field instance has a non-identity LOAD/SAVE (all others are
`ID_TRANS`).

## Multiplication

The product of a, b ∈ K is the carry-less (GF(2)[x]) product a·b, of degree
≤ 2b − 2, reduced modulo x^b + MOD.

### Generic fields (implementation 0)

`GFMul` (`src/int_utils.h:250`) is shift-and-add with reduction per step:

```
mul(a, c):              # a, c < 2^b
    r = 0
    for k in 0..b-1:
        if bit k of c: r ^= a
        a = mul2(a)
    return r

mul2(a):                # multiply by x
    t = a >> (b-1)
    a = (a << 1) & (2^b - 1)
    return a ^ (MOD if t else 0)
```

For loops that multiply many values by one scalar a, the generic field
builds a `Multiplier`: lookup tables of the GF(2)-linear map v ↦ a·v, one
table per chunk of 3–6 input bits (`RecLinTrans`/`LinTrans`,
`src/lintrans.h`; chunk sizes per field are the `DynTable` typedefs, for
example 4,4,4,4,3,3,3 for b = 25). The table for the chunk starting at bit
o holds all XOR combinations of a·x^o, a·x^(o+1), …. Building it costs
about 2^chunk XORs per chunk, which is why the decoder's tier choices
distinguish "cheap" and "table" multipliers (chapter 7).

### CLMUL fields (implementations 1 and 2)

`MulWithClMulReduce` / `MulTrinomial` (`src/fields/clmul_common_impl.h:173`
and `:218`): one 64×64 → 128-bit carry-less multiply (x86 `PCLMULQDQ`,
ARMv8 `PMULL`) of the two elements, then a reduction recipe. The
`Multiplier` just stores the scalar.

## Reduction recipes (C47)

Both recipes reduce a 128-bit carry-less product P of degree ≤ 2b − 2.
Write P_lo for P mod x^b and "fold" for the identity
x^b·h ≡ h·MOD (mod x^b + MOD). All shifts below are logical; "lane shift"
means a shift of each 64-bit half separately (`_mm_srli_epi64`), "128-bit
shift" a shift of the full value (built from lane and byte shifts).

### `ReduceWithClMul` (implementation 1), `src/fields/clmul_common_impl.h:123`

Two folds, each one carry-less multiply by MOD:

```
reduce(P):
    h1 = P >> b                 # degree <= b-2
    r1 = clmul(h1, MOD)         # degree <= b-2+deg(MOD)
    h2 = r1 >> b                # degree <= deg(MOD)-2
    r2 = clmul(h2, MOD)         # degree <= 2*deg(MOD)-2 < b
    return (P ^ r1 ^ r2) & (2^b - 1)
```

Correctness requires 2·deg(MOD) − 2 < b, true for every instance (largest
case b = 62, deg 29: 56 < 62). The code specializes the shifts:

| case | h1 | h2 |
|---|---|---|
| b ≤ 32 | lane shift (P fits in the low lane) | lane shift |
| b = 64 | high 64-bit word of P (`Clmul01`: hi(P) × MOD) | high word of r1 |
| b ∈ {40, 48, 56} | 128-bit shift by b/8 bytes | 128-bit byte shift |
| other b > 32 | 128-bit shift by b | lane shift if MOD >> (66 − b) == 0 (that is, deg(MOD) < 66 − b, so r1 fits in 64 bits), else 128-bit shift |

**ARM, b = 64, MOD < 256.** On AArch64 the 64-bit field (MOD = 27) folds
the high word in general registers instead of two dependent PMULLs:

```
lo, hi = P
o   = (hi ⊗ MOD) >> 64          # bits of the 72-bit product above 2^64, o < 2^7
return lo ^ ((hi ⊗ MOD) mod 2^64) ^ ((o ⊗ MOD) mod 2^64)
```

where ⊗ is carry-less multiplication by the 8-bit constant, written out
as shifts (`ShiftMul8`, `FoldOverflow`, unrolled by hand: compilers left
the constant-trip loops rolled at -O2, five times slower). The second fold
cannot overflow since o ⊗ MOD < 2^15. Measured on Cortex-A76: Mul 3.0 vs
4.0 ns throughput, 6.4 vs 10.0 ns latency; x86-64 keeps PCLMULQDQ (the
shift form is 2× slower on Haswell).

### `ReduceTrinomial` (implementation 2), `src/fields/clmul_common_impl.h:179`

For x^b + x^POS + 1, MOD = 1 + x^POS, so a fold is h ^ (h << POS), no
multiply. Requires POS < b and 2·POS − 2 < b (a `static_assert`), so two
folds suffice:

```
reduce_tri(P):
    h1 = P >> b                         # degree <= b-2
    r1 = h1 ^ (h1 << POS)               # see the b > 32 cases below
    if POS == 1: return (P ^ r1) & mask # r1 has degree <= b-1: done
    h2 = r1 >> b                        # degree <= POS-2
    r2 = h2 ^ (h2 << POS)               # degree <= 2*POS-2 < b
    return (P ^ r1 ^ r2) & mask
```

- b ≤ 32: lane shifts throughout.
- b > 32: h1 by 128-bit shift. If b + POS ≤ 65, h1 (degree ≤ b − 2)
  shifted by POS still fits a 64-bit lane, and the lane-shift form above
  is used. Otherwise (b = 58, POS = 19 and b = 62, POS = 29) the first fold
  is a carry-less multiply, r1 = clmul(h1, 1 + 2^POS), with h2 by 128-bit
  shift, and the second fold is the shift form.

### Linearity

Both recipes are GF(2)-linear in P (C39f): reduce(P ⊕ Q) =
reduce(P) ⊕ reduce(Q), and an XOR of products of reduced elements still has
degree ≤ 2b − 2. This is what makes lazy reduction valid (below).

## Lazy (unreduced) accumulation

CLMUL fields expose (`src/fields/clmul_common_impl.h:256`):

| operation | meaning |
|---|---|
| `Wide` | 128-bit unreduced value |
| `MulWide(a, b)` | carry-less product, not reduced |
| `WideXor(u, v)` | XOR |
| `WideFrom(a)` | a reduced element as a Wide (Reduce(WideFrom(a)) = a) |
| `Reduce(w)` | the reduction recipe |
| `LAZY_REDUCE = true` | |

A sum Σ a_i·b_i is computed as Reduce(⊕ MulWide(a_i, b_i)): one carry-less
multiply per term instead of three, and the accumulator lives in a
register. The discipline: only XOR together products of reduced elements
and `WideFrom` of reduced elements; reduce once when the sum is complete.
The accumulator never overflows (degree ≤ 2b − 2 < 128).

Generic fields expose trivial versions (`Wide` = element, `MulWide` =
`Mul`, `Reduce` = identity) with `LAZY_REDUCE = false`. Every hot loop in
the decoder has two code paths: a lazy one for CLMUL fields and a
`Multiplier`-table one for generic fields (CLAUDE.md: "change and test
both"). They compute the same values.

The decoder also distinguishes fields by multiplier cost:
`TraceModCheapScalarMul` is true when `sizeof(Multiplier) ≤ 2·sizeof(Elem)`,
that is, for the CLMUL fields (both implementations), false for generic
fields (`src/sketch_impl.h:387`).

## Squaring

Squaring is GF(2)-linear: (Σ a_k x^k)^2 = Σ a_k x^(2k) mod M. All fields
have a table `SQR_TABLE_b` with the images of each bit, x^(2k) mod M,
applied through chunked lookup tables (`RecLinTrans`, chunks of 5–6 bits).
CLMUL fields with b > 32 instead compute `Mul(a, a)` (one carry-less
product: 64-bit field 2.0 vs 3.4 ns on x86-64, 3.7 vs 6.1 ns on A76;
`src/fields/clmul_common_impl.h:269`). CLMUL fields also carry tables
`SQR2/4/8/16_TABLE_b` for a ↦ a^(2^2), a^(2^4), …, used by inversion.

## Qrt: solving y^2 + y = a

`Qrt(a)` returns y with y^2 + y = a when a solution exists, i.e. when
Tr(a) = 0. The two solutions are y and y + 1; minisketch returns the one
with bit 0 (the coefficient of x^0) equal to 0. `Qrt` is a GF(2)-linear map
given by `QRT_TABLE_b` (the image of each bit; generated by `qrt_table` in
`doc/gen_params.sage`, which clears bit 0 of every image, so every output
has bit 0 clear). For inputs with Tr(a) = 1 the output is meaningless; every
caller either checks y^2 + y = a or knows Tr(a) = 0.

A port may compute Qrt any way that returns the bit-0-clear solution. For
odd b the half-trace H(a) = Σ_{i=0}^{(b−1)/2} a^(2^(2i)) is a solution
(clear bit 0 afterwards with y ^= y & 1); for any b, Gaussian elimination on
the GF(2)-linear map y ↦ y^2 + y works (`solve_artin_schreier` in
`tests/pyminisketch.py`, then clear bit 0). The choice of solution affects
only the Cantor basis (chapter 3) and hence output order.

## Inversion

`Inv(a)` for a ≠ 0 (Inv(0) returns 0 in both implementations; the decoder
never relies on that):

- Generic: binary extended Euclid on the bit polynomials
  (`InvExtGCD`, `src/int_utils.h:254`).
- CLMUL: an addition-chain ladder for a^(2^b − 2) (`InvLadder`,
  `src/int_utils.h:280`): with e = b − 1 it forms
  x_(2^k) = a^(2^(2^k) − 1) for 2^k ≤ e using the SQR2^k tables, combines
  them by the bits of e into a^(2^e − 1), and squares once.

Any correct inversion conforms. The decoder avoids inversions in hot loops
(batch inversion, inversion-free Euclid steps; chapters 3 and 5).

## Field interface used by the rest of this spec

| operation | |
|---|---|
| `Bits()` | b |
| `Mul(a, c)`, `Sqr(a)`, `Inv(a)`, `Qrt(a)` | as above |
| `Mul2(a)` | a·x (x the field's generator, internal representation) |
| `Multiplier(a)` | callable v ↦ a·v (table or wrapper) |
| `Wide`, `MulWide`, `WideXor`, `WideFrom`, `Reduce`, `LAZY_REDUCE` | lazy reduction |
| `FromUint64`, `ToUint64` | mask to b bits and LOAD; SAVE |
| `Serialize`, `Deserialize` | b bits of SAVE(v) / LOAD of b bits (chapter 2) |
| `FromSeed(seed)` | nonzero element from a seed (chapter 2) |
