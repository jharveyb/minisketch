# 2. Sketches

This chapter is normative except where marked. Sources: `src/sketch_impl.h`
(`SketchImpl`, `AddToOddSyndromes*`, from line 2716), `src/minisketch.cpp`
(C API, implementation selection), `src/int_utils.h` (bit I/O, SipHash),
`src/false_positives.h` (capacity functions), `include/minisketch.h` (API
contract).

## Contents

A sketch has a field size b (2..64), an implementation number (0, 1 or 2,
see below), a capacity c ≥ 0, and c field elements

```
u[j] = S_(2j+1) = Σ_{m ∈ M} m^(2j+1),   j = 0 .. c-1
```

the odd power sums of the set M ⊆ K \ {0} it represents. The even power
sums are implied: S_(2k) = S_k^2 (C4). The c odd sums of a set with at
most c elements determine it (C3), and sketches are linear: the sketch of
A Δ B is the XOR of the sketches of A and B (C1).

A sketch also holds, outside its contents, a **basis** element (a nonzero
field element derived from a seed, used only by the root finder) and an
add-batch size. Neither affects serialization or the decoded set.

## Add

`minisketch_add_uint64(v)` (`AddToOddSyndromes`, `src/sketch_impl.h:2716`):

```
add(v):
    m = FromUint64(v)            # low b bits (then LOAD, chapter 1)
    sq = m^2
    d = m
    for j in 0 .. c-1:
        u[j] ^= d                # d = m^(2j+1)
        d = d * sq
```

If the low b bits of v are 0, m = 0 and the loop XORs zeros: a no-op.
Adding the same element twice removes it.

`minisketch_add_uint64s(elements, count)` gives exactly the same result as
calling `add` on each element in order. Internally (informative) it takes K
elements at a time (K = the batch size: 1, 2, 4 or 8, default 4, set with
`minisketch_set_batch_size`, which returns 1 when it accepts the value and
0, changing nothing, for other values) and walks u once per group, with K independent multiplication
chains; leftover elements are added one at a time.

## Merge

`minisketch_merge(s, t)`:

- If s and t differ in b or in implementation number, return 0 and leave s
  unchanged.
- Otherwise set c_s = min(c_s, c_t), u_s[j] ^= u_t[j] for j < c_s, and
  return the new c_s.

For equal capacities, merging equals XOR of the two serializations. A
merge with a capacity-0 sketch truncates s to capacity 0 and returns 0,
which the API cannot tell apart from failure; unlike a failed merge, s
*is* modified.

`minisketch_clone(s)` creates a sketch with the same b, implementation and
capacity and merges s into it. The clone gets a fresh random seed (the
basis is not copied), so it may decode in a different order.

## Serialization

The serialized size is ⌈b·c / 8⌉ bytes. The bytes are the little-endian
encoding of the integer

```
V = Σ_{j<c} SAVE(u[j]) · 2^(b·j)
```

that is, element j occupies bits b·j .. b·j + b − 1 of the byte string,
least significant bit first, bit k of the string being bit (k mod 8) of
byte ⌊k/8⌋. Padding bits in the last byte are 0. SAVE is the identity
except for implementation 1 at b = 58 (chapter 1), so serializations are
the same for every implementation.

```
serialize(u):                   # BitWriter, src/int_utils.h:55
    out = bytes(ceil(b*c/8)); bitpos = 0
    for j in 0 .. c-1:
        v = SAVE(u[j])
        for k in 0 .. b-1:
            if bit k of v: out[bitpos/8] |= 1 << (bitpos % 8)
            bitpos += 1
    return out

deserialize(bytes):             # BitReader, src/int_utils.h:104
    for j in 0 .. c-1:
        u[j] = LOAD(bits b*j .. b*j+b-1 of bytes, as an integer)
```

Deserialization never fails and does not check the padding bits; any byte
string of the right length is a valid sketch (it may not decode). The
caller must pass exactly `minisketch_serialized_size` bytes.

Example (`tests/golden_vectors.json`): b = 8, c = 4, set {31, 66, 150}
serializes to `cb32df3f`. The file has 28 such vectors for
b ∈ {2, 3, 8, 11, 16, 32, 64} and c ∈ {1, 2, 4, 8}.

## Seeds and the basis

Decoding splits polynomials with parameters derived from the sketch's basis
element (chapter 5). The basis changes only the order of decoded elements
and the work done, never the set or the failure (C25: any GF(2)-independent
split parameters separate all roots).

- At creation: basis = `FromSeed(r)` with r a 64-bit value from
  `std::random_device` (fuzzing builds with `MINISKETCH_FUZZ_DETERMINISTIC`
  use r = `0x6d696e69736b6574`).
- `minisketch_set_seed(seed)`: if seed = 2^64 − 1 (−1), basis = 1; else
  basis = `FromSeed(seed)`.

`FromSeed` (`src/fields/generic_common_impl.h:60`,
`src/fields/clmul_common_impl.h:278`):

```
FromSeed(seed):
    k0 = 0x496e744669656c64      # "IntField": implementation 0
         0x434c4d554c466c64      # "CLMULFld": implementations 1 and 2
    count = b << 32
    repeat:
        r = SipHash24(k0, seed, count) & (2^b - 1)
        count += 1
    until r != 0
    return LOAD(r)               # LOAD only for CLMUL fields (identity except b = 58)
```

`SipHash24(k0, k1, m)` (`src/int_utils.h:35`) is standard SipHash-2-4 with
key (k0, k1) on the 8-byte little-endian encoding of the 64-bit message m:

```
v0 = k0 ^ 0x736f6d6570736575;  v1 = k1 ^ 0x646f72616e646f6d
v2 = k0 ^ 0x6c7967656e657261;  v3 = k1 ^ 0x7465646279746573
v3 ^= m;  SipRound x2;  v0 ^= m
v3 ^= 0x0800000000000000;  SipRound x2;  v0 ^= 0x0800000000000000
v2 ^= 0xff;  SipRound x4
return v0 ^ v1 ^ v2 ^ v3

SipRound:
    v0 += v1; v1 = rotl(v1,13); v1 ^= v0; v0 = rotl(v0,32)
    v2 += v3; v3 = rotl(v3,16); v3 ^= v2
    v0 += v3; v3 = rotl(v3,21); v3 ^= v0
    v2 += v1; v1 = rotl(v1,17); v1 ^= v2; v2 = rotl(v2,32)
```

Implementations 0 and 1/2 therefore get different bases from the same seed,
and so may return decoded elements in different orders.

## Implementations and their selection

| number | name | arithmetic | exists for b |
|---|---|---|---|
| 0 | generic | shift-and-add, Multiplier tables | all 2..64 |
| 1 | CLMUL | carry-less multiply + `ReduceWithClMul` | MOD ≠ 3 (table in chapter 1) |
| 2 | CLMUL_TRI | carry-less multiply + `ReduceTrinomial` | trinomial moduli |

All three give identical serializations and decoded sets for the same
inputs.

- `minisketch_implementation_max()` is 2 if the library was compiled with
  carry-less multiply support (`HAVE_CLMUL`), else 0.
- `minisketch_implementation_supported(b, i)`: false unless b is compiled in
  and i ≤ max; then true iff constructing the sketch succeeds, i.e. the
  field instance exists for (b, i) and, for i ∈ {1, 2}, the CPU passes the
  runtime check.
- `minisketch_create(b, i, c)` returns NULL when unsupported (or on
  allocation failure); a new sketch has u = 0.

Runtime check (`EnableClmul`, `src/minisketch.cpp:80`): on x86, CPUID
leaf 1, ECX bit 1 (PCLMULQDQ); on AArch64 Linux `getauxval(AT_HWCAP) &
HWCAP_PMULL` (32-bit ARM: `AT_HWCAP2 & HWCAP2_PMULL`); FreeBSD
`elf_aux_info(AT_HWCAP)` with `HWCAP_PMULL`; Windows
`IsProcessorFeaturePresent(PF_ARM_V8_CRYPTO_INSTRUCTIONS_AVAILABLE)`;
Apple: always true; other targets: trust the compile-time gate.

## Decode (API)

`minisketch_decode(s, max_elements, out)` returns the number of decoded
elements (written to `out`, at most `max_elements`, as `ToUint64` values)
or −1 on failure. Chapter 6 defines which sets are returned.

## Capacity functions

`minisketch_compute_capacity(b, n, fpbits)` returns the capacity c such
that a sketch decoded with `max_elements = n` produces a wrong answer for a
random overfull sketch with probability at most about 2^−fpbits (the model
is C40). `minisketch_compute_max_elements(b, c, fpbits)` is its inverse.
Both are exact integer functions (`src/false_positives.h`), which a port
must reproduce bit for bit.

```
ComputeCapacity(bits, n, fpbits):
    if bits == 0: return 0
    if n > 0xffffffff: return n
    base = BaseFPBits(bits, n)
    if base >= fpbits: return n
    return n + ceil((fpbits - base) / bits)       # (fpbits - base + bits - 1) / bits

ComputeMaxElements(bits, c, fpbits):
    if bits == 0: return 0
    if c > 0xffffffff: return c
    n = c
    loop:
        cap = ComputeCapacity(bits, n, fpbits)    # always >= c
        if cap <= c: return n
        adjust = cap - c
        if n < adjust: return 0
        n -= adjust
```

`ComputeMaxElements` returns max{n ≤ c : ComputeCapacity(n) ≤ c}, or 0
(C10): lowering n by N lowers the capacity by at most N, so the step never
skips a solution.

`BaseFPBits(b, c)` is meant to be ⌊log2(2^(b·c) / Σ_{k≤c} C(2^b − 1, k))⌋
= b·c − ⌈log2 Σ_{k≤c} C(2^b − 1, k)⌉, the false-positive bits of a sketch
whose capacity equals its element count. It equals that value whenever the
value is below 256, and both are ≥ 256 otherwise (C9, checked for
b = 2..64, c ≤ 128 by `tests/differential.py check_capacity_functions`
against the big-integer `base_fpbits` in `tests/pyminisketch.py`). So
`compute_capacity` is exact for fpbits ≤ 256, the documented range.

```
BaseFPBits(bits, c):                              # src/false_positives.h:49
    if c == 0: return 0
    ret = 0
    if bits < 32 and c >= 2^bits:
        ret = bits * (c - 2^bits + 1)
        c = 2^bits - 1
    ret += Log2Factorial(c)
    switch bits:
      2:  return ret + (c <= 2 ? 0 : 1)
      3:  return ret + (c <= 2 ? 0 : (0x2a5 >> 2*(c-3)) & 3)
      4:  return ret + (c <= 3 ? 0 : (0xb6d91a449 >> 3*(c-4)) & 7)
      5:  return ret + (c <= 4 ? 0 : ADD5[c-5])
      6:  return ret + (c <= 4 ? 0 : c > 54 ? 25 : ADD6[c-5])
      7:  return ret + (c <= 4 ? 0 : c > 57 ? 21 : ADD7[c-5])
      8:  return ret + (c <= 9 ? 0 : c > 56 ? 10 : ADD8[c-10])
      9:  return ret + (c <= 11 ? 0 : c > 54 ? 5 : ADD9[c-12])
      10: return ret + (c <= 21 ? 0 : c > 50 ? 2 : (0x1a6665545555041 >> 2*(c-22)) & 3)
      11: return ret + (c <= 21 ? 0 : c > 45 ? 1 : (0x5b3dc1 >> (c-22)) & 1)
      12: return ret + (c <= 21 ? 0 : c > 57 ? 0 : (0xe65522041 >> (c-22)) & 1)
      13: return ret + (c <= 27 ? 0 : c > 55 ? 0 : (0x8904081 >> (c-28)) & 1)
      14: return ret + (c <= 47 ? 0 : c > 48 ? 0 : 1)
      default: return ret
```

(For b ≤ 4, c ≤ 2^b − 1 after the clamp, so the shifts stay in range.) The
correction tables, from `src/false_positives.h` (generated by
`doc/gen_basefpbits.sage`):

```
ADD5 = 1,1,1,1,2,2,2,3,4,4,5,5,6,7,8,8,9,10,10,10,11,11,11,12,12,12,12
ADD6 = 1,0,0,0,1,1,1,2,2,2,2,3,3,4,4,4,5,6,6,6,7,8,8,10,10,11,12,12,13,14,15,15,
       16,17,18,18,19,20,20,21,21,22,22,23,23,23,24,24,24,24
ADD7 = 1,0,0,0,0,1,0,1,1,1,1,1,1,2,2,2,2,3,3,3,4,4,4,5,5,5,6,6,7,7,8,7,8,9,9,9,
       10,11,11,12,12,13,13,15,15,15,16,17,17,18,19,20,20
ADD8 = 1,0,1,0,0,0,0,1,1,1,1,1,2,1,1,2,2,2,3,3,3,3,3,3,4,4,3,4,4,5,4,5,5,5,6,6,
       6,6,7,7,7,8,8,8,8,9,9
ADD9 = 1,0,0,0,0,0,0,1,0,1,1,1,0,1,1,1,2,1,1,1,2,2,2,2,1,2,2,2,2,2,3,2,3,3,3,3,
       4,3,3,4,4,4,4
```

`Log2Factorial(x)` for 1 ≤ x < 2^32 (`src/false_positives.h:19`) equals
⌊log2 x!⌋ for x ≤ 57 and never exceeds it (C8). With
T = 0,4,9,13,18,22,26,30,34,37,41,45,48,52,55,58,62,65,68,71,74,77,80,82,
85,88,90,93,96,98,101,103 (T[i] = ⌊106·log2(1 + i/32)⌋):

```
Log2Factorial(x):                      # x is a uint32, x >= 1
    bits = bit_width(x)                # floor(log2 x) + 1
    i = ((x << (32 - bits)) mod 2^32) >> 26 & 31   # the 5 bits after the leading 1
    l2_106 = 106 * (bits - 1) + T[i]   # underestimate of floor(106 log2 x)
    return (418079 * (2x + 1) * l2_106 - 127870026 * x + 117504694
            + 88632748 * (x < 3 ? 1 : 0)) / 88632748     # 64-bit unsigned arithmetic
```

All arithmetic is unsigned 64-bit and wraps: evaluated left to right the
numerator goes "negative" for small x (x = 1 subtracts 127870026 before the
additions bring it back), so a port with checked arithmetic must use
wrapping operations or add the positive terms first. The final numerator
stays below 2^64 for every x < 2^32. `BaseFPBits`, `ComputeCapacity` and
`ComputeMaxElements` likewise work in unsigned 64-bit integers throughout.

The documented example: b = 32, fpbits = 16, n = 8 gives capacity 9;
`compute_max_elements(32, 9, 16)` returns 9.
