# The decode optimizations from the coding-theory side

The 2026-10 decode optimizations described for readers who know BCH codes.
[decode_explained.md](decode_explained.md) covers the same ideas without coding
theory, and [math_decode.md](math_decode.md) in the style of [math.md](math.md).

## The code

PinSketch is a binary narrow-sense BCH code of length 2ᵇ − 1 and designed
distance 2c+1. Set elements are the **error locators** Xᵢ, and the sketch is
the odd syndromes S₁, S₃, …, S₂c₋₁. Decoding is bounded-distance decoding with
ν ≤ t = c errors: find the error-locator polynomial Λ(x) = Π(1 − Xᵢx), then its
roots. Chien search is impossible over a 2⁶⁴-element field, so roots are found
with Berlekamp's trace algorithm (BTA).

## Simplified BM for binary codes

This is Berlekamp's binary-BCH simplification: S₂ⱼ = Sⱼ², so every second
discrepancy vanishes and the iteration needs only t steps instead of 2t. The
textbook treatment usually assumes a genuine error pattern. The decoder relies
on it for arbitrary Frobenius-structured syndrome vectors, and an exhaustive
oracle comparison on small fields confirms it.

## Early termination as miscorrection detection

With ν ≪ t, the locator is determined after about 2ν syndromes, at which point
the remaining syndromes are pure redundancy. A bounded-distance decoder must
still check that all 2t syndromes satisfy the key equation. Otherwise a
received word outside every decoding sphere could be miscorrected, and here
low-weight codewords make that adversarially constructible.

The check is an LFSR-consistency test, i.e. Λ applied as an FIR filter to the
syndrome stream. Two improvements:

- **It reduces to the odd syndromes.** Writing C = ΛS gives
  ΛC = C² + Λ²·S_odd. The odd syndromes uⱼ = S₂ⱼ₊₁ = Σ Xᵢ·(Xᵢ²)ʲ are themselves
  a *nonbinary* syndrome sequence, as for a Reed–Solomon code, with locators
  Xᵢ², error values Xᵢ and locator polynomial Λ with each coefficient squared.
  Past the convergence point, checking that sequence is equivalent to the full
  key-equation check.
- **It's done blockwise.** Direct dot products for small ν (below 64),
  transposed Karatsuba (fast FIR filtering) for medium ν (64 to 511), the
  additive FFT for large ν (512 and up).

The cost goes from O(t·ν) sequential work to O(ν²) plus O(t·ν^0.585) parallel
work.

## BTA with a shared Frobenius orbit

BTA splits Λ by gcd(Λ, Tr(βx) mod Λ). The traces for the first log₂m + 2
splitting levels all come from one computation of the orbit x^(2ⁱ) mod Λ at the
root: every factor divides Λ, so reduction commutes down the factor tree. One
more squaring tests whether x^(2ᵐ) ≡ x mod Λ. That is exactly the
decoding-failure test "Λ has ν distinct roots in GF(2ᵐ)", done before any
splitting instead of being discovered deep in the tree.

## Tower-basis split parameters

BTA only needs the β's to be GF(2)-linearly independent. Taking β = γ·xˡ·zʲ
with z generating GF(2⁸) makes the first eight traces combinations of only 8
Frobenius-residue buckets.

## Arithmetic

Delayed modular reduction with carry-less multiplication is standard in
GF(2ᵐ) hardware and software. Fraction-free (pseudo-remainder) Euclid, also
standard, avoids per-step inversions in the gcd.

## What did not change

None of this changes the decoding radius or the failure-detection behaviour.
It is the same bounded-distance decoder, just doing less work.

See [optimizations.md](optimizations.md) for the derivations, measurements and
negative results, and the footnote in [math.md](math.md) for how PinSketch maps
onto BCH decoding.
