# The decode optimizations from the coding-theory side

The 2026-10 decode optimizations described for readers who know BCH codes,
including the root-finding techniques for large error counts (from
"Evaluation-domain Montgomery reduction" on).
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

## Large error counts: evaluation-domain Montgomery reduction

For ν in the hundreds or more, the Frobenius orbit, the node remainders and
the gcds of BTA dominate. Three changes address them; the first replaces the
modular squarings of the orbit.

- Additive FFTs over subspaces (Gao–Mateer / Lin–Chung–Han, the tools of
  binary Reed–Solomon erasure decoding) evaluate on an affine subspace C.
  That makes F[x]/(Z_C) ≅ F^N, the CRT at the points of C.
- With R = Z_C, Montgomery reduction mod Λ becomes pointwise:
  q = a·Λ⁻¹ on C.
- The quotient (a + qΛ)/Z_C is recovered on C from derivatives:
  r(w) = h′(w)/Z_C′(w). This is the same "residue at a simple pole" formula
  as Forney's algorithm and Lagrange interpolation over a vanishing
  polynomial.
- For a subspace (linearized) polynomial Z_C′ is constant (here 1), and for
  squared inputs the a′ term vanishes, as in the binary-BCH simplification.
- Net cost per Frobenius step: one IFFT and one FFT on N ≥ deg Λ points.

## Trace codes and the GF(16) tower

- Tr₆₄→₁ = Tr₄→₁ ∘ Tr₆₄→₄. The relative trace T = Tr₆₄→₄(βx) mod g
  classifies the error locators by 16 values.
- The binary splits at 4 consecutive levels are the Tr₄→₁(zᵏT) for a
  GF(2)-basis {zᵏ} of GF(16). (For large Λ the tower uses GF(16) rather than
  GF(2⁸): carrying a relative trace down the tree is cheaper than the extra
  squarings a larger subfield costs per node.)
- T is carried down the factorization tree as a residue mod each node.
  Moving to a node's factor uses CRT on the halves of the evaluation coset:
  Z_C = Z_C₀·Z_C₁ with Z_C₀ ≡ 1 on C₁.

## Syndrome splitting instead of gcds

- The trace residue t is an **idempotent** of F[x]/(Λ): it is the CRT
  indicator of the error subset A.
- The syndromes of the sub-pattern A are a correlation of t with the
  syndromes: S_k^(A) = Σⱼ tⱼS_{j+k}. This is the generating-function
  identity P(y) = rev(t·Λ′ mod Λ)/rev(Λ).
- Because error values in a binary code are 1, the masked sequence keeps
  S₂ₖ = Sₖ². So the half-length (odd-syndrome) BM applies, and it returns
  Λ_A for about deg Λ/2 errors.
- Euclid on (Λ, t), the Sugiyama form, would run the full half-gcd from
  deg Λ down to deg Λ/2 without being able to exploit the binary structure.
  For a non-binary (Reed–Solomon-like) code the masked syndromes would be
  weighted, and you'd lose the factor of 2.
- Λ_B = Λ/Λ_A is computed in the evaluation domain. The exact-degree check
  certifies Λ_A | Λ, which protects against too few syndromes or
  inconsistent syndromes, for example from an overfull sketch.

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
