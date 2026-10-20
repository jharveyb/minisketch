# Fast root finding for large sketches: the coding-theory view

This note restates the root-finding optimizations of the 2026-10 "coset
Montgomery / power-sum" series (see [optimizations.md](optimizations.md)
for costs and measurements) in coding-theory terms. An introduction without
that background is [root-finding-intro.md](root-finding-intro.md).

## Setting

PinSketch is a narrow-sense **binary BCH code** of length 2ᵇ − 1 with
designed distance 2c + 1:

- The set elements are the **error locators** X_i.
- The sketch is the syndrome vector S_k = Σᵢ X_iᵏ.
- Decoding is error-locator computation (BM, i.e. solving the key
  equation) followed by root finding.

A Chien search over 2⁶⁴ positions is impossible. So root finding is
equal-degree factorization, here the **Berlekamp trace algorithm** (the
characteristic-2 Cantor–Zassenhaus variant): gcd(Λ, Tr(βx) mod Λ) for a
GF(2)-linearly independent sequence of β.

## Evaluation-domain Montgomery reduction

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

## Trace codes and the subfield tower

- Tr₆₄→₁ = Tr₄→₁ ∘ Tr₆₄→₄. The relative trace T = Tr₆₄→₄(βx) mod g
  classifies the error locators by 16 values.
- The binary splits at 4 consecutive levels are the Tr₄→₁(zᵏT) for a
  GF(2)-basis {zᵏ} of GF(16).
- T is carried down the factorization tree as a residue mod each node.
  Moving to a node's factor uses CRT on the halves of the evaluation coset:
  Z_C = Z_C₀·Z_C₁ with Z_C₀ ≡ 1 on C₁.

## Syndrome splitting

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
