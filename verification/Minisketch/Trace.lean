/-
Copyright (c) 2026 The minisketch developers
Distributed under the MIT software license, see the accompanying
file LICENSE or http://www.opensource.org/licenses/mit-license.php.
-/
import Mathlib

/-!
# Trace splitting

Root finding uses the Berlekamp trace algorithm: the trace `Tr(βx) mod g` takes
values in `{0, 1}` at the roots of `g`, so its gcd with `g` splits off the roots
with trace 0. This file proves claims C20–C23, C26 and C28 of
`doc/spec/claims.md` for a field `K` with `2^b` elements.
-/

namespace Minisketch.Trace

open Polynomial

variable {K : Type*} [Field K] [CharP K 2]

/-- The absolute trace `a + a^2 + a^4 + … + a^(2^(b-1))`. -/
def tr (b : ℕ) (a : K) : K := ∑ i ∈ Finset.range b, a ^ 2 ^ i

theorem frob_add (a c : K) (n : ℕ) : (a + c) ^ 2 ^ n = a ^ 2 ^ n + c ^ 2 ^ n :=
  add_pow_char_pow a c 2 n

/-- `(Σ_{i<n} y^(2^i))^2 + Σ_{i<n} y^(2^i) = y^(2^n) + y` in any commutative ring
of characteristic 2. -/
theorem sum_sq_add {R : Type*} [CommRing R] [CharP R 2] (y : R) (n : ℕ) :
    (∑ i ∈ Finset.range n, y ^ 2 ^ i) ^ 2 + ∑ i ∈ Finset.range n, y ^ 2 ^ i = y ^ 2 ^ n + y := by
  induction n with
  | zero => simp [CharTwo.add_self_eq_zero]
  | succ n ih =>
    simp only [Finset.sum_range_succ] at ih ⊢
    rw [CharTwo.add_sq, add_add_add_comm, ih, ← pow_mul, ← pow_succ]
    ring_nf
    rw [show y ^ 2 ^ n * 2 = y ^ 2 ^ n + y ^ 2 ^ n by ring, CharTwo.add_self_eq_zero]
    ring

/-- Tr is additive. -/
theorem tr_add (b : ℕ) (a c : K) : tr b (a + c) = tr b a + tr b c := by
  simp only [tr, frob_add, Finset.sum_add_distrib]

/-- `Tr(a)^2 + Tr(a) = a^(2^b) + a`; in particular `Tr(a)^2 = Tr(a)` when
`a^(2^b) = a`. -/
theorem tr_sq_add (b : ℕ) (a : K) : tr b a ^ 2 + tr b a = a ^ 2 ^ b + a :=
  sum_sq_add a b

/-- C20: `Tr(a) ∈ {0, 1}` when `a^(2^b) = a` (every element of `GF(2^b)`). -/
theorem tr_mem (b : ℕ) (a : K) (ha : a ^ 2 ^ b = a) : tr b a = 0 ∨ tr b a = 1 := by
  have h := tr_sq_add b a
  rw [ha, CharTwo.add_self_eq_zero] at h
  have : tr b a * (tr b a - 1) = 0 := by
    rw [mul_sub, mul_one, ← sq, sub_eq_add_neg, CharTwo.neg_eq]; exact h
  rcases mul_eq_zero.mp this with h0 | h1
  · exact Or.inl h0
  · exact Or.inr (sub_eq_zero.mp h1)

/-- C28 (and the bucket trick of C27): for `z` in the subfield `GF(2^s)`
(`z^(2^s) = z`) and `b = s t`, `Tr(z a) = Σ_{r < s} (z T(a))^(2^r)` with the
relative trace `T(a) = Σ_{i < t} a^(2^(s i))` into `GF(2^s)`. So one relative
trace serves the `s` split parameters `z^k β` at the cost of `s - 1` squarings. -/
theorem tr_tower (s t : ℕ) (z a : K) (hz : z ^ 2 ^ s = z) :
    tr (s * t) (z * a) =
      ∑ r ∈ Finset.range s, (z * ∑ i ∈ Finset.range t, a ^ 2 ^ (s * i)) ^ 2 ^ r := by
  have hzk : ∀ i r, z ^ 2 ^ (s * i + r) = z ^ 2 ^ r := by
    intro i r
    induction i with
    | zero => simp
    | succ i ih =>
      rw [show s * (i + 1) + r = s + (s * i + r) by ring, pow_add, pow_mul, hz, ih]
  induction t with
  | zero => simp [tr]
  | succ t ih =>
    rw [tr, show s * (t + 1) = s * t + s by ring, Finset.sum_range_add, ← tr, ih,
      ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl fun r _ => ?_
    rw [Finset.sum_range_succ, mul_add, frob_add, mul_pow, mul_pow, mul_pow, hzk t r,
      ← pow_mul, ← pow_add]

/-! ### The trace as a polynomial -/

/-- `Tr(βX) = Σ_{i < b} (βX)^(2^i)`. -/
noncomputable def trPoly (b : ℕ) (β : K) : K[X] := ∑ i ∈ Finset.range b, (C β * X) ^ 2 ^ i

omit [CharP K 2] in
theorem trPoly_eval (b : ℕ) (β α : K) : (trPoly b β).eval α = tr b (β * α) := by
  simp [trPoly, tr, eval_finsetSum]

omit [CharP K 2] in
/-- C26, the orbit form: `Tr(βX) = Σ_i β^(2^i) X^(2^i)`, so its residue mod `f` is
`Σ_i β^(2^i) X_i` with `X_i = X^(2^i) mod f`. -/
theorem trPoly_orbit (b : ℕ) (β : K) (f : K[X]) :
    trPoly b β %ₘ f = ∑ i ∈ Finset.range b, C (β ^ 2 ^ i) * (X ^ 2 ^ i %ₘ f) := by
  rw [trPoly]
  simp_rw [mul_pow, ← C_pow]
  rw [← modByMonicHom_apply, map_sum]
  refine Finset.sum_congr rfl fun i _ => ?_
  rw [modByMonicHom_apply, ← smul_eq_C_mul, ← smul_eq_C_mul, smul_modByMonic]

omit [CharP K 2] in
/-- C26: reducing mod `f` and then mod a factor `g` is reducing mod `g`. -/
theorem modByMonic_of_dvd (h f g : K[X]) (hg : g.Monic) (hgf : g ∣ f) :
    (h %ₘ f) %ₘ g = h %ₘ g := by
  rw [modByMonic_eq_sub_mul_div h f, sub_modByMonic,
    (modByMonic_eq_zero_iff_dvd hg).mpr (dvd_mul_of_dvd_left hgf _), sub_zero]

/-- C21, the fail-fast test: `T^2 + T = β (X^(2^b) + X)` for `T = Tr(βX)`, given
`β^(2^b) = β`. -/
theorem trPoly_sq_add (b : ℕ) (β : K) (hβ : β ^ 2 ^ b = β) :
    trPoly b β ^ 2 + trPoly b β = C β * (X ^ 2 ^ b + X) := by
  have key := sum_sq_add (C β * X) b
  rw [trPoly, key, mul_pow, ← C_pow, hβ, mul_add]

/-- C23: a root `α` of `g` in any extension `L` at which `T ≡ Tr(βX) (mod g)`
vanishes lies in `K`'s Frobenius fixed field: `α^(2^b) = α`. So the gcd side of
a split only holds roots in `GF(2^b)`, never an irreducible factor of higher
degree. -/
theorem gcd_side_fixed (b : ℕ) (β : K) (hβ : β ^ 2 ^ b = β) (hβ0 : β ≠ 0)
    {L : Type*} [Field L] [Algebra K L] (g T : K[X]) (hgT : g ∣ trPoly b β - T) (α : L)
    (hg : aeval α g = 0) (hT : aeval α T = 0) : α ^ 2 ^ b = α := by
  have hTr : aeval α (trPoly b β) = 0 := by
    obtain ⟨c, hc⟩ := hgT
    have := congrArg (aeval α) hc
    rw [map_sub, map_mul, hg, zero_mul, hT, sub_zero] at this
    exact this
  have h := congrArg (aeval α) (trPoly_sq_add b β hβ)
  rw [map_add, map_pow, hTr, map_mul, aeval_C, map_add, map_pow, aeval_X] at h
  simp only [ne_eq, OfNat.ofNat_ne_zero, not_false_eq_true, zero_pow, add_zero] at h
  have hβL : algebraMap K L β ≠ 0 := (map_ne_zero_iff _ (algebraMap K L).injective).mpr hβ0
  have hsum : α ^ 2 ^ b + α = 0 := by
    rcases mul_eq_zero.mp h.symm with h1 | h1
    · exact absurd h1 hβL
    · exact h1
  have : CharP L 2 := charP_of_injective_algebraMap (algebraMap K L).injective 2
  rw [← sub_eq_zero, sub_eq_add_neg, CharTwo.neg_eq]
  exact hsum

omit [CharP K 2] in
/-- C20: with `T ≡ Tr(βX) (mod g)` and `g = Π_{α ∈ R} (X - α)`, the roots with
trace 0 divide `T` and the roots with trace 1 are coprime to it. So
`gcd(g, T) = Π_{α ∈ R, Tr(βα) = 0} (X - α)`. -/
theorem trace_split [DecidableEq K] (b : ℕ) (β : K) (R : Finset K)
    (T : K[X])
    (hT : (∏ α ∈ R, (X - C α)) ∣ trPoly b β - T) :
    (∏ α ∈ R.filter (fun α => tr b (β * α) = 0), (X - C α)) ∣ T ∧
      IsCoprime (∏ α ∈ R.filter (fun α => tr b (β * α) ≠ 0), (X - C α)) T := by
  have hval : ∀ α ∈ R, T.eval α = tr b (β * α) := by
    intro α hα
    have hdiv : (X - C α) ∣ trPoly b β - T :=
      (Finset.dvd_prod_of_mem (fun α => X - C α) hα).trans hT
    rw [dvd_iff_isRoot, IsRoot, eval_sub, trPoly_eval, sub_eq_zero] at hdiv
    exact hdiv.symm
  constructor
  · apply Finset.prod_dvd_of_coprime
    · intro w _ w' _ hne
      exact pairwise_coprime_X_sub_C Function.injective_id hne
    · intro α hα
      rw [Finset.mem_filter] at hα
      rw [dvd_iff_isRoot, IsRoot, hval α hα.1, hα.2]
  · apply IsCoprime.prod_left
    intro α hα
    rw [Finset.mem_filter] at hα
    rw [(irreducible_X_sub_C α).coprime_iff_not_dvd, dvd_iff_isRoot, IsRoot, hval α hα.1]
    exact hα.2

/-! ### The full-split test -/

omit [CharP K 2] in
/-- `X^q - X = Π_{a ∈ K} (X - a)` for a finite field with `q` elements. -/
theorem X_pow_card_sub_X_eq_prod [Fintype K] :
    (X ^ Fintype.card K - X : K[X]) = ∏ a : K, (X - C a) := by
  have h1 : 1 < Fintype.card K := Fintype.one_lt_card
  have hdeg := FiniteField.X_pow_card_sub_X_natDegree_eq K h1
  have hmon : (X ^ Fintype.card K - X : K[X]).Monic := by
    apply (monic_X_pow _).sub_of_left
    rw [degree_X, degree_X_pow]
    exact_mod_cast h1
  rw [← prod_multiset_X_sub_C_of_monic_of_roots_card_eq hmon
    (by rw [FiniteField.roots_X_pow_card_sub_X, hdeg]; rfl),
    FiniteField.roots_X_pow_card_sub_X]
  rfl

/-- C22: a monic `f` divides `X^(2^b) - X` (equivalently `X^(2^b) ≡ X mod f`) iff it
is a product of distinct linear factors `Π_{α ∈ R} (X - α)`. -/
theorem dvd_frobenius_iff [Fintype K] [DecidableEq K] (b : ℕ) (hq : Fintype.card K = 2 ^ b)
    (f : K[X]) (hf : f.Monic) :
    f ∣ X ^ 2 ^ b - X ↔ ∃ R : Finset K, f = ∏ α ∈ R, (X - C α) := by
  rw [← hq, X_pow_card_sub_X_eq_prod]
  constructor
  · intro hdvd
    have hne : (∏ a : K, (X - C a)) ≠ 0 :=
      Finset.prod_ne_zero_iff.mpr fun a _ => X_sub_C_ne_zero a
    have hsplit : (∏ a : K, (X - C a)).Splits := by
      rw [splits_iff_card_roots, roots_prod_X_sub_C, natDegree_prod_of_monic _ _
        fun a _ => monic_X_sub_C a]
      simp
    have hfs : f.Splits := hsplit.of_dvd hne hdvd
    have hsep : f.Separable := by
      refine Separable.of_dvd ?_ hdvd
      rw [← X_pow_card_sub_X_eq_prod, hq]
      exact galois_poly_separable 2 (2 ^ b) (dvd_pow_self 2 (by
        rintro rfl; simp at hq; exact absurd hq (by have := Fintype.one_lt_card (α := K); omega)))
    refine ⟨f.roots.toFinset, ?_⟩
    rw [Finset.prod_eq_multiset_prod, Multiset.toFinset_val, (nodup_roots hsep).dedup,
      prod_multiset_X_sub_C_of_monic_of_roots_card_eq hf (splits_iff_card_roots.mp hfs)]
  · rintro ⟨R, rfl⟩
    exact Finset.prod_dvd_prod_of_subset R Finset.univ _ (Finset.subset_univ R)

omit [CharP K 2] in
/-- The test as the code runs it: `X^(2^b) mod f = X mod f`. -/
theorem frobenius_mod_iff (b : ℕ) (f : K[X]) (hf : f.Monic) :
    X ^ 2 ^ b %ₘ f = X %ₘ f ↔ f ∣ X ^ 2 ^ b - X := by
  rw [← modByMonic_eq_zero_iff_dvd hf, sub_modByMonic, sub_eq_zero]

end Minisketch.Trace
