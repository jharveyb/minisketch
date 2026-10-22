/-
Copyright (c) 2026 The minisketch developers
Distributed under the MIT software license, see the accompanying
file LICENSE or http://www.opensource.org/licenses/mit-license.php.
-/
import Mathlib

/-!
# Sketches as power sums

A sketch of a finite set `M` of nonzero field elements stores the power sums
`s_k = Σ_{m ∈ M} m^k`. These are the foundation claims C1, C3 and C4 of
`doc/spec/claims.md`.
-/

namespace Minisketch.Sketch

variable {K : Type*} [Field K]

/-- The `k`-th power sum of a finite set. -/
def psum (M : Finset K) (k : ℕ) : K := ∑ m ∈ M, m ^ k

/-- C1: sketches are linear. Adding the sketches of `A` and `B` gives the sketch
of their symmetric difference, because the common elements cancel in
characteristic 2. -/
theorem psum_symmDiff [DecidableEq K] [CharP K 2] (A B : Finset K) (k : ℕ) :
    psum A k + psum B k = psum (symmDiff A B) k := by
  unfold psum
  have hsub : A ∩ B ⊆ A ∪ B := Finset.inter_subset_left.trans Finset.subset_union_left
  rw [symmDiff_eq_sup_sdiff_inf, ← Finset.sum_union_inter (s₁ := A),
    ← Finset.sum_sdiff hsub (f := fun m => m ^ k)]
  simp only [Finset.sup_eq_union, Finset.inf_eq_inter]
  rw [add_assoc, CharTwo.add_self_eq_zero, add_zero]

/-- C4: in characteristic 2, `s_{2k} = s_k^2`. -/
theorem psum_two_mul [CharP K 2] (M : Finset K) (k : ℕ) :
    psum M (2 * k) = psum M k ^ 2 := by
  unfold psum
  rw [CharTwo.sum_sq]
  refine Finset.sum_congr rfl fun m _ => ?_
  ring

/-- C3, core form: a nonempty set `D` of nonzero elements has a nonzero power
sum `s_k` with `1 ≤ k ≤ |D|` (Vandermonde). -/
theorem exists_psum_ne_zero (D : Finset K) (hD : D.Nonempty) (h0 : (0 : K) ∉ D) :
    ∃ k, 1 ≤ k ∧ k ≤ D.card ∧ psum D k ≠ 0 := by
  by_contra hall
  push Not at hall
  set n := D.card with hn
  let e : Fin n ≃ D := D.equivFin.symm
  let f : Fin n → K := fun j => (e j : K)
  have hf : Function.Injective f := fun a b h => e.injective (Subtype.ext h)
  -- The transposed Vandermonde system (Σ_m m^i · m = s_{i+1} = 0) forces f = 0.
  have hw : (Matrix.vandermonde f).transpose.mulVec f = 0 := by
    funext i
    simp only [Matrix.mulVec, dotProduct, Matrix.transpose_apply, Matrix.vandermonde_apply,
      Pi.zero_apply]
    have hsum : ∑ j, f j ^ (i : ℕ) * f j = psum D (i + 1) := by
      unfold psum
      rw [← Finset.sum_coe_sort D, ← e.sum_comp]
      exact Finset.sum_congr rfl fun j _ => by rw [pow_succ]
    rw [hsum]
    exact hall (i + 1) (by omega) (by omega)
  have hdet : (Matrix.vandermonde f).transpose.det ≠ 0 := by
    rw [Matrix.det_transpose]
    exact Matrix.det_vandermonde_ne_zero_iff.mpr hf
  have hz := Matrix.eq_zero_of_mulVec_eq_zero hdet hw
  have hpos : 0 < n := Finset.card_pos.mpr hD
  have h0' : f ⟨0, hpos⟩ = 0 := congrFun hz _
  exact h0 (h0' ▸ (e ⟨0, hpos⟩).2)

/-- Weighted form of C3: if `Σ_{α ∈ D} w(α) α^k = 0` for `k = 1..|D|`, with `D`
a set of nonzero elements, then `w` vanishes on `D`. -/
theorem weight_eq_zero_of_sum_pow (D : Finset K) (h0 : (0 : K) ∉ D) (w : K → K)
    (h : ∀ k, 1 ≤ k → k ≤ D.card → ∑ α ∈ D, w α * α ^ k = 0) : ∀ α ∈ D, w α = 0 := by
  set n := D.card
  let e : Fin n ≃ D := D.equivFin.symm
  let f : Fin n → K := fun j => (e j : K)
  have hf : Function.Injective f := fun a b h => e.injective (Subtype.ext h)
  have hw : (Matrix.vandermonde f).transpose.mulVec (fun j => w (f j) * f j) = 0 := by
    funext i
    simp only [Matrix.mulVec, dotProduct, Matrix.transpose_apply, Matrix.vandermonde_apply,
      Pi.zero_apply]
    have hsum : ∑ j, f j ^ (i : ℕ) * (w (f j) * f j) = ∑ α ∈ D, w α * α ^ ((i : ℕ) + 1) := by
      rw [← Finset.sum_coe_sort D, ← e.sum_comp]
      exact Finset.sum_congr rfl fun j _ => by rw [pow_succ]; ring
    rw [hsum]
    exact h (i + 1) (by omega) (by omega)
  have hdet : (Matrix.vandermonde f).transpose.det ≠ 0 := by
    rw [Matrix.det_transpose]
    exact Matrix.det_vandermonde_ne_zero_iff.mpr hf
  have hz := Matrix.eq_zero_of_mulVec_eq_zero hdet hw
  intro α hα
  have hj := congrFun hz (e.symm ⟨α, hα⟩)
  simp only [f, Equiv.apply_symm_apply, Pi.zero_apply] at hj
  rcases mul_eq_zero.mp hj with h1 | h1
  · exact h1
  · exact absurd (h1 ▸ hα) h0

/-- C3: the power sums `s_1..s_{2n}` determine a set of at most `n` nonzero
elements. -/
theorem eq_of_psum_eq [DecidableEq K] [CharP K 2] {n : ℕ} (M M' : Finset K)
    (hM : M.card ≤ n) (hM' : M'.card ≤ n) (h0 : (0 : K) ∉ M) (h0' : (0 : K) ∉ M')
    (h : ∀ k, 1 ≤ k → k ≤ 2 * n → psum M k = psum M' k) : M = M' := by
  set D := symmDiff M M'
  by_contra hne
  have hD : D.Nonempty := by
    rw [Finset.nonempty_iff_ne_empty]
    intro he
    exact hne (symmDiff_eq_bot.mp he)
  have hD0 : (0 : K) ∉ D := by
    intro h0D
    rcases Finset.mem_symmDiff.mp h0D with ⟨h1, _⟩ | ⟨h1, _⟩
    · exact h0 h1
    · exact h0' h1
  have hcard : D.card ≤ 2 * n := by
    calc D.card ≤ (M ∪ M').card := Finset.card_le_card Finset.symmDiff_subset_union
      _ ≤ M.card + M'.card := Finset.card_union_le _ _
      _ ≤ 2 * n := by omega
  obtain ⟨k, hk1, hkD, hk⟩ := exists_psum_ne_zero D hD hD0
  apply hk
  rw [← psum_symmDiff, h k hk1 (hkD.trans hcard), CharTwo.add_self_eq_zero]

/-- With C4 the odd power sums `s_1, s_3, …, s_{2n-1}` (the stored sketch)
determine all of `s_1..s_{2n}`. -/
theorem psum_eq_of_odd [CharP K 2] {n : ℕ} (M M' : Finset K)
    (h : ∀ j, j < n → psum M (2 * j + 1) = psum M' (2 * j + 1)) :
    ∀ k, 1 ≤ k → k ≤ 2 * n → psum M k = psum M' k := by
  intro k
  induction k using Nat.strong_induction_on with
  | _ k ih =>
    intro hk1 hk2
    rcases Nat.even_or_odd' k with ⟨j, rfl | rfl⟩
    · rw [psum_two_mul, psum_two_mul, ih j (by omega) (by omega) (by omega)]
    · exact h j (by omega)

/-- C3 as the decoder uses it: two sets of at most `n` nonzero elements with the
same odd power sums `s_1, s_3, …, s_{2n-1}` are equal. -/
theorem eq_of_odd_psum_eq [DecidableEq K] [CharP K 2] {n : ℕ} (M M' : Finset K)
    (hM : M.card ≤ n) (hM' : M'.card ≤ n) (h0 : (0 : K) ∉ M) (h0' : (0 : K) ∉ M')
    (h : ∀ j, j < n → psum M (2 * j + 1) = psum M' (2 * j + 1)) : M = M' :=
  eq_of_psum_eq M M' hM hM' h0 h0' (psum_eq_of_odd M M' h)

open Polynomial in
theorem reverse_one_sub (a : K) (ha : a ≠ 0) : (1 - C a * X : K[X]).reverse = X - C a := by
  have hd : (1 - C a * X : K[X]).natDegree = 1 := by
    rw [natDegree_sub_eq_right_of_natDegree_lt] <;> rw [natDegree_C_mul_X a ha]; simp
  ext i
  rw [reverse, hd, coeff_reflect]
  rcases i with _ | _ | i
  · rw [revAt_le (by norm_num)]; simp [coeff_X, coeff_C, coeff_one]
  · rw [revAt_le (by norm_num)]; simp [coeff_X, coeff_C, coeff_one]
  · rw [revAt_eq_self_of_lt (by omega)]
    simp [coeff_X, coeff_one]

open Polynomial in
/-- C6: for nonzero elements, reversing the locator `Π_{m ∈ M} (1 - m x)` gives
`Π_{m ∈ M} (x - m)`, whose roots are the elements themselves. -/
theorem reverse_locator [DecidableEq K] (M : Finset K) (h0 : (0 : K) ∉ M) :
    (∏ m ∈ M, (1 - C m * X)).reverse = ∏ m ∈ M, (X - C m) := by
  induction M using Finset.induction_on with
  | empty => simp [reverse]
  | insert a M ha ih =>
    rw [Finset.prod_insert ha, Finset.prod_insert ha, reverse_mul_of_domain,
      ih fun h => h0 (Finset.mem_insert_of_mem h),
      reverse_one_sub a fun h => h0 (h ▸ Finset.mem_insert_self a M)]

end Minisketch.Sketch