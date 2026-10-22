/-
Copyright (c) 2026 The minisketch developers
Distributed under the MIT software license, see the accompanying
file LICENSE or http://www.opensource.org/licenses/mit-license.php.
-/
import Minisketch.Sketch

/-!
# Power-sum splits

`PowerSumSplit` (`src/sketch_impl.h`) replaces the gcd that splits a node `g`
of the root-finding tree. It computes the power sums of one side directly from
the node's own power sums, and recovers that side's polynomial with
Berlekamp–Massey. This file proves the algebra behind it: claim C35 of
`doc/spec/claims.md`.

`R` is the set of roots of the node polynomial `g`; `t` is the split trace,
which takes values in `{0, 1}` on `R`; `A = {α ∈ R : t(α) = 1}`.
-/

namespace Minisketch.PowerSumSplit

open Polynomial Minisketch.Sketch

variable {K : Type*} [Field K]

/-- Applying a polynomial `f = Σ f_i X^i` as a filter to the power sums of `A`
gives `Σ_{α ∈ A} α^m f(α)`. -/
theorem filter_psum (A : Finset K) (f : K[X]) (m : ℕ) :
    ∑ i ∈ Finset.range (f.natDegree + 1), f.coeff i * psum A (m + i) =
      ∑ α ∈ A, α ^ m * f.eval α := by
  unfold psum
  simp_rw [eval_eq_sum_range, Finset.mul_sum]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun α _ => Finset.sum_congr rfl fun i _ => ?_
  rw [pow_add]
  ring

/-- C35, the split formula: `P_k = Σ_j t_j Q_{j+k}`. The power sums of the side
`A` are a correlation of the trace's coefficients with the node's power sums. -/
theorem psum_side [DecidableEq K] (R : Finset K) (t : K[X])
    (h01 : ∀ α ∈ R, t.eval α = 0 ∨ t.eval α = 1) (k : ℕ) :
    psum (R.filter fun α => t.eval α = 1) k =
      ∑ j ∈ Finset.range (t.natDegree + 1), t.coeff j * psum R (k + j) := by
  rw [filter_psum, psum, Finset.sum_filter]
  refine Finset.sum_congr rfl fun α hα => ?_
  rcases h01 α hα with h | h
  · rw [ite_eq_right (by rw [h]; exact zero_ne_one), h, mul_zero]
  · rw [ite_eq_left h, h, mul_one]

/-- C35, uniqueness of the recovered side: let `A, B ⊆ R` with `0 ∉ R` and
`|R| ≤ 2M`. If `f_B = Π_{β ∈ B} (X - β)` generates the power sums of `A` (the
LFSR condition on `P_1..P_{2M}`, written with the reversed connection
polynomial) and `|B| ≤ |A|` (BM returns a shortest LFSR, claim C5, and `Λ_A`
is one of length `|A|`), then `B = A`. So the side that passes `CosetDivide`
inherits correct power sums. -/
theorem side_unique [DecidableEq K] (R A B : Finset K) (hA : A ⊆ R) (hB : B ⊆ R)
    (h0 : (0 : K) ∉ R) (M : ℕ) (hR : R.card ≤ 2 * M) (hBA : B.card ≤ A.card)
    (hgen : ∀ m, 1 ≤ m → m + B.card ≤ 2 * M →
      ∑ i ∈ Finset.range ((∏ β ∈ B, (X - C β)).natDegree + 1),
        (∏ β ∈ B, (X - C β)).coeff i * psum A (m + i) = 0) :
    A = B := by
  set f := ∏ β ∈ B, (X - C β)
  have hfB : ∀ α, f.eval α = 0 ↔ α ∈ B := by
    intro α
    simp [f, eval_prod, Finset.prod_eq_zero_iff, sub_eq_zero]
  -- Σ_{α ∈ A \ B} α^m f(α) = 0 for m = 1..2M - |B|, which covers |A \ B| values.
  have hsub : ∀ α ∈ A \ B, f.eval α = 0 := by
    refine weight_eq_zero_of_sum_pow (A \ B) (fun h => h0 (hA (Finset.mem_sdiff.mp h).1))
      (fun α => f.eval α) fun m hm1 hm2 => ?_
    have hcard : (A \ B).card + B.card ≤ 2 * M := by
      have : (A \ B).card ≤ (R \ B).card := Finset.card_le_card (Finset.sdiff_subset_sdiff hA le_rfl)
      rw [Finset.card_sdiff_of_subset hB] at this
      have := Finset.card_le_card hB
      omega
    have h := hgen m hm1 (by omega)
    rw [filter_psum, ← Finset.sum_sdiff (Finset.inter_subset_left : A ∩ B ⊆ A)] at h
    have hint : ∑ α ∈ A ∩ B, α ^ m * f.eval α = 0 :=
      Finset.sum_eq_zero fun α hα => by
        rw [(hfB α).mpr (Finset.mem_inter.mp hα).2, mul_zero]
    rw [hint, add_zero, show A \ (A ∩ B) = A \ B by ext; simp] at h
    rw [← h]
    exact Finset.sum_congr rfl fun α _ => mul_comm _ _
  have hAB : A ⊆ B := by
    intro α hα
    by_contra hnot
    exact absurd ((hfB α).mp (hsub α (Finset.mem_sdiff.mpr ⟨hα, hnot⟩))) hnot
  exact Finset.eq_of_subset_of_card_le hAB hBA

/-- C35, empty side: a BM result of length 1 (`f = 1`, all power sums of `A`
vanish) means `A = ∅`. -/
theorem side_empty [DecidableEq K] (R A : Finset K) (hA : A ⊆ R) (h0 : (0 : K) ∉ R)
    (M : ℕ) (hR : R.card ≤ 2 * M) (hz : ∀ k, 1 ≤ k → k ≤ 2 * M → psum A k = 0) : A = ∅ := by
  by_contra hne
  obtain ⟨k, hk1, hk2, hk⟩ := exists_psum_ne_zero A (Finset.nonempty_iff_ne_empty.mpr hne)
    (fun h => h0 (hA h))
  exact hk (hz k hk1 (hk2.trans ((Finset.card_le_card hA).trans hR)))

/-- C35, the complement: the other side's power sums are `Q_k - P_k`. -/
theorem psum_complement [DecidableEq K] (R : Finset K) (p : K → Prop) [DecidablePred p]
    (k : ℕ) : psum (R.filter fun α => ¬p α) k = psum R k - psum (R.filter p) k := by
  unfold psum
  rw [eq_sub_iff_add_eq, add_comm, Finset.sum_filter_add_sum_filter_not]

end Minisketch.PowerSumSplit
