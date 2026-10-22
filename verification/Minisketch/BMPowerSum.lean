/-
Copyright (c) 2026 The minisketch developers
Distributed under the MIT software license, see the accompanying
file LICENSE or http://www.opensource.org/licenses/mit-license.php.
-/
import Minisketch.BM
import Minisketch.PowerSumSplit

/-!
# BM on power sums

The locator `Π_{α ∈ A} (1 - α x)` generates the power sums of `A` with length
`|A|`, so by Massey's theorem (C5) BM's result on them is no longer than `|A|`.
This is the minimality hypothesis of `PowerSumSplit.side_unique` (C35).
-/

namespace Minisketch.BMPowerSum

open Polynomial Minisketch.Sketch Minisketch.BM Minisketch.PowerSumSplit

variable {K : Type*} [Field K] [DecidableEq K]

/-- Power sums as a sequence, with the unused `s_0` set to 1 (BM's convention). -/
def seq (A : Finset K) : ℕ → K := Function.update (psum A) 0 1

/-- The locator's coefficients `c_i = f_A.coeff (|A| - i)`, `f_A = Π (X - α)`. -/
noncomputable def locator (A : Finset K) : ℕ → K :=
  fun i => (∏ α ∈ A, (X - C α)).coeff (A.card - i)

omit [DecidableEq K] in
theorem locator_zero (A : Finset K) : locator A 0 = 1 := by
  have hm : (∏ α ∈ A, (X - C α)).Monic := monic_prod_of_monic _ _ fun α _ => monic_X_sub_C α
  have hd : (∏ α ∈ A, (X - C α)).natDegree = A.card := by
    rw [natDegree_prod_of_monic _ _ fun α _ => monic_X_sub_C α]; simp
  simp only [locator, Nat.sub_zero]
  rw [← hd]; exact hm.leadingCoeff

omit [DecidableEq K] in
/-- The locator generates the power sums of `A` with length `|A|`, for every
number of terms. -/
theorem locator_gen (A : Finset K) (n : ℕ) : Gen (seq A) (locator A) A.card n := by
  intro k hk1 _
  set f := ∏ α ∈ A, (X - C α)
  have hd : f.natDegree = A.card := by
    rw [natDegree_prod_of_monic _ _ fun α _ => monic_X_sub_C α]; simp
  unfold disc locator
  -- Reindex i ↦ |A| - i and use Σ_i f_i psum (m + i) = Σ_α α^m f(α), m = k - |A|.
  rw [← Finset.sum_range_reflect]
  have h := filter_psum A f (k - A.card)
  rw [hd] at h
  have hz : ∑ α ∈ A, α ^ (k - A.card) * f.eval α = 0 :=
    Finset.sum_eq_zero fun α hα => by
      rw [eval_prod, Finset.prod_eq_zero hα (by simp), mul_zero]
  rw [hz] at h
  rw [← h]
  refine Finset.sum_congr rfl fun i hi => ?_
  rw [Finset.mem_range] at hi
  rw [show A.card - (A.card + 1 - 1 - i) = i by omega, seq, Function.update_of_ne (by omega),
    show k - (A.card + 1 - 1 - i) = k - A.card + i by omega]

/-- C35's minimality hypothesis: BM on the power sums of `A` returns a length
of at most `|A|`. -/
theorem bm_len_le_card (A : Finset K) (n : ℕ) : (bm (seq A) n).L ≤ A.card :=
  (bm_minimal (seq A) (by simp [seq]) n).2.2 _ _ (locator_zero A) (locator_gen A n)

end Minisketch.BMPowerSum
