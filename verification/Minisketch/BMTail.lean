/-
Copyright (c) 2026 The minisketch developers
Distributed under the MIT software license, see the accompanying
file LICENSE or http://www.opensource.org/licenses/mit-license.php.
-/
import Mathlib

/-!
# The odd-only tail check of Berlekamp–Massey

`BerlekampMasseyOdd` (`src/sketch_impl.h`) stops once its connection
polynomial `Λ` has been stable for a while, and then checks the remaining
syndromes in bulk with `BerlekampMasseyVerifyTail`. That check only looks at
the stored odd syndromes `u_j = s_{2j+1}` and the recurrence with `Λ²`. This file
proves that it is exactly the full LFSR check: claims C13, C14 and C15 of
`doc/spec/claims.md`.

Everything is over power series in characteristic 2. `S = Σ_{k ≥ 1} s_k X^k` is
the syndrome series of a Frobenius-structured sequence (`s_{2k} = s_k^2` for
`k ≥ 1`), and `C = Λ · S`. For `k > L = deg Λ`, `C_k` is the LFSR discrepancy
`s_k + Σ_{i=1}^{L} Λ_i s_{k-i}`.
-/

namespace Minisketch.BMTail

open PowerSeries

variable {K : Type*} [Field K] [CharP K 2]

/-- Squaring a power series in characteristic 2 squares the coefficients and
doubles the exponents. -/
theorem coeff_sq (S : K⟦X⟧) (n : ℕ) :
    coeff n (S ^ 2) = if Even n then coeff (n / 2) S ^ 2 else 0 := by
  rw [sq, coeff_mul, ← Finset.sum_filter_add_sum_filter_not _ (fun p : ℕ × ℕ => p.1 = p.2)]
  have hoff : ∑ p ∈ (Finset.HasAntidiagonal.antidiagonal n).filter (fun p : ℕ × ℕ => ¬p.1 = p.2),
      coeff p.1 S * coeff p.2 S = 0 := by
    refine Finset.sum_involution (fun p _ => p.swap) ?_ ?_ ?_ ?_
    · intro p _
      rw [Prod.fst_swap, Prod.snd_swap, mul_comm, CharTwo.add_self_eq_zero]
    · intro p hp _ h
      rw [Finset.mem_filter] at hp
      exact hp.2 (by simpa using (congrArg Prod.fst h).symm)
    · intro p hp
      rw [Finset.mem_filter, Finset.HasAntidiagonal.mem_antidiagonal] at hp ⊢
      exact ⟨by rw [Prod.fst_swap, Prod.snd_swap, add_comm]; exact hp.1,
        fun h => hp.2 h.symm⟩
    · intro p _
      exact Prod.swap_swap p
  rw [hoff, add_zero]
  split_ifs with he
  · obtain ⟨m, rfl⟩ := he
    have hdiag : (Finset.HasAntidiagonal.antidiagonal (m + m)).filter (fun p : ℕ × ℕ => p.1 = p.2) = {(m, m)} := by
      ext ⟨a, b⟩
      simp only [Finset.mem_filter, Finset.HasAntidiagonal.mem_antidiagonal, Finset.mem_singleton, Prod.mk.injEq]
      omega
    rw [hdiag, Finset.sum_singleton, show (m + m) / 2 = m by omega, sq]
  · have hdiag : (Finset.HasAntidiagonal.antidiagonal n).filter (fun p : ℕ × ℕ => p.1 = p.2) = ∅ := by
      ext ⟨a, b⟩
      simp only [Finset.mem_filter, Finset.HasAntidiagonal.mem_antidiagonal, Finset.notMem_empty, iff_false,
        not_and]
      intro hab h
      exact he ⟨a, by omega⟩
    rw [hdiag, Finset.sum_empty]

/-- The syndrome series `Σ_{k ≥ 1} s_k X^k`. -/
noncomputable def series (s : ℕ → K) : K⟦X⟧ := mk fun k => if k = 0 then 0 else s k

/-- Its odd part `Σ_j s_{2j+1} X^{2j+1}`, the stored sketch. -/
noncomputable def oddPart (s : ℕ → K) : K⟦X⟧ := mk fun k => if Odd k then s k else 0

/-- A sequence is Frobenius-structured if `s_{2k} = s_k^2` for `k ≥ 1`. -/
def Frobenius (s : ℕ → K) : Prop := ∀ k, 1 ≤ k → s (2 * k) = s k ^ 2

/-- C13, first half: a Frobenius-structured series is its odd part plus its
own square. -/
theorem series_eq (s : ℕ → K) (hs : Frobenius s) :
    series s = oddPart s + series s ^ 2 := by
  ext n
  rw [map_add, coeff_sq, series, oddPart, coeff_mk, coeff_mk, coeff_mk]
  rcases Nat.even_or_odd' n with ⟨m, rfl | rfl⟩
  · have hno : ¬Odd (2 * m) := Nat.not_odd_iff_even.mpr (even_two_mul m)
    rw [ite_eq_right hno, ite_eq_left (even_two_mul m), zero_add, show 2 * m / 2 = m by omega]
    rcases Nat.eq_zero_or_pos m with rfl | hm
    · simp
    · rw [ite_eq_right (by omega), ite_eq_right (by omega), hs m hm]
  · have ho : Odd (2 * m + 1) := odd_two_mul_add_one m
    rw [ite_eq_right (by omega), ite_eq_left ho, ite_eq_right (Nat.not_even_iff_odd.mpr ho), add_zero]

/-- C13: for `C = Λ · S`, `Λ · C = C² + Λ² · S_odd`. -/
theorem tail_identity (s : ℕ → K) (hs : Frobenius s) (Λ : K⟦X⟧) :
    Λ * (Λ * series s) = (Λ * series s) ^ 2 + Λ ^ 2 * oddPart s := by
  conv_lhs => rw [series_eq s hs]
  ring

/-- `Λ² · S_odd` has no even coefficients: `Λ²` only has even exponents. -/
theorem coeff_sq_mul_odd_even (s : ℕ → K) (Λ : K⟦X⟧) (j : ℕ) :
    coeff (2 * j) (Λ ^ 2 * oddPart s) = 0 := by
  rw [coeff_mul]
  refine Finset.sum_eq_zero fun p hp => ?_
  rw [Finset.HasAntidiagonal.mem_antidiagonal] at hp
  rw [coeff_sq, oddPart, coeff_mk]
  by_cases he : Even p.1
  · have : ¬Odd p.2 := by
      intro ho
      obtain ⟨a, ha⟩ := he
      obtain ⟨b, hb⟩ := ho
      omega
    rw [ite_eq_right this, mul_zero]
  · rw [ite_eq_right he, zero_mul]

/-- C14: `C_{2j} = C_j² + Σ_{i ≥ 1} Λ_i C_{2j-i}` when `Λ_0 = 1`. -/
theorem even_coeff (s : ℕ → K) (hs : Frobenius s) (Λ : K⟦X⟧) (h0 : constantCoeff Λ = 1)
    (j : ℕ) :
    coeff (2 * j) (Λ * series s) = coeff j (Λ * series s) ^ 2 +
      ∑ p ∈ (Finset.HasAntidiagonal.antidiagonal (2 * j)).filter (fun p : ℕ × ℕ => p.1 ≠ 0),
        coeff p.1 Λ * coeff p.2 (Λ * series s) := by
  have h := congrArg (coeff (2 * j)) (tail_identity s hs Λ)
  rw [map_add, coeff_sq_mul_odd_even, add_zero, coeff_sq, ite_eq_left (even_two_mul j),
    show 2 * j / 2 = j by omega, coeff_mul,
    ← Finset.sum_filter_add_sum_filter_not _ (fun p : ℕ × ℕ => p.1 ≠ 0)] at h
  have hz : ∑ p ∈ (Finset.HasAntidiagonal.antidiagonal (2 * j)).filter (fun p : ℕ × ℕ => ¬p.1 ≠ 0),
      coeff p.1 Λ * coeff p.2 (Λ * series s) = coeff (2 * j) (Λ * series s) := by
    have : (Finset.HasAntidiagonal.antidiagonal (2 * j)).filter (fun p : ℕ × ℕ => ¬p.1 ≠ 0) = {(0, 2 * j)} := by
      ext ⟨a, b⟩
      simp only [Finset.mem_filter, Finset.HasAntidiagonal.mem_antidiagonal, Finset.mem_singleton, Prod.mk.injEq,
        ne_eq, not_not]
      omega
    rw [this, Finset.sum_singleton, coeff_zero_eq_constantCoeff_apply, h0, one_mul]
  rw [hz] at h
  rw [← h, add_comm (∑ p ∈ _, _) (coeff (2 * j) _), add_assoc, CharTwo.add_self_eq_zero, add_zero]

/-- C14 at odd positions: `C_{2j+1} + Σ_{i ≥ 1} Λ_i C_{2j+1-i}` is the odd
recurrence residual `[Λ² S_odd]_{2j+1} = u_j + Σ_{a ≥ 1} Λ_a² u_{j-a}`. -/
theorem odd_coeff (s : ℕ → K) (hs : Frobenius s) (Λ : K⟦X⟧) (h0 : constantCoeff Λ = 1)
    (j : ℕ) :
    coeff (2 * j + 1) (Λ * series s) +
      ∑ p ∈ (Finset.HasAntidiagonal.antidiagonal (2 * j + 1)).filter (fun p : ℕ × ℕ => p.1 ≠ 0),
        coeff p.1 Λ * coeff p.2 (Λ * series s) =
      coeff (2 * j + 1) (Λ ^ 2 * oddPart s) := by
  have h := congrArg (coeff (2 * j + 1)) (tail_identity s hs Λ)
  rw [map_add, coeff_sq, ite_eq_right (Nat.not_even_iff_odd.mpr (odd_two_mul_add_one j)), zero_add,
    coeff_mul, ← Finset.sum_filter_add_sum_filter_not _ (fun p : ℕ × ℕ => p.1 ≠ 0)] at h
  have hz : ∑ p ∈ (Finset.HasAntidiagonal.antidiagonal (2 * j + 1)).filter (fun p : ℕ × ℕ => ¬p.1 ≠ 0),
      coeff p.1 Λ * coeff p.2 (Λ * series s) = coeff (2 * j + 1) (Λ * series s) := by
    have : (Finset.HasAntidiagonal.antidiagonal (2 * j + 1)).filter (fun p : ℕ × ℕ => ¬p.1 ≠ 0) =
        {(0, 2 * j + 1)} := by
      ext ⟨a, b⟩
      simp only [Finset.mem_filter, Finset.HasAntidiagonal.mem_antidiagonal, Finset.mem_singleton, Prod.mk.injEq,
        ne_eq, not_not]
      omega
    rw [this, Finset.sum_singleton, coeff_zero_eq_constantCoeff_apply, h0, one_mul]
  rw [hz] at h
  rw [← h, add_comm]

omit [CharP K 2] in
/-- The sum `Σ_{i ≥ 1} Λ_i C_{k-i}` vanishes when `Λ` has degree at most `L` and
`C` vanishes on `(L, k)`, provided `k > 2L`. -/
theorem shifted_sum_zero (Λ C : K⟦X⟧) (L k : ℕ) (hΛ : ∀ i, L < i → coeff i Λ = 0)
    (hk : 2 * L < k) (hC : ∀ m, L < m → m < k → coeff m C = 0) :
    ∑ p ∈ (Finset.HasAntidiagonal.antidiagonal k).filter (fun p : ℕ × ℕ => p.1 ≠ 0),
      coeff p.1 Λ * coeff p.2 C = 0 := by
  refine Finset.sum_eq_zero fun p hp => ?_
  rw [Finset.mem_filter, Finset.HasAntidiagonal.mem_antidiagonal] at hp
  by_cases hL : L < p.1
  · rw [hΛ _ hL, zero_mul]
  · rw [hC _ (by omega) (by omega), mul_zero]

/-- C15: the odd-only tail check is exact. Let `Λ` have `Λ_0 = 1` and degree at
most `L`, and suppose the discrepancies `C_k` already vanish for `L < k ≤ P`,
with `P ≥ 2L`. Then they vanish for all `L < k ≤ N` iff the odd recurrence
residuals `[Λ² S_odd]_k` vanish at every odd `k` in `(P, N]`. -/
theorem tail_check_iff (s : ℕ → K) (hs : Frobenius s) (Λ : K⟦X⟧) (h0 : constantCoeff Λ = 1)
    (L P N : ℕ) (hΛ : ∀ i, L < i → coeff i Λ = 0) (hP : 2 * L ≤ P)
    (hpre : ∀ k, L < k → k ≤ P → coeff k (Λ * series s) = 0) :
    (∀ k, P < k → k ≤ N → coeff k (Λ * series s) = 0) ↔
      (∀ j, P < 2 * j + 1 → 2 * j + 1 ≤ N → coeff (2 * j + 1) (Λ ^ 2 * oddPart s) = 0) := by
  constructor
  · intro htail j hj1 hj2
    rw [← odd_coeff s hs Λ h0 j, htail _ hj1 hj2, zero_add]
    refine shifted_sum_zero Λ _ L _ hΛ (by omega) fun m hm1 hm2 => ?_
    by_cases hmP : m ≤ P
    · exact hpre m hm1 hmP
    · exact htail m (by omega) (by omega)
  · intro hodd k
    induction k using Nat.strong_induction_on with
    | _ k ih =>
      intro hk1 hk2
      have hbelow : ∀ m, L < m → m < k → coeff m (Λ * series s) = 0 := by
        intro m hm1 hm2
        by_cases hmP : m ≤ P
        · exact hpre m hm1 hmP
        · exact ih m hm2 (by omega) (by omega)
      have hsum := shifted_sum_zero Λ (Λ * series s) L k hΛ (by omega) hbelow
      rcases Nat.even_or_odd' k with ⟨j, rfl | rfl⟩
      · rw [even_coeff s hs Λ h0 j, hsum, add_zero, hbelow j (by omega) (by omega), sq,
          mul_zero]
      · have := odd_coeff s hs Λ h0 j
        rw [hsum, add_zero, hodd j hk1 hk2] at this
        exact this

end Minisketch.BMTail
