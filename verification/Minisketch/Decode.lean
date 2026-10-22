/-
Copyright (c) 2026 The minisketch developers
Distributed under the MIT software license, see the accompanying
file LICENSE or http://www.opensource.org/licenses/mit-license.php.
-/
import Minisketch.BMEven
import Minisketch.Sketch

/-!
# A successful decode re-encodes

Claims C2 and C37 of `doc/spec/claims.md`. If Berlekamp–Massey's connection
polynomial on a Frobenius-structured syndrome sequence factors as
`Π_{α ∈ R} (1 - α x)` over nonzero `α`, then the power sums of `R` are the input
syndromes (so the decoded set re-encodes to the sketch), and `|R|` is BM's
length (so the reversed polynomial, whose roots are returned, has no root 0).
-/

namespace Minisketch.Decode

open PowerSeries Minisketch.BM Minisketch.BMTail Minisketch.BMEven Minisketch.Sketch

variable {K : Type*} [Field K]

/-- The Euler operator `X d/dX`: `Σ k F_k X^k`. -/
noncomputable def euler (F : K⟦X⟧) : K⟦X⟧ := mk fun k => (k : K) * coeff k F

theorem euler_mul (F G : K⟦X⟧) : euler (F * G) = euler F * G + F * euler G := by
  ext k
  simp only [euler, coeff_mk, map_add, coeff_mul, Finset.mul_sum, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun p hp => ?_
  rw [Finset.HasAntidiagonal.mem_antidiagonal] at hp
  rw [← hp]
  push_cast
  ring

/-- The locator series `Π_{α ∈ R} (1 - α X)`. -/
noncomputable def loc (R : Finset K) : K⟦X⟧ := ∏ α ∈ R, (1 - C α * X)

/-- `Σ_{k ≥ 1} β^k X^k`. -/
noncomputable def geom (β : K) : K⟦X⟧ := mk fun k => if k = 0 then 0 else β ^ k

theorem one_sub_mul_geom (β : K) : (1 - C β * X) * geom β = C β * X := by
  ext k
  rcases k with _ | k
  · simp [geom]
  · rw [sub_mul, one_mul, map_sub, mul_assoc, coeff_C_mul, coeff_succ_X_mul, coeff_C_mul,
      coeff_X]
    simp only [geom, coeff_mk]
    rcases k with _ | k
    · simp
    · simp [pow_succ]
      ring

/-- C2: `Λ_R · S_R = -X Λ_R'` for any finite set `R` (the power-sum series
`S_R = Σ_{k ≥ 1} s_k X^k`). -/
theorem loc_mul_series [DecidableEq K] (R : Finset K) :
    loc R * series (psum R) = -euler (loc R) := by
  induction R using Finset.induction_on with
  | empty =>
    ext k
    simp [loc, euler, series, psum]
    intro hk; simp [hk]
  | insert β R hβ ih =>
    have hs : series (psum (insert β R)) = series (psum R) + geom β := by
      ext k
      simp only [series, geom, coeff_mk, map_add, psum, Finset.sum_insert hβ]
      split_ifs <;> ring
    rw [loc, Finset.prod_insert hβ, ← loc, hs, euler_mul]
    have he : euler (1 - C β * X : K⟦X⟧) = -(C β * X) := by
      ext k
      rcases k with _ | k
      · simp [euler]
      · rcases k with _ | k <;> simp [euler, coeff_X]
    rw [he]
    calc (1 - C β * X) * loc R * (series (psum R) + geom β)
        = (1 - C β * X) * (loc R * series (psum R)) + loc R * ((1 - C β * X) * geom β) := by ring
      _ = _ := by rw [ih, one_sub_mul_geom]; ring

/-- In characteristic 2, `-X F'` is the odd part. -/
theorem neg_euler [CharP K 2] (c : ℕ → K) : -euler (ps c) = oddPart c := by
  ext k
  simp only [euler, ps, oddPart, coeff_mk, map_neg, CharTwo.neg_eq]
  rcases Nat.even_or_odd k with he | ho
  · have hk : (k : K) = 0 := (CharP.cast_eq_zero_iff K 2 k).mpr (even_iff_two_dvd.mp he)
    rw [ite_eq_right (Nat.not_odd_iff_even.mpr he), hk, zero_mul]
  · rw [ite_eq_left ho]
    have : (k : K) = 1 := by
      obtain ⟨a, rfl⟩ := ho
      push_cast
      rw [show (2 : K) = 0 from CharTwo.two_eq_zero, zero_mul, zero_add]
    rw [this, one_mul]

/-- If `Λ_0 = 1` and `[Λ F]_j = 0` for `j ≤ N`, then `F_j = 0` for `j ≤ N`. -/
theorem coeff_eq_zero_of_mul (Λ F : K⟦X⟧) (h0 : constantCoeff Λ = 1) (N : ℕ)
    (h : ∀ j, j ≤ N → coeff j (Λ * F) = 0) : ∀ j, j ≤ N → coeff j F = 0 := by
  intro j
  induction j using Nat.strong_induction_on with
  | _ j ih =>
    intro hj
    have hm := h j hj
    rw [coeff_mul, Finset.sum_eq_single (0, j)] at hm
    · rwa [coeff_zero_eq_constantCoeff_apply, h0, one_mul] at hm
    · intro p hp hne
      rw [Finset.HasAntidiagonal.mem_antidiagonal] at hp
      have : p.2 < j := by
        rcases p with ⟨a, b⟩
        simp only [ne_eq, Prod.mk.injEq] at hne hp ⊢
        omega
      rw [ih p.2 this (by omega), mul_zero]
    · intro hn
      exact absurd (Finset.HasAntidiagonal.mem_antidiagonal.mpr (by simp)) hn

/-- `Π_{α ∈ R} (1 - α X)` has degree `|R|`, with top coefficient `Π (-α)`. -/
theorem loc_coeff [DecidableEq K] (R : Finset K) :
    (∀ i, R.card < i → coeff i (loc R) = 0) ∧ coeff R.card (loc R) = ∏ α ∈ R, (-α) := by
  induction R using Finset.induction_on with
  | empty =>
    refine ⟨fun i hi => ?_, by simp [loc]⟩
    simp [loc, coeff_one]; omega
  | insert β R hβ ih =>
    obtain ⟨hz, htop⟩ := ih
    have hstep : ∀ i, coeff (i + 1) (loc (insert β R)) = coeff (i + 1) (loc R) - β * coeff i (loc R) := by
      intro i
      rw [loc, Finset.prod_insert hβ, ← loc, sub_mul, one_mul, map_sub, mul_assoc, coeff_C_mul,
        coeff_succ_X_mul]
    rw [Finset.card_insert_of_notMem hβ, Finset.prod_insert hβ]
    refine ⟨fun i hi => ?_, ?_⟩
    · obtain ⟨i, rfl⟩ : ∃ j, i = j + 1 := ⟨i - 1, by omega⟩
      rw [hstep, hz _ (by omega), hz _ (by omega), mul_zero, sub_zero]
    · rw [hstep, hz _ (by omega), htop]
      ring

/-- C37: if BM's connection polynomial after `N` terms of a Frobenius-structured
sequence factors as `Π_{α ∈ R} (1 - α x)` with `0 ∉ R`, then `R`'s power sums
are `s_1..s_N` (the decoded set re-encodes to the sketch) and `|R|` is BM's
length, so the reversed polynomial `x^L Λ(1/x) = Π (x - α)` has no root `0`. -/
theorem decode_reencodes [CharP K 2] [DecidableEq K] (s : ℕ → K) (hs : Frobenius s)
    (hs0 : s 0 = 1) (N : ℕ) (R : Finset K) (h0R : (0 : K) ∉ R)
    (hfac : ps (bm s N).C = loc R) :
    (∀ k, 1 ≤ k → k ≤ N → s k = psum R k) ∧ R.card = (bm s N).L := by
  have hinv := inv_bm s hs0 N
  have hkey := (inv2_bm s hs hs0 N).keyC
  have h0 : constantCoeff (ps (bm s N).C) = 1 := by simp [ps, hinv.C0]
  have hC2 : ps (bm s N).C * series (psum R) = oddPart (bm s N).C := by
    rw [hfac, loc_mul_series, ← hfac, neg_euler]
  have hz := coeff_eq_zero_of_mul (ps (bm s N).C) (series s - series (psum R)) h0 N
    fun j hj => by rw [mul_sub, map_sub, hkey j hj, hC2, sub_self]
  obtain ⟨hlz, hltop⟩ := loc_coeff R
  have hCcoeff : ∀ i, (bm s N).C i = coeff i (loc R) := fun i => by
    rw [← hfac, ps, coeff_mk]
  refine ⟨fun k hk1 hk2 => ?_, ?_⟩
  · have := hz k hk2
    rw [map_sub, sub_eq_zero, series, series, coeff_mk, coeff_mk, ite_eq_right (by omega),
      ite_eq_right (by omega)] at this
    exact this
  · have hsupp : ∀ i, R.card < i → (bm s N).C i = 0 := fun i hi => by rw [hCcoeff, hlz i hi]
    have hgen : Gen s (bm s N).C R.card N := by
      intro k hk1 hk2
      rw [disc_eq_coeff s _ _ _ hsupp hk1, hkey k hk2, oddPart_coeff, hsupp k hk1]
      simp
    have hle := (bm_minimal s hs0 N).2.2 _ _ hinv.C0 hgen
    have hge : R.card ≤ (bm s N).L := by
      by_contra hlt
      push Not at hlt
      have h1 := hinv.Csupp R.card hlt
      rw [hCcoeff, hltop] at h1
      exact Finset.prod_ne_zero_iff.mpr (fun α hα => neg_ne_zero.mpr
        (fun h => h0R (h ▸ hα))) h1
    omega

/-- A product of two series supported on `[0, a]` and `[0, b]` is supported on
`[0, a + b]`. -/
theorem coeff_mul_supp (F G : K⟦X⟧) (a b : ℕ) (hF : ∀ i, a < i → coeff i F = 0)
    (hG : ∀ i, b < i → coeff i G = 0) (k : ℕ) (hk : a + b < k) : coeff k (F * G) = 0 := by
  rw [coeff_mul]
  refine Finset.sum_eq_zero fun p hp => ?_
  rw [Finset.HasAntidiagonal.mem_antidiagonal] at hp
  by_cases ha : a < p.1
  · rw [hF _ ha, zero_mul]
  · rw [hG _ (by omega), mul_zero]

/-- C16: for a set `R`, the stored odd syndromes `u_j = s_{2j+1}` satisfy the
recurrence with `Λ_R^2`: `Λ_R^2 · S_odd` is a polynomial of degree at most
`2|R|`. (From C13 and C2: `Λ² S_odd = Λ C + C²` with `C = Λ S` the odd part of
`Λ`, both of degree at most `|R|`.) -/
theorem odd_recurrence [CharP K 2] [DecidableEq K] (R : Finset K) (k : ℕ) (hk : 2 * R.card < k) :
    coeff k (loc R ^ 2 * oddPart (psum R)) = 0 := by
  have hs : Frobenius (psum R) := fun k _ => psum_two_mul R k
  have hloc : loc R = ps fun i => coeff i (loc R) := by ext i; simp [ps]
  have hC : loc R * series (psum R) = oddPart (fun i => coeff i (loc R)) := by
    rw [loc_mul_series, hloc, neg_euler, ← hloc]
  have hid := tail_identity (psum R) hs (loc R)
  rw [hC] at hid
  obtain ⟨hz, _⟩ := loc_coeff R
  have hodd : ∀ i, R.card < i → coeff i (oddPart fun i => coeff i (loc R)) = 0 := by
    intro i hi; rw [oddPart_coeff, hz i hi]; simp
  have : CharP K⟦X⟧ 2 := charP_of_injective_algebraMap (algebraMap K K⟦X⟧).injective 2
  have e : loc R ^ 2 * oddPart (psum R) =
      loc R * oddPart (fun i => coeff i (loc R)) + oddPart (fun i => coeff i (loc R)) ^ 2 := by
    rw [hid, add_comm, ← add_assoc, CharTwo.add_self_eq_zero, zero_add]
  rw [e, map_add, coeff_mul_supp _ _ R.card R.card hz hodd k (by omega), sq,
    coeff_mul_supp _ _ R.card R.card hodd hodd k (by omega), add_zero]

end Minisketch.Decode
