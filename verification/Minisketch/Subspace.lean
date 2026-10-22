/-
Copyright (c) 2026 The minisketch developers
Distributed under the MIT software license, see the accompanying
file LICENSE or http://www.opensource.org/licenses/mit-license.php.
-/
import Minisketch.Montgomery

/-!
# Subspace polynomials and evaluation cosets

The additive FFT (`src/additive_fft.h`) evaluates on cosets of the subspaces
`W_m = span(β_0, …, β_{m-1})` of a Cantor chain `β_0 = 1`, `β_{j+1}^2 + β_{j+1} = β_j`.
Their vanishing polynomials are the subspace polynomials `S^m`, with
`S^0 = X` and `S^{m+1} = (S^m)^2 + S^m`. This file proves claim C29 of
`doc/spec/claims.md` (a coset's vanishing polynomial is `S^m + c`, so its
derivative is 1) and the coset-halving fact `Z_{C0} = 1` on `C1` used by C31,
and combines them with `Montgomery.lean`.
-/

namespace Minisketch.Subspace

open Polynomial Minisketch.Montgomery

variable {K : Type*} [Field K]

/-- The subspace polynomial `S^m`. -/
noncomputable def S : ℕ → K[X]
  | 0 => X
  | m + 1 => S m ^ 2 + S m

theorem S_succ (m : ℕ) : (S (m + 1) : K[X]) = S m ^ 2 + S m := rfl

theorem S_monic (m : ℕ) : (S m : K[X]).Monic ∧ (S m : K[X]).natDegree = 2 ^ m := by
  induction m with
  | zero => exact ⟨monic_X, by simp [S]⟩
  | succ m ih =>
    obtain ⟨hm, hd⟩ := ih
    have hsq : (S m ^ 2 : K[X]).natDegree = 2 ^ (m + 1) := by
      rw [hm.natDegree_pow, hd, pow_succ, mul_comm]
    have hlt : (S m : K[X]).degree < (S m ^ 2 : K[X]).degree := by
      rw [degree_eq_natDegree hm.ne_zero, degree_eq_natDegree (hm.pow 2).ne_zero, hsq, hd]
      exact_mod_cast Nat.pow_lt_pow_right (by norm_num) (Nat.lt_succ_self m)
    refine ⟨(hm.pow 2).add_of_left hlt, ?_⟩
    rw [S_succ, natDegree_add_eq_left_of_degree_lt hlt, hsq]

variable [CharP K 2]

/-- `S^m` is additive (a linearized polynomial). -/
theorem S_eval_add (m : ℕ) (a b : K) :
    (S m).eval (a + b) = (S m).eval a + (S m).eval b := by
  induction m with
  | zero => simp [S]
  | succ m ih =>
    simp only [S_succ, eval_add, eval_pow, ih]
    rw [CharTwo.add_sq]
    ring

/-- C29: `(S^m)' = 1`. -/
theorem S_derivative (m : ℕ) : derivative (S m : K[X]) = 1 := by
  induction m with
  | zero => simp [S]
  | succ m ih => rw [S_succ, derivative_add, derivative_sq_eq_zero, ih, zero_add]

theorem S_eval_zero (m : ℕ) : (S m).eval (0 : K) = 0 := by
  have h := S_eval_add m (0 : K) 0
  rw [add_zero] at h
  have : (S m).eval (0 : K) + (S m).eval 0 = 0 := CharTwo.add_self_eq_zero _
  rw [← h] at this
  simpa using this

omit [CharP K 2] in
/-- `S^{a+b} = S^a ∘ S^b` as functions. -/
theorem S_eval_comp (a b : ℕ) (x : K) : (S (a + b)).eval x = (S a).eval ((S b).eval x) := by
  induction a with
  | zero => simp [S]
  | succ a ih =>
    rw [show a + 1 + b = (a + b) + 1 by omega, S_succ, S_succ, eval_add, eval_pow, eval_add,
      eval_pow, ih]

section Cantor

variable (β : ℕ → K) (hβ0 : β 0 = 1) (hβ : ∀ j, β (j + 1) ^ 2 + β (j + 1) = β j)
include hβ

omit [CharP K 2] in
/-- On a Cantor chain, `S^m(β_{j+m}) = β_j`. -/
theorem S_cantor (m j : ℕ) : (S m).eval (β (j + m)) = β j := by
  induction m generalizing j with
  | zero => simp [S]
  | succ m ih =>
    rw [S_succ, eval_add, eval_pow, show j + (m + 1) = (j + 1) + m by omega, ih (j + 1), hβ]

include hβ0 in
/-- `S^m(β_j) = 0` for `j < m`: the chain's first `m` elements span the roots. -/
theorem S_cantor_zero (m j : ℕ) (hj : j < m) : (S m).eval (β j) = 0 := by
  obtain ⟨k, rfl⟩ : ∃ k, m = (k + 1) + j := ⟨m - j - 1, by omega⟩
  rw [S_eval_comp, show β j = β (0 + j) by simp, S_cantor β hβ j 0, hβ0, S_eval_comp k 1,
    S_succ, S]
  simp [CharTwo.add_self_eq_zero, S_eval_zero]

include hβ0 in
/-- `S^m` is constant on `η + span(β_0..β_{m-1})`. -/
theorem S_coset (m : ℕ) (η : K) (T : Finset ℕ) (hT : ∀ j ∈ T, j < m) :
    (S m).eval (η + ∑ j ∈ T, β j) = (S m).eval η := by
  rw [S_eval_add]
  suffices h : (S m).eval (∑ j ∈ T, β j) = 0 by rw [h, add_zero]
  induction T using Finset.induction_on with
  | empty => simp [S_eval_zero]
  | insert j T hj ih =>
    rw [Finset.sum_insert hj, S_eval_add, S_cantor_zero β hβ0 hβ m j (hT j (by simp)),
      ih fun i hi => hT i (by simp [hi]), add_zero]

/-- The other half of a coset: `S^m(η + β_m + w) = S^m(η) + 1` for `w` in `W_m`. -/
theorem S_coset_shift (m : ℕ) (η : K) (T : Finset ℕ) (hT : ∀ j ∈ T, j < m) (hβ0' : β 0 = 1) :
    (S m).eval (η + β m + ∑ j ∈ T, β j) = (S m).eval η + 1 := by
  rw [S_coset β hβ0' hβ m (η + β m) T hT, S_eval_add, show β m = β (0 + m) by simp,
    S_cantor β hβ m 0, hβ0']

include hβ0 in
/-- A Cantor chain is linearly independent over GF(2): a nonempty subset sum is
nonzero (apply `S^t` for the largest index `t`). -/
theorem cantor_sum_ne_zero (T : Finset ℕ) (hT : T.Nonempty) : ∑ j ∈ T, β j ≠ 0 := by
  set t := T.max' hT
  have ht : t ∈ T := T.max'_mem hT
  intro h0
  have h := congrArg (S t).eval h0
  rw [← Finset.add_sum_erase T _ ht, S_eval_add, S_eval_zero,
    show β t = β (0 + t) by simp, S_cantor β hβ t 0, hβ0] at h
  have hz : (S t).eval (∑ j ∈ T.erase t, β j) = 0 := by
    have := S_coset β hβ0 hβ t 0 (T.erase t) fun j hj =>
      lt_of_le_of_ne (T.le_max' j (Finset.mem_of_mem_erase hj)) (Finset.ne_of_mem_erase hj)
    rwa [zero_add, S_eval_zero] at this
  rw [hz, add_zero] at h
  exact one_ne_zero h

/-- The evaluation coset `η + W_m`, `W_m = span(β_0..β_{m-1})`. -/
noncomputable def coset [DecidableEq K] (η : K) (m : ℕ) : Finset K :=
  (Finset.range m).powerset.image fun T => η + ∑ j ∈ T, β j

include hβ0 in
theorem coset_card [DecidableEq K] (η : K) (m : ℕ) : (coset β η m).card = 2 ^ m := by
  rw [coset, Finset.card_image_of_injOn, Finset.card_powerset, Finset.card_range]
  intro T _ T' _ h
  simp only [add_right_inj] at h
  by_contra hne
  have hsd : (symmDiff T T').Nonempty := by
    rw [Finset.nonempty_iff_ne_empty]
    exact fun he => hne (symmDiff_eq_bot.mp he)
  apply cantor_sum_ne_zero β hβ0 hβ _ hsd
  rw [symmDiff_eq_sup_sdiff_inf]
  simp only [Finset.sup_eq_union, Finset.inf_eq_inter]
  have hsub : T ∩ T' ⊆ T ∪ T' := Finset.inter_subset_left.trans Finset.subset_union_left
  have e1 := Finset.sum_sdiff hsub (f := β)
  have e2 := Finset.sum_union_inter (s₁ := T) (s₂ := T') (f := β)
  rw [h, CharTwo.add_self_eq_zero] at e2
  calc ∑ j ∈ (T ∪ T') \ (T ∩ T'), β j
      = (∑ j ∈ (T ∪ T') \ (T ∩ T'), β j + ∑ j ∈ T ∩ T', β j) + ∑ j ∈ T ∩ T', β j := by
        rw [add_assoc, CharTwo.add_self_eq_zero, add_zero]
    _ = 0 := by rw [e1, e2]

include hβ0 in
theorem coset_S [DecidableEq K] (η : K) (m : ℕ) :
    ∀ w ∈ coset β η m, (S m).eval w = (S m).eval η := by
  intro w hw
  obtain ⟨T, hT, rfl⟩ := Finset.mem_image.mp hw
  exact S_coset β hβ0 hβ m η T fun j hj => Finset.mem_range.mp (Finset.mem_powerset.mp hT hj)

end Cantor

/-- C29: if `S^m` takes the value `S^m(η)` at all `2^m` points of `C`, the
vanishing polynomial of `C` is `S^m + S^m(η)`, so `Z_C' = 1`. -/
theorem vanishing_eq (m : ℕ) (η : K) (C : Finset K) (hC : C.card = 2 ^ m)
    (hval : ∀ w ∈ C, (S m).eval w = (S m).eval η) :
    vanishing C = S m + Polynomial.C ((S m).eval η) := by
  obtain ⟨hm, hd⟩ := S_monic (K := K) m
  set R := S m + Polynomial.C ((S m).eval η)
  have hlt : (Polynomial.C ((S m).eval η)).degree < (S m : K[X]).degree :=
    degree_C_le.trans_lt (by
      rw [degree_eq_natDegree hm.ne_zero, hd]
      exact_mod_cast pow_pos two_pos m)
  have hR : R.Monic := hm.add_of_left hlt
  have hRdeg : R.natDegree = 2 ^ m := by rw [natDegree_add_eq_left_of_degree_lt hlt, hd]
  have hdvd : vanishing C ∣ R := by
    apply Finset.prod_dvd_of_coprime
    · intro w _ w' _ hne
      exact pairwise_coprime_X_sub_C Function.injective_id hne
    · intro w hw
      rw [dvd_iff_isRoot, IsRoot, eval_add, eval_C, hval w hw, CharTwo.add_self_eq_zero]
  exact (eq_of_monic_of_dvd_of_natDegree_le (vanishing_monic C) hR hdvd
    (by rw [hRdeg, vanishing_natDegree, hC])).symm

/-- C29: the vanishing polynomial of such a coset has derivative 1. -/
theorem vanishing_derivative (m : ℕ) (η : K) (C : Finset K) (hC : C.card = 2 ^ m)
    (hval : ∀ w ∈ C, (S m).eval w = (S m).eval η) :
    derivative (vanishing C) = 1 := by
  rw [vanishing_eq m η C hC hval, derivative_add, derivative_C, add_zero, S_derivative]

/-- Coset halving (used by C31): if `S^m` is `S^m(η)` on the `2^m` points of
`C0` and `S^m(η) + 1` on `C1`, then `Z_{C0} = 1` on `C1`. -/
theorem vanishing_eq_one (m : ℕ) (η : K) (C0 C1 : Finset K) (hC0 : C0.card = 2 ^ m)
    (h0 : ∀ w ∈ C0, (S m).eval w = (S m).eval η)
    (h1 : ∀ w ∈ C1, (S m).eval w = (S m).eval η + 1) :
    ∀ w ∈ C1, (vanishing C0).eval w = 1 := by
  intro w hw
  rw [vanishing_eq m η C0 hC0 h0, eval_add, eval_C, h1 w hw, add_assoc,
    add_comm (1 : K), ← add_assoc, CharTwo.add_self_eq_zero, zero_add]

/-- C30 with C29 discharged: one Montgomery squaring step on a coset `C` of
`2^m` points on which `S^m` is constant. -/
theorem reduce_coset [DecidableEq K] (m : ℕ) (η : K) (C : Finset K) (hC : C.card = 2 ^ m)
    (hval : ∀ w ∈ C, (S m).eval w = (S m).eval η) (s f q : K[X])
    (hf : ∀ w ∈ C, f.eval w ≠ 0) (hq : ∀ w ∈ C, q.eval w = (s ^ 2).eval w / f.eval w) :
    ∃ r : K[X], r * vanishing C = s ^ 2 + q * f ∧ f ∣ r * vanishing C - s ^ 2 ∧
      (∀ w ∈ C, r.eval w = (derivative q).eval w * f.eval w + q.eval w * (derivative f).eval w) ∧
      ((s ^ 2 + q * f).natDegree < 2 * C.card → r.natDegree < C.card) := by
  obtain ⟨r, h1, h2, h3, h4⟩ := reduce_spec C s f q hf hq
  refine ⟨r, h1, h2, fun w hw => h3 w hw ?_, h4⟩
  rw [vanishing_derivative m η C hC hval, eval_one]

/-- C31 with the coset facts discharged: carrying a residue from `C = C0 ⊔ C1`
to the child coset `C1`. -/
theorem carry_coset [DecidableEq K] (m : ℕ) (η : K) (C0 C1 : Finset K) (hC0 : C0.card = 2 ^ m)
    (h0 : ∀ w ∈ C0, (S m).eval w = (S m).eval η)
    (h1 : ∀ w ∈ C1, (S m).eval w = (S m).eval η + 1) (v g q : K[X])
    (hg : ∀ w ∈ C0, g.eval w ≠ 0) (hq : ∀ w ∈ C0, q.eval w = v.eval w / g.eval w) :
    ∃ r : K[X], g ∣ r * vanishing C0 - v ∧ ∀ w ∈ C1, r.eval w = (v + q * g).eval w :=
  carry_to_factor C0 C1 v g q hg hq (vanishing_eq_one m η C0 C1 hC0 h0 h1)

/-- C29 on the actual evaluation cosets: `Z_{η + W_m}' = 1`. -/
theorem cantor_vanishing_derivative [DecidableEq K] (β : ℕ → K) (hβ0 : β 0 = 1)
    (hβ : ∀ j, β (j + 1) ^ 2 + β (j + 1) = β j) (η : K) (m : ℕ) :
    derivative (vanishing (coset β η m)) = 1 :=
  vanishing_derivative m η _ (coset_card β hβ0 hβ η m) (coset_S β hβ0 hβ η m)

/-- C34: the novel basis `X_T = Π_{j ∈ T} S^j` has derivative
`Σ_{j ∈ T} X_{T \ {j}}`, since every `(S^j)' = 1`. With `T` the bit set of `k`
this is `X_k' = Σ_{j ∈ bits(k)} X_{k - 2^j}`. -/
theorem novel_derivative (T : Finset ℕ) :
    derivative (∏ j ∈ T, (S j : K[X])) = ∑ j ∈ T, ∏ i ∈ T.erase j, (S i : K[X]) := by
  rw [derivative_prod_finset]
  refine Finset.sum_congr rfl fun j _ => ?_
  rw [S_derivative, mul_one]

end Minisketch.Subspace
