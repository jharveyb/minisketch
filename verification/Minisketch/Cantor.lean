/-
Copyright (c) 2026 The minisketch developers
Distributed under the MIT software license, see the accompanying
file LICENSE or http://www.opensource.org/licenses/mit-license.php.
-/
import Minisketch.Subspace
import Minisketch.Trace

/-!
# The Cantor basis exists

Claim C38 of `doc/spec/claims.md`: in `GF(2^b)` with `b` a power of two there is
a chain `β_0 = 1`, `β_{j+1}^2 + β_{j+1} = β_j` for `j + 1 < b`, which the
additive FFT (`src/additive_fft.h`) builds with the field's `Qrt`.
-/

namespace Minisketch.Cantor

open Polynomial Minisketch.Subspace Minisketch.Trace

variable {K : Type*} [Field K] [CharP K 2]

/-- `Tr(y^2 + y) = 0` in `GF(2^b)`. -/
theorem tr_artin_schreier (b : ℕ) (y : K) (hy : y ^ 2 ^ b = y) : tr b (y ^ 2 + y) = 0 := by
  have h := tr_sq_add b y
  rw [hy, CharTwo.add_self_eq_zero] at h
  -- tr(y^2) = tr(y)^2 by Frobenius additivity.
  have hsq : tr b (y ^ 2) = tr b y ^ 2 := by
    simp only [tr]
    rw [CharTwo.sum_sq]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [← pow_mul, ← pow_mul, mul_comm]
  rw [tr_add, hsq, h]

/-- In `GF(2^b)`, `y^2 + y = c` is solvable iff `Tr(c) = 0`. -/
theorem artin_schreier [Fintype K] [DecidableEq K] (b : ℕ) (hb : 1 ≤ b)
    (hq : Fintype.card K = 2 ^ b) (c : K) (hc : tr b c = 0) : ∃ y, y ^ 2 + y = c := by
  have hfix : ∀ y : K, y ^ 2 ^ b = y := fun y => by rw [← hq]; exact FiniteField.pow_card y
  set φ : K → K := fun y => y ^ 2 + y
  -- φ(y) = φ(y') iff y' ∈ {y, y + 1}, so the image has |K| / 2 elements.
  have hfib : ∀ y y', φ y = φ y' → y' = y ∨ y' = y + 1 := by
    intro y y' h
    have : (y' + y) * (y' + y + 1) = 0 := by
      have e : (y' + y) * (y' + y + 1) = φ y' + φ y := by
        simp only [φ]; linear_combination (y' * y) * CharTwo.two_eq_zero (R := K)
      rw [e, h, CharTwo.add_self_eq_zero]
    rcases mul_eq_zero.mp this with h1 | h1
    · left; rw [← sub_eq_zero, sub_eq_add_neg, CharTwo.neg_eq]; exact h1
    · right; rw [← sub_eq_zero, sub_eq_add_neg, CharTwo.neg_eq, ← add_assoc]; exact h1
  set img := Finset.univ.image φ
  set ker := Finset.univ.filter fun x : K => tr b x = 0
  have hsub : img ⊆ ker := by
    intro x hx
    obtain ⟨y, _, rfl⟩ := Finset.mem_image.mp hx
    simp only [ker, Finset.mem_filter, Finset.mem_univ, true_and]
    exact tr_artin_schreier b y (hfix y)
  have hlow : Fintype.card K ≤ 2 * img.card := by
    rw [← Finset.card_univ]
    refine Finset.card_le_mul_card_image _ 2 fun a ha => ?_
    obtain ⟨y, _, rfl⟩ := Finset.mem_image.mp ha
    calc (Finset.univ.filter fun x => φ x = φ y).card ≤ ({y, y + 1} : Finset K).card :=
          Finset.card_le_card fun x hx => by
            simp only [Finset.mem_filter, Finset.mem_univ, true_and] at hx
            rcases hfib y x hx.symm with rfl | rfl <;> simp
      _ ≤ 2 := Finset.card_le_two
  obtain ⟨hdeg, hne0⟩ := trPoly_one_natDegree (K := K) b hb
  have hhigh : ker.card ≤ 2 ^ (b - 1) := by
    rw [← hdeg]
    exact Polynomial.card_le_degree_of_subset_roots fun x hx => by
      rw [Polynomial.mem_roots hne0, IsRoot, trPoly_eval, one_mul]
      simpa [ker] using hx
  have hpow : 2 ^ b = 2 * 2 ^ (b - 1) := by
    rw [← pow_succ']; congr 1; omega
  have heq : img = ker := Finset.eq_of_subset_of_card_le hsub (by
    have := Finset.card_le_card hsub
    omega)
  have hcmem : c ∈ img := by rw [heq]; simpa [ker] using hc
  obtain ⟨y, _, hy⟩ := Finset.mem_image.mp hcmem
  exact ⟨y, hy⟩

/-- `S^(2^ℓ)(x) = x^(2^(2^ℓ)) + x`. -/
theorem S_two_pow (ℓ : ℕ) (x : K) : (S (2 ^ ℓ)).eval x = x ^ 2 ^ 2 ^ ℓ + x := by
  induction ℓ generalizing x with
  | zero => simp [S]
  | succ ℓ ih =>
    rw [show 2 ^ (ℓ + 1) = 2 ^ ℓ + 2 ^ ℓ by ring, S_eval_comp, ih, ih, frob_add, ← pow_mul,
      ← pow_add, show 2 ^ ℓ + 2 ^ ℓ = 2 ^ (ℓ + 1) by ring, add_assoc, ← add_assoc (x ^ 2 ^ 2 ^ ℓ),
      CharTwo.add_self_eq_zero, zero_add]

/-- Frobenius powers distribute over sums. -/
theorem sum_pow_two_pow {ι : Type*} (s : Finset ι) (f : ι → K) (n : ℕ) :
    (∑ i ∈ s, f i) ^ 2 ^ n = ∑ i ∈ s, f i ^ 2 ^ n := by
  classical
  induction s using Finset.induction_on with
  | empty => simp
  | insert a s ha ih => rw [Finset.sum_insert ha, Finset.sum_insert ha, frob_add, ih]

/-- `S^(2^ℓ - 1) = Tr` on `GF(2^(2^ℓ))` (all binomials `C(2^ℓ - 1, j)` are odd). -/
theorem S_tr (ℓ : ℕ) (x : K) : (S (2 ^ ℓ - 1)).eval x = tr (2 ^ ℓ) x := by
  induction ℓ generalizing x with
  | zero => simp [S, tr]
  | succ ℓ ih =>
    have h1 : 1 ≤ 2 ^ ℓ := Nat.one_le_two_pow
    rw [show 2 ^ (ℓ + 1) - 1 = 2 ^ ℓ + (2 ^ ℓ - 1) by rw [pow_succ]; omega, S_eval_comp, ih,
      S_two_pow, tr, sum_pow_two_pow, tr, show 2 ^ (ℓ + 1) = 2 ^ ℓ + 2 ^ ℓ by ring,
      Finset.sum_range_add, add_comm]
    congr 1
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [← pow_mul, ← pow_add, add_comm]

omit [CharP K 2] in
/-- A partial Cantor chain satisfies `S^m(β_{j+m}) = β_j` within its length. -/
theorem S_chain (β : ℕ → K) (n : ℕ) (hβ : ∀ j, j + 1 < n → β (j + 1) ^ 2 + β (j + 1) = β j) :
    ∀ m j, j + m < n → (S m).eval (β (j + m)) = β j := by
  intro m
  induction m with
  | zero => intro j _; simp [S]
  | succ m ih =>
    intro j hj
    rw [S_succ, eval_add, eval_pow, show j + (m + 1) = (j + 1) + m by omega, ih (j + 1) (by omega),
      hβ j (by omega)]

theorem S_one (k : ℕ) (hk : 1 ≤ k) : (S k).eval (1 : K) = 0 := by
  obtain ⟨k, rfl⟩ : ∃ k', k = k' + 1 := ⟨k - 1, by omega⟩
  rw [show k + 1 = k + 1 by rfl, S_eval_comp k 1, S_succ, S]
  simp [CharTwo.add_self_eq_zero, S_eval_zero]

/-- C38: in `GF(2^b)` with `b = 2^ℓ` there is a Cantor chain of length `b`:
`β_0 = 1` and `β_{j+1}^2 + β_{j+1} = β_j` for `j + 1 < b`. -/
theorem cantor_basis_exists [Fintype K] [DecidableEq K] (ℓ : ℕ)
    (hq : Fintype.card K = 2 ^ 2 ^ ℓ) :
    ∃ β : ℕ → K, β 0 = 1 ∧ ∀ j, j + 1 < 2 ^ ℓ → β (j + 1) ^ 2 + β (j + 1) = β j := by
  set b := 2 ^ ℓ
  have hb : 1 ≤ b := Nat.one_le_two_pow
  suffices h : ∀ n, n ≤ b → ∃ β : ℕ → K, β 0 = 1 ∧ ∀ j, j + 1 < n → β (j + 1) ^ 2 + β (j + 1) = β j
    from h b le_rfl
  intro n
  induction n with
  | zero => exact fun _ => ⟨fun _ => 1, rfl, fun j hj => by omega⟩
  | succ n ih =>
    intro hn
    obtain ⟨β, h0, hβ⟩ := ih (by omega)
    rcases Nat.eq_zero_or_pos n with rfl | hpos
    · exact ⟨β, h0, fun j hj => by omega⟩
    -- Extend by a root of y^2 + y = β_{n-1}; its trace is S^{b-n}(β_0) = S^{b-n}(1) = 0.
    have htr : tr b (β (n - 1)) = 0 := by
      rw [← S_tr, show b - 1 = (b - n) + (n - 1) by omega, S_eval_comp,
        show β (n - 1) = β (0 + (n - 1)) by simp, S_chain β n hβ (n - 1) 0 (by omega), h0,
        S_one _ (by omega)]
    obtain ⟨y, hy⟩ := artin_schreier b hb hq _ htr
    refine ⟨Function.update β n y, by rw [Function.update_of_ne (by omega), h0], fun j hj => ?_⟩
    rcases Nat.lt_or_ge (j + 1) n with hj' | hj'
    · rw [Function.update_of_ne (by omega), Function.update_of_ne (by omega), hβ j hj']
    · obtain rfl : j = n - 1 := by omega
      rw [show n - 1 + 1 = n by omega, Function.update_self, Function.update_of_ne (by omega), hy]

end Minisketch.Cantor
