/-
Copyright (c) 2026 The minisketch developers
Distributed under the MIT software license, see the accompanying
file LICENSE or http://www.opensource.org/licenses/mit-license.php.
-/
import Mathlib

/-!
# Two-fold modular reduction

The CLMUL/PMULL fields (`ReduceWithClMul`, `ReduceTrinomial` in
`src/fields/clmul_common_impl.h`) reduce a carry-less product modulo
`x^B + MOD` with two folds `hi · x^B ≡ -hi · MOD`. This file proves the
algebraic half of claim C47 of `doc/spec/claims.md`: two folds give the
remainder exactly when `2 · deg MOD - 2 < B`. (The lane-width bounds of the
SIMD implementation are a separate, bit-level matter.)
-/

namespace Minisketch.Reduce

open Polynomial

variable {K : Type*} [Field K]

/-- One fold: `p mod x^B - (p div x^B) · m`. -/
noncomputable def fold (B : ℕ) (m p : K[X]) : K[X] := p %ₘ X ^ B - (p /ₘ X ^ B) * m

theorem fold_add (B : ℕ) (m p : K[X]) : fold B m p + (X ^ B + m) * (p /ₘ X ^ B) = p := by
  have := modByMonic_add_div p (X ^ B : K[X])
  rw [fold]; linear_combination this

theorem fold_natDegree (B : ℕ) (hB : 1 ≤ B) (m p : K[X]) :
    (fold B m p).natDegree ≤ max (B - 1) (p.natDegree - B + m.natDegree) := by
  have h1 : (p %ₘ X ^ B).natDegree ≤ B - 1 := by
    have := natDegree_modByMonic_lt p (monic_X_pow (R := K) B)
      (by intro h; have := congrArg natDegree h; simp at this; omega)
    rw [natDegree_X_pow] at this; omega
  have h2 : ((p /ₘ X ^ B) * m).natDegree ≤ p.natDegree - B + m.natDegree := by
    refine natDegree_mul_le.trans ?_
    rw [natDegree_divByMonic p (monic_X_pow B), natDegree_X_pow]
  exact (natDegree_sub_le _ _).trans (max_le_max h1 h2)

/-- C47 (algebra): for `f = x^B + m` with `deg m = t < B` and `2t - 2 < B`,
two folds of any `p` with `deg p ≤ 2B - 2` (a product of two reduced
elements) give exactly `p mod f`. -/
theorem two_folds (B : ℕ) (hB : 1 ≤ B) (m p : K[X]) (ht : m.natDegree < B)
    (ht2 : 2 * m.natDegree < B + 2) (hp : p.natDegree ≤ 2 * B - 2) :
    fold B m (fold B m p) = p %ₘ (X ^ B + m) := by
  set t := m.natDegree
  have hlt : m.degree < (X ^ B : K[X]).degree := by
    rw [degree_X_pow]
    exact (degree_le_natDegree).trans_lt (by exact_mod_cast ht)
  have hf : (X ^ B + m).Monic := (monic_X_pow B).add_of_left hlt
  have hfdeg : (X ^ B + m).degree = (B : WithBot ℕ) := by
    rw [degree_add_eq_left_of_degree_lt hlt, degree_X_pow]
  have d1 := fold_natDegree B hB m p
  have d2 := fold_natDegree B hB m (fold B m p)
  have hd2 : (fold B m (fold B m p)).natDegree ≤ B - 1 := by
    refine d2.trans (max_le le_rfl ?_)
    have : (fold B m p).natDegree ≤ max (B - 1) (B - 2 + t) := d1.trans (by omega)
    omega
  refine ((div_modByMonic_unique (f := p) ((p /ₘ X ^ B) + (fold B m p /ₘ X ^ B))
    (fold B m (fold B m p)) hf ⟨?_, ?_⟩).2).symm
  · have e1 := fold_add B m p
    have e2 := fold_add B m (fold B m p)
    linear_combination e1 + e2
  · rw [hfdeg]
    exact (degree_le_natDegree).trans_lt (by exact_mod_cast (by omega : _ < B))

end Minisketch.Reduce
