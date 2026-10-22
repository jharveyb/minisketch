/-
Copyright (c) 2026 The minisketch developers
Distributed under the MIT software license, see the accompanying
file LICENSE or http://www.opensource.org/licenses/mit-license.php.
-/
import Mathlib

/-!
# Coset Montgomery reduction

Root finding for large root polynomials keeps residues mod `f` as values on an
evaluation coset `C`, with Montgomery radix `R = Z_C = Π_{w ∈ C} (x - w)`
(`MontCoset` in `src/sketch_impl.h`). These are claims C30, C32 and C33 of
`doc/spec/claims.md`. The facts about `Z_C` itself (monic, `Z_C' = 1` for a
coset of an additive subgroup) are C29 and live with the additive FFT; here
they are hypotheses.
-/

namespace Minisketch.Montgomery

open Polynomial

variable {K : Type*} [Field K]

/-- C30, Montgomery step: choosing `q = a / f` pointwise on `C` makes
`a + q f` vanish on `C` (characteristic 2), so it is divisible by `Z_C`. -/
theorem vanishing_dvd [CharP K 2] [DecidableEq K] (C : Finset K) (a f q : K[X])
    (hf : ∀ w ∈ C, f.eval w ≠ 0) (hq : ∀ w ∈ C, q.eval w = a.eval w / f.eval w) :
    (∏ w ∈ C, (X - Polynomial.C w)) ∣ a + q * f := by
  apply Finset.prod_dvd_of_coprime
  · intro w _ w' _ hne
    exact pairwise_coprime_X_sub_C (Function.injective_id) hne
  · intro w hw
    rw [dvd_iff_isRoot, IsRoot, eval_add, eval_mul, hq w hw, div_mul_cancel₀ _ (hf w hw),
      CharTwo.add_self_eq_zero]

/-- C30, Montgomery congruence: if `r · Z_C = a + q f` then `r · Z_C ≡ a (mod f)`,
i.e. `r` is `a · Z_C⁻¹ mod f` whenever `Z_C` is a unit mod `f`. -/
theorem reduce_modEq (a f q r R : K[X]) (h : r * R = a + q * f) :
    f ∣ r * R - a :=
  ⟨q, by rw [h]; ring⟩

/-- C30, the derivative trick: from `r · R = a + q f`, at a point `w` where
`R(w) = 0` and `R'(w) = 1`, `r(w) = a'(w) + q'(w) f(w) + q(w) f'(w)`. Everything
except `q'` is a pointwise product of known values. -/
theorem reduce_eval (a f q r R : K[X]) (h : r * R = a + q * f) (w : K)
    (hR : R.eval w = 0) (hR' : (derivative R).eval w = 1) :
    r.eval w = (derivative a).eval w + (derivative q).eval w * f.eval w +
      q.eval w * (derivative f).eval w := by
  have hd := congrArg (fun p => (derivative p).eval w) h
  simp only [derivative_mul, derivative_add, eval_add, eval_mul, hR, hR'] at hd
  linear_combination hd

/-- C30 for squares: in characteristic 2 the derivative of a square vanishes, so
for `a = s^2` the `a'` term drops out. -/
theorem derivative_sq_eq_zero [CharP K 2] (s : K[X]) : derivative (s ^ 2) = 0 := by
  rw [Polynomial.derivative_sq, CharTwo.two_eq_zero, map_zero, zero_mul, zero_mul]

/-- C32: if `Z` is coprime to `g`, then `g` has the same common divisors with
`T · Z mod g` as with `T`. So the gcd split of a node can use a trace left in
Montgomery form. -/
theorem dvd_mod_mul_unit_iff (g T Z d : K[X]) (hZ : IsCoprime Z g) (hd : d ∣ g) :
    d ∣ (T * Z) % g ↔ d ∣ T := by
  have hdz : IsCoprime d Z := (hZ.symm.of_isCoprime_of_dvd_left hd)
  have hmod : (T * Z) % g = T * Z - g * (T * Z / g) := by
    rw [EuclideanDomain.mod_eq_sub_mul_div]
  rw [hmod, dvd_sub_left (dvd_mul_of_dvd_left hd _)]
  exact ⟨fun h => hdz.dvd_of_dvd_mul_right h, fun h => dvd_mul_of_dvd_left h _⟩

/-- C33, `CosetDivide`: if `g` and `h` are monic, `deg g ≤ |C|`, and `q` agrees
with `g / h` at every point of `C`, then `q` monic of degree `deg g - deg h`
already forces `q · h = g`. The difference `q h - g` has degree below `deg g`
and vanishes on `|C| ≥ deg g` points. -/
theorem coset_divide_exact (C : Finset K) (g h q : K[X]) (hg : g.Monic) (hh : h.Monic)
    (hq : q.Monic) (hdeg : q.natDegree + h.natDegree = g.natDegree)
    (hN : g.natDegree ≤ C.card) (hval : ∀ w ∈ C, q.eval w * h.eval w = g.eval w) :
    q * h = g := by
  have hqh : (q * h).Monic := hq.mul hh
  have hqhdeg : (q * h).natDegree = g.natDegree := by rw [hq.natDegree_mul hh, hdeg]
  have hdg : (q * h).degree = g.degree := by
    rw [degree_eq_natDegree hqh.ne_zero, degree_eq_natDegree hg.ne_zero, hqhdeg]
  apply Polynomial.eq_of_degree_sub_lt_of_eval_finset_eq C
  · calc (q * h - g).degree < (q * h).degree :=
          degree_sub_lt_left hdg hqh.ne_zero (by rw [hqh.leadingCoeff, hg.leadingCoeff])
      _ = g.degree := hdg
      _ ≤ C.card := by
          rw [degree_eq_natDegree hg.ne_zero]
          exact_mod_cast hN
  · intro w hw
    rw [eval_mul, hval w hw]

/-- The vanishing polynomial of a finite set of points. -/
noncomputable def vanishing (C : Finset K) : K[X] := ∏ w ∈ C, (X - Polynomial.C w)

theorem vanishing_monic (C : Finset K) : (vanishing C).Monic :=
  monic_prod_of_monic _ _ fun w _ => monic_X_sub_C w

theorem vanishing_natDegree (C : Finset K) : (vanishing C).natDegree = C.card := by
  rw [vanishing, natDegree_prod_of_monic _ _ fun w _ => monic_X_sub_C w]
  simp

/-- C30, combined: one Montgomery squaring step on the coset `C`. With `q = a/f`
on `C`, there is `r` with `r · Z_C = a + q f`; it represents `a · Z_C⁻¹ mod f`;
its values on `C` are `q' f + q f'` when `a` is a square and `Z_C' = 1` on `C`
(claim C29); and `deg r < |C|` when `deg a, deg (q f) < 2|C|`. -/
theorem reduce_spec [CharP K 2] [DecidableEq K] (C : Finset K) (s f q : K[X])
    (hf : ∀ w ∈ C, f.eval w ≠ 0) (hq : ∀ w ∈ C, q.eval w = (s ^ 2).eval w / f.eval w) :
    ∃ r : K[X], r * vanishing C = s ^ 2 + q * f ∧ f ∣ r * vanishing C - s ^ 2 ∧
      (∀ w ∈ C, (derivative (vanishing C)).eval w = 1 →
        r.eval w = (derivative q).eval w * f.eval w + q.eval w * (derivative f).eval w) ∧
      ((s ^ 2 + q * f).natDegree < 2 * C.card → r.natDegree < C.card) := by
  obtain ⟨r, hr⟩ := vanishing_dvd C (s ^ 2) f q hf hq
  have hr' : r * vanishing C = s ^ 2 + q * f := by rw [hr, vanishing]; ring
  refine ⟨r, hr', reduce_modEq _ _ _ _ _ hr', fun w hw hd => ?_, fun hdeg => ?_⟩
  · have hz : (vanishing C).eval w = 0 := by
      rw [vanishing, eval_prod]
      exact Finset.prod_eq_zero hw (by simp)
    have := reduce_eval (s ^ 2) f q r (vanishing C) hr' w hz hd
    rw [derivative_sq_eq_zero, eval_zero, zero_add] at this
    exact this
  · by_cases hr0 : r = 0
    · rw [hr0, natDegree_zero]
      omega
    · rw [← hr', natDegree_mul hr0 (vanishing_monic C).ne_zero, vanishing_natDegree] at hdeg
      omega

/-- C31, carry to a child: if `Z_C = Z_{C0} · Z_{C1}` and `Z_{C0} = 1` on `C1`
(claim C29), the Montgomery step on `C0` (with `q = v/g` there) yields `r` with
`r · Z_{C0} ≡ v (mod g)`, i.e. `r = v · Z_{C0}⁻¹ mod g`, whose values on `C1`
are just those of `v + q g`. -/
theorem carry_to_factor [CharP K 2] [DecidableEq K] (C0 C1 : Finset K) (v g q : K[X])
    (hg : ∀ w ∈ C0, g.eval w ≠ 0) (hq : ∀ w ∈ C0, q.eval w = v.eval w / g.eval w)
    (h1 : ∀ w ∈ C1, (vanishing C0).eval w = 1) :
    ∃ r : K[X], g ∣ r * vanishing C0 - v ∧ ∀ w ∈ C1, r.eval w = (v + q * g).eval w := by
  obtain ⟨r, hr⟩ := vanishing_dvd C0 v g q hg hq
  have hr' : r * vanishing C0 = v + q * g := by rw [hr, vanishing]; ring
  refine ⟨r, reduce_modEq _ _ _ _ _ hr', fun w hw => ?_⟩
  rw [← hr', eval_mul, h1 w hw, mul_one]

end Minisketch.Montgomery
