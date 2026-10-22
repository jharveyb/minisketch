/-
Copyright (c) 2026 The minisketch developers
Distributed under the MIT software license, see the accompanying
file LICENSE or http://www.opensource.org/licenses/mit-license.php.
-/
import Minisketch.Trace
import Minisketch.Cantor

/-!
# Tower split parameters are independent

Claim C27 of `doc/spec/claims.md`: the split parameters
`β_(sl+j) = basis · x^l · z^j` (`l < b/s`, `j < s`; `TowerBetas` in
`src/sketch_impl.h`) are GF(2)-linearly independent, where `x` generates
`GF(2^b)` and `z = SubfieldGenerator(field, s)` generates `GF(2^s)`.

Both facts are used through Frobenius periods: `x^(2^k) ≠ x` for `0 < k < b`
(`x` is a root of the irreducible modulus, C7), and `z^(2^s) = z` but
`z^(2^k) ≠ z` for `0 < k < s`, which `SubfieldGenerator`'s acceptance test
implies (`subfield_generator_period`). The proof only counts roots: a
dependency would give a nonzero polynomial with more distinct roots (Frobenius
conjugates) than its degree.
-/

namespace Minisketch.Tower

open Polynomial Minisketch.Trace

variable {K : Type*} [Field K] [CharP K 2]

/-- Frobenius powers `y ↦ y^(2^k)` are injective. -/
theorem frob_inj (k : ℕ) {a c : K} (h : a ^ 2 ^ k = c ^ 2 ^ k) : a = c := by
  have := (iterateFrobenius_inj K 2 k)
  exact this (by simpa [iterateFrobenius_def] using h)

omit [CharP K 2] in
theorem pow_two_pow_add (y : K) (a c : ℕ) : y ^ 2 ^ (a + c) = (y ^ 2 ^ a) ^ 2 ^ c := by
  rw [← pow_mul, ← pow_add]

omit [Field K] [CharP K 2] in
theorem mul_lt_of_lt_div (n d : ℕ) (hd : 0 < d) (i : ℕ) (hi : i < n / d) : d * i < n :=
  calc d * i < d * (i + 1) := by nlinarith
    _ ≤ d * (n / d) := Nat.mul_le_mul_left _ hi
    _ ≤ n := Nat.mul_div_le n d

/-- Distinct Frobenius conjugates: if `y` has period at least `n`
(`y^(2^k) ≠ y` for `0 < k < n`), then `y^(2^(d i))`, `i < n / d`, are distinct. -/
theorem conj_injOn (y : K) (n d : ℕ) (hd : 0 < d) (hy : ∀ k, 0 < k → k < n → y ^ 2 ^ k ≠ y) :
    Set.InjOn (fun i => y ^ 2 ^ (d * i)) (Finset.range (n / d) : Set ℕ) := by
  intro i hi i' hi' h
  simp only [Finset.coe_range, Set.mem_Iio] at hi hi'
  have hdi := mul_lt_of_lt_div n d hd i hi
  have hdi' := mul_lt_of_lt_div n d hd i' hi'
  by_contra hne
  rcases Nat.lt_or_gt_of_ne hne with hlt | hlt
  · apply hy (d * i' - d * i) (by
      have := Nat.mul_lt_mul_of_pos_left hlt hd; omega) (by omega)
    apply frob_inj (d * i)
    simp only at h
    rw [← pow_two_pow_add, show d * i' - d * i + d * i = d * i' by
      have := Nat.mul_lt_mul_of_pos_left hlt hd; omega, ← h]
  · apply hy (d * i - d * i') (by
      have := Nat.mul_lt_mul_of_pos_left hlt hd; omega) (by omega)
    apply frob_inj (d * i')
    simp only at h
    rw [← pow_two_pow_add, show d * i - d * i' + d * i' = d * i by
      have := Nat.mul_lt_mul_of_pos_left hlt hd; omega, h]

omit [CharP K 2] in
/-- A nonzero polynomial of degree `< n` cannot vanish on `n` distinct points. -/
theorem no_roots {P : K[X]} (hP : P ≠ 0) (n : ℕ) (hdeg : P.natDegree < n) (r : ℕ → K)
    (hinj : Set.InjOn r (Finset.range n : Set ℕ)) (hr : ∀ i < n, P.eval (r i) = 0) : False := by
  classical
  have hcard : ((Finset.range n).image r).card = n := by
    rw [Finset.card_image_of_injOn hinj, Finset.card_range]
  have := Polynomial.card_le_degree_of_subset_roots (p := P) (Z := (Finset.range n).image r)
    (fun w hw => by
      obtain ⟨i, hi, rfl⟩ := Finset.mem_image.mp (Finset.mem_val.mp hw)
      rw [Polynomial.mem_roots hP]; exact hr i (Finset.mem_range.mp hi))
  omega

/-- Step (A): `1, z, …, z^(s-1)` are GF(2)-independent when `z` has period `s`. -/
theorem z_indep (s : ℕ) (hs : 0 < s) (z : K) (hz : ∀ k, 0 < k → k < s → z ^ 2 ^ k ≠ z)
    (S : Finset ℕ) (hS : S ⊆ Finset.range s) (h0 : ∑ j ∈ S, z ^ j = 0) : S = ∅ := by
  classical
  by_contra hne
  obtain ⟨j0, hj0⟩ := Finset.nonempty_iff_ne_empty.mpr hne
  set P : K[X] := ∑ j ∈ S, X ^ j
  have hcoeff : P.coeff j0 = 1 := by
    simp [P, finsetSum_coeff, coeff_X_pow, hj0]
  have hP : P ≠ 0 := fun h => by rw [h, coeff_zero] at hcoeff; exact zero_ne_one hcoeff
  have hdeg : P.natDegree < s := by
    refine lt_of_le_of_lt (natDegree_sum_le_of_forall_le _ _ (n := s - 1) fun j hj => ?_) (by omega)
    rw [natDegree_X_pow]; have := Finset.mem_range.mp (hS hj); omega
  refine no_roots hP s hdeg (fun k => z ^ 2 ^ (1 * k)) ?_ fun k _ => ?_
  · have := conj_injOn z s 1 one_pos hz
    simpa using this
  · simp only [P, eval_finsetSum, eval_pow, eval_X, one_mul]
    rw [show (∑ j ∈ S, (z ^ 2 ^ k) ^ j) = (∑ j ∈ S, z ^ j) ^ 2 ^ k by
      rw [Minisketch.Cantor.sum_pow_two_pow]
      refine Finset.sum_congr rfl fun j _ => by rw [← pow_mul, ← pow_mul, mul_comm]]
    rw [h0, zero_pow (by positivity)]

omit [CharP K 2] in
theorem pow_fixed_iter (c : K) (s : ℕ) (hc : c ^ 2 ^ s = c) (i : ℕ) : c ^ 2 ^ (s * i) = c := by
  induction i with
  | zero => simp
  | succ i ih => rw [show s * (i + 1) = s * i + s by ring, pow_two_pow_add, ih, hc]

/-- Step (B): `1, x, …, x^(b/s - 1)` are independent over `GF(2^s)` when `x` has
period `b`: coefficients `c_l` fixed by `y ↦ y^(2^s)` with `Σ c_l x^l = 0` vanish. -/
theorem x_indep (b s : ℕ) (hs : 0 < s) (x : K) (hx : ∀ k, 0 < k → k < b → x ^ 2 ^ k ≠ x)
    (c : ℕ → K) (hc : ∀ l, c l ^ 2 ^ s = c l)
    (h0 : ∑ l ∈ Finset.range (b / s), c l * x ^ l = 0) : ∀ l < b / s, c l = 0 := by
  by_contra hne
  push Not at hne
  obtain ⟨l0, hl0, hcl0⟩ := hne
  set P : K[X] := ∑ l ∈ Finset.range (b / s), C (c l) * X ^ l
  have hcoeff : P.coeff l0 = c l0 := by
    simp [P, finsetSum_coeff, coeff_C_mul, coeff_X_pow, hl0]
  have hP : P ≠ 0 := fun h => by rw [h, coeff_zero] at hcoeff; exact hcl0 hcoeff.symm
  have hdeg : P.natDegree < b / s := by
    refine lt_of_le_of_lt (natDegree_sum_le_of_forall_le _ _ (n := b / s - 1) fun l hl => ?_)
      (by omega)
    refine (natDegree_C_mul_le _ _).trans ?_
    rw [natDegree_X_pow]; have := Finset.mem_range.mp hl; omega
  refine no_roots hP (b / s) hdeg (fun i => x ^ 2 ^ (s * i)) (conj_injOn x b s hs hx) fun i _ => ?_
  simp only [P, eval_finsetSum, eval_mul, eval_C, eval_pow, eval_X]
  rw [show (∑ l ∈ Finset.range (b / s), c l * (x ^ 2 ^ (s * i)) ^ l) =
      (∑ l ∈ Finset.range (b / s), c l * x ^ l) ^ 2 ^ (s * i) by
    rw [Minisketch.Cantor.sum_pow_two_pow]
    refine Finset.sum_congr rfl fun l _ => ?_
    rw [mul_pow, pow_fixed_iter (c l) s (hc l) i, ← pow_mul, ← pow_mul, mul_comm (2 ^ (s * i))]]
  rw [h0, zero_pow (by positivity)]

/-- C27 core: if `x` has period `b` and `z` has period `s` (`s ∣ b`), the
products `x^l z^j`, `l < b/s`, `j < s`, are GF(2)-independent: no nonempty set
of them sums to zero. -/
theorem tower_indep (b s : ℕ) (hs : 0 < s) (x z : K)
    (hx : ∀ k, 0 < k → k < b → x ^ 2 ^ k ≠ x) (hzs : z ^ 2 ^ s = z)
    (hz : ∀ k, 0 < k → k < s → z ^ 2 ^ k ≠ z) (T : Finset (ℕ × ℕ))
    (hT : T ⊆ Finset.range (b / s) ×ˢ Finset.range s)
    (h0 : ∑ p ∈ T, x ^ p.1 * z ^ p.2 = 0) : T = ∅ := by
  classical
  set S : ℕ → Finset ℕ := fun l => (Finset.range s).filter fun j => (l, j) ∈ T
  set c : ℕ → K := fun l => ∑ j ∈ S l, z ^ j
  have hsum : ∑ l ∈ Finset.range (b / s), c l * x ^ l = 0 := by
    rw [← h0]
    have e1 : ∑ p ∈ T, x ^ p.1 * z ^ p.2 = ∑ p ∈ Finset.range (b / s) ×ˢ Finset.range s,
        (if p ∈ T then x ^ p.1 * z ^ p.2 else 0) := by
      rw [← Finset.sum_filter, Finset.filter_mem_eq_inter, Finset.inter_eq_right.mpr hT]
    rw [e1, Finset.sum_product]
    refine Finset.sum_congr rfl fun l _ => ?_
    simp only [c, S, Finset.sum_mul, Finset.sum_filter]
    refine Finset.sum_congr rfl fun j _ => ?_
    split_ifs <;> ring
  have hc : ∀ l, c l ^ 2 ^ s = c l := fun l => by
    simp only [c]
    rw [Minisketch.Cantor.sum_pow_two_pow]
    refine Finset.sum_congr rfl fun j _ => ?_
    rw [← pow_mul, mul_comm, pow_mul, hzs]
  have hcz := x_indep b s hs x hx c hc hsum
  ext p
  simp only [Finset.notMem_empty, iff_false]
  intro hp
  have hp' := hT hp
  rw [Finset.mem_product, Finset.mem_range, Finset.mem_range] at hp'
  have hS := z_indep s hs z hz (S p.1) (Finset.filter_subset _ _) (hcz p.1 hp'.1)
  have : p.2 ∈ S p.1 := by simp [S, hp'.2, hp]
  rw [hS] at this
  exact Finset.notMem_empty _ this

theorem sum_symmDiff {ι' : Type*} [DecidableEq ι'] (A B : Finset ι') (f : ι' → K) :
    ∑ i ∈ symmDiff A B, f i = ∑ i ∈ A, f i + ∑ i ∈ B, f i := by
  rw [symmDiff_eq_sup_sdiff_inf]
  simp only [Finset.sup_eq_union, Finset.inf_eq_inter]
  have hsub : A ∩ B ⊆ A ∪ B := Finset.inter_subset_left.trans Finset.subset_union_left
  have e1 := Finset.sum_sdiff hsub (f := f)
  have e2 := Finset.sum_union_inter (s₁ := A) (s₂ := B) (f := f)
  calc ∑ i ∈ (A ∪ B) \ (A ∩ B), f i
      = (∑ i ∈ (A ∪ B) \ (A ∩ B), f i + ∑ i ∈ A ∩ B, f i) + ∑ i ∈ A ∩ B, f i := by
        rw [add_assoc, CharTwo.add_self_eq_zero, add_zero]
    _ = _ := by rw [e1, e2]

/-- C27: the tower split parameters `β_i = basis · x^(i / s) · z^(i % s)`,
`i < b`, have distinct subset sums (they are GF(2)-independent), which is the
hypothesis of `Trace.separate` (C25). -/
theorem tower_betas_indep (b s : ℕ) (hs : 0 < s) (hsb : s ∣ b) (basis x z : K) (hbasis : basis ≠ 0)
    (hx : ∀ k, 0 < k → k < b → x ^ 2 ^ k ≠ x) (hzs : z ^ 2 ^ s = z)
    (hz : ∀ k, 0 < k → k < s → z ^ 2 ^ k ≠ z) :
    Set.InjOn (fun T : Finset ℕ => ∑ i ∈ T, basis * x ^ (i / s) * z ^ (i % s))
      ((Finset.range b).powerset : Set (Finset ℕ)) := by
  classical
  intro T hT T' hT' h
  simp only [Finset.coe_powerset, Set.mem_preimage, Set.mem_powerset_iff, Finset.coe_subset] at hT hT'
  set U := symmDiff T T'
  have hU0 : ∑ i ∈ U, basis * x ^ (i / s) * z ^ (i % s) = 0 := by
    rw [sum_symmDiff]; simp only at h; rw [h, CharTwo.add_self_eq_zero]
  have hUb : U ⊆ Finset.range b := fun i hi => by
    rcases Finset.mem_symmDiff.mp hi with ⟨h1, _⟩ | ⟨h1, _⟩
    · exact hT h1
    · exact hT' h1
  set e : ℕ → ℕ × ℕ := fun i => (i / s, i % s)
  have he : Set.InjOn e (U : Set ℕ) := fun i _ j _ hij => by
    simp only [e, Prod.mk.injEq] at hij
    rw [← Nat.div_add_mod i s, ← Nat.div_add_mod j s, hij.1, hij.2]
  have himg : ∑ p ∈ U.image e, x ^ p.1 * z ^ p.2 = 0 := by
    rw [Finset.sum_image he]
    have : basis * ∑ i ∈ U, x ^ (i / s) * z ^ (i % s) = 0 := by
      rw [Finset.mul_sum, ← hU0]
      exact Finset.sum_congr rfl fun i _ => by ring
    exact (mul_eq_zero.mp this).resolve_left hbasis
  have hsub : U.image e ⊆ Finset.range (b / s) ×ˢ Finset.range s := by
    intro p hp
    obtain ⟨i, hi, rfl⟩ := Finset.mem_image.mp hp
    have hib := Finset.mem_range.mp (hUb hi)
    simp only [e, Finset.mem_product, Finset.mem_range]
    exact ⟨Nat.div_lt_div_of_lt_of_dvd hsb hib, Nat.mod_lt i hs⟩
  have hempty := tower_indep b s hs x z hx hzs hz _ hsub himg
  have hU : U = ∅ := Finset.image_eq_empty.mp hempty
  exact symmDiff_eq_bot.mp hU

/-! ### `SubfieldGenerator` gives an element of period `s` -/

section Period

variable (w : K)

/-- `a` is a Frobenius period of `w`. -/
def Fixed (a : ℕ) : Prop := w ^ 2 ^ a = w

omit [CharP K 2] in
theorem fixed_add {a c : ℕ} (ha : Fixed w a) (hc : Fixed w c) : Fixed w (a + c) := by
  unfold Fixed at *; rw [pow_two_pow_add, ha, hc]

omit [CharP K 2] in
theorem fixed_mul {a : ℕ} (ha : Fixed w a) (m : ℕ) : Fixed w (a * m) := by
  induction m with
  | zero => simp [Fixed]
  | succ m ih => rw [Nat.mul_succ]; exact fixed_add w ih ha

theorem fixed_sub {a c : ℕ} (ha : Fixed w a) (hc : Fixed w c) (h : c ≤ a) : Fixed w (a - c) := by
  unfold Fixed at *
  apply frob_inj c
  rw [← pow_two_pow_add, Nat.sub_add_cancel h, ha, hc]

theorem fixed_gcd {a c : ℕ} (ha : Fixed w a) (hc : Fixed w c) : Fixed w (Nat.gcd a c) := by
  induction a, c using Nat.gcd.induction with
  | H0 n => simpa using hc
  | H1 m n hm ih =>
    rw [Nat.gcd_rec]
    refine ih ?_ ha
    rw [Nat.mod_eq_sub_mul_div]
    exact fixed_sub w hc (fixed_mul w ha _) (Nat.mul_div_le n m)

end Period

/-- The relative trace `w = Σ_{l < b/s} y^(2^(s l))` lies in `GF(2^s)`. -/
theorem rel_trace_fixed (b s : ℕ) (hsb : s ∣ b) (y : K) (hy : y ^ 2 ^ b = y) :
    (∑ l ∈ Finset.range (b / s), y ^ 2 ^ (s * l)) ^ 2 ^ s =
      ∑ l ∈ Finset.range (b / s), y ^ 2 ^ (s * l) := by
  rw [Minisketch.Cantor.sum_pow_two_pow]
  have hf : ∀ l, (y ^ 2 ^ (s * l)) ^ 2 ^ s = y ^ 2 ^ (s * (l + 1)) := fun l => by
    rw [← pow_two_pow_add, show s * l + s = s * (l + 1) by ring]
  simp_rw [hf]
  have h1 := Finset.sum_range_succ' (fun l => y ^ 2 ^ (s * l)) (b / s)
  have h2 := Finset.sum_range_succ (fun l => y ^ 2 ^ (s * l)) (b / s)
  rw [Nat.mul_div_cancel' hsb, hy] at h2
  simp only [mul_zero, pow_zero, pow_one] at h1
  have := h1.symm.trans h2
  exact add_right_cancel this

/-- C27's hypothesis on `z`: `SubfieldGenerator(field, s)` (with `s = 2^k`)
returns a relative trace `w ∈ GF(2^s)` with `w^(2^(s/2)) ≠ w`; such a `w` has
period exactly `s`. -/
theorem subfield_generator_period (s k : ℕ) (hsk : s = 2 ^ k) (hk : 1 ≤ k) (w : K)
    (hws : w ^ 2 ^ s = w) (hw : w ^ 2 ^ (s / 2) ≠ w) :
    ∀ k', 0 < k' → k' < s → w ^ 2 ^ k' ≠ w := by
  intro k' hk1 hk2 hfix
  have hg := fixed_gcd w (a := k') (c := s) hfix hws
  have hdvd : Nat.gcd k' s ∣ 2 ^ k := hsk ▸ Nat.gcd_dvd_right k' s
  obtain ⟨j, hj, hgj⟩ := (Nat.dvd_prime_pow Nat.prime_two).mp hdvd
  have hgle : Nat.gcd k' s ≤ k' := Nat.gcd_le_left _ hk1
  have hjk : j < k := by
    by_contra hge
    have : j = k := by omega
    rw [this, ← hsk] at hgj
    omega
  have hhalf : s / 2 = 2 ^ (k - 1) := by
    rw [hsk, show k = (k - 1) + 1 by omega, pow_succ, Nat.mul_div_cancel _ two_pos]; simp
  apply hw
  rw [hhalf]
  have : 2 ^ (k - 1) = 2 ^ j * 2 ^ (k - 1 - j) := by rw [← pow_add]; congr 1; omega
  rw [this]
  exact fixed_mul w (hgj ▸ hg) _

/-- C27 with C25, as the code runs it: in `GF(2^b)`, with `x` of period `b`
(the field generator), `z` the relative trace `Σ_{l<b/s} y^(2^(s l))` accepted
by `SubfieldGenerator` (`z^(2^(s/2)) ≠ z`, `s = 2^k` dividing `b`) and any
nonzero `basis`, the split parameters `β_i = basis · x^(i/s) · z^(i%s)`
separate every pair of distinct roots. -/
theorem tower_separates [Fintype K] [DecidableEq K] (b : ℕ) (hb : 1 ≤ b)
    (hq : Fintype.card K = 2 ^ b) (s k : ℕ) (hsk : s = 2 ^ k) (hk : 1 ≤ k) (hsb : s ∣ b)
    (basis x y : K) (hbasis : basis ≠ 0) (hx : ∀ k, 0 < k → k < b → x ^ 2 ^ k ≠ x)
    (hw : (∑ l ∈ Finset.range (b / s), y ^ 2 ^ (s * l)) ^ 2 ^ (s / 2) ≠
      ∑ l ∈ Finset.range (b / s), y ^ 2 ^ (s * l))
    (α α' : K) (hne : α ≠ α') :
    ∃ i, i < b ∧ tr b (basis * x ^ (i / s) * (∑ l ∈ Finset.range (b / s), y ^ 2 ^ (s * l)) ^ (i % s) * α) ≠
      tr b (basis * x ^ (i / s) * (∑ l ∈ Finset.range (b / s), y ^ 2 ^ (s * l)) ^ (i % s) * α') := by
  have hs : 0 < s := by rw [hsk]; positivity
  have hy : y ^ 2 ^ b = y := by rw [← hq]; exact FiniteField.pow_card y
  have hzs := rel_trace_fixed b s hsb y hy
  have hz := subfield_generator_period s k hsk hk _ hzs hw
  exact separate b hb hq _ (tower_betas_indep b s hs hsb basis x _ hbasis hx hzs hz) α α' hne

end Minisketch.Tower
