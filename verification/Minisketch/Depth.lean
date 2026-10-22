/-
Copyright (c) 2026 The minisketch developers
Distributed under the MIT software license, see the accompanying
file LICENSE or http://www.opensource.org/licenses/mit-license.php.
-/
import Minisketch.Trace

/-!
# Depth bound of the root-finding tree

Each split of the trace algorithm sorts a node's roots by one trace bit
`Tr(β_i α)`. A node reached after `k` splits therefore only holds roots with
`Tr(β_i α) = c_i` for `i < k` and fixed bits `c_i`. This file proves claim C24
of `doc/spec/claims.md`: in `GF(2^b)` with independent split parameters, at most
`2^(b-k)` field elements satisfy those `k` constraints, which is what the C++
safety check `(deg - 1) >> (B - depth) == 0` in `RecFindRoots` asserts.

The proof avoids linear algebra: by C25 (`Trace.separate`) the full bit vector
`(Tr(β_i α))_{i < b}` determines `α`. On the constrained set the first `k` bits
are fixed, so the remaining `b - k` bits already determine `α`, giving an
injection into `Fin (b - k) → Bool`.
-/

namespace Minisketch.Depth

open Minisketch.Trace

variable {K : Type*} [Field K] [CharP K 2]

/-- C24: in `GF(2^b)` with GF(2)-independent split parameters `β_0..β_{b-1}`
(distinct subset sums, as in `Trace.separate`), the elements `α` with
`Tr(β_i α) = c_i` for all `i < k` number at most `2^(b-k)`. No `k ≤ b` is
needed: for `k ≥ b` the bound is `1`. -/
theorem depth_bound [Fintype K] [DecidableEq K] (b : ℕ) (hq : Fintype.card K = 2 ^ b)
    (β : ℕ → K)
    (hind : Set.InjOn (fun T : Finset ℕ => ∑ i ∈ T, β i) ((Finset.range b).powerset : Set (Finset ℕ)))
    (k : ℕ) (c : ℕ → K) :
    (Finset.univ.filter fun α : K => ∀ i < k, tr b (β i * α) = c i).card ≤ 2 ^ (b - k) := by
  -- `b ≥ 1`, since a field has at least two elements.
  have hb : 1 ≤ b := by
    rcases Nat.eq_zero_or_pos b with rfl | h
    · have := Fintype.one_lt_card (α := K); rw [hq] at this; simp at this
    · exact h
  -- Every trace bit is 0 or 1 (C20), as every element satisfies `a^(2^b) = a`.
  have hmem : ∀ a : K, tr b a = 0 ∨ tr b a = 1 := fun a =>
    tr_mem b a (by rw [← hq]; exact FiniteField.pow_card a)
  -- The bits `k, k+1, …, b-1` as a Boolean vector.
  let g : K → (Fin (b - k) → Bool) := fun α j => decide (tr b (β (k + j) * α) = 1)
  -- `g` is injective on the constrained set.
  have hinj : Set.InjOn g (Finset.univ.filter fun α : K => ∀ i < k, tr b (β i * α) = c i) := by
    intro α hα α' hα' hg
    simp only [Finset.coe_filter, Finset.mem_univ, true_and, Set.mem_ofPred_eq] at hα hα'
    by_contra hne
    obtain ⟨i, hib, hi⟩ := separate b hb hq β hind α α' hne
    apply hi
    rcases Nat.lt_or_ge i k with hik | hik
    · -- A fixed bit: both equal `c_i`.
      rw [hα i hik, hα' i hik]
    · -- A free bit: `g` records it, and a bit is determined by whether it is 1.
      have hj := congrFun hg ⟨i - k, by omega⟩
      simp only [g, Nat.add_sub_cancel' hik, decide_eq_decide] at hj
      rcases hmem (β i * α) with h0 | h1 <;> rcases hmem (β i * α') with h0' | h1'
      · rw [h0, h0']
      · exact absurd (hj.mpr h1') (by rw [h0]; exact zero_ne_one)
      · exact absurd (hj.mp h1) (by rw [h0']; exact zero_ne_one)
      · rw [h1, h1']
  calc _ ≤ (Finset.univ : Finset (Fin (b - k) → Bool)).card :=
        Finset.card_le_card_of_injOn g (fun _ _ => Finset.mem_coe.mpr (Finset.mem_univ _)) hinj
    _ = 2 ^ (b - k) := by simp

/-- C24, as the C++ check uses it: a node's root set `R`, all of whose roots
satisfy the first `k` trace constraints, has `|R| ≤ 2^(b-k)`, i.e.
`(|R| - 1) >> (b - k) = 0` for nonempty `R`. -/
theorem depth_bound_roots [Fintype K] [DecidableEq K] (b : ℕ) (hq : Fintype.card K = 2 ^ b)
    (β : ℕ → K)
    (hind : Set.InjOn (fun T : Finset ℕ => ∑ i ∈ T, β i) ((Finset.range b).powerset : Set (Finset ℕ)))
    (k : ℕ) (c : ℕ → K) (R : Finset K) (hR : ∀ α ∈ R, ∀ i < k, tr b (β i * α) = c i) :
    R.card ≤ 2 ^ (b - k) :=
  (Finset.card_le_card fun α hα => by simpa using hR α hα).trans
    (depth_bound b hq β hind k c)

end Minisketch.Depth
