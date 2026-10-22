/-
Copyright (c) 2026 The minisketch developers
Distributed under the MIT software license, see the accompanying
file LICENSE or http://www.opensource.org/licenses/mit-license.php.
-/
import Mathlib

/-!
# Berlekamp–Massey

An executable model of textbook Berlekamp–Massey (as in `berlekamp_massey` in
`tests/pyminisketch.py` and `BerlekampMasseyRef` in `src/test_refimpl.h`) and
Massey's theorem that it returns a shortest LFSR: claim C5 of
`doc/spec/claims.md`.

Sequences are `s : ℕ → K` with terms `s 1, s 2, …`; `s 0` is never read by an
LFSR condition. A connection polynomial is a coefficient function `c` with
`c 0 = 1`; it generates `s_1..s_n` with length `L` if every `k` with `L < k ≤ n`
has zero discrepancy `Σ_{i ≤ L} c_i s_{k-i}`.
-/

namespace Minisketch.BM

variable {K : Type*} [Field K]

/-- The discrepancy of `c` (length `L`) at position `k`. -/
def disc (s c : ℕ → K) (L k : ℕ) : K := ∑ i ∈ Finset.range (L + 1), c i * s (k - i)

/-- `c` with length `L` generates `s_1..s_n`. -/
def Gen (s c : ℕ → K) (L n : ℕ) : Prop := ∀ k, L < k → k ≤ n → disc s c L k = 0

theorem disc_extend (s c : ℕ → K) (L L2 k : ℕ) (hc : ∀ i, L < i → c i = 0) (h : L ≤ L2) :
    disc s c L2 k = disc s c L k := by
  unfold disc
  rw [← Finset.sum_range_add_sum_Ico _ (by omega : L + 1 ≤ L2 + 1)]
  rw [Finset.sum_eq_zero (s := Finset.Ico (L + 1) (L2 + 1)) fun i hi => by
    rw [Finset.mem_Ico] at hi; rw [hc i (by omega), zero_mul], add_zero]

theorem disc_succ' (s c : ℕ → K) (L k : ℕ) (hc0 : c 0 = 1) :
    disc s c L k = s k + ∑ j ∈ Finset.range L, c (j + 1) * s (k - (j + 1)) := by
  rw [disc, Finset.sum_range_succ', hc0, one_mul, Nat.sub_zero, add_comm]

/-- Massey's lemma: if `c` (length `L < n + 1`) generates `s_1..s_n` but not
`s_{n+1}`, every LFSR generating `s_1..s_{n+1}` has length at least `n + 1 - L`. -/
theorem massey_lemma (s c c' : ℕ → K) (L L' n : ℕ) (hc0 : c 0 = 1) (hc'0 : c' 0 = 1)
    (hg : Gen s c L n) (hd : disc s c L (n + 1) ≠ 0) (hg' : Gen s c' L' (n + 1)) :
    n + 1 ≤ L + L' := by
  by_contra hlt
  push Not at hlt
  apply hd
  -- Each s_{n+1-i}, i ≤ L, is predicted by c' (its position exceeds L').
  have key : ∀ i ∈ Finset.range (L + 1), s (n + 1 - i) =
      -∑ j ∈ Finset.range L', c' (j + 1) * s (n - j - i) := by
    intro i hi
    rw [Finset.mem_range] at hi
    have h := hg' (n + 1 - i) (by omega) (by omega)
    rw [disc_succ' s c' L' _ hc'0] at h
    rw [eq_neg_iff_add_eq_zero, ← h]
    congr 1
    refine Finset.sum_congr rfl fun j hj => ?_
    rw [Finset.mem_range] at hj
    congr 2
    omega
  unfold disc
  rw [Finset.sum_congr rfl fun i hi => by rw [key i hi]]
  simp_rw [mul_neg, Finset.mul_sum, Finset.sum_neg_distrib]
  rw [Finset.sum_comm, neg_eq_zero]
  refine Finset.sum_eq_zero fun j hj => ?_
  rw [Finset.mem_range] at hj
  have h := hg (n - j) (by omega) (by omega)
  unfold disc at h
  calc ∑ i ∈ Finset.range (L + 1), c i * (c' (j + 1) * s (n - j - i))
      = c' (j + 1) * ∑ i ∈ Finset.range (L + 1), c i * s (n - j - i) := by
        rw [Finset.mul_sum]; exact Finset.sum_congr rfl fun i _ => by ring
    _ = 0 := by rw [h, mul_zero]

/-! ### The algorithm -/

variable [DecidableEq K]

/-- BM state: connection polynomial `C` of length `L`, and the last
connection polynomial `B` (length `LB`) before the most recent length change,
its failing discrepancy `b`, and the shift `m` since then. -/
structure St (K : Type*) where
  C : ℕ → K
  L : ℕ
  B : ℕ → K
  LB : ℕ
  b : K
  m : ℕ

/-- `C - (d/b) x^m B`. -/
def update (st : St K) (d : K) : ℕ → K :=
  fun i => st.C i - d / st.b * (if st.m ≤ i then st.B (i - st.m) else 0)

/-- One BM step, processing `s_{n+1}`. -/
def step (s : ℕ → K) (n : ℕ) (st : St K) : St K :=
  let d := disc s st.C st.L (n + 1)
  if d = 0 then { st with m := st.m + 1 }
  else if 2 * st.L ≤ n then
    { C := update st d, L := n + 1 - st.L, B := st.C, LB := st.L, b := d, m := 1 }
  else { st with C := update st d, m := st.m + 1 }

def init : St K :=
  { C := fun i => if i = 0 then 1 else 0, L := 0, B := fun i => if i = 0 then 1 else 0,
    LB := 0, b := 1, m := 1 }

/-- Textbook BM after `n` terms. -/
def bm (s : ℕ → K) : ℕ → St K
  | 0 => init
  | n + 1 => step s n (bm s n)

/-- The BM invariant after `n` terms. -/
structure Inv (s : ℕ → K) (n : ℕ) (st : St K) : Prop where
  C0 : st.C 0 = 1
  Csupp : ∀ i, st.L < i → st.C i = 0
  gen : Gen s st.C st.L n
  minimal : ∀ c' L', c' 0 = 1 → Gen s c' L' n → st.L ≤ L'
  B0 : st.B 0 = 1
  Bsupp : ∀ i, st.LB < i → st.B i = 0
  m_pos : 1 ≤ st.m
  m_le : st.m ≤ n + 1
  len : st.L + st.LB = n + 1 - st.m
  Bgen : Gen s st.B st.LB (n - st.m)
  Bdisc : disc s st.B st.LB (n + 1 - st.m) = st.b
  b_ne : st.b ≠ 0

theorem Gen.mono {s c : ℕ → K} {L n n' : ℕ} (h : Gen s c L n) (hn : n' ≤ n) : Gen s c L n' :=
  fun k hk1 hk2 => h k hk1 (hk2.trans hn)

/-- The shifted term of the update contributes `disc B` at `k - m`. -/
theorem disc_update (s : ℕ → K) (st : St K) (d : K) (L2 k : ℕ)
    (hB : ∀ i, st.LB < i → st.B i = 0) (hL2 : st.m + st.LB ≤ L2) (hk : st.m ≤ k) :
    disc s (update st d) L2 k =
      disc s st.C L2 k - d / st.b * disc s st.B st.LB (k - st.m) := by
  unfold disc update
  simp_rw [sub_mul, Finset.sum_sub_distrib, mul_assoc, ← Finset.mul_sum]
  congr 2
  rw [← Finset.sum_range_add_sum_Ico _ (by omega : st.m ≤ L2 + 1)]
  rw [Finset.sum_eq_zero fun i hi => by
    rw [Finset.mem_range] at hi; rw [ite_eq_right (by omega), zero_mul], zero_add]
  rw [Finset.sum_Ico_eq_sum_range]
  rw [← Finset.sum_range_add_sum_Ico _ (by omega : st.LB + 1 ≤ L2 + 1 - st.m)]
  rw [Finset.sum_eq_zero (s := Finset.Ico _ _) fun i hi => by
    rw [Finset.mem_Ico] at hi; rw [ite_eq_left (by omega), show st.m + i - st.m = i by omega,
      hB i (by omega), zero_mul], add_zero]
  refine Finset.sum_congr rfl fun i hi => ?_
  rw [ite_eq_left (by omega), show st.m + i - st.m = i by omega, show k - (st.m + i) = k - st.m - i by omega]

theorem inv_init (s : ℕ → K) (hs0 : s 0 = 1) : Inv s 0 (init : St K) where
  C0 := by simp [init]
  Csupp := fun i hi => by simp [init]; omega
  gen := fun k hk1 hk2 => by omega
  minimal := fun _ _ _ _ => Nat.zero_le _
  B0 := by simp [init]
  Bsupp := fun i hi => by simp [init]; omega
  m_pos := le_rfl
  m_le := le_rfl
  len := rfl
  Bgen := fun k hk1 hk2 => by
    have : (init : St K).m = 1 := rfl
    omega
  Bdisc := by simp [init, disc, hs0]
  b_ne := one_ne_zero

theorem update_zero (st : St K) (d : K) (hm : 1 ≤ st.m) : update st d 0 = st.C 0 := by
  simp only [update]
  rw [ite_eq_right (by omega), mul_zero, sub_zero]

theorem update_supp (st : St K) (d : K) (L2 : ℕ) (hC : ∀ i, st.L < i → st.C i = 0)
    (hB : ∀ i, st.LB < i → st.B i = 0) (h1 : st.L ≤ L2) (h2 : st.m + st.LB ≤ L2) :
    ∀ i, L2 < i → update st d i = 0 := by
  intro i hi
  simp only [update]
  rw [hC i (by omega), ite_eq_left (by omega), hB _ (by omega), mul_zero, sub_zero]

/-- The updated polynomial `C - (d/b) x^m B` generates `s_1..s_{n+1}` with any
length `L2 ≥ max(L, m + LB)`: it cancels the discrepancy `d` at `n + 1`. -/
theorem update_gen (s : ℕ → K) (n : ℕ) (st : St K) (h : Inv s n st) (L2 : ℕ)
    (h1 : st.L ≤ L2) (h2 : st.m + st.LB ≤ L2) :
    Gen s (update st (disc s st.C st.L (n + 1))) L2 (n + 1) := by
  intro k hk1 hk2
  rw [disc_update s st _ L2 k h.Bsupp h2 (by omega), disc_extend s st.C st.L L2 k h.Csupp h1]
  rcases Nat.lt_or_ge k (n + 1) with hk | hk
  · rw [h.gen k (by omega) (by omega), h.Bgen (k - st.m) (by omega) (by omega)]
    simp
  · obtain rfl : k = n + 1 := by omega
    rw [h.Bdisc]
    field_simp [h.b_ne]
    ring

theorem inv_step (s : ℕ → K) (n : ℕ) (st : St K) (h : Inv s n st) :
    Inv s (n + 1) (step s n st) := by
  have hlen := h.len
  have hm1 := h.m_pos
  have hm2 := h.m_le
  by_cases hd : disc s st.C st.L (n + 1) = 0
  · -- zero discrepancy: only the shift grows
    have e : step s n st = { st with m := st.m + 1 } := by simp [step, hd]
    rw [e]
    exact {
      C0 := h.C0, Csupp := h.Csupp,
      gen := fun k hk1 hk2 => by
        rcases Nat.lt_or_ge k (n + 1) with hk | hk
        · exact h.gen k hk1 (by omega)
        · obtain rfl : k = n + 1 := by omega
          exact hd
      minimal := fun c' L' hc' hg' => h.minimal c' L' hc' (hg'.mono (by omega)),
      B0 := h.B0, Bsupp := h.Bsupp,
      m_pos := by simp, m_le := by simp; omega,
      len := by simp; omega,
      Bgen := by simpa [show n + 1 - (st.m + 1) = n - st.m by omega] using h.Bgen,
      Bdisc := by simpa [show n + 1 + 1 - (st.m + 1) = n + 1 - st.m by omega] using h.Bdisc,
      b_ne := h.b_ne }
  by_cases h2 : 2 * st.L ≤ n
  · -- length change: L becomes n + 1 - L
    have e : step s n st = St.mk (update st (disc s st.C st.L (n + 1))) (n + 1 - st.L)
        st.C st.L (disc s st.C st.L (n + 1)) 1 := by
      simp [step, hd, h2]
    rw [e]
    exact {
      C0 := by simp only; rw [update_zero st _ hm1, h.C0],
      Csupp := update_supp st _ _ h.Csupp h.Bsupp (by dsimp only; omega) (by dsimp only; omega),
      gen := update_gen s n st h _ (by dsimp only; omega) (by dsimp only; omega),
      minimal := fun c' L' hc' hg' => by
        have := massey_lemma s st.C c' st.L L' n h.C0 hc' h.gen hd hg'
        simp only; omega,
      B0 := h.C0, Bsupp := h.Csupp,
      m_pos := le_rfl, m_le := by dsimp only; omega,
      len := by simp only; omega,
      Bgen := by simpa using h.gen,
      Bdisc := by simp,
      b_ne := hd }
  · -- no length change
    have e : step s n st = { st with C := update st (disc s st.C st.L (n + 1)), m := st.m + 1 } := by
      simp [step, hd, h2]
    rw [e]
    exact {
      C0 := by simp only; rw [update_zero st _ hm1, h.C0],
      Csupp := update_supp st _ _ h.Csupp h.Bsupp le_rfl (by dsimp only; omega),
      gen := update_gen s n st h _ le_rfl (by dsimp only; omega),
      minimal := fun c' L' hc' hg' => h.minimal c' L' hc' (hg'.mono (by omega)),
      B0 := h.B0, Bsupp := h.Bsupp,
      m_pos := by simp, m_le := by simp; omega,
      len := by simp; omega,
      Bgen := by simpa [show n + 1 - (st.m + 1) = n - st.m by omega] using h.Bgen,
      Bdisc := by simpa [show n + 1 + 1 - (st.m + 1) = n + 1 - st.m by omega] using h.Bdisc,
      b_ne := h.b_ne }

theorem inv_bm (s : ℕ → K) (hs0 : s 0 = 1) (n : ℕ) : Inv s n (bm s n) := by
  induction n with
  | zero => exact inv_init s hs0
  | succ n ih => exact inv_step s n _ ih

/-- C5, Massey's theorem: BM's connection polynomial after `n` terms generates
`s_1..s_n`, and no LFSR generating them is shorter. -/
theorem bm_minimal (s : ℕ → K) (hs0 : s 0 = 1) (n : ℕ) :
    (bm s n).C 0 = 1 ∧ Gen s (bm s n).C (bm s n).L n ∧
      ∀ c' L', c' 0 = 1 → Gen s c' L' n → (bm s n).L ≤ L' :=
  let h := inv_bm s hs0 n
  ⟨h.C0, h.gen, h.minimal⟩

end Minisketch.BM
