/-
Copyright (c) 2026 The minisketch developers
Distributed under the MIT software license, see the accompanying
file LICENSE or http://www.opensource.org/licenses/mit-license.php.
-/
import Minisketch.BM
import Minisketch.BMTail

/-!
# Even discrepancies vanish

Claim C12 of `doc/spec/claims.md`: on any Frobenius-structured sequence
(`s_{2k} = s_k^2`), textbook Berlekamp–Massey has zero discrepancy at every even
step, so `BerlekampMasseyOdd` may skip them. The proof carries the binary key
equation `Λ S ≡ Λ_odd (mod X^{n+1})` (`Λ_odd = X Λ'` is the odd part of `Λ`)
through the BM steps; at an even position `2k` it forces the discrepancy to be
`[Λ_odd · Λ_even]_{2k}`, which is zero because an odd series times an even one
has only odd exponents.
-/

namespace Minisketch.BMEven

open PowerSeries Minisketch.BM Minisketch.BMTail

variable {K : Type*} [Field K] [CharP K 2]

/-- The series of a coefficient function. -/
noncomputable def ps (c : ℕ → K) : K⟦X⟧ := mk c

/-- The key equation up to `n`: `[Λ S]_j = [Λ_odd]_j` for `j ≤ n`. -/
def Key (s c : ℕ → K) (n : ℕ) : Prop :=
  ∀ j, j ≤ n → coeff j (ps c * series s) = coeff j (oddPart c)

omit [CharP K 2] in
/-- An odd series times an even one has no even coefficients. -/
theorem coeff_odd_mul_even (O E : K⟦X⟧) (hO : ∀ j, ¬Odd j → coeff j O = 0)
    (hE : ∀ j, Odd j → coeff j E = 0) (k : ℕ) : coeff (2 * k) (O * E) = 0 := by
  rw [coeff_mul]
  refine Finset.sum_eq_zero fun p hp => ?_
  rw [Finset.HasAntidiagonal.mem_antidiagonal] at hp
  by_cases ho : Odd p.1
  · have : Odd p.2 := by
      obtain ⟨a, ha⟩ := ho
      exact ⟨k - a - 1, by omega⟩
    rw [hE _ this, mul_zero]
  · rw [hO _ ho, zero_mul]

/-- C12, the even step: if `Λ` (with `Λ_0 = 1`) satisfies the key equation up
to `2k - 1`, then `[Λ S]_{2k} = 0`. -/
theorem even_step (s c : ℕ → K) (hs : Frobenius s) (hc0 : c 0 = 1) (k : ℕ) (hk : 1 ≤ k)
    (hkey : Key s c (2 * k - 1)) : coeff (2 * k) (ps c * series s) = 0 := by
  set T := ps c * series s
  set O := oddPart c
  have h0 : constantCoeff (ps c) = 1 := by simp [ps, hc0]
  -- [Λ T]_{2k} = T_{2k} + [Λ O]_{2k}, since T and O agree below 2k and O_{2k} = 0.
  have hO2k : coeff (2 * k) O = 0 := by
    simp [O, oddPart, coeff_mk]
  have hdiff : coeff (2 * k) (ps c * (T - O)) = coeff (2 * k) T := by
    rw [coeff_mul, Finset.sum_eq_single (0, 2 * k)]
    · rw [coeff_zero_eq_constantCoeff_apply, h0, one_mul, map_sub, hO2k, sub_zero]
    · intro p hp hne
      rw [Finset.HasAntidiagonal.mem_antidiagonal] at hp
      have : p.2 < 2 * k := by
        rcases p with ⟨a, b⟩
        simp only [ne_eq, Prod.mk.injEq] at hne hp ⊢
        omega
      rw [map_sub, hkey p.2 (by omega), sub_self, mul_zero]
    · intro h
      exact absurd (Finset.HasAntidiagonal.mem_antidiagonal.mpr (by simp)) h
  -- [Λ T]_{2k} = T_{2k} + [Λ O]_{2k}.
  have hlhs : coeff (2 * k) (ps c * T) = coeff (2 * k) T + coeff (2 * k) (ps c * O) := by
    rw [← hdiff, ← map_add, mul_sub, sub_add_cancel]
  -- [T^2]_{2k} = T_k^2 = O_k^2 = [O^2]_{2k}.
  have hsq : coeff (2 * k) (T ^ 2) = coeff (2 * k) (O ^ 2) := by
    rw [coeff_sq, coeff_sq, ite_eq_left (even_two_mul k), ite_eq_left (even_two_mul k),
      show 2 * k / 2 = k by omega, hkey k (by omega)]
  have hid := congrArg (coeff (2 * k)) (tail_identity s hs (ps c))
  rw [hlhs, map_add, coeff_sq_mul_odd_even, add_zero, hsq] at hid
  -- So T_{2k} = [O^2 + Λ O]_{2k} = [O (O + Λ)]_{2k}, and O + Λ is even.
  have hT : coeff (2 * k) T = coeff (2 * k) (O * (O + ps c)) := by
    rw [mul_add, map_add, ← sq, ← hid, mul_comm O (ps c), add_assoc, CharTwo.add_self_eq_zero,
      add_zero]
  rw [hT]
  apply coeff_odd_mul_even
  · intro j hj
    simp [O, oddPart, coeff_mk, hj]
  · intro j hj
    simp only [O, oddPart, ps, map_add, coeff_mk, ite_eq_left hj, CharTwo.add_self_eq_zero]

/-! ### Carrying the key equation through BM -/

omit [CharP K 2] in
/-- For `k > L` (and `c` supported on `[0, L]`), BM's discrepancy is the series
coefficient `[Λ S]_k`. -/
theorem disc_eq_coeff (s c : ℕ → K) (L k : ℕ) (hc : ∀ i, L < i → c i = 0) (hk : L < k) :
    disc s c L k = coeff k (ps c * series s) := by
  rw [ps, series, coeff_mul, Finset.Nat.sum_antidiagonal_eq_sum_range_succ_mk]
  simp only [coeff_mk]
  rw [← Finset.sum_range_add_sum_Ico _ (by omega : L + 1 ≤ k + 1)]
  rw [Finset.sum_eq_zero (s := Finset.Ico (L + 1) (k + 1)) fun i hi => by
    rw [Finset.mem_Ico] at hi; rw [hc i (by omega), zero_mul], add_zero]
  unfold disc
  refine Finset.sum_congr rfl fun i hi => ?_
  rw [Finset.mem_range] at hi
  rw [ite_eq_right (by omega)]

omit [CharP K 2] in
theorem coeff_update (st : St K) (d : K) (F : K⟦X⟧) (j : ℕ) :
    coeff j (ps (update st d) * F) = coeff j (ps st.C * F) -
      d / st.b * (if st.m ≤ j then coeff (j - st.m) (ps st.B * F) else 0) := by
  have hps : ps (update st d) = ps st.C - (d / st.b) • (X ^ st.m * ps st.B) := by
    ext i
    simp only [ps, update, coeff_mk, map_sub, coeff_smul, coeff_X_pow_mul', smul_eq_mul]
  rw [hps, sub_mul, map_sub, smul_mul_assoc, coeff_smul, mul_assoc, coeff_X_pow_mul', smul_eq_mul]

variable [DecidableEq K]

/-- The extra BM invariant on Frobenius-structured sequences: the key equation
for `C`, and for `B` (at its failure position, which is odd) unless `B` is still
the initial one. -/
structure Inv2 (s : ℕ → K) (n : ℕ) (st : St K) : Prop where
  keyC : Key s st.C n
  Bcase : (st.C = (init : St K).C ∧ st.B = (init : St K).B ∧ st.L = 0 ∧ st.LB = 0 ∧ st.b = 1) ∨
    (1 ≤ st.L ∧ Odd (n + 1 - st.m) ∧ Key s st.B (n - st.m) ∧
      coeff (n + 1 - st.m) (ps st.B * series s) = st.b)

omit [CharP K 2] in
theorem coeff_zero_mul_series (c s : ℕ → K) : coeff 0 (ps c * series s) = 0 := by
  simp [ps, series, coeff_mul]

theorem inv2_init (s : ℕ → K) : Inv2 s 0 (init : St K) where
  keyC := fun j hj => by
    obtain rfl : j = 0 := by omega
    rw [coeff_zero_mul_series]
    simp [oddPart]
  Bcase := Or.inl ⟨rfl, rfl, rfl, rfl, rfl⟩

/-- C12 at one step: an even position has zero discrepancy. -/
theorem even_disc (s : ℕ → K) (hs : Frobenius s) (n : ℕ) (st : St K) (h : Inv s n st)
    (h2 : Inv2 s n st) (he : Even (n + 1)) : disc s st.C st.L (n + 1) = 0 := by
  obtain ⟨k, hk⟩ := he
  have hlen := h.len
  have hm := h.m_pos
  rw [disc_eq_coeff s _ _ _ h.Csupp (by omega), show n + 1 = 2 * k by omega]
  exact even_step s st.C hs h.C0 k (by omega) fun j hj => h2.keyC j (by omega)

theorem oddPart_coeff (c : ℕ → K) (j : ℕ) : coeff j (oddPart c) = if Odd j then c j else 0 := by
  simp [oddPart, coeff_mk]

theorem inv2_step (s : ℕ → K) (hs : Frobenius s) (n : ℕ) (st : St K) (h : Inv s n st)
    (h2 : Inv2 s n st) : Inv2 s (n + 1) (step s n st) := by
  have hlen := h.len
  have hm1 := h.m_pos
  have hm2 := h.m_le
  set d := disc s st.C st.L (n + 1) with hd_def
  have hdC : coeff (n + 1) (ps st.C * series s) = d := (disc_eq_coeff s _ _ _ h.Csupp (by omega)).symm
  have hC1 : st.C (n + 1) = 0 := h.Csupp _ (by omega)
  by_cases hd : d = 0
  · have e : step s n st = { st with m := st.m + 1 } := by simp [step, ← hd_def, hd]
    rw [e]
    refine ⟨fun j hj => ?_, ?_⟩
    · rcases Nat.lt_or_ge j (n + 1) with hj' | hj'
      · exact h2.keyC j (by omega)
      · obtain rfl : j = n + 1 := by omega
        simp only
        rw [hdC, hd, oddPart_coeff, hC1]
        simp
    · rcases h2.Bcase with hb | ⟨hL, hodd, hkey, hb⟩
      · exact Or.inl hb
      · refine Or.inr ⟨hL, ?_, ?_, ?_⟩
        · simpa [show n + 1 + 1 - (st.m + 1) = n + 1 - st.m by omega] using hodd
        · simpa [show n + 1 - (st.m + 1) = n - st.m by omega] using hkey
        · simpa [show n + 1 + 1 - (st.m + 1) = n + 1 - st.m by omega] using hb
  -- A nonzero discrepancy only happens at an odd position.
  have hodd1 : Odd (n + 1) := by
    rcases Nat.even_or_odd (n + 1) with he | ho
    · exact absurd (even_disc s hs n st h h2 he) hd
    · exact ho
  -- The key equation for the updated C.
  have hkey' : Key s (update st d) (n + 1) := by
    intro j hj
    rw [coeff_update, oddPart_coeff]
    simp only [update]
    rcases h2.Bcase with ⟨hC, hB, hL, hLB, hb⟩ | ⟨hL, hoddB, hkeyB, hbB⟩
    · have hmn : st.m = n + 1 := by omega
      rcases Nat.lt_or_ge j (n + 1) with hj' | hj'
      · have hmj : ¬st.m ≤ j := by omega
        simp only [hmj, ↓reduceIte, mul_zero, sub_zero]
        rw [h2.keyC j (by omega), oddPart_coeff]
      · obtain rfl : j = n + 1 := by omega
        have hmj : st.m ≤ n + 1 := by omega
        simp only [hmj, hodd1, ↓reduceIte]
        rw [hdC, hmn, Nat.sub_self, coeff_zero_mul_series, hC1, hB, hb]
        simp [init, CharTwo.neg_eq]
    · have hmeven : Even st.m := by
        obtain ⟨a, ha⟩ := hodd1
        obtain ⟨b', hb'⟩ := hoddB
        exact ⟨a - b', by omega⟩
      rcases Nat.lt_or_ge j (n + 1) with hj' | hj'
      · rw [h2.keyC j (by omega), oddPart_coeff]
        by_cases hmj : st.m ≤ j
        · rw [hkeyB (j - st.m) (by omega), oddPart_coeff]
          by_cases hoj : Odd j
          · have : Odd (j - st.m) := by
              obtain ⟨a, ha⟩ := hoj
              obtain ⟨c', hc'⟩ := hmeven
              exact ⟨a - c', by omega⟩
            simp only [hmj, hoj, this, ↓reduceIte]
          · have : ¬Odd (j - st.m) := by
              intro ho
              apply hoj
              obtain ⟨a, ha⟩ := ho
              obtain ⟨c', hc'⟩ := hmeven
              exact ⟨a + c', by omega⟩
            simp only [hmj, hoj, this, ↓reduceIte, mul_zero, sub_zero]
        · simp only [hmj, ↓reduceIte, mul_zero, sub_zero]
      · obtain rfl : j = n + 1 := by omega
        have hmj : st.m ≤ n + 1 := by omega
        simp only [hmj, hodd1, ↓reduceIte]
        rw [hdC, hbB, div_mul_cancel₀ _ h.b_ne, sub_self, hC1, h.Bsupp _ (by omega), mul_zero,
          sub_zero]
  by_cases h2L : 2 * st.L ≤ n
  · have e : step s n st = St.mk (update st d) (n + 1 - st.L) st.C st.L d 1 := by
      simp [step, ← hd_def, hd, h2L]
    rw [e]
    refine ⟨hkey', Or.inr ⟨by dsimp only; omega, ?_, ?_, ?_⟩⟩
    · simpa using hodd1
    · simpa using h2.keyC
    · simpa using hdC
  · have e : step s n st = { st with C := update st d, m := st.m + 1 } := by
      simp [step, ← hd_def, hd, h2L]
    rw [e]
    refine ⟨hkey', ?_⟩
    rcases h2.Bcase with ⟨_, _, hL, _, _⟩ | ⟨hL, hodd, hkey, hb⟩
    · omega
    · refine Or.inr ⟨hL, ?_, ?_, ?_⟩
      · simpa [show n + 1 + 1 - (st.m + 1) = n + 1 - st.m by omega] using hodd
      · simpa [show n + 1 - (st.m + 1) = n - st.m by omega] using hkey
      · simpa [show n + 1 + 1 - (st.m + 1) = n + 1 - st.m by omega] using hb

theorem inv2_bm (s : ℕ → K) (hs : Frobenius s) (hs0 : s 0 = 1) (n : ℕ) : Inv2 s n (bm s n) := by
  induction n with
  | zero => exact inv2_init s
  | succ n ih => exact inv2_step s hs n _ (inv_bm s hs0 n) ih

/-- C12: on a Frobenius-structured sequence, textbook BM's discrepancy at every
even position is zero, so skipping the even steps changes nothing but the shift
counter. -/
theorem even_discrepancy_zero (s : ℕ → K) (hs : Frobenius s) (hs0 : s 0 = 1) (n : ℕ)
    (he : Even (n + 1)) : disc s (bm s n).C (bm s n).L (n + 1) = 0 :=
  even_disc s hs n _ (inv_bm s hs0 n) (inv2_bm s hs hs0 n) he

/-! ### Early termination (C48) -/

omit [CharP K 2] in
/-- If every discrepancy of the state after `P` terms vanishes on `(P, P + t]`,
BM ends at `P + t` with the same connection polynomial and length. -/
theorem bm_stable (s : ℕ → K) (P t : ℕ)
    (hz : ∀ k, P < k → k ≤ P + t → disc s (bm s P).C (bm s P).L k = 0) :
    (bm s (P + t)).C = (bm s P).C ∧ (bm s (P + t)).L = (bm s P).L := by
  induction t with
  | zero => exact ⟨rfl, rfl⟩
  | succ t ih =>
    obtain ⟨hC, hL⟩ := ih fun k hk1 hk2 => hz k hk1 (by omega)
    have hd : disc s (bm s (P + t)).C (bm s (P + t)).L (P + t + 1) = 0 := by
      rw [hC, hL]; exact hz _ (by omega) (by omega)
    have e : bm s (P + (t + 1)) = { bm s (P + t) with m := (bm s (P + t)).m + 1 } := by
      show step s (P + t) (bm s (P + t)) = _
      simp [step, hd]
    rw [e]
    exact ⟨hC, hL⟩

/-- C48: early termination is exact. After `P ≥ 2L` terms (`L` the current
length), if the odd syndromes satisfy the `Λ²` recurrence at every odd position
in `(P, N]` (the bulk tail check, C15), textbook BM on all `N` terms returns the
same connection polynomial and length as the state after `P` terms. With C12
(even steps are no-ops) this is `BerlekampMasseyOdd`'s result. -/
theorem bm_odd_early_eq_bm (s : ℕ → K) (hs : Frobenius s) (hs0 : s 0 = 1) (P N : ℕ)
    (hPN : P ≤ N) (hP : 2 * (bm s P).L ≤ P)
    (hcheck : ∀ j, P < 2 * j + 1 → 2 * j + 1 ≤ N →
      coeff (2 * j + 1) (ps (bm s P).C ^ 2 * oddPart s) = 0) :
    (bm s N).C = (bm s P).C ∧ (bm s N).L = (bm s P).L := by
  have h := inv_bm s hs0 P
  have h0 : constantCoeff (ps (bm s P).C) = 1 := by simp [ps, h.C0]
  have hΛ : ∀ i, (bm s P).L < i → coeff i (ps (bm s P).C) = 0 := fun i hi => by
    simp [ps, h.Csupp i hi]
  have hpre : ∀ k, (bm s P).L < k → k ≤ P → coeff k (ps (bm s P).C * series s) = 0 :=
    fun k hk1 hk2 => by rw [← disc_eq_coeff s _ _ _ h.Csupp hk1]; exact h.gen k hk1 hk2
  have htail := (tail_check_iff s hs (ps (bm s P).C) h0 (bm s P).L P N hΛ hP hpre).mpr hcheck
  obtain ⟨t, rfl⟩ : ∃ t, N = P + t := ⟨N - P, by omega⟩
  exact bm_stable s P t fun k hk1 hk2 => by
    rw [disc_eq_coeff s _ _ _ h.Csupp (by omega)]; exact htail k hk1 hk2

end Minisketch.BMEven
