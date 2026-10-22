/-
Copyright (c) 2026 The minisketch developers
Distributed under the MIT software license, see the accompanying
file LICENSE or http://www.opensource.org/licenses/mit-license.php.
-/
import Mathlib

/-!
# Polynomial arithmetic identities

Claim C17 of `doc/spec/claims.md`: the transposed Karatsuba middle product used
by the BM tail check and the power-sum split (`LazyMiddleProduct` in
`src/sketch_impl.h`, `middle_product_kara` in `tests/pyminisketch.py`).
-/

namespace Minisketch.PolyArith

variable {K : Type*} [Field K] [CharP K 2]

/-- The middle product `MP(c, w)_i = Σ_{j < n} c_j w_{i+j}`. -/
def mp (c w : ℕ → K) (n i : ℕ) : K := ∑ j ∈ Finset.range n, c j * w (i + j)

/-- C17, low half: with `c = c0 + x^h c1` and windows `A = w`, `B = w[h..]`,
`MP(c, w)_i = α_i + MP(c0, A + B)_i` for `i < h`, `α = MP(c0 + c1, B)`. -/
theorem mp_low (c w : ℕ → K) (h i : ℕ) :
    mp c w (2 * h) i = mp (fun j => c j + c (h + j)) (fun j => w (h + j)) h i +
      mp c (fun j => w j + w (h + j)) h i := by
  unfold mp
  rw [two_mul, Finset.sum_range_add]
  simp only [add_mul, mul_add, Finset.sum_add_distrib]
  have e : ∀ j, w (h + (i + j)) = w (i + (h + j)) := fun j => by congr 1; ring
  have e2 : ∀ j, w (i + (h + j)) = w (h + (i + j)) := fun j => by congr 1; ring
  simp only [e2]
  linear_combination (-∑ j ∈ Finset.range h, c j * w (h + (i + j))) * CharTwo.two_eq_zero (R := K)

/-- C17, high half: `MP(c, w)_{h+i} = α_i + MP(c1, B + D)_i`, `D = w[2h..]`. -/
theorem mp_high (c w : ℕ → K) (h i : ℕ) :
    mp c w (2 * h) (h + i) = mp (fun j => c j + c (h + j)) (fun j => w (h + j)) h i +
      mp (fun j => c (h + j)) (fun j => w (h + j) + w (2 * h + j)) h i := by
  unfold mp
  rw [two_mul, Finset.sum_range_add]
  simp only [add_mul, mul_add, Finset.sum_add_distrib]
  have e1 : ∀ j, w (h + i + j) = w (h + (i + j)) := fun j => by congr 1; ring
  have e2 : ∀ j, w (h + (i + (h + j))) = w (h + h + (i + j)) := fun j => by congr 1; ring
  simp only [e1, e2]
  linear_combination (-∑ j ∈ Finset.range h, c (h + j) * w (h + (i + j))) * CharTwo.two_eq_zero (R := K)

end Minisketch.PolyArith
