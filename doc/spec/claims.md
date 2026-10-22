# Claim registry

Every mathematical claim the library's correctness or the docs rely on, with
where it is stated, how it is checked, and how far that check goes. This is the
index that ties together the Lean proofs (`verification/`), the Python oracle
(`tests/pyminisketch.py`), the C++ tests and the spec (`doc/spec/`).

`verification/check_claims.py` parses this file: every `lean:` name must exist
in the Lean project, and a claim marked `proved` must depend on no `sorry` and
no axiom beyond Lean's standard three (`propext`, `Classical.choice`,
`Quot.sound`), plus `Lean.ofReduceBool` where `native_decide: allowed`.

Field meanings:

- **kind**: `textbook` (published, standard), `novel` (derived for this
  library, or a new combination), `heuristic` (a model or measurement, not a
  theorem).
- **status**: `proved` (Lean, no sorry), `stated` (Lean statement present,
  proof pending or `sorry`), `tested` (checked by tests only), `open` (neither
  yet), `not-a-theorem` (heuristic; the evidence is listed instead).
- **depends**: claims a Lean proof takes as hypotheses. A `proved` claim is
  only unconditional once everything it depends on is proved too.
- **lean / oracle / tests / spec**: the Lean theorem, the pyminisketch function,
  the C++ test group or fuzz target, the spec section. `-` means none yet;
  `(planned)` marks a name that does not exist yet.

Notation: K = GF(2^b) with b the field size, char K = 2. S(x) = Σ_k s_k x^k is
the power-sum generating function of a set M ⊆ K\{0}, s_k = Σ_{m∈M} m^k. Λ(x) =
Π_{m∈M}(1 − m x) is the locator, f = rev(Λ) = Π(x − m). Tr is the absolute
trace K → GF(2).

## Foundations

### C1: sketches are linear
- statement: s_k(A) + s_k(B) = s_k(A Δ B) for all k; so XOR of sketches is the sketch of the symmetric difference.
- source: doc/math.md:30
- kind: textbook
- lean: Minisketch.Sketch.psum_symmDiff
- oracle: Minisketch.merge
- tests: TestSyndromes; fuzz roundtrip
- spec: 02-sketch
- status: proved

### C2: generating-function identity
- statement: as power series, Λ_R(x)·S_R(x) = −x·Λ_R′(x) for any finite set R (Λ_R = Π(1 − αx), S_R = Σ_{k≥1} s_k x^k); in characteristic 2 this is the odd part of Λ_R. In particular Λ·S is a polynomial of degree ≤ |R|.
- source: doc/math.md:40
- kind: textbook
- lean: Minisketch.Decode.loc_mul_series, Minisketch.Decode.neg_euler
- oracle: -
- tests: TestBerlekampMassey
- spec: 02-sketch
- status: proved

### C3: 2n power sums determine a set of at most n elements
- statement: a nonempty D ⊆ K\{0} with |D| ≤ 2n has s_k(D) ≠ 0 for some 1 ≤ k ≤ 2n (Vandermonde). With C1: if |M|, |M'| ≤ n and s_k(M) = s_k(M') for k = 1..2n then M = M'. With C4 the odd k < 2n suffice.
- source: doc/math.md:9, doc/math.md:54
- kind: textbook
- lean: Minisketch.Sketch.eq_of_odd_psum_eq
- oracle: -
- tests: TestExhaustive (count check)
- spec: 02-sketch
- status: proved

### C4: Frobenius halving
- statement: s_{2k} = s_k^2.
- source: doc/math.md:91, doc/decode_explained.md:34
- kind: textbook
- lean: Minisketch.Sketch.psum_two_mul
- oracle: Minisketch.decode (reconstructs even sums)
- tests: TestSyndromes
- spec: 02-sketch
- status: proved

### C5: BM finds the shortest LFSR
- statement: textbook Berlekamp–Massey (Massey 1969) after n terms returns a connection polynomial C (C_0 = 1) of length L that generates s_1..s_n, and every LFSR generating s_1..s_n has length ≥ L. Key step (Massey's lemma): if a length-L LFSR generates s_1..s_n but not s_{n+1}, every LFSR generating s_1..s_{n+1} has length ≥ n + 1 − L. With N = 2c and a true set of size ≤ c this gives Λ (with C3).
- source: doc/math.md:102
- kind: textbook
- lean: Minisketch.BM.bm_minimal, Minisketch.BM.massey_lemma
- oracle: berlekamp_massey
- tests: TestBerlekampMassey
- spec: 04-bm
- status: proved

### C6: elements are the roots of rev(Λ)
- statement: the roots of f = x^n Λ(1/x) are exactly the m ∈ M.
- source: doc/math.md:106
- kind: textbook
- lean: Minisketch.Sketch.reverse_locator_roots
- oracle: Minisketch.decode
- tests: TestFindRoots
- spec: 05-roots
- status: open

### C7: field moduli
- statement: each of the 63 moduli x^b + MOD (b = 2..64) in doc/moduli.md is irreducible over GF(2), and is the minimal-weight, lexicographically first such polynomial of degree b.
- source: doc/moduli.md, doc/gen_params.sage:25
- kind: textbook (computation)
- lean: Minisketch.Moduli.all_irreducible
- native_decide: allowed
- oracle: gf2x_irreducible
- tests: pyminisketch TestModuli (Rabin test for b = 2..64, plus no lower-weight and no lexicographically smaller same-weight irreducible); TestFieldOps (bits ≤ 8, via all inverses)
- spec: 01-field
- status: tested

### C8: Log2Factorial
- statement: Log2Factorial(x) = floor(log2 x!) for x ≤ 57 and ≤ floor(log2 x!) for all x < 2^32.
- source: src/false_positives.h:18
- kind: novel (constants)
- lean: Minisketch.FalsePositive.log2factorial_le
- oracle: base_fpbits
- tests: tests/differential.py check_capacity_functions
- spec: 02-sketch
- status: tested

### C9: BaseFPBits
- statement: BaseFPBits(b, c) = b·c − ceil(log2 Σ_{k≤c} C(2^b − 1, k)) whenever that value is < 256, and both are ≥ 256 otherwise.
- source: src/false_positives.h:46, doc/gen_basefpbits.sage
- kind: novel (tables)
- lean: -
- oracle: base_fpbits
- tests: tests/differential.py check_capacity_functions (b = 2..64, c ≤ 128, through compute_capacity)
- spec: 02-sketch
- status: tested

### C10: ComputeMaxElements terminates and inverts ComputeCapacity
- statement: decreasing max_elements by N decreases ComputeCapacity by at most N; the loop returns max{m ≤ c : ComputeCapacity(m) ≤ c}, or 0.
- source: src/false_positives.h:102
- kind: novel
- lean: Minisketch.FalsePositive.max_elements_spec
- oracle: compute_max_elements
- tests: tests/differential.py check_capacity_functions
- spec: 02-sketch
- status: tested

### C11: exactly C(2^b − 1, i) sketches decode to i elements
- statement: the map {M ⊆ K\{0} : |M| ≤ c} → K^c, M ↦ (s_1, s_3, …, s_{2c−1}) is injective, and the decoder inverts it on its image and fails elsewhere.
- source: src/test.cpp:118
- kind: textbook
- lean: Minisketch.Decode.decode_spec
- oracle: Minisketch.decode
- tests: TestExhaustive
- spec: 06-decode
- status: tested

## Berlekamp–Massey

### C12: even discrepancies vanish
- statement: for any sequence with s_{2k} = s_k^2 (not only true sketches), the textbook BM discrepancy at every even position is 0, so BM on the odd steps only gives the same result as full BM. Proof: BM's iterates satisfy the binary key equation Λ·S ≡ Λ_odd (mod x^{n+1}) (Λ_odd = xΛ′, its odd part), carried through every update because length changes only happen at odd positions, so the shift m is even; at an even position 2k it makes the discrepancy [Λ_odd·Λ_even]_{2k}, which is 0 (odd times even series has only odd exponents).
- source: doc/optimizations.md:453, doc/decode_explained.md:38, src/sketch_impl.h:2594
- kind: textbook (Berlekamp 1968), relied on for arbitrary inputs
- lean: Minisketch.BMEven.even_discrepancy_zero
- oracle: bm_odd
- tests: TestBerlekampMasseyExhaustive (GF(4..16)), TestBerlekampMasseyVsRef; pyminisketch TestOracle.test_bm_odd_exhaustive (every vector and max_degree: GF(4) cap 5 … GF(256) cap 2)
- spec: 04-bm
- status: proved

### C13: tail identity
- statement: for a sequence with s_{2k} = s_k^2, S = S_odd + S^2 as power series (the even part is Σ s_k^2 x^{2k} = S(x)^2). Hence for C = Λ·S: Λ·C = C^2 + Λ^2·S_odd.
- source: doc/optimizations.md:474, src/sketch_impl.h:2606
- kind: novel
- lean: Minisketch.BMTail.tail_identity
- oracle: verify_tail_odd
- tests: TestBerlekampMasseyVsRef
- spec: 04-bm
- status: proved

### C14: even-position recurrence
- statement: C_{2j} = C_j^2 + Σ_{i≥1} Λ_i C_{2j−i}; hence if C_k = 0 for L < k < 2j then C_{2j} = C_j^2. (The docs state the second form as an identity; it needs that hypothesis.)
- source: doc/optimizations.md:467, src/sketch_impl.h:2611
- kind: novel
- lean: Minisketch.BMTail.even_coeff
- oracle: verify_tail_odd
- tests: TestBerlekampMasseyVsRef
- spec: 04-bm
- status: proved

### C15: odd-only tail check is exact
- statement: let Λ have Λ_0 = 1 and degree ≤ L, and let the discrepancies C_k = [Λ·S]_k vanish for L < k ≤ P with P ≥ 2L. Then C_k = 0 for all L < k ≤ N iff u_j + Σ_{a≥1} Λ_a^2 u_{j−a} = 0 (u_j = s_{2j+1}) at every odd 2j+1 in (P, N]. So checking the stored odd syndromes against Λ^2 is exactly the full LFSR check.
- source: doc/optimizations.md:465, src/sketch_impl.h:2450, src/sketch_impl.h:2606
- kind: novel
- lean: Minisketch.BMTail.tail_check_iff
- oracle: bm_odd
- tests: TestBerlekampMasseyVsRef (all tail tiers)
- spec: 04-bm
- status: proved

### C48: early termination returns the textbook BM result
- statement: BerlekampMasseyOdd with early termination and the bulk tail check returns exactly what textbook BM over all 2c syndromes returns, including the max_degree failure. Follows from C12 and C15: once the check passes, every remaining discrepancy is 0, and BM makes no update on a zero discrepancy; a failing check resumes BM where it stopped.
- source: doc/optimizations.md:465, src/sketch_impl.h:2623
- kind: novel
- lean: Minisketch.BMEven.bm_odd_early_eq_bm, Minisketch.BMEven.bm_stable
- depends: C12, C15
- oracle: bm_odd
- tests: TestBerlekampMasseyExhaustive (GF(4..16)), TestBerlekampMasseyVsRef; pyminisketch TestOracle.test_bm_odd_exhaustive (every vector and max_degree: GF(4) cap 5 … GF(256) cap 2)
- spec: 04-bm
- status: proved

### C16: odd syndromes are a Λ^2 recurrence
- statement: u_j = s_{2j+1} = Σ_m m·(m^2)^j, so u satisfies the LFSR with connection polynomial Λ^{(2)} (each coefficient squared) whose roots are the m^2.
- source: src/sketch_impl.h:2450, doc/decode_coding_theory.md:38
- kind: textbook
- lean: Minisketch.Decode.odd_recurrence
- oracle: verify_tail_odd
- tests: TestBerlekampMasseyVsRef
- spec: 04-bm
- status: proved

### C17: transposed Karatsuba middle product
- statement: for c = c0 + x^h c1 and w = w0 + x^h w1 + x^{2h} w2, MP(c, w) = (α + MP(c0, w0 + w1), α + MP(c1, w1 + w2)) with α = MP(c0 + c1, w1), in characteristic 2.
- source: doc/optimizations.md:489, src/sketch_impl.h:2405
- kind: textbook (transposition principle)
- lean: Minisketch.PolyArith.mp_low, Minisketch.PolyArith.mp_high
- oracle: middle_product_kara
- tests: via VerifyTail and PowerSumSplit only
- spec: 03-polyarith
- status: proved

### C18: BM cost
- statement: BM goes from O(c·n) sequential work to O(n^2) sequential plus O(c·n), O(c·n^0.585) or O(c log n) parallel tail checking, by tier.
- source: doc/optimizations.md:18, doc/optimizations.md:501
- kind: heuristic (complexity analysis plus measurement)
- lean: -
- oracle: -
- tests: bench
- spec: 04-bm
- status: not-a-theorem

### C19: probabilistic tail verification
- statement: a random-evaluation check accepts a wrong tail with probability ≤ N/2^b per round (Schwartz–Zippel). Not implemented: it would change the exactness guarantee. The remark that exact verification is "Ω(c log n)-ish" is not a known lower bound.
- source: doc/optimizations.md:690
- kind: textbook (bound); the lower-bound remark is unjustified
- lean: -
- oracle: -
- tests: -
- spec: -
- status: not-a-theorem

## Root finding

### C20: trace splitting
- statement: Tr is GF(2)-linear with values in {0,1}; for f square-free and split, T = Tr(βx) mod f gives gcd(f, T) = Π_{α: Tr(βα)=0}(x − α).
- source: doc/math_decode.md:19, doc/decode_explained.md:81
- kind: textbook
- lean: Minisketch.Trace.trace_split, Minisketch.Trace.tr_mem
- oracle: poly_find_roots
- tests: TestFindRoots
- spec: 05-roots
- status: proved

### C21: fail-fast trace test
- statement: T + T^2 ≡ β(x + x^{2^b}) mod g for T = Tr(βx) mod g; if g is square-free and split this is 0.
- source: src/sketch_impl.h:2193
- kind: textbook
- lean: Minisketch.Trace.trPoly_sq_add
- oracle: -
- tests: TestFindRoots
- spec: 05-roots
- status: proved

### C22: full-split test
- statement: x^{2^b} ≡ x mod f iff f is square-free and splits into distinct linear factors over K. (Decodability additionally needs BM success, deg f ≤ capacity and f(0) ≠ 0, see C37.)
- source: doc/optimizations.md:535, src/sketch_impl.h:1769, doc/math_decode.md:25
- kind: textbook
- lean: Minisketch.Trace.dvd_frobenius_iff, Minisketch.Trace.frobenius_mod_iff
- oracle: poly_frobeniusmod
- tests: TestFindRoots (irreducible quadratic, double root)
- spec: 05-roots
- status: proved

### C23: the trace-0 side never holds an irreducible factor
- statement: in the non-fully-factorizable path, at a root α ∉ K of g, t + t^2 = β(α + α^{2^b}) ≠ 0, so t(α) ∉ {0,1}; hence gcd(g, t) has only roots in K and recursing on it cannot fail.
- source: src/sketch_impl.h:2338
- kind: textbook
- lean: Minisketch.Trace.gcd_side_fixed
- oracle: -
- tests: TestFindRoots
- spec: 05-roots
- status: proved

### C24: depth bound
- statement: a node reached after k splits with linearly independent β_1..β_k has its roots in an affine subspace of dimension b − k, so at most 2^{b−k} roots.
- source: src/sketch_impl.h:2239
- kind: textbook
- lean: Minisketch.Trace.depth_bound
- oracle: -
- tests: TestFindRootsDeepTree
- spec: 05-roots
- status: open

### C25: independent β's separate all roots
- statement: in GF(2^b), if β_0..β_{b−1} are GF(2)-linearly independent (distinct subset sums), then for α ≠ α′ some Tr(β_i α) ≠ Tr(β_i α′). Proof: the subset sums are then all of K, so otherwise Tr(γ(α − α′)) = 0 for every γ, but Tr is a nonzero polynomial of degree 2^{b−1} < |K|.
- source: src/sketch_impl.h:1758
- kind: textbook
- lean: Minisketch.Trace.separate
- oracle: -
- tests: TestFindRootsOrbit (rank check)
- spec: 05-roots
- status: proved

### C26: shared Frobenius orbit
- statement: for g | f, (h mod f) mod g = h mod g; and Tr(βx) mod f = Σ_i β^{2^i}·X_i with X_i = x^{2^i} mod f. So one orbit at the root gives every node's trace.
- source: doc/optimizations.md:512, doc/math_decode.md:21
- kind: textbook
- lean: Minisketch.Trace.trPoly_orbit, Minisketch.Trace.modByMonic_of_dvd
- oracle: frobenius_orbit
- tests: TestFindRootsOrbit
- spec: 05-roots
- status: proved

### C27: tower split parameters
- statement: for 8 | b and z generating GF(2^8), {γ·x^l·z^j : l < b/8, j < 8} is a GF(2)-basis of K (x a generator of K over GF(2^8)); since z^{2^8} = z, Tr(γ z^k x) = Σ_{r<8} (z^k)^{2^r}·Z_r with 8 bucket sums Z_r. Same with GF(16) and 4 | b. The bucket identity is C28 (proved); the basis property is checked by rank for b ∈ {8, 16, 24, 32, 64}.
- source: doc/optimizations.md:524, src/sketch_impl.h:1758
- kind: novel (application)
- lean: Minisketch.Trace.tower_basis
- oracle: tower_betas
- tests: TestFindRootsOrbit (rank check, RefTrace); pyminisketch TestOracle.test_orbit_and_betas
- spec: 05-roots
- status: open

### C28: relative-trace transitivity
- statement: Tr_{K/GF(2)} = Tr_{GF(16)/GF(2)} ∘ Tr_{K/GF(16)}; so Tr(β z^k x) = Tr_{GF(16)/GF(2)}(z^k·T) with T = Tr_{K/GF(16)}(βx), and one carried T serves 4 depths.
- source: doc/optimizations.md:789, doc/math_decode.md:52
- kind: textbook
- lean: Minisketch.Trace.tr_tower
- oracle: -
- tests: TestFindRootsMontgomery (RefTrace)
- spec: 05-roots
- status: proved

### C29: coset vanishing polynomial
- statement: for a Cantor chain (β_0 = 1, β_{j+1}^2 + β_{j+1} = β_j) the subspace polynomials S^0 = x, S^{m+1} = (S^m)^2 + S^m are additive, monic of degree 2^m, with (S^m)' = 1, S^m(β_{j+m}) = β_j and S^m(β_j) = 0 for j < m. The chain is GF(2)-independent, so the coset C = η + span(β_0..β_{m−1}) has 2^m points, Z_C = S^m + S^m(η) and Z_C' = 1. Halving: Z_{C0} = 1 on C1 = C0 + β_{m−1}.
- source: src/sketch_impl.h:1319, doc/optimizations.md:740, src/additive_fft.h:28
- kind: textbook (Cantor)
- lean: Minisketch.Subspace.cantor_vanishing_derivative
- oracle: cantor_basis
- tests: TestAdditiveFFT
- spec: 03-polyarith
- status: proved

### C30: Montgomery squaring by the derivative
- statement: let N = |C|, R = Z_C, deg a < 2N, gcd(f, R) = 1, q ≡ a/f on C (deg q < N) and h = a + q f. Then R | h, r = h/R satisfies r ≡ a·R^{-1} mod f, deg r < max(N, deg f), and r(w) = a'(w) + q'(w) f(w) + q(w) f'(w) for w ∈ C. If a is a square, a' = 0.
- source: doc/optimizations.md:743, src/sketch_impl.h:1325, doc/math_decode.md:39
- kind: novel
- lean: Minisketch.Subspace.reduce_coset
- depends: C29
- oracle: mont_square
- tests: TestFindRootsMontgomery
- spec: 05-roots
- status: proved

### C31: coset halving
- statement: C = C0 ⊔ C1 (cosets of W_{m−1}), Z_C = Z_{C0}·Z_{C1}, and Z_{C0} ≡ 1 on C1. So a value vector of v·Z_C mod g on C converts to v·Z_{C1} mod g on C1 with one pointwise division and one size-N/2 interpolation.
- source: src/sketch_impl.h:1339, doc/optimizations.md:780
- kind: novel
- lean: Minisketch.Subspace.carry_coset
- depends: C29
- oracle: mont_carry_to_factor
- tests: TestFindRootsMontgomery
- spec: 05-roots
- status: proved

### C32: unit multiples do not change the gcd
- statement: if gcd(Z_C, g) = 1 then gcd(g, T·Z_C mod g) = gcd(g, T).
- source: src/sketch_impl.h:1596, doc/optimizations.md:810
- kind: novel (application)
- lean: Minisketch.Montgomery.dvd_mod_mul_unit_iff
- oracle: -
- tests: TestFindRootsMontgomery
- spec: 05-roots
- status: proved

### C33: CosetDivide is exact
- statement: if g, h are monic, deg g ≤ N, and the polynomial q interpolated from g(w)/h(w) on C (|C| = N, h nonzero on C) is monic of degree deg g − deg h, then q·h = g.
- source: src/sketch_impl.h:1968, doc/optimizations.md:815
- kind: novel
- lean: Minisketch.Montgomery.coset_divide_exact
- oracle: coset_divide
- tests: TestFindRootsMontgomery
- spec: 05-roots
- status: proved

### C34: novel-basis derivative
- statement: with X_k = Π_{j∈bits(k)} S^j and (S^j)' = 1, X_k' = Σ_{j∈bits(k)} X_{k − 2^j}.
- source: doc/optimizations.md:767, src/additive_fft.h:312
- kind: novel (application of Lin–Chung–Han)
- lean: Minisketch.Subspace.novel_derivative
- oracle: novel_derivative
- tests: TestAdditiveFFT
- spec: 03-polyarith
- status: proved

### C35: power-sum split
- statement: for t with t(α) ∈ {0,1} at every root α of g (root set R, 0 ∉ R) and A = {α : t(α) = 1}: (a) P_k = Σ_α t(α)α^k = Σ_j t_j Q_{j+k}, so P_{2k} = P_k^2; (b) if |R| ≤ 2M and BM on P_1..P_{2M} returns a polynomial whose reversal f_B = Π_{β∈B}(x − β) divides g (so B ⊆ R; this is what CosetDivide certifies), then B = A, so the children inherit correct power sums. The proof needs BM's minimality (C5: |B| ≤ |A|, since Λ_A also generates P): divisibility alone does not exclude a superset of A. Then Σ_{α∈A∖B} α^m f_B(α) = 0 for m = 1..2M−|B| ≥ |A∖B| forces A ⊆ B (weighted Vandermonde). (c) A BM result of length 1 means A = ∅. (d) The other side has power sums Q_k − P_k.
- source: src/sketch_impl.h:1994, doc/optimizations.md:822, doc/math_decode.md:59
- kind: novel
- lean: Minisketch.PowerSumSplit.side_unique, Minisketch.BMPowerSum.bm_len_le_card
- depends: C5, C48
- oracle: power_sum_split
- tests: TestFindRootsMontgomery (intact and corrupted sums)
- spec: 05-roots
- status: proved

### C36: Montgomery full-split test
- statement: with r the Montgomery reduction of X̃_b − x̃ (so r ≡ X^{2^b} − X mod f) and deg r < deg f ≤ N = |C|: r vanishes on C iff f | X^{2^b} − X, i.e. iff x^{2^b} ≡ x mod f (C22).
- source: src/sketch_impl.h:1691
- kind: novel
- lean: Minisketch.Montgomery.split_test
- oracle: montgomery_orbit
- tests: TestFindRootsMontgomery; pyminisketch TestOracle.test_montgomery_orbit
- spec: 05-roots
- status: proved

### C37: a successful decode re-encodes, and 0 is never a root
- statement: if BM's connection polynomial on a Frobenius-structured syndrome sequence s_1..s_N factors as Π_{α∈R}(1 − αx) over nonzero α, then s_k = Σ_{α∈R} α^k for k ≤ N (the decoded set re-encodes to the sketch) and |R| = L, so the reversed polynomial f = x^L Λ(1/x) has no root 0. Proof: BM's key equation Λ·S ≡ Λ_odd (from the C12 invariant) and C2 for R give Λ·(S − S_R) ≡ 0, and Λ is invertible; then Λ generates s with length |R|, so L ≤ |R| by C5, and deg Λ = |R| ≤ L.
- source: implicit; src/test.cpp:98 and fuzz decode test it
- kind: textbook (binary BCH), not written down anywhere in the docs
- lean: Minisketch.Decode.decode_reencodes
- depends: C2, C5, C12
- oracle: Minisketch.decode
- tests: TestExhaustive, TestRandomized, fuzz decode; pyminisketch TestOracle.test_decode_fast_exhaustive (re-encoding and no zero element, every sketch of the tiny fields)
- spec: 06-decode
- status: proved

### C38: the Cantor basis exists for b = 2^ℓ
- statement: for b a power of two there are β_0 = 1 and β_j with β_j^2 + β_j = β_{j−1} for 1 ≤ j < b, linearly independent. Reason: y^2 + y = β_{j−1} is solvable iff Tr(β_{j−1}) = 0, and with S^{b−1} = Tr (b = 2^ℓ) one gets Tr(β_{j−1}) = S^{b−j}(1) = 0 for j < b.
- source: src/additive_fft.h:21, doc/optimizations.md:158
- kind: textbook (Cantor 1989)
- lean: Minisketch.AdditiveFFT.cantor_basis_exists
- oracle: cantor_basis
- tests: TestAdditiveFFT
- spec: 03-polyarith
- status: open

### C39: polynomial arithmetic primitives
- statement: (a) S^i(x) = Σ_{j: C(i,j) odd} x^{2^j}; (b) the LCH additive FFT evaluates at η_0..η_{2^m−1} in index order; (c) char-2 Karatsuba; (d) Newton inverse series g' = f g^2 mod x^{2m}; (e) the inversion-free pseudo-remainder equals β^2·(a mod b); (f) carry-less reduction is GF(2)-linear, so Σ reduce(a_i b_i) = reduce(Σ a_i b_i); (g) the half-gcd matrix depends only on the top 2·d − 1 coefficients.
- source: src/additive_fft.h, src/sketch_impl.h:447, src/sketch_impl.h:636, src/sketch_impl.h:266, src/sketch_impl.h:923, doc/optimizations.md:356
- kind: textbook (a–d, f, g), novel application (e)
- lean: Minisketch.PolyArith.primitives
- oracle: additive_fft
- tests: TestAdditiveFFT, TestTraceModReducers, TestSubquadraticGCD
- spec: 03-polyarith
- status: tested

### C47: CLMUL reduction recipes
- statement: ReduceWithClMul (two folds, needs 2·deg MOD − 2 < b), ReduceTrinomial (needs b + POS ≤ 65 on the lane-shift path and 2·POS − 2 < b), and the ARM shift-fold for b = 64, MOD < 256 compute the remainder of a 2b − 1 bit carry-less product mod x^b + MOD.
- source: src/fields/clmul_common_impl.h:123, src/fields/clmul_common_impl.h:179
- kind: novel (implementation)
- lean: Minisketch.Field.reduce_correct
- native_decide: allowed
- oracle: GF2Ops
- tests: TestFieldOps (exhaustive for b ≤ 8)
- spec: 01-field
- status: tested

## Heuristic and analytic statements

### C40: false-positive model
- statement: a sketch of an overfull set behaves like a uniformly random vector, so P(false decode) ≈ #decodable/2^{bc} and each extra capacity unit buys b bits.
- source: doc/math.md, src/false_positives.h
- kind: heuristic
- lean: -
- oracle: base_fpbits
- tests: -
- spec: 02-sketch
- status: not-a-theorem

### C41: split retry probability
- statement: a split at a node of degree d fails with probability about 2^{1−d} for a random basis. A fixed basis with adversarial roots defeats this (see TestFindRootsDeepTree).
- source: doc/optimizations.md:51, doc/optimizations.md:966
- kind: heuristic
- lean: -
- oracle: -
- tests: TestFindRootsDeepTree
- spec: 05-roots
- status: not-a-theorem

### C42: no closed-form root formula for degree ≥ 5
- statement: as written, this misapplies Abel–Ruffini: over a finite field every polynomial is solvable by radicals (Galois groups are cyclic). The intended point is that no cheap closed form exists.
- source: doc/optimizations.md:975
- kind: heuristic (wording to fix)
- lean: -
- oracle: -
- tests: -
- spec: -
- status: not-a-theorem

### C43: no cheaper route around the Frobenius chain
- statement: "circular by the transposition principle"; the alternatives to the squaring chain all cost more.
- source: doc/optimizations.md:548
- kind: heuristic
- lean: -
- oracle: -
- tests: -
- spec: -
- status: not-a-theorem

### C44: Toom-3 field size
- statement: Toom-3 needs 5 evaluation points including ∞, which GF(4) already has, so b ≥ 2 suffices (the doc says b ≥ 3).
- source: doc/optimizations.md:972
- kind: heuristic (wording to fix)
- lean: -
- oracle: -
- tests: -
- spec: -
- status: not-a-theorem

### C45: cost and crossover numbers
- statement: all transform counts, cost ratios and MINISKETCH_* cutoffs are measurements.
- source: doc/optimizations.md
- kind: heuristic
- lean: -
- oracle: -
- tests: bench, tools/bench_ab.sh
- spec: 07-tuning
- status: not-a-theorem

### C46: decode is exactly the textbook decoder
- statement: for every input (including adversarial and overfull sketches), decode returns the same set, or the same failure, as textbook BM + Berlekamp trace root finding. Output order is not specified.
- source: doc/optimizations.md:445, doc/decode_explained.md:255, doc/decode_coding_theory.md:120
- kind: novel (conjunction of C12–C15, C20–C37)
- lean: Minisketch.Decode.decode_eq_textbook
- oracle: decode_fast
- tests: pyminisketch TestOracle.test_decode_fast* (exhaustive for GF(4..32) at capacity ≤ 3 and GF(256) at 2, random to capacity 30, a forced Montgomery/power-sum path at 300 elements), tests/differential.py (capacity ≤ 16), TestExhaustive, fuzz decode
- spec: 06-decode
- status: tested
