# Machine-checked proofs

Lean 4 + Mathlib proofs of the mathematical claims listed in
[`doc/spec/claims.md`](../doc/spec/claims.md). The registry says which claims
are proved, which only have a statement, and which are covered by tests only.

## Setup

`elan` is in the flake devshell (`nix develop`); it installs the Lean version
pinned in `lean-toolchain`. Mathlib is pinned in `lakefile.lean` and
`lake-manifest.json`; fetch its compiled cache instead of building it:

```
cd verification
lake exe cache get
lake build
```

## Checking the registry

```
python3 verification/check_claims.py           # needs a built project
python3 verification/check_claims.py --no-lean # registry and oracle names only
```

It fails if a claim marked `proved` depends on `sorryAx` or on any axiom other
than `propext`, `Classical.choice` and `Quot.sound` (plus `Lean.ofReduceBool`
for claims marked `native_decide: allowed`), if a listed Lean theorem does not
exist, or if a listed oracle function is missing from `tests/pyminisketch.py`.

## Layout

One module per area, matching the sections of the registry:

| Module | Claims |
|---|---|
| `Minisketch/Sketch.lean` | power sums: linearity, Frobenius halving, uniqueness (C1, C3, C4) |
| `Minisketch/BM.lean` | executable textbook BM model; Massey's lemma; BM returns a shortest LFSR (C5) |
| `Minisketch/BMPowerSum.lean` | a set's locator generates its power sums, so BM on them has length ≤ the set size (C35's minimality hypothesis) |
| `Minisketch/BMEven.lean` | even BM discrepancies vanish on Frobenius sequences (C12); early termination with the odd-only tail check returns textbook BM's result (C48) |
| `Minisketch/Decode.lean` | Λ_R S_R = −xΛ_R′ (C2); odd-syndrome Λ² recurrence (C16); a decode BM accepts re-encodes to the input and has no zero root (C37) |
| `Minisketch/BMTail.lean` | BM tail identity, even/odd discrepancy recurrences, exactness of the odd-only tail check (C13–C15) |
| `Minisketch/Trace.lean` | trace splitting, fail-fast test, gcd side stays in GF(2^b), shared orbit, full-split test, relative trace, separation by independent split parameters (C20–C23, C25, C26, C28) |
| `Minisketch/Depth.lean` | depth bound: a node after k trace splits has at most 2^(b−k) roots (C24) |
| `Minisketch/Reduce.lean` | two folds reduce a carry-less product exactly when 2·deg MOD − 2 < B (C47) |
| `Minisketch/Cantor.lean` | y² + y = c solvable iff Tr(c) = 0; S^(2^ℓ−1) = Tr; the Cantor chain exists for b = 2^ℓ (C38) |
| `Minisketch/Tower.lean` | tower split parameters basis·x^l·z^j are independent; SubfieldGenerator's z has period s; they separate all roots (C27) |
| `Minisketch/PolyArith.lean` | transposed Karatsuba middle product (C17) |
| `Minisketch/PowerSumSplit.lean` | power-sum split formula, uniqueness of the BM-recovered side, empty side, complement (C35) |
| `Minisketch/Montgomery.lean` | coset Montgomery step, carry to a child, Montgomery-form gcd, CosetDivide, coset full-split test (C30–C33, C36) |
| `Minisketch/Subspace.lean` | subspace polynomials on a Cantor chain, coset vanishing polynomials with derivative 1, coset halving, novel-basis derivative; C30/C31 on real cosets (C29, C34) |

Theorems are stated over an arbitrary field `K` of characteristic 2
(`[Field K] [CharP K 2]`), adding `[Fintype K]` with `card K = 2^b` where the
field size matters. None of this proves that the C++ code implements these
algorithms; that link is made by the differential tests against the Python
oracle (`tests/differential.py`).
