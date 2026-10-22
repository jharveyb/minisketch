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
| `Minisketch/Montgomery.lean` | coset Montgomery step, carry to a child, Montgomery-form gcd, CosetDivide (C30–C33) |

Theorems are stated over an arbitrary field `K` of characteristic 2
(`[Field K] [CharP K 2]`), adding `[Fintype K]` with `card K = 2^b` where the
field size matters. None of this proves that the C++ code implements these
algorithms; that link is made by the differential tests against the Python
oracle (`tests/differential.py`).
