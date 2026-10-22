import Lake
open Lake DSL

-- Machine-checked proofs of the claims in doc/spec/claims.md.
package minisketch where
  leanOptions := #[⟨`autoImplicit, false⟩]

require mathlib from git
  "https://github.com/leanprover-community/mathlib4" @ "v4.35.0"

@[default_target]
lean_lib Minisketch
