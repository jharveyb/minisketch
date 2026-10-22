# 6. Decode

This chapter defines what decode returns (normative) and how libminisketch
computes it (informative). Source: `SketchImpl::Decode` in
`src/sketch_impl.h` and `minisketch_decode` in `src/minisketch.cpp`.

## Conformance: the textbook decoder

For a sketch with odd syndromes u_0..u_(c−1) and a limit
max_elements = n_max, the result is defined by this procedure
(`Minisketch.decode` in `tests/pyminisketch.py`):

```
textbook_decode(u, n_max):
    S_(2j+1) = u_j;  S_(2k) = S_k^2            for k = 1..c
    Λ = BerlekampMassey(S_1 .. S_(2c))        # textbook, chapter 4
    L = |Λ| - 1                               # LFSR length
    if L == 0: return {}                      # empty set
    if L > n_max: return FAIL
    if Λ_L == 0: return FAIL                  # see below
    f = rev(Λ) = x^L Λ(1/x)                   # monic, roots are the elements (C6)
    if f is not square-free with all roots in K: return FAIL
    return the L roots of f
```

A conforming decoder returns the same set or FAIL for every input,
including vectors that are not sketches of any set of at most c elements.
The order of the returned elements is unspecified.

Remarks on the definition:

- **Λ_L = 0.** The Python textbook decoder does not test this; it goes on
  to find the roots of f, which has the root 0 (f(0) = Λ_L). C37 says that
  whenever BM's Λ (of any length L) factors as Π(1 − αx) over nonzero α,
  then |R| = L, so f = rev(Λ) has no root 0, and the roots re-encode to
  the input. So a Λ with Λ_L = 0 never leads to a successful decode, the
  Python decoder fails too, and the explicit test changes nothing; it
  guarantees that 0 is never returned. C37 is proved in Lean for every L,
  including L > c (possible when max_elements > c). A port must still
  return FAIL when Λ_L = 0.
- **Uniqueness.** When the input is the sketch of a set of at most c
  elements, the result is that set (C3, C5, C6), provided n_max ≥ its
  size. Exactly C(2^b − 1, i) of the 2^(bc) input vectors decode to i
  elements (C11).
- **`max_elements` = 0** decodes only the empty set (L = 0).
- **Capacity 0**: BM on no syndromes returns Λ = 1; decode returns the
  empty set.

## The libminisketch pipeline

```
Decode(max_count):                                        # src/sketch_impl.h:2912 (Decode)
    Λ = BerlekampMasseyOdd(u, n_syn = 2c, max_degree = max_count)    # chapter 4
    if FAIL: return -1                                    # length would exceed max_count,
                                                          # or Λ_L == 0
    if |Λ| == 1: return 0                                 # empty set
    if |Λ| > 1 + max_count: return -1                     # unreachable (BM enforces it)
    f = reverse(Λ)
    roots = FindRoots(f, basis, field, odd_power_sums = u)            # chapter 5
    if roots empty: return -1
    write ToUint64(root) for each root, in FindRoots order
    return |roots|
```

`max_count` is min(`max_elements`, 2c, INT_MAX − 1): BM's connection
polynomial never has more than 2c + 1 terms, so any larger limit means no
limit. (Before this was fixed, `max_elements` was cast to a C `int`, and at
INT_MAX the comparison overflowed, while values ≥ 2^31 wrapped negative;
both made every non-empty decode fail. A port must not reproduce that.)
The API documents `output` as an array of `max_elements` entries; the
decoder writes at most max_count values because L ≤ max_count.

## Call graph

```
Decode
├─ BerlekampMasseyOdd(u, 2c, max_count)                       ch. 4
│   ├─ odd-step BM iteration (even steps skipped, C12)
│   └─ BerlekampMasseyVerifyTail (early termination, C15/C48)
│        tiers: lazy direct | lazy Karatsuba (LazyMiddleProduct) | table direct | FFT | Karatsuba blocks
└─ FindRoots(rev Λ, basis, u)                                 ch. 5
    ├─ deg < 8: RecFindRoots(ff = false), plain betas, no orbit
    └─ deg >= 8: ComputeTraceOrbit
    │   ├─ deg >= 224, b ∈ {16,32,64}, CLMUL: MontgomeryOrbit
    │   │     ├─ OK → GF(16) betas, root trace, Montgomery carry; psum if deg >= 256
    │   │     ├─ NOT_SPLIT → FAIL
    │   │     └─ NOT_APPLICABLE (no root-free coset, or m + 2 > b) → standard orbit
    │   └─ standard orbit: TraceMod squarings, GF(2^8)-tower or plain betas,
    │         traces for the first min(b, 2 + floor(log2 b)) depths; X_b != x → FAIL
    └─ RecFindRoots(ff = true, orbit, betas, carry, psum)
        ├─ deg 1: root; deg 2: Qrt (FAIL if a_1 = 0 or no solution)
        ├─ trace: carried T_l (MontTowerTrace | TowerTraceSmall)
        │         | orbit trace mod g (FastDivMod) | SmallTrace | TraceMod.trace
        ├─ fail-fast test (small path only, C21)
        ├─ split: PowerSumSplit (deg >= 256, Montgomery, sums known)
        │         ├─ SPLIT → children with power sums
        │         ├─ TRIVIAL → next depth
        │         └─ GCD_FALLBACK → GCD
        │         | GCD(t, g) (lazy, inversion-free Euclid)
        ├─ quotient: from PowerSumSplit | CosetDivide (Montgomery, deg >= 64) | FastDivMod
        ├─ carry to factors: MontCarryToFactor (deg > 16) | normal form (3..16)
        ├─ RecFindRoots(t = 1 factor, ff)        may FAIL
        └─ RecFindRoots(t = 0 factor, ff = true) never fails (C23)
```

## Failure conditions

Decode returns −1 in exactly these cases, each of which is a FAIL of the
textbook decoder:

1. **BM length exceeds max_count.** `BerlekampMasseyOdd` fails before a
   length change to L > max_count (equivalent to the textbook test
   L > n_max, since L never decreases; C48 makes the early-terminated
   result equal to textbook BM's).
2. **Λ_L = 0.** The final check in `BerlekampMasseyOdd`. By C37 rev(Λ)
   does not split into distinct linear factors then, so the root search
   would fail anyway (see the remark above).
3. **f not square-free or not split** (C22). For deg f ≥ 8 this is decided
   once at the root by X_b = x (standard orbit) or its Montgomery form
   (C36). For 3 ≤ deg f ≤ 7 it is detected during splitting: the fail-fast
   test t + t^2 ≡ 0 mod g (C21), a quadratic leaf with a_1 = 0 (double
   root) or without a root in K, or the depth bound (C24). For deg f = 2
   the quadratic case alone; deg f = 1 always succeeds.
4. **Root count mismatch** (defensive; unreachable when the above are
   implemented correctly).

A zero root is impossible (C37): f(0) = Λ_L ≠ 0 after check 2.

Every other branch in chapter 5 is a fallback with the same result, not a
failure: Montgomery not applicable → standard orbit; a power-sum split that
cannot be certified → gcd; a trivial split → next depth.

## Output order

The order is implementation-defined and should not be relied on. In
libminisketch it is the depth-first order of chapter 5 ("Output order"),
which depends on:

- the basis, which is random per sketch unless `minisketch_set_seed` was
  called, so two decodes of equal sketches can differ in order unless both
  sketches have the same seed;
- the implementation number: implementations 0 and 1/2 derive different
  bases from the same seed (different SipHash keys), and only
  implementation 1 at b ∈ {16, 32, 64} has the Montgomery path;
- the root degree: below 8 the plain betas, from 8 the GF(2^8) tower when
  8 | b, from 224 (Montgomery) the GF(16) tower;
- the field's internal representation (implementation 1 at b = 58), and
  the Cantor basis and coset chosen on the Montgomery path.

For a fixed library build, sketch contents, seed and implementation, the
order is deterministic. Tests compare decoded sets sorted
(`tests/differential.py` sorts both sides).

## Exactness argument (summary)

C46 states that decode equals the textbook decoder. It decomposes into:

- BM: even steps skipped (C12), odd-only tail check exact (C13–C15),
  early termination returns textbook BM's result (C48).
- Root finding: trace splitting (C20), shared orbit (C26), full-split test
  (C22, C36), independent tower betas (C25, C27, C28), Montgomery reduction
  and carrying (C29–C31), unit multiples in gcds (C32), certified coset
  division (C33), power-sum splits (C35), the never-failing second
  recursion (C23), the depth bound (C24).
- No zero root (C37).

See [claims.md](claims.md) for the status of each (most are proved in Lean
under `verification/`; C11 and C46 are tested, exhaustively on small
fields).
