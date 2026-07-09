# Testing minisketch

This document describes the test layers, how to run and reproduce them, and how
to measure test quality. The overall goal is enough confidence in the tests
that the implementation of individual decode stages (Berlekamp-Massey, root
finding, TraceMod reduction, ...) can be swapped out and validated quickly.

## Test taxonomy

| Layer | Binary / entry point | What it checks | Oracle |
| --- | --- | --- | --- |
| Exhaustive | `test-noverify` / `test-verify` (`src/test.cpp`, `TestExhaustive`) | every sketch for small bits×capacity | cross-implementation agreement, re-encode identity, decodable-count = C(2^bits−1, i) |
| Randomized end-to-end | same binaries (`TestRandomized`) | random sets through the public API | cross-implementation agreement, serialize roundtrip, exact recovery within capacity |
| Parameter functions | same binaries (`TestComputeFunctions`) | `minisketch_compute_capacity` / `max_elements` | monotonicity, bounds, roundtripping |
| Unit / property | `unit-tests` (`src/unit_tests*.cpp`, Boost.Test) | internal algorithm stages directly | independent bit-level GF(2^b) reference (`src/test_refimpl.h`), naive polynomial references, algebraic identities |
| Fuzzing | `fuzz` (`src/fuzz/`, libFuzzer) | attacker-controlled serializations, data-driven sets, internal poly ops | same invariants, coverage-guided, under ASan/UBSan + `MINISKETCH_VERIFY` |
| Differential (Python) | `tests/differential.py` | C++ library vs. pure-Python reimplementation | independent implementation (`tests/pyminisketch.py`) |
| Golden vectors | `tests/golden_vectors.json` | serialization stability | fixed known-good bytes, asserted by both C++ and Python suites |

The verify builds (`test-verify`, `unit-tests`, `fuzz`) compile with
`MINISKETCH_VERIFY`, which turns the internal `CHECK_SAFE` invariant
assertions into aborts.

## Running

```sh
# Everything that runs under ctest (library tests, unit tests, Python suites):
cmake -B build -DMINISKETCH_BUILD_TESTS=ON
cmake --build build -j && ctest --test-dir build

# Fuzzing (requires clang; instruments library with fuzzer,ASan,UBSan):
cmake -B build-fuzz -DMINISKETCH_BUILD_FUZZ=ON -DCMAKE_CXX_COMPILER=clang++
cmake --build build-fuzz -j
FUZZ=decode ./build-fuzz/bin/fuzz -max_len=512 build-fuzz/corpus/decode src/fuzz/corpus/decode
```

Fuzz targets: `decode`, `roundtrip`, `poly_ops` (run the binary without `FUZZ`
set to list them). Each also runs as a short smoke test under plain `ctest` in
a fuzz build.

## Reproducing failures

Every randomized test prints its seed at startup:

```text
Test seed: 0x00000000deadbeef (reproduce with --seed=... or MINISKETCH_TEST_SEED)
```

- `test-noverify` / `test-verify`: `./test-noverify [complexity] --seed=0x...`
- `unit-tests` (Boost.Test consumes its own arguments; pass ours after `--`):
  `./unit-tests -- --seed=0x...`, or `MINISKETCH_TEST_SEED=0x... ./unit-tests`.
  Run a single case with `./unit-tests -t generic_field_11`.
- Fuzz crashes: libFuzzer writes a `crash-*` file;
  `FUZZ=<target> ./fuzz crash-...` replays it, and
  `-minimize_crash=1 -exact_artifact_path=min crash-...` minimizes it.
- Python suites print and accept `MINISKETCH_TEST_SEED` as well.

All test randomness is a pure function of the seed (SipHash counter mode; no
`std::random_device`, no implementation-defined `std::uniform_int_distribution`),
so seeds reproduce across platforms.

## Swapping a decode-stage implementation

The unit layer is written so that testing a replacement stage is one
differential property. For example, to validate a new fixed-modulus reducer
for TraceMod (as in the `fast_tracemod_reducers` branch):

1. In `src/unit_tests_impl.h`, add a check that the new code path equals the
   existing naive reference (`ut::PolyReduceRef`-based TraceMod reference in
   `TestPolyOps`) — or, simplest, nothing at all if the new code replaces
   `TraceMod`/`PolyMod` internals: the existing `TraceMod` property already
   pins the result exactly, for every instantiated field, on seeded random
   moduli/params, at ~256 iterations per field per run.
2. Add the same check to `src/fuzz/poly_ops.cpp` so it also runs
   coverage-guided with sanitizers (this catches memory bugs in new vector
   arithmetic that fixed-size random tests miss).
3. Run `ctest`, then a fuzz session:
   `tools/fuzz_stats.sh 300` (5 minutes per target).
4. Compare `tools/coverage.sh` output before/after to confirm the new code is
   actually exercised (watch the per-file lines column for `sketch_impl.h`).

## Metrics

### Coverage

```sh
tools/coverage.sh            # writes build-coverage/coverage-report.txt + HTML
```

Uses clang source-based coverage over the ctest suite (library tests, unit
tests, Python differential tests via the shared library are *not* included —
they exercise a separately-built library). Baseline as of 2026-07 (complexity
4, all field sizes, clmul enabled):

| File | Lines | Branches |
| --- | --- | --- |
| `src/sketch_impl.h` | 94.3% | 83.1% |
| `src/minisketch.cpp` | 87.7% | 91.0% |
| `src/int_utils.h` | 92.0% | 81.8% |
| TOTAL | 73.9% | 92.1% |

(The generic/clmul `*_Nbytes.cpp` files show ~40-55% line coverage because
only representative field sizes are unit-tested and each file's
per-field-size constructor switch is data-heavy; the end-to-end tests do
construct every field size.)

### Fuzzing

```sh
tools/fuzz_stats.sh 60       # 60 seconds per target
```

Reports libFuzzer edge coverage (`cov:`), features (`ft:`), corpus size and
executions per target, keeping a persistent working corpus under
`build-fuzz/corpus/`. Reference numbers after ~5 minutes per target on one
core (2026-07): decode cov≈28500, roundtrip cov≈29200, poly_ops cov≈1000
(the poly_ops target only instantiates two generic fields, hence the smaller
denominator), no failures across ~500k executions.

A small seed corpus is committed under `src/fuzz/corpus/<target>/` so short CI
runs start from meaningful inputs; long local runs grow the working corpus.

### Test power (mutation spot-checks)

There is no automated mutation-testing setup, but planted-bug checks are the
quality bar for this suite: flipping an operator or an index in
`sketch_impl.h` (e.g. the discrepancy XOR in `BerlekampMassey`) should be
caught by the `unit-tests` binary in under a second. When adding new algorithm
code, plant a bug and confirm a test notices before trusting the green run.

Results so far (2026-07):

- master: flipping the discrepancy XOR to OR in `BerlekampMassey` — caught by
  23 of 24 unit-test cases.
- `fast_tracemod_reducers` cherry-picked onto this suite (see the
  `tracemod-check` branch, which also carries the test adaptation to the new
  `TraceMod<F>` class API and added properties for
  `TraceModPolyMulFull/Low`, `TraceModMulBySquareLow`, `TraceModInvSeries`,
  and the reciprocal reducer path): 5 planted mutations — dropped coefficient
  squaring in `SquareReduce`, a stale shift value in the square-table build,
  off-by-one quotient and remainder lengths in `ReciprocalReduce`, and a
  dropped Karatsuba middle-term correction — were each caught by `unit-tests`
  within seconds, at the failing property, while the full suite passes on the
  unmutated branch.
