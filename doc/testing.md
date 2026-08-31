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

Fuzz targets: `decode`, `roundtrip`, `poly_ops`, and `reconcile` (a two-party
set-reconciliation scenario adapted from Bitcoin Core's
`src/test/fuzz/minisketch.cpp`, with 32- or 64-bit elements). Run the binary
without `FUZZ` set to list them. Each also runs as a short smoke test under
plain `ctest` in a fuzz build.

## Capacity and field-size coverage

Decode cost grows quadratically with capacity (a full 64-bit decode measures
~60ms at capacity 256, ~0.8s at 1024, uninstrumented), so each layer trades
depth for iteration count deliberately:

| Layer | Capacity / degree | Field sizes |
| --- | --- | --- |
| `TestExhaustive` | bits×capacity ≤ ~40 | all |
| `TestRandomized` | 0–8 / 0–128 / 0–4096 tiers, plus a guaranteed 512–1024 tier | all; large tier: 32, 64 |
| unit tests, reference-checked properties | degree ≤ 12 (quadratic naive references) | 2–8 exhaustive, 11/16/27/32/64 |
| unit tests, `TestDecodeStagesLarge` (reference-free) | 256 per field; **1024** for the 64-bit field | same |
| fuzz `decode`/`roundtrip` | log-uniform 1–1024 | 2–64 |
| fuzz `reconcile` | log-uniform 1–1024 | 32, 64 |
| fuzz `poly_ops` | degree ≤ 24 (referenced) + log-uniform ≤ 256 root sets (reference-free) | 11, 32 |
| `differential.py` | ≤ 16 (pure-Python decode cost) | 2–16, 32, 64 |

The 64-bit field — the practically relevant element size — is covered at
every layer, including full-capacity-1024 decode in both the unit suite (both
generic and clmul implementations) and the guaranteed large tier of the
end-to-end tests.

Fuzz targets are fully deterministic: fuzz builds define
`MINISKETCH_FUZZ_DETERMINISTIC`, which replaces the `std::random_device`
basis draw in sketch construction with a fixed value (production builds are
unaffected), and the targets set an input-derived basis via `SetSeed()`.
Running the same input twice produces identical behavior and coverage, so
crashes always replay.

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
5. Check test power over the change: plant one or two design-level bugs by
   hand (see "Test power" below) and run `MULL_DIFF_REF=<base-ref>
   tools/mutation.sh` to sweep operator-level mutants over just the diff
   (see "Mutation testing").
6. If the change is performance-motivated, measure it with
   `tools/bench_ab.sh <base-ref> .` (see "Performance measurement") rather
   than ad-hoc timing.

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

### Corpus management

There are two corpora per target with different lifecycles:

- **Working corpus** (`build-fuzz/corpus/<target>/`): libFuzzer writes every
  coverage-increasing input here *as it finds it*, so it grows during any
  fuzzing session and persists across runs (including interrupted ones). It
  accumulates unminimized; delete it freely to start over.
- **Committed seed corpus** (`src/fuzz/corpus/<target>/`): a
  coverage-minimized seed set (a few hundred to a few thousand byte-sized
  files per target), so fresh checkouts start from meaningful inputs. It is
  *not* updated automatically. Update it with:

  ```sh
  tools/update_seed_corpus.sh   # additive; then review and commit
  ```

  This is additive and monotone: libFuzzer's `-merge=1` keeps all existing
  seeds and appends only working-corpus inputs that add coverage features
  the seeds lack, so committed coverage never decreases and running the
  script after a too-short session is harmless. The set grows only while a
  target still has undiscovered features.

  After changing a target's input layout (adding/removing `Consume*`
  calls), old seeds still run but no longer decode to the cases they were
  selected for; rebuild that target's set from scratch with
  `RESET=1 FUZZ_TARGETS=<target> tools/update_seed_corpus.sh` (guarded: a
  reset that would replay less coverage than the committed set is refused).

- **Slow corpus tier** (`src/fuzz/corpus-slow/<target>/`, decode and
  roundtrip only): decode cost is quadratic in the `capacity` parameter, so
  high-capacity inputs cost seconds each and dominate startup replay — the
  full decode corpus takes ~10 minutes to replay once, versus ~35s for the
  low-capacity majority. The committed corpus is therefore split at capacity
  128 (the same worst-case boundary `decode.cpp` uses to prune
  implementations): inputs with capacity ≤ 128 stay in `corpus/` (the fast
  tier, always loaded) and those above go to `corpus-slow/`. Measured on the
  2026-07 decode corpus, the slow tier is ~32% of the inputs and ~94% of the
  replay time but adds only ~0.2% edge coverage (it does add ~4% more
  *features* — the deep-recursion hit-count buckets that matter when
  optimizing root finding).

  The tiering applies to *both* corpora. `tools/fuzz_stats.sh` loads only the
  fast tier by default (so startup stays ~35s); before each run it also moves
  any high-capacity inputs out of the working corpus into
  `build-fuzz/corpus-slow/<target>/`, so the load stays fast even as fuzzing
  discovers deep-capacity inputs. Pass `SLOW=1` to load the slow tier as well.
  `update_seed_corpus.sh` keeps the committed tiers correct: a full run merges
  both working tiers over the union of both committed tiers and re-partitions
  by capacity (via `build-fuzz/bin/corpus-tier`), so coverage is never lost;
  `FAST=1` merges only the fast working corpus into the committed fast tier and
  skips the slow tier entirely, trading completeness for speed (~2 min vs
  ~10 min for decode). See the example campaigns below.

libFuzzer also writes informational `slow-unit-*` artifacts (inputs taking
more than ~10s; expected for adversarial capacity-1024 decodes) to
`build-fuzz/artifacts/<target>/`. They are not failures and are safe to
delete — but they make good worst-case reproducers, e.g. to benchmark a
decode-stage optimization against a known-bad input.

Reproducers worth keeping are committed under
`src/fuzz/slow-units/<target>/`. Nothing replays that directory
automatically, deliberately so: unlike seed-corpus files, each of these
costs seconds to *minutes*, so they must never be picked up by corpus
merges, coverage replays, or CI. In particular, do not move them into
`src/fuzz/corpus/`. Replay one explicitly with:

```sh
FUZZ=decode ./build-fuzz/bin/fuzz src/fuzz/slow-units/decode/<file>
```

Current inventory (decode replay wall-times in the sanitized fuzz build,
measured 2026-07-10; all three are adversarial capacity-~1024 decodes, the
expected worst-case family):

| file                   | size | decode replay | parameters                  |
| ---------------------- | ---- | ------------- | --------------------------- |
| `slow-unit-0865be…39b` | 13 B | ~20 s         | bits=54, capacity=1021      |
| `slow-unit-bf232e…943` | 11 B | ~10 s         | bits=56, capacity=1024      |
| `slow-unit-d2f8ac…b59` | 11 B | ~12 s         | bits=62, capacity=1024, exercises the `SetSeed(-1)` fixed-basis path (~10x slower root-finding than a hashed basis on this input) |

### Example campaigns

**Fast-only campaign** — repeated short interactive passes, then promote just
the low-capacity coverage. Each `fuzz_stats.sh` run loads only the fast tier,
so startup stays ~35s and the whole time budget is spent fuzzing:

```sh
FUZZ_TARGETS=decode tools/fuzz_stats.sh 120   # ~35s load + 120s fuzz
FUZZ_TARGETS=decode tools/fuzz_stats.sh 120   # repeat as often as you like;
FUZZ_TARGETS=decode tools/fuzz_stats.sh 120   # the working corpus persists

FAST=1 FUZZ_TARGETS=decode tools/update_seed_corpus.sh   # ~2 min, fast tier only
git add src/fuzz/corpus/decode
git commit -m "fuzz: grow decode fast seed tier"
```

`FAST=1` merges only the fast working corpus into the committed fast tier and
never replays the slow tier. Any high-capacity inputs the runs discovered
accumulate in `build-fuzz/corpus-slow/decode/` and are *not* promoted here (the
run reports how many were left behind); they wait for a full update. So the
only committed change is `src/fuzz/corpus/`.

**Mixed campaign** — fast passes for broad coverage plus occasional slow passes
that also exercise the high-capacity tier, then a full update that folds in
both:

```sh
FUZZ_TARGETS=decode tools/fuzz_stats.sh 120            # fast passes
SLOW=1 FUZZ_TARGETS=decode tools/fuzz_stats.sh 300     # slow pass: loads
                                                        # corpus-slow too, so
                                                        # startup is minutes

tools/update_seed_corpus.sh                            # full: ~10 min for decode
git add src/fuzz/corpus/decode src/fuzz/corpus-slow/decode
git commit -m "fuzz: grow decode seed corpus (both tiers)"
```

A full update (no `FAST`) merges both working tiers over both committed tiers
and re-partitions by capacity, so new worst-case inputs land in
`src/fuzz/corpus-slow/` — hence both directories are committed. Run this
occasionally to fold in the high-capacity finds that fast-only updates skip.

### Test power (mutation spot-checks)

Planted-bug checks are the design-level quality bar for this suite: flipping
an operator or an index in `sketch_impl.h` (e.g. the discrepancy XOR in
`BerlekampMassey`) should be caught by the `unit-tests` binary in under a
second. When adding new algorithm code, plant a bug and confirm a test
notices before trusting the green run. Hand-planted mutations remain the
right tool for *targeted* checks — "a bug of this specific shape must be
caught at this specific layer" (e.g. a bug that roundtrip identities cannot
see must fall to a direct-evaluation property); the automated sweep below
does not replace that judgment.

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

### Mutation testing (Mull)

`tools/mutation.sh` sweeps operator-level mutants with
[Mull](https://mull-project.com/) (installed as the `mull-21` package; the
version suffix must match the clang that compiles the mutated build). It
compiles `unit-tests` with Mull's IR plugin — which embeds every mutant into
the binary at compile time, scoped by `mull.yml` at the repo root — then
`mull-runner` executes the binary once per mutant and reports each as killed
(a test failed), survived (all tests passed), or timed out (counted as
killed: mutations that hang loops).

```
tools/mutation.sh                          # full sweep of the mull.yml scope
MULL_DIFF_REF=master tools/mutation.sh     # only mutants on lines changed vs master
CTEST_FIELDS="11;32;64" tools/mutation.sh  # thorough: add the 64-bit field cases
```

Operational notes:

- The killer binary is built with a reduced field list (default `11;32`) so
  a full test run takes seconds. The default omits 64 deliberately: Mull's
  per-mutation-point runtime guards slow the *generic* 64-bit cases ~80×
  (the whole suite goes from ~1 s to ~190 s), and while killed mutants die
  at the first failing check, the warmup and every surviving mutant pay the
  full run.
- Scope and mutators live in `mull.yml`: the algorithm headers
  (`sketch_impl.h`, `int_utils.h`, `additive_fft.h` where present), with
  `cxx_default` plus the bitwise and logical groups (the XOR-heavy GF(2^m)
  code makes those the most meaningful operators). Widening `includePaths`
  is the lever for broader campaigns.
- Surviving mutants are findings about the *tests*, not necessarily bugs:
  triage each as (a) a missing assertion worth adding, (b) semantically
  equivalent to the original (e.g. mutating a value that is only an upper
  bound), or (c) unreachable under the reduced field list. Record notable
  survivors here rather than chasing a score of 100%.

Baseline (2026-08, `test-additions`, default scope and fields): **199
mutants, 152 killed, 47 survived (76% kill rate), 12 minutes** on 8
workers. Survivor triage, as an example of the categories above:

- ~6 `reserve()` capacity hints (`BerlekampMassey`, `FindRoots`) —
  equivalent, results cannot change.
- The GCD operand-order swap condition — equivalent, the algorithm is
  order-insensitive.
- The schoolbook `PolyMod`/`DivMod` inner-loop bound (`<` → `<=`) — the
  extra iteration writes one element past the remainder, which only a
  sanitizer can see: invisible to the plain unit binary but covered by the
  ASan `poly_ops` fuzz target. A good reminder that the layers back each
  other up.
- Genuine minor gaps: the `minisketch_decode` API error-path guards and the
  non-factorizable fast-abort counter in `RecFindRoots` are not exercised
  under the reduced field list.

The `MULL_DIFF_REF` mode is the intended per-change workflow (step 5 of
"Swapping a decode-stage implementation"): it mutates only the diff, so a
run completes in minutes even though a full-scope campaign takes an hour or
more.

## Performance measurement

Performance claims in this repo have been burned twice by measurement
artifacts, so both benchmarking and profiling are scripted with the
protocol built in.

### Benchmark A/B protocol

`tools/bench_ab.sh <ref-A> [ref-B=HEAD]` builds `bench` from two git refs
('.' = the current working tree) with identical flags (`-g -O2`), measures
every configuration in alternating-order rounds, and reports per side the
best time plus the min→max spread of the per-round values:

```
tools/bench_ab.sh master .                        # current tree vs master
CONFIGS="4096 4096 3 64" RUNS=7 tools/bench_ab.sh <ref> .
```

Protocol rules, learned the hard way:

- **Quiet machine.** No builds, indexing, or other load during measurement:
  concurrent compilation has inflated short rows ~2.5×. The script builds
  everything before the first measurement; do not edit files while it runs.
- **Alternate the order.** Same-order sweeps show order-correlated
  thermal/frequency drift that reads as a consistent fake delta. The script
  alternates sides every round.
- **A delta must beat the spread.** The noise floor on desktop hardware is
  a few percent; rows whose spread exceeds `NOISE_PCT` (default 3%) are
  flagged. Sub-millisecond rows additionally show bimodal scheduling
  jitter — only warm min-of-many runs see through it.
- **Suspect code layout before the algorithm.** A consistent delta that the
  changed code cannot explain (e.g. it never executes in that
  configuration) can be pure inlining/code-placement effects on the hot
  paths. The control that settles it: build ref-B with only the changed
  call sites reverted — unused templates are not instantiated, so a control
  matching ref-A attributes the delta to layout. This is how the FastGCD
  NOINLINE regression was found; keeping large cold entry points `NOINLINE`
  avoids the class.

### Causal profiling (Coz)

`tools/coz.sh [syndromes] [errors] [bits] [loops]` builds `bench` with
progress points (`-DMINISKETCH_COZ=ON`) and runs its profiling mode (5th
argument: decode every state repeatedly, one `COZ_PROGRESS` per decode,
highest available implementation only) under the
[Coz causal profiler](https://github.com/plasma-umass/coz):

```
sudo sysctl kernel.perf_event_paranoid=1   # once per boot; coz needs perf events
tools/coz.sh 1024 1024 64 50
(cd build-coz && coz plot)                 # or load the .jsonl in /usr/share/coz/viewer
```

Coz runs virtual-speedup experiments and answers "if this source line were
p% faster, how much faster would whole decodes complete?" — the number that
actually decides whether an optimization or cutoff change is worth
pursuing. Decode is single-threaded, so the contention effects causal
profiling is famous for (downward slopes) cannot appear here; the profile
is a measured per-line Amdahl chart, and its value over a conventional
profiler is that it reports *end-to-end impact* rather than time share.
Interpretation: steep upward slope = optimize this line; flat = don't
bother regardless of its time share. Pick sizes where decodes complete a
few times per second (syndromes 512–2048); each experiment needs several
progress-point visits, so very slow configurations converge slowly (use
more loops, or `COZ_ARGS="--end-to-end"`).

Besides the per-decode throughput point in bench, `MINISKETCH_COZ` builds
compile **decode-stage markers** into the library
(`MINISKETCH_COZ_BEGIN/END` in `util.h`, no-ops otherwise): latency pairs
around syndrome reconstruction, Berlekamp-Massey, and root finding, and
inside root finding around the trace computation, gcd, and factor
division. `tools/coz_summary.py` (run automatically by coz.sh) turns these
into a stage-share table plus the top causal lines. Two reading caveats:
stage shares are built from one active/inactive snapshot per experiment,
so their granularity is ~1/#experiments; and `-O2` inlining collapses much
of decode into few source lines (the `minisketch_decode` API line
"winning" means "the whole decode" and carries no information) — the
stage markers are the reliable axis, and line entries within
`sketch_impl.h` refine them.

Baseline example (this branch, pre-optimization algorithms, 64-bit,
512/512 and 1024/1024): root finding is ~100% of decode wall time, its
trace computation 94–100%, factor division ~6% at 512, and
syndromes/BM/gcd all sample at 0% — with `sketch_impl.h`'s `PolyMod` call
inside the `TraceMod` loop measuring ~94% program speedup at 100% line
speedup. That reproduces, by measurement, the cost model that motivated
the TraceMod-reducer and additive-FFT work on the perf branches: the
trace-loop reduction is the only thing worth optimizing at these sizes on
this branch.
