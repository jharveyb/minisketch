# Notes for Claude Code

## Commits and pull requests

- Never add a `Claude-Session:` trailer, a claude.ai session link, or any
  other session/conversation URL to commit messages, PR descriptions, or
  code. This repository is public. (`.claude/settings.json` also disables
  it via `attribution.sessionUrl: false`.)
- A `Co-Authored-By: Claude ... <noreply@anthropic.com>` trailer is fine.

## Build and test

- `cmake -B build -DCMAKE_BUILD_TYPE=Release -DMINISKETCH_BUILD_BENCHMARK=ON && cmake --build build` - library, `bench`, `test-verify`/`test-noverify`.
- Boost.Test `unit-tests` only builds inside `nix develop` (Boost comes from the flake devshell).
- `build/bin/test-verify [complexity]` - legacy randomized suite over all fields and implementations; `unit-tests` - property tests (BM vs textbook oracle, FindRoots, poly ops).
- `-DMINISKETCH_FIELDS="32;64"` speeds builds, but the fuzz build needs field 11 and `tests/differential.py` needs field 2: build those with all fields.
- Fuzzing: `cmake -B build-fuzz -DMINISKETCH_BUILD_FUZZ=ON -DCMAKE_CXX_COMPILER=clang++`; libFuzzer writes new inputs into the *first* corpus dir, so pass a scratch dir first, then `src/fuzz/corpus/<target>`.
- `MULL_DIFF_REF=<ref> tools/mutation.sh <build-dir>` - Mull on changed lines only; the build dir must be repo-relative.

## Benchmarking

- Pin first: `sudo tools/bench_host.sh setup` ... `sudo tools/bench_host.sh restore`, and run with `BENCH_PREFIX="taskset -c <cpu it suggests>"`. Unpinned numbers (turbo, ondemand/schedutil) are not comparable.
- `tools/bench_ab.sh <refA> <refB>` - interleaved A/B; `IMPL=1` = CLMUL only, `DATA_DIR=<dir>` keeps per-round TSV for plotting; `tools/bench_progression.sh <host-tag>` runs the branch-progression plan.
- `-DMINISKETCH_STAGE_TIMING` - per-stage decode timings to stderr from the COZ markers; its timer overhead inflates tiny per-node calls, so confirm small wins with a real A/B.
- Algorithm cutoffs in `src/sketch_impl.h` (`MINISKETCH_*` macros) were tuned on an i7-4770 and a Pi 5 (Cortex-A76); retune with `CXXFLAGS_B=-DMINISKETCH_<NAME>=<v> tools/bench_ab.sh . .`.

## Code gotchas (src/sketch_impl.h)

- Most hot loops have two paths: `F::LAZY_REDUCE` (CLMUL/PMULL, unreduced `Wide` accumulation) and Multiplier tables (generic fields). Change and test both.
- `BerlekampMassey` requires Frobenius-structured syndromes (S_2k = S_k^2); decode calls `BerlekampMasseyOdd` on the sketch's odd syndromes.
- Don't pass class-scope `static const` members to `std::min`/`std::max` (ODR-use): it links at -O2 but fails in unoptimized builds such as fuzz.
- See `doc/optimizations.md` for the decode cost model, measured results and negative results before trying a new optimization.
