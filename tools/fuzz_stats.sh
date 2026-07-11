#!/usr/bin/env bash
# Run each fuzz target for a time budget and report corpus/coverage statistics.
#
# Usage: tools/fuzz_stats.sh [seconds-per-target] [build-dir]
#
# Uses (and grows) a persistent working corpus under <build-dir>/corpus/,
# seeded from the committed corpus in src/fuzz/corpus/. Reports libFuzzer's
# edge-coverage (cov:), feature (ft:) and corpus-size counters per target.
#
# Parallelism: set FORK to a worker count (e.g. FORK=$(nproc)) to fuzz each
# target with libFuzzer's fork mode, spreading mutation across cores. Note
# the seconds-per-target budget bounds only the *fuzzing* phase, not the
# initial replay of the working corpus: once a target's corpus grows to many
# thousands of inputs (especially the slow high-capacity decode/roundtrip
# ones), that startup replay dominates wall time and the time budget has
# little visible effect. Fork mode parallelizes fuzzing but the startup still
# processes the whole corpus, so it helps long sessions more than quick snaps.
#
# Environment:
#   FUZZ_TARGETS  space-separated targets to run (default: all)
#   FUZZ_ARGS     extra libFuzzer arguments (e.g. -max_len=1024)
#   FUZZ_MAX_LEN  max input length (default 4096)
#   FORK          worker count for libFuzzer fork mode (default: unset = 1 process)
#   SLOW          1 = also load the high-capacity slow corpus tier
#                 (src/fuzz/corpus-slow/); default loads only the fast tier

set -euo pipefail

cd "$(dirname "$0")/.."
SECONDS_PER_TARGET="${1:-60}"
BUILD_DIR="${2:-build-fuzz}"
TARGETS="${FUZZ_TARGETS:-decode roundtrip poly_ops reconcile}"

# Fork-mode flags: ignore timeouts/OOMs (slow high-capacity inputs are
# expected, not bugs) but NOT crashes, so a real FUZZ_CHECK failure still
# stops the run and is reported below.
FORK_ARGS=()
if [ -n "${FORK:-}" ]; then
    FORK_ARGS=(-fork="$FORK" -ignore_timeouts=1 -ignore_ooms=1)
fi

if [ ! -x "$BUILD_DIR/bin/fuzz" ]; then
    echo "error: $BUILD_DIR/bin/fuzz not found." >&2
    echo "Configure with: cmake -B $BUILD_DIR -DMINISKETCH_BUILD_FUZZ=ON -DCMAKE_CXX_COMPILER=clang++ && cmake --build $BUILD_DIR -j" >&2
    exit 1
fi

printf '%-12s %10s %10s %10s %12s\n' "target" "cov" "ft" "corpus" "execs"
for target in $TARGETS; do
    workdir="$BUILD_DIR/corpus/$target"
    artifacts="$BUILD_DIR/artifacts/$target"
    mkdir -p "$workdir" "$artifacts"
    # Seed the working corpus from the committed fast tier. The slow tier
    # (high-capacity decode/roundtrip inputs) is loaded only when SLOW=1: it
    # is ~a third of the corpus but dominates startup replay while adding a
    # fraction of a percent of edge coverage, so skipping it keeps snapshot
    # runs interactive. See tools/update_seed_corpus.sh and doc/testing.md.
    for src in "src/fuzz/corpus/$target" \
               ${SLOW:+"src/fuzz/corpus-slow/$target"}; do
        if [ -d "$src" ]; then
            cp -n "$src"/* "$workdir/" 2> /dev/null || true
        fi
    done
    log="$BUILD_DIR/fuzz-$target.log"
    rc=0
    FUZZ="$target" "$BUILD_DIR/bin/fuzz" \
        -max_total_time="$SECONDS_PER_TARGET" -max_len="${FUZZ_MAX_LEN:-4096}" -print_final_stats=1 \
        -artifact_prefix="$artifacts/" \
        "${FORK_ARGS[@]}" ${FUZZ_ARGS:-} "$workdir" > "$log" 2>&1 || rc=$?
    if [ "$rc" -ne 0 ]; then
        if grep -q "run interrupted" "$log"; then
            echo "Interrupted during target $target (no failure); see $log" >&2
        else
            echo "FAILURE in target $target; see $log (crash input saved under $artifacts/)" >&2
        fi
        exit "$rc"
    fi
    # Tolerate missing matches (|| true): under `set -e` an empty grep would
    # otherwise abort the script. Fork mode does not emit the -print_final_stats
    # line, so fall back to the last iteration counter (#N:) for execs.
    cov=$(grep -o 'cov: [0-9]*' "$log" | tail -1 | cut -d' ' -f2 || true)
    ft=$(grep -o 'ft: [0-9]*' "$log" | tail -1 | cut -d' ' -f2 || true)
    execs=$(grep -o 'stat::number_of_executed_units: *[0-9]*' "$log" | grep -o '[0-9]*$' | tail -1 || true)
    [ -n "$execs" ] || execs=$(grep -oE '^#[0-9]+' "$log" | tr -d '#' | tail -1 || true)
    printf '%-12s %10s %10s %10s %12s\n' "$target" "${cov:-?}" "${ft:-?}" "$(ls "$workdir" | wc -l)" "${execs:-?}"
done
echo
echo "Working corpora are kept in $BUILD_DIR/corpus/<target>/ and reused across runs."
echo "Crash and slow-unit artifacts (informational) go to $BUILD_DIR/artifacts/<target>/."
