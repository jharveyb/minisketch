#!/usr/bin/env bash
# Run each fuzz target for a time budget and report corpus/coverage statistics.
#
# Usage: tools/fuzz_stats.sh [seconds-per-target] [build-dir]
#
# Uses (and grows) a persistent working corpus under <build-dir>/corpus/,
# seeded from the committed corpus in src/fuzz/corpus/. Reports libFuzzer's
# edge-coverage (cov:), feature (ft:) and corpus-size counters per target.
#
# Environment:
#   FUZZ_TARGETS  space-separated targets to run (default: all)
#   FUZZ_ARGS     extra libFuzzer arguments (e.g. -max_len=1024)

set -euo pipefail

cd "$(dirname "$0")/.."
SECONDS_PER_TARGET="${1:-60}"
BUILD_DIR="${2:-build-fuzz}"
TARGETS="${FUZZ_TARGETS:-decode roundtrip poly_ops reconcile}"

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
    if [ -d "src/fuzz/corpus/$target" ]; then
        cp -n src/fuzz/corpus/"$target"/* "$workdir/" 2> /dev/null || true
    fi
    log="$BUILD_DIR/fuzz-$target.log"
    rc=0
    FUZZ="$target" "$BUILD_DIR/bin/fuzz" \
        -max_total_time="$SECONDS_PER_TARGET" -max_len="${FUZZ_MAX_LEN:-4096}" -print_final_stats=1 \
        -artifact_prefix="$artifacts/" \
        ${FUZZ_ARGS:-} "$workdir" > "$log" 2>&1 || rc=$?
    if [ "$rc" -ne 0 ]; then
        if grep -q "run interrupted" "$log"; then
            echo "Interrupted during target $target (no failure); see $log" >&2
        else
            echo "FAILURE in target $target; see $log (crash input saved under $artifacts/)" >&2
        fi
        exit "$rc"
    fi
    cov=$(grep -o 'cov: [0-9]*' "$log" | tail -1 | cut -d' ' -f2)
    ft=$(grep -o 'ft: [0-9]*' "$log" | tail -1 | cut -d' ' -f2)
    execs=$(grep -o 'stat::number_of_executed_units: *[0-9]*' "$log" | grep -o '[0-9]*$')
    printf '%-12s %10s %10s %10s %12s\n' "$target" "${cov:-?}" "${ft:-?}" "$(ls "$workdir" | wc -l)" "${execs:-?}"
done
echo
echo "Working corpora are kept in $BUILD_DIR/corpus/<target>/ and reused across runs."
echo "Crash and slow-unit artifacts (informational) go to $BUILD_DIR/artifacts/<target>/."
