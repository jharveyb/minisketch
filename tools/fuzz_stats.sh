#!/usr/bin/env bash
# Run each fuzz target for a time budget and report corpus/coverage statistics.
#
# Usage: tools/fuzz_stats.sh [seconds-per-target] [build-dir]
#
# Uses (and grows) a persistent working corpus under <build-dir>/corpus/,
# seeded from the committed corpus in src/fuzz/corpus/. Reports libFuzzer's
# edge-coverage (cov:), feature (ft:) and corpus-size counters per target.
#
# Corpus tiering: the seconds-per-target budget bounds only the *fuzzing*
# phase, not the initial replay of the working corpus, so a corpus full of
# slow high-capacity decode/roundtrip inputs would make even a short run take
# minutes. To keep runs interactive, the working corpus is tiered by capacity
# just like the committed one: before each run, high-capacity (> 128) inputs
# are moved out of <build-dir>/corpus/<target>/ into <build-dir>/corpus-slow/,
# which is loaded only when SLOW=1. libFuzzer writes newly-discovered inputs
# into the fast dir, so any deep-capacity finds are re-tiered on the next run.
#
# Parallelism: set FORK to a worker count (e.g. FORK=$(nproc)) to fuzz each
# target with libFuzzer's fork mode, spreading mutation across cores.
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

# decode/roundtrip inputs vary hugely in decode cost (quadratic in capacity);
# above this capacity an input joins the slow tier. Matches update_seed_corpus.sh
# and the worst-case pruning boundary in src/fuzz/decode.cpp.
SLOW_CAP_THRESHOLD=128
is_tiered() { case "$1" in decode|roundtrip) return 0;; *) return 1;; esac; }

if [ ! -x "$BUILD_DIR/bin/fuzz" ]; then
    echo "error: $BUILD_DIR/bin/fuzz not found." >&2
    echo "Configure with: cmake -B $BUILD_DIR -DMINISKETCH_BUILD_FUZZ=ON -DCMAKE_CXX_COMPILER=clang++ && cmake --build $BUILD_DIR -j" >&2
    exit 1
fi

printf '%-12s %10s %10s %10s %12s\n' "target" "cov" "ft" "corpus" "execs"
for target in $TARGETS; do
    fast_work="$BUILD_DIR/corpus/$target"
    slow_work="$BUILD_DIR/corpus-slow/$target"
    artifacts="$BUILD_DIR/artifacts/$target"
    mkdir -p "$fast_work" "$artifacts"

    # Seed the fast working corpus from the committed fast tier (always); with
    # SLOW=1 also seed the slow working corpus from the committed slow tier.
    cp -n "src/fuzz/corpus/$target"/* "$fast_work/" 2> /dev/null || true
    if [ -n "${SLOW:-}" ]; then
        mkdir -p "$slow_work"
        cp -n "src/fuzz/corpus-slow/$target"/* "$slow_work/" 2> /dev/null || true
    fi

    # Keep the fast working corpus fast: move any high-capacity inputs (seeded,
    # or discovered by earlier fuzzing sessions) out to the slow working corpus.
    # Without this the working corpus re-accumulates the slow inputs that the
    # committed tiering keeps out of the fast seed set, and load time balloons.
    if is_tiered "$target"; then
        if [ -x "$BUILD_DIR/bin/corpus-tier" ]; then
            mkdir -p "$slow_work"
            while IFS=$'\t' read -r cap name; do
                if [ "$cap" -gt "$SLOW_CAP_THRESHOLD" ]; then
                    mv "$fast_work/$name" "$slow_work/" 2> /dev/null || true
                fi
            done < <("$BUILD_DIR/bin/corpus-tier" "$target" "$fast_work")
        else
            echo "warning: $BUILD_DIR/bin/corpus-tier missing; $target load may be slow" >&2
            echo "  build it with: cmake --build $BUILD_DIR --target corpus-tier" >&2
        fi
    fi

    # Corpora to load and fuzz into. libFuzzer writes new inputs to the first
    # dir, so discoveries land in the fast dir and get re-tiered next run.
    corpora=("$fast_work")
    [ -n "${SLOW:-}" ] && corpora+=("$slow_work")

    log="$BUILD_DIR/fuzz-$target.log"
    rc=0
    FUZZ="$target" "$BUILD_DIR/bin/fuzz" \
        -max_total_time="$SECONDS_PER_TARGET" -max_len="${FUZZ_MAX_LEN:-4096}" -print_final_stats=1 \
        -artifact_prefix="$artifacts/" \
        "${FORK_ARGS[@]}" ${FUZZ_ARGS:-} "${corpora[@]}" > "$log" 2>&1 || rc=$?
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
    # Corpus count = the loaded set: fast tier by default, both tiers on SLOW=1.
    corpus_count=$(ls "$fast_work" 2>/dev/null | wc -l)
    [ -n "${SLOW:-}" ] && corpus_count=$((corpus_count + $(ls "$slow_work" 2>/dev/null | wc -l)))
    printf '%-12s %10s %10s %10s %12s\n' "$target" "${cov:-?}" "${ft:-?}" "$corpus_count" "${execs:-?}"
done
echo
echo "Working corpora are kept in $BUILD_DIR/corpus/<target>/ (fast) and"
echo "$BUILD_DIR/corpus-slow/<target>/ (slow tier, loaded with SLOW=1), reused across runs."
echo "Crash and slow-unit artifacts (informational) go to $BUILD_DIR/artifacts/<target>/."
