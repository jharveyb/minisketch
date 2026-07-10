#!/usr/bin/env bash
# Refresh the committed seed corpus (src/fuzz/corpus/<target>/) from the
# working corpora that fuzzing sessions grow under <build-dir>/corpus/.
#
# Usage: tools/update_seed_corpus.sh [build-dir]
#
# For each target: minimize (committed seeds + working corpus) with
# libFuzzer's -merge=1 to the subset preserving all observed coverage, then
# commit a small size-stratified sample of it (smallest / median / largest
# files) as the new seeds. Review and git add the result.
#
# Run this after long fuzzing sessions, and always after changing a target's
# input layout (added/removed Consume* calls) — old seeds still run but no
# longer mean what they did.
#
# Environment:
#   FUZZ_TARGETS      space-separated targets (default: all)
#   SEEDS_PER_TARGET  committed seeds per target (default: 24, rounded down
#                     to a multiple of 3)

set -euo pipefail

cd "$(dirname "$0")/.."
BUILD_DIR="${1:-build-fuzz}"
TARGETS="${FUZZ_TARGETS:-decode roundtrip poly_ops reconcile}"
PER_STRATUM=$(( ${SEEDS_PER_TARGET:-24} / 3 ))

if [ ! -x "$BUILD_DIR/bin/fuzz" ]; then
    echo "error: $BUILD_DIR/bin/fuzz not found." >&2
    exit 1
fi

for target in $TARGETS; do
    work="$BUILD_DIR/corpus/$target"
    seeds="src/fuzz/corpus/$target"
    if [ ! -d "$work" ] || [ -z "$(ls -A "$work" 2>/dev/null)" ]; then
        echo "$target: no working corpus in $work; skipping (run tools/fuzz_stats.sh first)"
        continue
    fi
    tmp=$(mktemp -d)
    FUZZ="$target" "$BUILD_DIR/bin/fuzz" -merge=1 "$tmp" "$seeds" "$work" > /dev/null 2>&1
    files=($(ls -S "$tmp")); n=${#files[@]}
    if [ "$n" -lt $(( 3 * PER_STRATUM )) ]; then
        echo "$target: minimized corpus has only $n files; keeping all of them"
        rm -f "$seeds"/*
        cp "$tmp"/* "$seeds/"
    else
        rm -f "$seeds"/*
        for (( i = 0; i < PER_STRATUM; ++i )); do
            cp "$tmp/${files[$i]}" "$seeds/"                       # largest
            cp "$tmp/${files[$(( n / 2 - PER_STRATUM / 2 + i ))]}" "$seeds/"  # median
            cp "$tmp/${files[$(( n - 1 - i ))]}" "$seeds/"         # smallest
        done
    fi
    rm -rf "$tmp"
    replay=$(FUZZ="$target" "$BUILD_DIR/bin/fuzz" -runs=0 "$seeds" 2>&1 | grep -o 'cov: [0-9]*' | tail -1)
    echo "$target: minimized $n -> $(ls "$seeds" | wc -l) seeds, replay $replay"
done
echo
echo "Review with 'git status' and commit the updated src/fuzz/corpus/."
