#!/usr/bin/env bash
# Update the committed seed corpus (src/fuzz/corpus/<target>/) from the
# working corpora that fuzzing sessions grow under <build-dir>/corpus/.
#
# Usage: tools/update_seed_corpus.sh [build-dir]
#
# Default (additive): libFuzzer's -merge=1 into the seed directory keeps all
# existing seeds and appends only those working-corpus inputs that add
# coverage features the seeds don't already have. This is monotone (committed
# coverage never decreases) and self-deduplicating; the set grows only while
# a target still has undiscovered features.
#
# RESET=1: rebuild a target's seed set from scratch (minimize the working
# corpus alone and replace). Use after changing a target's input layout
# (added/removed Consume* calls), when old seeds no longer decode to the
# cases they were selected for. A coverage guard refuses a reset that would
# replay less coverage than the committed set.
#
# Corpus tiering: decode and roundtrip inputs vary enormously in decode cost
# (quadratic in the `capacity` parameter), so their committed corpus is split
# into a fast tier (capacity <= 128, src/fuzz/corpus/) that is always loaded
# and a slow tier (capacity > 128, src/fuzz/corpus-slow/) loaded only on
# demand. This script keeps both tiers correct: it merges over their union so
# coverage is never lost, then re-partitions the result by capacity using the
# corpus-tier helper. poly_ops and reconcile are small and untiered.
#
# Environment:
#   FUZZ_TARGETS  space-separated targets (default: all)
#   RESET         1 = replace instead of append (guarded)
#
# Merges can take many minutes (every seed and working-corpus input is
# replayed once). Progress is observable via the per-target control file
# printed at the start of each merge; an interrupted merge resumes from
# that file on the next run.

set -euo pipefail

cd "$(dirname "$0")/.."
BUILD_DIR="${1:-build-fuzz}"
TARGETS="${FUZZ_TARGETS:-decode roundtrip poly_ops reconcile}"
# Capacity above which a decode/roundtrip input joins the slow tier. Matches
# the worst-case pruning boundary in src/fuzz/decode.cpp.
SLOW_CAP_THRESHOLD=128

if [ ! -x "$BUILD_DIR/bin/fuzz" ]; then
    echo "error: $BUILD_DIR/bin/fuzz not found." >&2
    exit 1
fi

is_tiered() { case "$1" in decode|roundtrip) return 0;; *) return 1;; esac; }

count() { ls "$1" 2>/dev/null | wc -l; }

copy_into() { # copy_into <dstdir> <srcdir>: copy srcdir's files if it exists
    [ -d "$2" ] && cp "$2"/* "$1/" 2>/dev/null || true
}

replay_cov() { # replay_cov <target> <dir>: edge coverage of replaying a corpus
    mkdir -p "$BUILD_DIR/artifacts/$1"
    FUZZ="$1" "$BUILD_DIR/bin/fuzz" -runs=0 -artifact_prefix="$BUILD_DIR/artifacts/$1/" "$2" 2>&1 | grep -o 'cov: [0-9]*' | tail -1 | cut -d' ' -f2
}

retier() { # retier <target> <srcdir>: split srcdir into the two committed tiers
    local target="$1" src="$2" fast="src/fuzz/corpus/$1" slow="src/fuzz/corpus-slow/$1"
    if [ ! -x "$BUILD_DIR/bin/corpus-tier" ]; then
        echo "error: $BUILD_DIR/bin/corpus-tier not found (build it: cmake --build $BUILD_DIR --target corpus-tier)." >&2
        exit 1
    fi
    local tf ts; tf=$(mktemp -d); ts=$(mktemp -d)
    while IFS=$'\t' read -r cap name; do
        if [ "$cap" -gt "$SLOW_CAP_THRESHOLD" ]; then cp "$src/$name" "$ts/"; else cp "$src/$name" "$tf/"; fi
    done < <("$BUILD_DIR/bin/corpus-tier" "$target" "$src")
    mkdir -p "$fast" "$slow"
    rm -f "$fast"/* "$slow"/*
    copy_into "$fast" "$tf"; copy_into "$slow" "$ts"
    rm -rf "$tf" "$ts"
    echo "$target: fast=$(count "$fast") slow=$(count "$slow")"
}

for target in $TARGETS; do
    work="$BUILD_DIR/corpus/$target"
    seeds="src/fuzz/corpus/$target"
    if [ ! -d "$work" ] || [ -z "$(ls -A "$work" 2>/dev/null)" ]; then
        echo "$target: no working corpus in $work; skipping (run tools/fuzz_stats.sh first)"
        continue
    fi
    mkdir -p "$seeds"
    ctl="$BUILD_DIR/merge-$target.control"
    artifacts="$BUILD_DIR/artifacts/$target"
    mkdir -p "$artifacts"
    [ -f "$ctl" ] && echo "$target: resuming interrupted merge from $ctl"
    work_count=$(count "$work")

    # The existing seed set is the union of both tiers (just the fast dir for
    # untiered targets). Merge into a copy of it so the merge sees every
    # committed input as already-present and only appends genuinely new ones.
    existing=$(mktemp -d)
    copy_into "$existing" "$seeds"
    is_tiered "$target" && copy_into "$existing" "src/fuzz/corpus-slow/$target"
    old_count=$(count "$existing")
    old_cov=$(replay_cov "$target" "$existing"); old_cov=${old_cov:-0}

    if [ "${RESET:-0}" = "1" ]; then
        candidate=$(mktemp -d)
        echo "$target: RESET-merging $work_count inputs; watch progress with: grep -c ^STARTED $ctl"
        FUZZ="$target" "$BUILD_DIR/bin/fuzz" -merge=1 -merge_control_file="$ctl" -artifact_prefix="$artifacts/" "$candidate" "$work" > /dev/null 2>&1
        rm -f "$ctl"
        new_cov=$(replay_cov "$target" "$candidate"); new_cov=${new_cov:-0}
        if [ "$new_cov" -lt "$old_cov" ]; then
            echo "$target: RESET candidate replays cov $new_cov < committed $old_cov; keeping committed set (fuzz longer first)"
        elif is_tiered "$target"; then
            echo "$target: RESET replay cov $old_cov -> $new_cov"
            retier "$target" "$candidate"
        else
            rm -f "$seeds"/*; cp "$candidate"/* "$seeds/"
            echo "$target: RESET $old_count -> $(count "$seeds") seeds, replay cov $old_cov -> $new_cov"
        fi
        rm -rf "$candidate"
    else
        echo "$target: merging $((old_count + work_count)) inputs; watch progress with: grep -c ^STARTED $ctl"
        FUZZ="$target" "$BUILD_DIR/bin/fuzz" -merge=1 -merge_control_file="$ctl" -artifact_prefix="$artifacts/" "$existing" "$work" > /dev/null 2>&1
        rm -f "$ctl"
        new_cov=$(replay_cov "$target" "$existing"); new_cov=${new_cov:-0}
        new_count=$(count "$existing")
        echo "$target: $old_count -> $new_count seeds (+$((new_count - old_count))), replay cov $old_cov -> $new_cov"
        if is_tiered "$target"; then
            retier "$target" "$existing"
        else
            rm -f "$seeds"/*; cp "$existing"/* "$seeds/"
        fi
    fi
    rm -rf "$existing"
done
echo
echo "Review with 'git status' and commit the updated src/fuzz/corpus/ and src/fuzz/corpus-slow/."
