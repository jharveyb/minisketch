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

if [ ! -x "$BUILD_DIR/bin/fuzz" ]; then
    echo "error: $BUILD_DIR/bin/fuzz not found." >&2
    exit 1
fi

replay_cov() { # replay_cov <target> <dir>: edge coverage of replaying a corpus
    mkdir -p "$BUILD_DIR/artifacts/$1"
    FUZZ="$1" "$BUILD_DIR/bin/fuzz" -runs=0 -artifact_prefix="$BUILD_DIR/artifacts/$1/" "$2" 2>&1 | grep -o 'cov: [0-9]*' | tail -1 | cut -d' ' -f2
}

for target in $TARGETS; do
    work="$BUILD_DIR/corpus/$target"
    seeds="src/fuzz/corpus/$target"
    if [ ! -d "$work" ] || [ -z "$(ls -A "$work" 2>/dev/null)" ]; then
        echo "$target: no working corpus in $work; skipping (run tools/fuzz_stats.sh first)"
        continue
    fi
    mkdir -p "$seeds"
    old_count=$(ls "$seeds" | wc -l)
    work_count=$(ls "$work" | wc -l)
    old_cov=$(replay_cov "$target" "$seeds"); old_cov=${old_cov:-0}

    ctl="$BUILD_DIR/merge-$target.control"
    artifacts="$BUILD_DIR/artifacts/$target"
    mkdir -p "$artifacts"
    [ -f "$ctl" ] && echo "$target: resuming interrupted merge from $ctl"

    if [ "${RESET:-0}" = "1" ]; then
        candidate=$(mktemp -d)
        echo "$target: RESET-merging $work_count inputs; watch progress with: grep -c ^STARTED $ctl"
        FUZZ="$target" "$BUILD_DIR/bin/fuzz" -merge=1 -merge_control_file="$ctl" -artifact_prefix="$artifacts/" "$candidate" "$work" > /dev/null 2>&1
        rm -f "$ctl"
        new_cov=$(replay_cov "$target" "$candidate"); new_cov=${new_cov:-0}
        if [ "$new_cov" -lt "$old_cov" ]; then
            echo "$target: RESET candidate replays cov $new_cov < committed $old_cov; keeping committed set (fuzz longer first)"
        else
            rm -f "$seeds"/*
            cp "$candidate"/* "$seeds/"
            echo "$target: RESET $old_count -> $(ls "$seeds" | wc -l) seeds, replay cov $old_cov -> $new_cov"
        fi
        rm -rf "$candidate"
    else
        echo "$target: merging $((old_count + work_count)) inputs; watch progress with: grep -c ^STARTED $ctl"
        FUZZ="$target" "$BUILD_DIR/bin/fuzz" -merge=1 -merge_control_file="$ctl" -artifact_prefix="$artifacts/" "$seeds" "$work" > /dev/null 2>&1
        rm -f "$ctl"
        new_cov=$(replay_cov "$target" "$seeds"); new_cov=${new_cov:-0}
        echo "$target: $old_count -> $(ls "$seeds" | wc -l) seeds (+$(( $(ls "$seeds" | wc -l) - old_count ))), replay cov $old_cov -> $new_cov"
    fi
done
echo
echo "Review with 'git status' and commit the updated src/fuzz/corpus/."
