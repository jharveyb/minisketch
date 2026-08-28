#!/usr/bin/env bash
# Causal profiling of sketch decode with Coz (https://github.com/plasma-umass/coz).
#
# Usage: tools/coz.sh [syndromes=1024] [errors=syndromes] [bits=64] [loops=50] [build-dir=build-coz]
#
# Coz runs virtual-speedup experiments: it measures how much speeding up each
# source line by p% would speed up the whole decode (the rate of the "decode"
# progress point in bench's profiling mode). Decode is single-threaded, so
# no contention effects (downward slopes) can appear; the profile is a
# measured per-line what-if chart — "optimize this line, gain that much" —
# which is exactly the number needed to decide whether an optimization or a
# cutoff change is worth pursuing.
#
# Configures a dedicated build directory with progress points enabled, builds
# bench, and runs it under coz in profiling mode (bench's 5th argument; only
# the highest available field implementation runs, so line samples are not
# smeared across implementations). Writes
#   <build-dir>/coz-s<syndromes>-e<errors>-b<bits>.jsonl
# View with `coz plot` from the profile's directory, or open the local viewer
# (/usr/share/coz/viewer/index.htm) and load the file.
#
# Pick sizes so that decodes complete a few times per second (syndromes
# 512-2048); a single experiment needs several progress-point visits, so very
# slow configurations (4096+) converge slowly - raise loops and let it run,
# or pass COZ_ARGS="--end-to-end" for one experiment per execution.
#
# Environment:
#   COZ_ARGS   extra arguments for `coz run` (e.g. --fixed-line file.h:123,
#              --fixed-speedup 50, --end-to-end)
#   SCOPE      --source-scope pattern ('%' wildcards; default %src% so all
#              library code is eligible for experiments)

set -euo pipefail

cd "$(dirname "$0")/.."
SYNDROMES="${1:-1024}"
ERRORS="${2:-$SYNDROMES}"
BITS="${3:-64}"
LOOPS="${4:-50}"
BUILD_DIR="${5:-build-coz}"
SCOPE="${SCOPE:-%src%}"

if ! command -v coz > /dev/null; then
    echo "error: coz not found; install the Coz profiler (e.g. apt install coz-profiler)" >&2
    exit 1
fi

PARANOID=$(cat /proc/sys/kernel/perf_event_paranoid)
if [ "$PARANOID" -gt 2 ]; then
    echo "error: kernel.perf_event_paranoid is $PARANOID; Coz needs perf events (<= 2)." >&2
    echo "Run once:  sudo sysctl kernel.perf_event_paranoid=1" >&2
    exit 1
fi

cmake -B "$BUILD_DIR" \
    -DMINISKETCH_BUILD_TESTS=OFF \
    -DMINISKETCH_BUILD_BENCHMARK=ON \
    -DMINISKETCH_COZ=ON \
    -DCMAKE_CXX_FLAGS="-g -O2 -fno-omit-frame-pointer"
cmake --build "$BUILD_DIR" --target bench -j "$(nproc)"

PROFILE="$BUILD_DIR/coz-s${SYNDROMES}-e${ERRORS}-b${BITS}.jsonl"
# shellcheck disable=SC2086  # COZ_ARGS is intentionally word-split
coz run --source-scope "$SCOPE" --output "$PROFILE" ${COZ_ARGS:-} --- \
    "$BUILD_DIR/bin/bench" "$SYNDROMES" "$ERRORS" 8 "$BITS" "$LOOPS"

echo
echo "Profile: $PROFILE"
echo "View:    (cd $(dirname "$PROFILE") && coz plot)   # or open /usr/share/coz/viewer/index.htm"
