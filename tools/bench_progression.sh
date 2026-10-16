#!/usr/bin/env bash
# Pinned A/B benchmarks across the decode-optimization branches, with all
# measurements kept for plotting.
#
# Usage (from a clone that has the branches below):
#   tools/bench_progression.sh <host-tag>        e.g. i7_4770, pi5_a76
#
# Pins the CPU with tools/bench_host.sh (needs passwordless sudo), runs
# each A/B with tools/bench_ab.sh (CLMUL implementation only), restores and
# verifies the CPU settings on exit, and writes to
# ~/ms-bench-results/<host-tag>/:
#   <host>_ab_<name>.txt         bench_ab.sh's printed output
#   <host>_ab_<name>/rounds.tsv  every per-round measurement (plot from this)
#   <host>_ab_<name>/summary.tsv, meta.txt
#   <host>_clock.log             pinned-core clock/temperature/throttling
#                                after each A/B
#   <host>_bench_host_*.txt      pinning setup/restore/verify output
#
# Environment: RESULTS (output root, default ~/ms-bench-results), ONLY (a
# '|'-separated list of A/B names to run, default all), FREQ_KHZ (passed to
# bench_host.sh).

set -uo pipefail

HOST="${1:?usage: tools/bench_progression.sh <host-tag>}"
cd "$(dirname "$0")/.."
R="${RESULTS:-$HOME/ms-bench-results}/$HOST"
STATE=/var/tmp/minisketch-bench-host.state
ONLY="${ONLY:-}"

if [ -e "$STATE" ]; then
    echo "error: $STATE exists (a pinning is active or stale); restore or remove it first" >&2
    exit 1
fi
mkdir -p "$R"

# Restore, then check every CPU against the recorded state and fix any
# governor / frequency limit that did not take (seen on intel_pstate after an
# interrupted restore) through sysfs directly.
restore_and_verify() {
    sudo tools/bench_host.sh restore 2>&1 | tee "$R/${HOST}_bench_host_restore.txt"
    local rec="$STATE.restored" line c gov min max d fixed=0
    [ -e "$rec" ] || rec="$STATE"
    [ -e "$rec" ] || { echo "verify: no recorded state found" | tee -a "$R/${HOST}_bench_host_restore.txt"; return; }
    while read -r line; do
        c=$(sed -n 's/^cpu\([0-9]*\) .*/\1/p' <<< "$line")
        [ -n "$c" ] || continue
        gov=$(sed -n 's/.*governor=\([^ ]*\).*/\1/p' <<< "$line")
        min=$(sed -n 's/.* min=\([0-9]*\).*/\1/p' <<< "$line")
        max=$(sed -n 's/.* max=\([0-9]*\).*/\1/p' <<< "$line")
        d=/sys/devices/system/cpu/cpu$c/cpufreq
        [ -n "$gov" ] && [ "$(cat $d/scaling_governor)" != "$gov" ] && { echo "$gov" | sudo tee $d/scaling_governor > /dev/null; fixed=$((fixed + 1)); }
        if [ -n "$max" ] && [ "$(cat $d/scaling_max_freq)" != "$max" ]; then echo "$max" | sudo tee $d/scaling_max_freq > /dev/null; fixed=$((fixed + 1)); fi
        if [ -n "$min" ] && [ "$(cat $d/scaling_min_freq)" != "$min" ]; then echo "$min" | sudo tee $d/scaling_min_freq > /dev/null; fixed=$((fixed + 1)); fi
    done < <(grep '^cpu' "$rec")
    {
        echo "verify: fixed $fixed setting(s) via sysfs"
        sudo tools/bench_host.sh status
    } 2>&1 | tee -a "$R/${HOST}_bench_host_restore.txt"
}

sudo FREQ_KHZ="${FREQ_KHZ:-}" tools/bench_host.sh setup 2>&1 | tee "$R/${HOST}_bench_host_setup.txt"
[ -e "$STATE" ] || { echo "error: pinning setup failed" >&2; exit 1; }
trap restore_and_verify EXIT
PIN=$(sed -n 's/.*BENCH_PREFIX="taskset -c \([0-9]*\)".*/\1/p' "$R/${HOST}_bench_host_setup.txt")
[ -n "$PIN" ] || { echo "error: no pin CPU in the setup output" >&2; exit 1; }

export FIELDS="32;64" RUNS=5 IMPL=1 BENCH_PREFIX="taskset -c $PIN"
# Full-capacity configs (the progression protocol).
FULL="128 128 7 64|512 512 7 64|1024 1024 5 64|2048 2048 3 64|4096 4096 3 64|1024 1024 5 32|4096 4096 3 32"
# High capacity, low-to-medium differences, plus 50/75% fill.
TARGET="1024 32 30 64|1024 64 20 64|1024 128 10 64|1024 256 7 64|2048 64 15 64|2048 256 5 64|4096 32 15 64|4096 64 10 64|4096 128 7 64|4096 256 5 64|8192 64 7 64|8192 256 4 64|1024 512 4 64|1024 768 3 64|4096 2048 2 64|4096 3072 2 64|1024 64 20 32|1024 256 7 32|4096 256 5 32"
# Root degrees around the CLMUL table/FFT reducer crossover (704 vs 1024).
CUTOFF="1024 512 4 64|1024 720 3 64|1024 768 3 64|1024 832 3 64|1024 896 3 64|1024 960 2 64|1024 1000 2 64|2048 1000 2 64|1024 768 3 32|1024 1000 2 32"

clock() {
    {
        echo "--- $(date -Is) after $1"
        grep MHz /proc/cpuinfo | sed -n "$((PIN + 1))p"
        echo "scaling_cur_freq $(cat /sys/devices/system/cpu/cpu$PIN/cpufreq/scaling_cur_freq)"
        echo "temp $(cat /sys/class/thermal/thermal_zone0/temp 2> /dev/null)"
        command -v vcgencmd > /dev/null && sudo vcgencmd get_throttled
        true
    } >> "$R/${HOST}_clock.log" 2>&1
}

ab() { # name refA refB configs [VAR=value...]
    local name="$1" a="$2" b="$3" cfg="$4"
    shift 4
    if [ -n "$ONLY" ] && [[ "|$ONLY|" != *"|$name|"* ]]; then return; fi
    echo "=== $(date -Is) $name: $a -> $b"
    env CONFIGS="$cfg" DATA_DIR="$R/${HOST}_ab_${name}" "$@" tools/bench_ab.sh "$a" "$b" 2>&1 | tee "$R/${HOST}_ab_${name}.txt"
    clock "$name"
}

clock start
ab 1-reducer test-additions fast-tracemod-reduce-allocs "$FULL"
ab 2-fft fast-tracemod-reduce-allocs additive-fft "$FULL"
ab 3-gcd additive-fft subquadratic-gcd "$FULL"
ab 4-cumulative test-additions subquadratic-gcd "$FULL"
ab 5-math subquadratic-gcd general_math_speedup "$FULL"
ab 6-cumulative-math test-additions general_math_speedup "$FULL"
ab 7-math-target subquadratic-gcd general_math_speedup "$TARGET"
ab 8-table-cutoff-1024 general_math_speedup general_math_speedup "$CUTOFF" CXXFLAGS_B=-DMINISKETCH_TRACEMOD_FFT_TABLE_CUTOFF=1024
echo "=== $(date -Is) all done"
