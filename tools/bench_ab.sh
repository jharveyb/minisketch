#!/usr/bin/env bash
# Interleaved A/B decode benchmark between two versions of the library.
#
# Usage: tools/bench_ab.sh <ref-A> [ref-B=HEAD]
#
# Builds bench from two git refs (use '.' for the current working tree,
# including uncommitted changes) with identical canonical flags (-g -O2),
# then measures each configuration in RUNS rounds with the side order
# alternating every round, and reports per side the best (min) time and the
# min->max spread of the per-round values. The protocol encodes the lessons
# from the decode-optimization series:
#
#  - Strictly no builds or other load during measurement: concurrent
#    compilation/indexing has been observed to inflate short rows ~2.5x.
#  - Alternating order decorrelates thermal/frequency drift; same-order
#    sweeps have shown consistent order-correlated bias.
#  - A conclusion requires the A/B delta to exceed the reported spread. On
#    this class of hardware the noise floor is a few percent; spreads above
#    NOISE_PCT are flagged. If a surprising delta survives, suspect code
#    layout before the algorithm: build a control at the new ref with the
#    changed call sites reverted (unused templates are not instantiated) -
#    a control matching the old ref attributes the delta to layout/inlining,
#    as happened with the FastGCD NOINLINE find.
#
# The bench binary needs the bits argument (present since the additive-FFT
# series); older refs need that commit cherry-picked to be measurable.
#
# Environment:
#   CONFIGS    '|'-separated "syndromes errors iters bits" rows
#              (default "1024 1024 7 64|2048 2048 5 64|4096 4096 3 64")
#   RUNS       rounds per configuration (default 5)
#   NOISE_PCT  spread percentage above which a row is flagged (default 3)

set -euo pipefail

cd "$(dirname "$0")/.."
REF_A="${1:?usage: tools/bench_ab.sh <ref-A> [ref-B=HEAD]}"
REF_B="${2:-HEAD}"
CONFIGS="${CONFIGS:-32 32 7 64| 128 128 7 64| 512 512 7 64| 1024 1024 7 64|2048 2048 5 64|4096 4096 3 64}"
RUNS="${RUNS:-5}"
NOISE_PCT="${NOISE_PCT:-3}"

ROOT="$(mktemp -d "${TMPDIR:-/tmp}/minisketch-bench-ab.XXXXXX")"
WORKTREES=()
cleanup() {
    for wt in "${WORKTREES[@]:-}"; do
        [ -n "$wt" ] && git worktree remove --force "$wt" 2> /dev/null || true
    done
    rm -rf "$ROOT"
}
trap cleanup EXIT

# Build bench for one side; prints the binary path. Configure/build output
# (including cmake's stderr summary and compiler warnings) goes to a log,
# shown only on failure.
build_side() {
    local ref="$1" dir="$ROOT/$2" src="." log="$ROOT/$2-build.log"
    if [ "$ref" != "." ]; then
        git worktree add --force --detach "$dir-src" "$ref" > /dev/null 2>&1
        WORKTREES+=("$dir-src")
        src="$dir-src"
    fi
    if ! {
        cmake -S "$src" -B "$dir" \
            -DMINISKETCH_BUILD_TESTS=OFF -DMINISKETCH_BUILD_BENCHMARK=ON \
            -DCMAKE_CXX_FLAGS="-g -O2" &&
        cmake --build "$dir" --target bench -j "$(nproc)"
    } > "$log" 2>&1; then
        cat "$log" >&2
        echo "error: build of $ref failed" >&2
        exit 1
    fi
    echo "$dir/bin/bench"
}

# A bench binary predating the bits argument prints its usage line to stdout
# (which the measurement loop redirects into the data file) and exits 1 -
# probe each binary up front so that fails loudly instead of silently.
probe_side() {
    local ref="$1" bin="$2" out
    if ! out=$("$bin" 16 16 1 64 2>&1); then
        printf '%s\n' "$out" >&2
        echo "error: bench at $ref does not accept the 4-argument form (bits)." >&2
        echo "Cherry-pick the bench bits-argument commit onto it first, e.g.:" >&2
        echo "  git worktree add /tmp/wt $ref && git -C /tmp/wt cherry-pick 722fed56" >&2
        echo "then benchmark the resulting commit." >&2
        exit 1
    fi
}

echo "Building A: $REF_A"
BIN_A=$(build_side "$REF_A" a)
echo "Building B: $REF_B"
BIN_B=$(build_side "$REF_B" b)
probe_side "$REF_A" "$BIN_A"
probe_side "$REF_B" "$BIN_B"
echo "Measuring: $RUNS alternating rounds per config; keep the machine idle."
echo

IFS='|' read -ra CONFIG_ROWS <<< "$CONFIGS"
DATA="$ROOT/data"
mkdir -p "$DATA"

for cfg in "${CONFIG_ROWS[@]}"; do
    read -r syn err it bits <<< "$cfg"
    tag="s${syn}-e${err}-b${bits}"
    for (( run = 1; run <= RUNS; run++ )); do
        if (( run % 2 )); then order="A B"; else order="B A"; fi
        for side in $order; do
            bin="$BIN_A"; [ "$side" = B ] && bin="$BIN_B"
            if ! "$bin" "$syn" "$err" "$it" "$bits" \
                >> "$DATA/${side}_${tag}.txt"; then
                echo "error: bench (side $side, config '$cfg') failed; its output:" >&2
                tail -5 "$DATA/${side}_${tag}.txt" >&2
                exit 1
            fi
        done
    done
done

printf '%-18s %-4s %-22s %-22s\n' "config" "side" "impl0 min (spread%)" "impl1 min (spread%)"
for cfg in "${CONFIG_ROWS[@]}"; do
    read -r syn err it bits <<< "$cfg"
    tag="s${syn}-e${err}-b${bits}"
    for side in A B; do
        awk -v b="$bits" -v tag="$tag" -v side="$side" -v noise="$NOISE_PCT" '
            $1 == "recover[ms]" && $2 == b {
                for (i = 3; i <= 4; i++) {
                    c = i - 2
                    if ($i == "-") continue
                    if (!seen[c] || $i < min[c]) min[c] = $i
                    if (!seen[c] || $i > max[c]) max[c] = $i
                    seen[c] = 1
                }
            }
            END {
                printf "%-18s %-4s", tag, side
                for (c = 1; c <= 2; c++) {
                    if (!seen[c]) { printf " %-22s", "-"; continue }
                    spread = 100 * (max[c] - min[c]) / min[c]
                    flag = spread > noise ? " NOISY" : ""
                    printf " %-22s", sprintf("%.3f (%.1f%%)%s", min[c], spread, flag)
                }
                printf "\n"
            }' "$DATA/${side}_${tag}.txt"
    done
done

echo
echo "min over $RUNS alternating rounds per side; spread = (max-min)/min of the"
echo "per-round values. Trust a delta only if it exceeds both sides' spreads."
