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
# Both bench CLIs are supported, detected per side: the flag CLI
# (--syndromes/--errors/--iters/--bits/--decode-only, long-format TSV output)
# and the older positional one ("syndromes errors iters bits", one wide row per
# field size). Positional refs need the bits argument (present since the
# additive-FFT series); older ones need that commit cherry-picked.
#
# Environment:
#   CONFIGS    '|'-separated "syndromes errors iters bits" rows
#              (default "1024 1024 7 64|2048 2048 5 64|4096 4096 3 64")
#   RUNS       rounds per configuration (default 5)
#   NOISE_PCT  spread percentage above which a row is flagged (default 3)
#   CXXFLAGS_A, CXXFLAGS_B
#              extra compiler flags for one side, appended to -g -O2. With
#              both refs '.', this A/Bs a build-time knob within one tree,
#              e.g. CXXFLAGS_B=-DMINISKETCH_HGCD_CUTOFF=2048 tools/bench_ab.sh . .
#   FIELDS     restrict the build to these field sizes (MINISKETCH_FIELDS,
#              ';'-separated, e.g. "32;64"); must cover every bits value in
#              CONFIGS. Cuts build time several-fold on slow machines.
#   BENCH_PREFIX
#              command prefixed to every bench run, e.g. "taskset -c 5" to
#              pin to one core on a big.LITTLE CPU (word-split)
#   DATA_DIR   if set, keep the measurements there for plotting: rounds.tsv
#              (every per-round bench row, as side/ref/commit/round followed
#              by bench's long-format columns), summary.tsv (the printed
#              table, one row per config, side and implementation) and
#              meta.txt (host, refs, commits, flags, CPU frequency state)

set -euo pipefail

cd "$(dirname "$0")/.."
REF_A="${1:?usage: tools/bench_ab.sh <ref-A> [ref-B=HEAD]}"
REF_B="${2:-HEAD}"
CONFIGS="${CONFIGS:-32 32 7 64| 128 128 7 64| 512 512 7 64| 1024 1024 7 64|2048 2048 5 64|4096 4096 3 64}"
RUNS="${RUNS:-5}"
NOISE_PCT="${NOISE_PCT:-3}"
CXXFLAGS_A="${CXXFLAGS_A:-}"
CXXFLAGS_B="${CXXFLAGS_B:-}"
FIELDS="${FIELDS:-}"
DATA_DIR="${DATA_DIR:-}"
read -ra BENCH_PREFIX_ARGS <<< "${BENCH_PREFIX:-}"

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
    local ref="$1" dir="$ROOT/$2" src="." log="$ROOT/$2-build.log" extra="$3"
    local fields_arg=()
    [ -n "$FIELDS" ] && fields_arg=("-DMINISKETCH_FIELDS=$FIELDS")
    if [ "$ref" != "." ]; then
        git worktree add --force --detach "$dir-src" "$ref" > /dev/null 2>&1
        WORKTREES+=("$dir-src")
        src="$dir-src"
    fi
    if ! {
        cmake -S "$src" -B "$dir" \
            -DMINISKETCH_BUILD_TESTS=OFF -DMINISKETCH_BUILD_BENCHMARK=ON \
            "${fields_arg[@]}" -DCMAKE_CXX_FLAGS="-g -O2 $extra" &&
        cmake --build "$dir" --target bench -j "$(nproc)"
    } > "$log" 2>&1; then
        cat "$log" >&2
        echo "error: build of $ref failed" >&2
        exit 1
    fi
    echo "$dir/bin/bench"
}

# Run one decode measurement with the CLI style detected for the binary.
run_bench() {
    local style="$1" bin="$2" syn="$3" err="$4" it="$5" bits="$6"
    if [ "$style" = flags ]; then
        "${BENCH_PREFIX_ARGS[@]}" "$bin" --syndromes "$syn" --errors "$err" --iters "$it" --bits "$bits" --decode-only
    else
        "${BENCH_PREFIX_ARGS[@]}" "$bin" "$syn" "$err" "$it" "$bits"
    fi
}

# Detect which CLI a binary speaks; prints "flags" or "positional". A binary
# rejecting both prints its usage line to stdout (which the measurement loop
# would redirect into the data file) and exits 1 - probe each binary up front
# so that fails loudly instead of silently.
probe_side() {
    local ref="$1" bin="$2" out style
    for style in flags positional; do
        if run_bench "$style" "$bin" 16 16 1 64 > /dev/null 2>&1; then
            echo "$style"
            return
        fi
    done
    out=$(run_bench positional "$bin" 16 16 1 64 2>&1 || true)
    printf '%s\n' "$out" >&2
    echo "error: bench at $ref accepts neither the flag CLI (--bits, --decode-only)" >&2
    echo "nor the 4-argument positional form (bits)." >&2
    echo "Cherry-pick the bench bits-argument commit onto it first, e.g.:" >&2
    echo "  git worktree add /tmp/wt $ref && git -C /tmp/wt cherry-pick 722fed56" >&2
    echo "then benchmark the resulting commit." >&2
    exit 1
}

echo "Building A: $REF_A${CXXFLAGS_A:+ ($CXXFLAGS_A)}"
BIN_A=$(build_side "$REF_A" a "$CXXFLAGS_A")
echo "Building B: $REF_B${CXXFLAGS_B:+ ($CXXFLAGS_B)}"
BIN_B=$(build_side "$REF_B" b "$CXXFLAGS_B")
STYLE_A=$(probe_side "$REF_A" "$BIN_A")
STYLE_B=$(probe_side "$REF_B" "$BIN_B")
COMMIT_A=$(git rev-parse --short "${REF_A/#./HEAD}")
COMMIT_B=$(git rev-parse --short "${REF_B/#./HEAD}")
[ "$REF_A" = . ] && COMMIT_A="$COMMIT_A+worktree"
[ "$REF_B" = . ] && COMMIT_B="$COMMIT_B+worktree"
if [ -n "$DATA_DIR" ]; then
    mkdir -p "$DATA_DIR"
    printf 'side\tref\tcommit\tround\tmetric\tbits\tcapacity\terrors\tdata_len\timplementation\tvalue\n' > "$DATA_DIR/rounds.tsv"
    {
        echo "host: $(uname -n) ($(uname -m)), $(date -Is)"
        echo "cpu: $(lscpu 2> /dev/null | sed -n 's/^Model name: *//p' | head -1)"
        echo "A: $REF_A ($COMMIT_A)${CXXFLAGS_A:+ CXXFLAGS_A=$CXXFLAGS_A}"
        echo "B: $REF_B ($COMMIT_B)${CXXFLAGS_B:+ CXXFLAGS_B=$CXXFLAGS_B}"
        echo "configs: $CONFIGS"
        echo "runs: $RUNS, fields: ${FIELDS:-all}, bench prefix: ${BENCH_PREFIX:-none}, flags: -g -O2"
        for f in /sys/devices/system/cpu/cpu[0-9]*/cpufreq; do
            [ -d "$f" ] || continue
            echo "$(basename "$(dirname "$f")"): governor $(cat "$f/scaling_governor") min $(cat "$f/scaling_min_freq") max $(cat "$f/scaling_max_freq") kHz"
        done
        [ -r /sys/devices/system/cpu/intel_pstate/no_turbo ] && echo "intel no_turbo: $(cat /sys/devices/system/cpu/intel_pstate/no_turbo)"
        [ -r /sys/devices/system/cpu/cpufreq/boost ] && echo "cpufreq boost: $(cat /sys/devices/system/cpu/cpufreq/boost)"
        true
    } > "$DATA_DIR/meta.txt"
fi
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
            bin="$BIN_A"; style="$STYLE_A"
            [ "$side" = B ] && { bin="$BIN_B"; style="$STYLE_B"; }
            if ! run_bench "$style" "$bin" "$syn" "$err" "$it" "$bits" \
                > "$DATA/round.txt"; then
                echo "error: bench (side $side, config '$cfg') failed; its output:" >&2
                tail -5 "$DATA/round.txt" >&2
                exit 1
            fi
            cat "$DATA/round.txt" >> "$DATA/${side}_${tag}.txt"
            if [ -n "$DATA_DIR" ]; then
                ref="$REF_A"; commit="$COMMIT_A"
                [ "$side" = B ] && { ref="$REF_B"; commit="$COMMIT_B"; }
                awk -v pre="$side\t$ref\t$commit\t$run" -v OFS='\t' '
                    $1 == "recover[ms]" && NF == 7 { print pre, $0 }
                    # Wide row (positional CLI): metric bits impl0 impl1
                    $1 == "recover[ms]" && NF != 7 {
                        print pre, $1, $2, "'"$syn"'", "'"$err"'", 0, "GENERIC", $3
                        if ($4 != "" && $4 != "-") print pre, $1, $2, "'"$syn"'", "'"$err"'", 0, "CLMUL", $4
                    }' "$DATA/round.txt" >> "$DATA_DIR/rounds.tsv"
            fi
        done
    done
done

[ -n "$DATA_DIR" ] && printf 'config\tside\tref\tcommit\timplementation\tmin_ms\tmax_ms\tspread_pct\n' > "$DATA_DIR/summary.tsv"
printf '%-18s %-4s %-22s %-22s\n' "config" "side" "impl0 min (spread%)" "impl1 min (spread%)"
for cfg in "${CONFIG_ROWS[@]}"; do
    read -r syn err it bits <<< "$cfg"
    tag="s${syn}-e${err}-b${bits}"
    for side in A B; do
        ref="$REF_A"; commit="$COMMIT_A"
        [ "$side" = B ] && { ref="$REF_B"; commit="$COMMIT_B"; }
        awk -v b="$bits" -v tag="$tag" -v side="$side" -v noise="$NOISE_PCT" \
            -v tsv="${DATA_DIR:+$DATA_DIR/summary.tsv}" -v ref="$ref" -v commit="$commit" '
            function add(c, v) {
                if (v == "-") return
                if (!seen[c] || v < min[c]) min[c] = v
                if (!seen[c] || v > max[c]) max[c] = v
                seen[c] = 1
            }
            # Long-format TSV (flag CLI): metric bits capacity errors data_len impl value
            $1 == "recover[ms]" && $2 == b && NF == 7 {
                if ($6 == "GENERIC") add(1, $7)
                else if ($6 == "CLMUL") add(2, $7)
            }
            # Wide row (positional CLI): metric bits impl0 impl1 ...
            $1 == "recover[ms]" && $2 == b && NF != 7 { add(1, $3); add(2, $4) }
            END {
                if (tsv != "") {
                    for (c = 1; c <= 2; c++) if (seen[c])
                        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.2f\n", tag, side, ref, commit, (c == 1 ? "GENERIC" : "CLMUL"), min[c], max[c], 100 * (max[c] - min[c]) / min[c] >> tsv
                }
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
