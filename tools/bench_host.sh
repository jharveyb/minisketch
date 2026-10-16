#!/usr/bin/env bash
# shellcheck disable=SC2034  # arrays are passed to apply() by name
# Prepare a Linux host for stable benchmarking with cpupower, and undo it.
#
# Usage: sudo tools/bench_host.sh setup     record current settings, then pin
#        sudo tools/bench_host.sh restore   restore the recorded settings
#             tools/bench_host.sh status    show current settings
#
# setup records the per-CPU governor, policy frequency limits, energy
# performance preference (EPP) and perf-bias (EPB), plus the boost state,
# to a human-readable state file (followed by the full `cpupower
# frequency-info` / `cpupower info` output for reference), then, as far as
# the platform supports each step:
#  - sets EPP to "performance" and EPB to 0 (max performance)
#  - disables turbo/boost
#  - sets the governor to "performance"
#  - pins every CPU to its own highest non-boost frequency (min = max):
#    each cluster of a big.LITTLE / hybrid CPU gets its own maximum
# and prints the fastest CPUs, i.e. where to pin the benchmark (taskset).
# Settings a platform lacks (EPP/EPB on ARM, boost on most SBCs) are
# skipped; commands that fail are reported as warnings, and the closing
# status shows what actually took effect. restore applies the recorded
# values in dependency order and renames the state file to <state>.restored.
#
# Boost: frequencies above the non-boost maximum (amd-pstate "nominal",
# Intel base clock) are opportunistic - firmware lowers them with
# temperature, power and active-core count - so pinning min = max there is
# a request, not a guarantee. The default therefore disables boost and pins
# the highest guaranteed clock. KEEP_BOOST=1 keeps boost on and pins the
# boost maximum (or FREQ_KHZ) instead, for CPUs whose non-boost clock is
# unrepresentatively low (e.g. 2.0 GHz on mobile Zen 5); verify the clock
# actually held during a run, e.g. `sudo cpupower monitor -m Mperf`.
#
# Requires cpupower (Debian: linux-cpupower) >= 6.9 for boost/EPP control
# (`cpupower set --turbo-boost/--epp`). cpupower 6.12 (and older) neither reads
# nor writes intel_pstate's turbo switch, so on intel_pstate systems the
# script reads and writes intel_pstate/no_turbo directly.
#
# Environment:
#   STATE     state file (default /var/tmp/minisketch-bench-host.state);
#             setup refuses to overwrite an existing one (a second setup
#             would record the pinned state as "original")
#   FREQ_KHZ  cap the pinned frequency (kHz), e.g. 2000000 for a Pi 5
#             without active cooling; snapped down to an available step
#   KEEP_BOOST  if set, leave boost enabled and pin the boost maximum (see
#             above)
#   DRY_RUN   if set, print the cpupower commands setup/restore would run
#             instead of running them (no root needed)

set -euo pipefail

STATE="${STATE:-/var/tmp/minisketch-bench-host.state}"
FREQ_KHZ="${FREQ_KHZ:-}"
DRY_RUN="${DRY_RUN:-}"
KEEP_BOOST="${KEEP_BOOST:-}"
NO_TURBO=/sys/devices/system/cpu/intel_pstate/no_turbo
WARNINGS=0

declare -a CPUS=()
declare -A GOV=() MIN=() MAX=() EPP=() EPB=()
BOOST=""

warn() { echo "warning: $*" >&2; WARNINGS=$((WARNINGS + 1)); }

# Run a cpupower command that changes settings (printed only in DRY_RUN).
run() {
    local out
    if [ -n "$DRY_RUN" ]; then
        echo "+ cpupower $*"
    elif ! out=$(cpupower "$@" 2>&1); then
        warn "cpupower $* failed: ${out//$'\n'/ }"
    fi
}

# "605.264 MHz" / "5.157895 GHz" (cpupower -n output) -> kHz
to_khz() {
    awk -v v="$1" -v u="$2" 'BEGIN {
        f = (u == "GHz") ? 1000000 : (u == "MHz") ? 1000 : (u == "Hz") ? 0.001 : 1
        printf "%.0f\n", v * f }'
}

list_cpus() {
    mapfile -t CPUS < <(cpupower -c all frequency-info -l 2> /dev/null |
        awk '/^analyzing CPU/ { sub(":", "", $3); print $3 }')
    [ "${#CPUS[@]}" -gt 0 ] || { echo "error: cpupower reports no cpufreq CPUs" >&2; exit 1; }
}

# Read the current settings of every CPU into GOV/MIN/MAX/EPP/EPB, BOOST.
read_current() {
    local c info line lo lou hi hiu
    list_cpus
    for c in "${CPUS[@]}"; do
        info=$(cpupower -c "$c" frequency-info -n 2> /dev/null)
        GOV[$c]=$(sed -n 's/.*The governor "\([^"]*\)".*/\1/p' <<< "$info")
        line=$(grep -o 'should be within .*' <<< "$info" || true)
        if [ -n "$line" ]; then
            read -r _ _ _ lo lou _ hi hiu <<< "$line"
            MIN[$c]=$(to_khz "$lo" "$lou")
            MAX[$c]=$(to_khz "$hi" "${hiu%.}")
        fi
        EPP[$c]=$(sed -n 's/.*energy performance preference: *//p' <<< "$info")
        EPB[$c]=$(cpupower -c "$c" info -b 2> /dev/null | sed -n 's/^perf-bias: *\([0-9]*\).*/\1/p')
        if [ "$c" = "${CPUS[0]}" ]; then
            if [ -r "$NO_TURBO" ]; then
                # cpupower 6.12 reports Intel boost as always active.
                BOOST=$((1 - $(cat "$NO_TURBO")))
            elif grep -q 'Supported: yes' <<< "$info"; then
                BOOST=$(grep -o 'Active: [a-z]*' <<< "$info" | head -1 | sed 's/Active: yes/1/;s/Active: no/0/')
            else
                BOOST=""
            fi
        fi
    done
}

# Enable (1) or disable (0) boost. intel_pstate's switch is written
# directly: cpupower 6.12 only knows the generic cpufreq/boost file,
# which intel_pstate does not provide.
set_boost() {
    if [ -e "$NO_TURBO" ]; then
        if [ -n "$DRY_RUN" ]; then
            echo "+ echo $((1 - $1)) > $NO_TURBO"
        elif ! echo $((1 - $1)) 2> /dev/null > "$NO_TURBO"; then
            warn "could not write $NO_TURBO (turbo disabled in firmware?)"
        fi
    else
        run set --turbo-boost "$1"
    fi
}

# Highest frequency CPU $1 may be pinned to: the hardware maximum (drivers
# lower it once boost is off), no higher than the amd-pstate nominal
# frequency (unless KEEP_BOOST) or FREQ_KHZ, snapped down to an available
# frequency step.
target_freq() {
    local c="$1" info best cap step lim line v u
    info=$(cpupower -c "$c" frequency-info -n 2> /dev/null)
    read -r _ best < <(cpupower -c "$c" frequency-info -l 2> /dev/null | tail -1)
    line=$(grep -o 'Nominal Frequency: [0-9.]* [kMG]*Hz' <<< "$info" || true)
    if [ -n "$line" ] && [ -z "$KEEP_BOOST" ]; then
        read -r _ _ v u <<< "$line"
        cap=$(to_khz "$v" "$u")
        [ "$cap" -lt "$best" ] && best=$cap
    fi
    [ -n "$FREQ_KHZ" ] && [ "$FREQ_KHZ" -lt "$best" ] && best=$FREQ_KHZ
    lim=$best
    line=$(sed -n 's/.*available frequency steps: *//p' <<< "$info")
    if [ -n "$line" ]; then
        best=0
        while read -r v u; do
            step=$(to_khz "$v" "$u")
            [ "$step" -le "$lim" ] && [ "$step" -gt "$best" ] && best=$step
        done < <(tr ',' '\n' <<< "$line")
        [ "$best" -gt 0 ] || best=$lim
    fi
    echo "$best"
}

# Apply one setting to all CPUs, one cpupower call per distinct value.
#   apply <array-name> <cpupower args before the value...>
apply() {
    local -n vals="$1"
    shift
    local c v
    declare -A groups=()
    for c in "${CPUS[@]}"; do
        v="${vals[$c]:-}"
        [ -n "$v" ] && groups[$v]+=",$c"
    done
    for v in "${!groups[@]}"; do
        run -c "${groups[$v]#,}" "$@" "$v"
    done
}

# One line per distinct setting combination, with the CPUs that share it.
status() {
    local c key
    read_current
    declare -A lines=()
    declare -a order=()
    for c in "${CPUS[@]}"; do
        key="governor ${GOV[$c]:-?}, min ${MIN[$c]:-?} kHz, max ${MAX[$c]:-?} kHz"
        [ -n "${EPP[$c]}" ] && key+=", EPP ${EPP[$c]}"
        [ -n "${EPB[$c]}" ] && key+=", EPB ${EPB[$c]}"
        [ -n "${lines[$key]+x}" ] || order+=("$key")
        lines[$key]+=",$c"
    done
    echo "boost: ${BOOST:-not supported}"
    for key in "${order[@]}"; do
        echo "cpus ${lines[$key]#,}: $key"
    done
    if command -v vcgencmd > /dev/null; then
        echo "vcgencmd: $(vcgencmd get_throttled) $(vcgencmd measure_temp)"
    fi
    return 0
}

need_root() {
    if [ -z "$DRY_RUN" ] && [ "$(id -u)" -ne 0 ]; then
        echo "error: '$1' needs root; run: sudo $0 $1 (or DRY_RUN=1 to preview)" >&2
        exit 1
    fi
}

setup() {
    need_root setup
    # A dry run records to a side file, so it never blocks the real setup.
    local out="$STATE"
    [ -n "$DRY_RUN" ] && out="$STATE.dry-run"
    if [ -e "$STATE" ]; then
        echo "error: $STATE exists (setup already ran?). Run 'restore' first," >&2
        echo "or delete it if it is stale (e.g. after a reboot)." >&2
        exit 1
    fi

    local c t
    read_current
    {
        echo "# minisketch bench_host.sh state, $(uname -n), $(date -Is)"
        echo "# restore with: sudo $0 restore"
        echo "boost=$BOOST"
        for c in "${CPUS[@]}"; do
            echo "cpu$c governor=${GOV[$c]} min=${MIN[$c]} max=${MAX[$c]} epp=${EPP[$c]} epb=${EPB[$c]}"
        done
        echo "#"
        echo "# --- cpupower frequency-info (reference only) ---"
        cpupower -c all frequency-info 2>&1 | sed 's/^/# /'
        echo "# --- cpupower info (reference only) ---"
        cpupower -c all info 2>&1 | sed 's/^/# /'
    } > "$out"
    echo "Recorded current settings to $out"

    # EPP/EPB before the governor: the pstate drivers reject EPP changes
    # while the governor is "performance".
    local -A perf=() zero=()
    for c in "${CPUS[@]}"; do
        [ -n "${EPP[$c]}" ] && perf[$c]=performance
        [ -n "${EPB[$c]}" ] && zero[$c]=0
    done
    apply perf set -e
    apply zero set -b
    if [ -n "$KEEP_BOOST" ]; then
        [ "$BOOST" = 0 ] && set_boost 1
    elif [ -n "$BOOST" ]; then
        set_boost 0
    fi

    local -A gov=()
    for c in "${CPUS[@]}"; do gov[$c]=performance; done
    apply gov frequency-set -g

    # Pin after boost is off: the hardware maximum is then non-boost. Each
    # CPU gets its own maximum; one cpupower call per distinct frequency
    # (i.e. per cluster / core type).
    local -A pins=()
    local top=0
    for c in "${CPUS[@]}"; do
        t=$(target_freq "$c")
        pins[$t]+=",$c"
        [ "$t" -gt "$top" ] && top=$t
    done
    for t in "${!pins[@]}"; do
        run -c "${pins[$t]#,}" frequency-set -d "$t" -u "$t"
    done
    # Suggest the fastest CPU: highest pinned clock, ties broken by the
    # recorded pre-setup maximum (the core type's boost clock, so a hybrid
    # CPU's big cores win even when all pin to the same nominal clock),
    # then the highest number (CPU 0 takes more interrupts).
    local best="" bmax=0 fastest=""
    for c in ${pins[$top]//,/ }; do
        if [ "${MAX[$c]:-0}" -ge "$bmax" ]; then
            [ "${MAX[$c]:-0}" -gt "$bmax" ] && fastest=""
            bmax="${MAX[$c]:-0}" best="$c" fastest+=",$c"
        fi
    done
    PIN_HINT="fastest CPUs (pinned at $top kHz): ${fastest#,} -> BENCH_PREFIX=\"taskset -c $best\""
    [ -n "$KEEP_BOOST" ] && warn "boost kept on: clocks above the non-boost maximum are not guaranteed to hold"

    echo
    echo "$PIN_HINT"
    [ -n "$DRY_RUN" ] && { echo "(dry run: nothing changed; preview state in $out)"; return 0; }
    status
    echo
    [ "$WARNINGS" -eq 0 ] && echo "Host pinned." ||
        echo "Host pinned with $WARNINGS warning(s); check the status above."
    echo "After benchmarking: sudo $0 restore"
}

check() { [ "$2" = "$3" ] || warn "$1 is '$2', recorded '$3'"; }

restore() {
    need_root restore
    [ -e "$STATE" ] || { echo "error: no state file at $STATE" >&2; exit 1; }

    local line name kv k v c
    local -A sgov=() smin=() smax=() sepp=() sepb=()
    local sboost=""
    CPUS=()
    while IFS= read -r line; do
        case "$line" in
            '#'* | '') continue ;;
            boost=*) sboost="${line#boost=}" ;;
            cpu*)
                read -r name kv <<< "$line"
                c="${name#cpu}"
                CPUS+=("$c")
                for k in $kv; do
                    v="${k#*=}"
                    case "$k" in
                        governor=*) sgov[$c]="$v" ;;
                        min=*) smin[$c]="$v" ;;
                        max=*) smax[$c]="$v" ;;
                        epp=*) sepp[$c]="$v" ;;
                        epb=*) sepb[$c]="$v" ;;
                    esac
                done
                ;;
        esac
    done < "$STATE"

    # Boost first (it bounds the frequency limits), then the governor (EPP
    # is locked while it is "performance"), then EPP/EPB, then the limits;
    # cpupower orders the min/max writes itself.
    [ -n "$sboost" ] && set_boost "$sboost"
    apply sgov frequency-set -g
    apply sepp set -e
    apply sepb set -b
    local -A limits=()
    for c in "${CPUS[@]}"; do
        if [ -n "${smin[$c]:-}" ] && [ -n "${smax[$c]:-}" ]; then
            limits["${smin[$c]} ${smax[$c]}"]+=",$c"
        fi
    done
    for k in "${!limits[@]}"; do
        run -c "${limits[$k]#,}" frequency-set -d "${k% *}" -u "${k#* }"
    done

    [ -n "$DRY_RUN" ] && { echo "(dry run: nothing changed)"; return 0; }

    read_current
    [ -z "$sboost" ] || [ "$BOOST" = "$sboost" ] || warn "boost is $BOOST, recorded $sboost"
    for c in "${CPUS[@]}"; do
        check "cpu$c governor" "${GOV[$c]:-}" "${sgov[$c]:-}"
        check "cpu$c min" "${MIN[$c]:-}" "${smin[$c]:-}"
        check "cpu$c max" "${MAX[$c]:-}" "${smax[$c]:-}"
        check "cpu$c EPP" "${EPP[$c]:-}" "${sepp[$c]:-}"
        check "cpu$c EPB" "${EPB[$c]:-}" "${sepb[$c]:-}"
    done
    mv "$STATE" "$STATE.restored"
    status
    echo
    [ "$WARNINGS" -eq 0 ] && echo "Restored (state kept as $STATE.restored)." ||
        echo "Restored with $WARNINGS warning(s) (state kept as $STATE.restored)."
}

command -v cpupower > /dev/null || { echo "error: cpupower not found (Debian: apt install linux-cpupower)" >&2; exit 1; }
case "${1:-}" in
    setup) setup ;;
    restore) restore ;;
    status) status ;;
    *) sed -n '3,7p' "$0" | sed 's/^# \{0,1\}//'; exit 1 ;;
esac
