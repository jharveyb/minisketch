#!/usr/bin/env bash
# Start, inspect and stop long benchmark jobs on a bench host without ever
# killing by name. Run it on the host (e.g. `ssh host 'cd ~/dir &&
# tools/bench_remote.sh start NAME -- ./driver.sh'`).
#
#   tools/bench_remote.sh start NAME -- CMD [ARGS...]
#       Runs CMD detached in its own session (setsid), so its process group
#       is its own. Records the group id and the leader's start time in
#       ~/.bench-jobs/NAME.job and its output in ~/.bench-jobs/NAME.log.
#       Refuses if a job NAME is still running.
#   tools/bench_remote.sh status [NAME]
#       Shows whether the job(s) are running.
#   tools/bench_remote.sh log NAME [LINES]
#       Tails the job's log (default 40 lines).
#   tools/bench_remote.sh stop NAME
#       Sends TERM to the recorded process group, waits up to 30 s, then
#       KILL to that group only if it is still there. Before signalling, it
#       checks that the group leader is still the recorded process (same
#       start time, owned by this user), so a reused PID is never hit.
#
# Never stop jobs with pkill/killall/pgrep patterns: `pgrep -f` matches the
# shell running it, and over Tailscale SSH that shell's process group belongs
# to tailscaled. A job's own EXIT trap (e.g. `bench_host.sh restore`) still
# runs on TERM.

set -euo pipefail

JOBDIR="${HOME}/.bench-jobs"
mkdir -p "$JOBDIR"

die() { echo "bench_remote.sh: $*" >&2; exit 1; }

usage() { sed -n '2,24p' "$0" | sed 's/^# \{0,1\}//'; exit 1; }

valid_name() { [[ "$1" =~ ^[A-Za-z0-9._-]+$ ]] || die "invalid job name '$1'"; }

# The start time (clock ticks since boot) of a process, or empty if gone.
start_time() { awk '{print $22}' "/proc/$1/stat" 2>/dev/null || true; }

# Is the recorded job still running: leader alive, same start time, ours?
job_alive() {
    local pgid="$1" stime="$2"
    [[ -n "$pgid" && "$pgid" -gt 1 ]] || return 1
    [[ "$(start_time "$pgid")" == "$stime" ]] || return 1
    [[ "$(stat -c %u "/proc/$pgid" 2>/dev/null)" == "$(id -u)" ]] || return 1
}

# Any process left in the group (the leader may have exited before its children)?
group_has_members() { pgrep -g "$1" >/dev/null 2>&1; }

load_job() {
    local f="$JOBDIR/$1.job"
    [[ -f "$f" ]] || die "no job '$1' (no $f)"
    # shellcheck disable=SC1090
    source "$f"
}

cmd_start() {
    local name="$1"; shift
    valid_name "$name"
    [[ "${1:-}" == "--" ]] || die "usage: start NAME -- CMD [ARGS...]"
    shift
    [[ $# -gt 0 ]] || die "no command given"
    if [[ -f "$JOBDIR/$name.job" ]]; then
        load_job "$name"
        if job_alive "$PGID" "$STIME"; then die "job '$name' is still running (pgid $PGID)"; fi
    fi
    local log="$JOBDIR/$name.log"
    setsid nohup "$@" > "$log" 2>&1 < /dev/null &
    local pid=$!
    # setsid makes the child a session and group leader: its PID is the PGID.
    sleep 0.2
    local stime; stime="$(start_time "$pid")"
    [[ -n "$stime" ]] || { echo "job exited immediately; log:" >&2; tail -n 20 "$log" >&2; exit 1; }
    {
        printf 'PGID=%q\n' "$pid"
        printf 'STIME=%q\n' "$stime"
        printf 'STARTED=%q\n' "$(date -Is)"
        printf 'CMD=%q\n' "$*"
        printf 'CWD=%q\n' "$PWD"
    } > "$JOBDIR/$name.job"
    echo "started '$name': pgid $pid, log $log"
}

cmd_status() {
    local f name
    if [[ $# -gt 0 ]]; then set -- "$JOBDIR/$1.job"; else set -- "$JOBDIR"/*.job; fi
    for f in "$@"; do
        [[ -f "$f" ]] || { echo "no jobs"; return; }
        name="$(basename "$f" .job)"
        # shellcheck disable=SC1090
        source "$f"
        if job_alive "$PGID" "$STIME"; then
            echo "$name: running (pgid $PGID, since $STARTED): $CMD"
        elif group_has_members "$PGID" && [[ -z "$(start_time "$PGID")" ]]; then
            echo "$name: leader exited, children of group $PGID still running: $CMD"
        else
            echo "$name: finished (started $STARTED): $CMD"
        fi
    done
}

cmd_log() {
    valid_name "$1"
    tail -n "${2:-40}" "$JOBDIR/$1.log"
}

cmd_stop() {
    local name="$1"
    valid_name "$name"
    load_job "$name"
    if ! job_alive "$PGID" "$STIME"; then
        # The leader is gone; only signal the group if its PID has not been
        # reused by an unrelated process (then start_time would be non-empty).
        if [[ -z "$(start_time "$PGID")" ]] && group_has_members "$PGID"; then
            echo "leader of '$name' exited; signalling remaining members of group $PGID"
        else
            echo "job '$name' is not running"; return
        fi
    fi
    ps -o pid,pgid,user,etime,args -g "$PGID" || true
    kill -TERM -- "-$PGID" 2>/dev/null || true
    for _ in $(seq 30); do
        group_has_members "$PGID" || { echo "stopped '$name'"; return; }
        sleep 1
    done
    echo "group $PGID still running after TERM; sending KILL to that group only"
    kill -KILL -- "-$PGID" 2>/dev/null || true
}

[[ $# -ge 1 ]] || usage
cmd="$1"; shift
case "$cmd" in
    start)  [[ $# -ge 1 ]] || usage; cmd_start "$@" ;;
    status) cmd_status "$@" ;;
    log)    [[ $# -ge 1 ]] || usage; cmd_log "$@" ;;
    stop)   [[ $# -eq 1 ]] || usage; cmd_stop "$1" ;;
    *)      usage ;;
esac
