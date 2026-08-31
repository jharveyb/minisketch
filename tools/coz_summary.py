#!/usr/bin/env python3
"""Summarize a coz profile.jsonl: stage time shares + top causal lines."""
#
# Usage: tools/coz_summary.py <profile.jsonl>
#
# Companion to tools/coz.sh (which runs it automatically): prints the
# decode-stage time shares (from the MINISKETCH_COZ_BEGIN/END latency pairs;
# each experiment contributes one active/inactive snapshot per stage, so
# shares have a granularity of about 1/#experiments) and the top lines by
# measured causal program speedup. Note that -O2 inlining collapses much of
# decode into few lines (the minisketch_decode API line "winning" means the
# whole decode, i.e. nothing) - the stage markers are the reliable axis;
# line-level entries inside sketch_impl.h refine them. Entries backed by a
# handful of experiments are noise; run longer for line-level detail
# (roughly 15 experiments accumulate per minute).
import json, sys, collections

path = sys.argv[1]
experiments = []  # (selected, speedup, duration, {name: record})
cur = None
for line in open(path):
    try:
        rec = json.loads(line)
    except json.JSONDecodeError:
        continue
    t = rec.get("type")
    if t == "experiment":
        cur = {"selected": rec["selected"], "speedup": rec["speedup"],
               "duration": rec["duration"], "points": {}}
        experiments.append(cur)
    elif t in ("latency_point", "throughput_point", "progress_point") and cur is not None:
        cur["points"][rec["name"]] = rec

# Stage shares: duration-weighted mean of `difference` over speedup-0
# experiments. For a stage that runs once per decode, the sampled
# difference (arrivals - departures, 0 or 1) equals "stage active now",
# so its mean is the fraction of wall time spent in the stage.
share_num = collections.defaultdict(float)
share_den = 0.0
for e in experiments:
    if e["speedup"] != 0:
        continue
    share_den += e["duration"]
    for name, rec in e["points"].items():
        if "difference" in rec:
            share_num[name] += rec["difference"] * e["duration"]
print(f"# {path}")
print(f"experiments: {len(experiments)} "
      f"(baseline: {sum(1 for e in experiments if e['speedup'] == 0)})")
if share_den:
    print("\nStage time shares (fraction of wall time inside each begin/end pair):")
    for name, num in sorted(share_num.items(), key=lambda kv: -kv[1]):
        print(f"  {name:28s} {100 * num / share_den:5.1f}%")

# Causal chart: program speedup per selected line, from the `decode`
# throughput point. Baseline period per line from its speedup-0 runs
# (falling back to the global baseline period).
def period(e):
    rec = e["points"].get("decode")
    if rec and rec.get("delta"):
        return e["duration"] / rec["delta"]
    return None

glob_base = [p for e in experiments if e["speedup"] == 0 if (p := period(e))]
glob_base = sum(glob_base) / len(glob_base) if glob_base else None
by_line = collections.defaultdict(list)
for e in experiments:
    by_line[e["selected"]].append(e)

rows = []
for sel, exps in by_line.items():
    base = [p for e in exps if e["speedup"] == 0 if (p := period(e))]
    base = sum(base) / len(base) if base else glob_base
    if not base:
        continue
    pts = collections.defaultdict(list)
    for e in exps:
        if e["speedup"] > 0 and (p := period(e)):
            pts[e["speedup"]].append(p)
    if not pts:
        continue
    # Report the mean program speedup at the largest measured line speedup.
    s = max(pts)
    ps = [1 - p / base for p in pts[s]]
    rows.append((sum(ps) / len(ps), s, len(exps), sel))

print("\nTop lines by measured program speedup"
      " (at their largest tested virtual speedup):")
print(f"  {'prog speedup':>12s} {'@line speedup':>13s} {'#exp':>5s}  line")
for ps, s, n, sel in sorted(rows, reverse=True)[:15]:
    sel = sel.replace("/home/jhb/2025/plebfi/minisketch/", "")
    print(f"  {100 * ps:11.1f}% {100 * s:12.0f}% {n:5d}  {sel}")
