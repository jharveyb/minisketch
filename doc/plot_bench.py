#!/usr/bin/env python3
"""Plot minisketch benchmark results from TSV produced by src/bench.

Dependencies: pip install pandas matplotlib

Generating the TSV data (pass one or more of these to the script):

    # Most plots -- sweep mode
    ./build/src/bench --sweep --sweep-errors --iters 10 > bench_sweep.tsv

    # Plot #1 (diff vs decode at fixed capacity) -- shell loop
    for c in 1024 2048 4096; do
      for e in 16 32 64 128 256 512 1024 2048 4096; do
        ./build/src/bench --syndromes $c --errors $e --iters 10
      done
    done > bench_diffs.tsv

    # Plot #4 (capacity vs decode at fixed errors) -- one run per errors value
    for e in 256 512 1024; do
      ./build/src/bench --sweep --errors $e --iters 10
    done > bench_fixed.tsv

Plotting:

    python3 doc/plot_bench.py bench_sweep.tsv bench_diffs.tsv bench_fixed.tsv \\
        --output-dir doc/generated/ \\
        --fixed-capacity 1024,2048,4096 --fixed-errors 256,512,1024
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, LogLocator, NullLocator, ScalarFormatter


IMPL_MARKERS = {
    "GENERIC": "o",
    "CLMUL": "x",
    "CLMUL_TRI": "+",
}

NUMERIC_COLS = ["bits", "capacity", "errors", "data_len", "value"]


def load_tsv(paths: list[Path]) -> pd.DataFrame:
    """Load and concatenate one or more bench TSVs, stripping repeated headers."""
    frames = [pd.read_csv(p, sep="\t") for p in paths]
    df = pd.concat(frames, ignore_index=True)
    df = df[df["metric"] != "metric"]
    for col in NUMERIC_COLS:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    return df.dropna(subset=["value"]).reset_index(drop=True)


def _format_time_ms(ms: float, _pos: int | None = None) -> str:
    """Render a millisecond value using the largest SI prefix that keeps
    it >= 1 (e.g. 0.001 -> '1 us', 1000 -> '1 s', 100 -> '100 ms')."""
    if ms <= 0:
        return ""
    seconds = ms * 1e-3
    for scale, unit in [(1.0, "s"), (1e-3, "ms"), (1e-6, "us"), (1e-9, "ns")]:
        if seconds >= scale * 0.999999:
            return f"{seconds / scale:g} {unit}"
    return f"{seconds * 1e12:g} ps"


def _style_loglog(ax, xlabel: str, ylabel: str, title: str) -> None:
    """X axis: log2 with plain power-of-2 labels (matches doc/plot_capacity.png
    and doc/plot_diff.png). Y axis: log10 with one evenly-spaced label per
    decade, rendered as SI-prefixed time (us / ms / s)."""
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.xaxis.set_major_formatter(ScalarFormatter())
    ax.yaxis.set_major_locator(LogLocator(base=10.0, subs=(1.0, 2.5, 5.0, 7.5), numticks=40))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_major_formatter(FuncFormatter(_format_time_ms))
    ax.grid(True, which="major", linestyle="--", alpha=0.4)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)


def _save(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    print(f"wrote {path}")
    plt.close(fig)


def _skip(name: str, reason: str) -> None:
    print(f"skipping {name}: {reason}")


def _color_palette(values: list) -> dict:
    """Map ordered values to colors from a sequential colormap. With one
    value, return a single mid-range color."""
    n = len(values)
    cmap = mpl.colormaps["viridis"]
    if n <= 1:
        return {v: cmap(0.5) for v in values}
    return {v: cmap(i / (n - 1) * 0.85) for i, v in enumerate(values)}


def _plot_grouped(sub: pd.DataFrame, *, x_col: str,
                  color_col: str, color_label: str,
                  xlabel: str, ylabel: str, title: str,
                  out: Path) -> None:
    """One line per (color_col, implementation) pair. Color encodes
    color_col, marker encodes implementation. The legend has two sections
    (color values + implementations) unless color_col == 'implementation',
    in which case both encodings collapse and a single legend is drawn."""
    color_values = sorted(sub[color_col].unique())
    impls = sorted(sub["implementation"].unique())
    palette = _color_palette(color_values)

    fig, ax = plt.subplots(figsize=(8, 5.5))
    for cv in color_values:
        for impl in impls:
            g = sub[(sub[color_col] == cv) & (sub["implementation"] == impl)]
            if g.empty:
                continue
            g = g.sort_values(x_col)
            marker = IMPL_MARKERS.get(impl, "o")
            ax.plot(g[x_col], g["value"],
                    color=palette[cv],
                    marker=marker,
                    markersize=7, markeredgewidth=1.5,
                    linewidth=1.4)

    _style_loglog(ax, xlabel, ylabel, title)
    _draw_legend(ax, color_values, color_label, palette, impls)
    _save(fig, out)


def _draw_legend(ax, color_values, color_label, palette, impls):
    if color_label == "implementation":
        # Color already encodes impl; merge marker + color into one handle per impl.
        ax.legend(handles=[Line2D([0], [0], color=palette[v],
                                  marker=IMPL_MARKERS.get(v, "o"),
                                  linewidth=2, markersize=7,
                                  markeredgewidth=1.5, label=str(v))
                           for v in color_values],
                  loc="best")
        return
    color_handles = [Line2D([0], [0], color=palette[v], linewidth=2,
                            label=f"{color_label}={int(v)}")
                     for v in color_values]
    impl_handles = [Line2D([0], [0], marker=IMPL_MARKERS.get(i, "o"),
                           color="black", linestyle="none",
                           markersize=8, markeredgewidth=1.5, label=i)
                    for i in impls]
    ax.legend(handles=color_handles + impl_handles, loc="best")


def _auto_pick(df: pd.DataFrame, col: str) -> int | None:
    """Pick the largest value of `col` with >=2 rows per implementation."""
    if df.empty:
        return None
    by_val = df.groupby([col, "implementation"]).size().unstack(fill_value=0)
    good = by_val[(by_val >= 2).all(axis=1)].index
    if len(good) == 0:
        return None
    return int(max(good))


def _resolve_list(values: list[int] | None, sub: pd.DataFrame,
                  col: str) -> list[int] | None:
    """Honor an explicit list, otherwise fall back to a single auto-picked value."""
    if values:
        return values
    picked = _auto_pick(sub, col)
    return [picked] if picked is not None else None


def plot_diff_vs_decode(df: pd.DataFrame, *, capacities: list[int] | None,
                        bits: int, out: Path,
                        title_suffix: str = "") -> None:
    """Plot #1: decode time vs differences.

    One line per (capacity, implementation): color = capacity, marker = impl."""
    sub = df[(df["metric"] == "recover[ms]") & (df["bits"] == bits)]
    capacities = _resolve_list(capacities, sub, "capacity")
    if capacities is None:
        _skip(out.name,
              f"no recover[ms] rows at bits={bits} with multiple errors per capacity")
        return
    sub = sub[sub["capacity"].isin(capacities)]
    if sub.empty or sub["errors"].nunique() < 2:
        _skip(out.name,
              f"need recover[ms] at bits={bits}, capacity in {capacities} with varying errors "
              f"(try: for c in {' '.join(map(str, capacities))}; do "
              f"for e in 16 32 64 128 256; do ./bench --syndromes $c --errors $e; done; done)")
        return
    _plot_grouped(sub, x_col="errors",
                  color_col="capacity", color_label="capacity",
                  xlabel="Differences", ylabel="Decode time",
                  title=f"Decode time vs differences (bits={bits}){title_suffix}",
                  out=out)


def plot_capacity_vs_decode_50pct(df: pd.DataFrame, *, bits: int,
                                  out: Path,
                                  title_suffix: str = "") -> None:
    """Plot #2: decode time vs capacity, at errors = capacity/2.

    Single dimension; color and marker both encode implementation."""
    sub = df[(df["metric"] == "recover[ms]") & (df["bits"] == bits)].copy()
    sub = sub[sub["errors"] * 2 == sub["capacity"]]
    if sub.empty or sub["capacity"].nunique() < 2:
        _skip(out.name,
              f"need recover[ms] at bits={bits} with errors == capacity/2 "
              f"(try: ./bench --sweep --sweep-errors)")
        return
    _plot_grouped(sub, x_col="capacity",
                  color_col="implementation", color_label="implementation",
                  xlabel="Capacity", ylabel="Decode time",
                  title=f"Decode time vs capacity (bits={bits}, 50% differences){title_suffix}",
                  out=out)


def plot_capacity_vs_create(df: pd.DataFrame, *, bits: int,
                            errors: list[int] | None, out: Path,
                            title_suffix: str = "") -> None:
    """Plot #3: total creation time vs capacity, at one or more fixed errors values.

    One line per (errors, implementation): color = errors, marker = impl."""
    sub = df[(df["metric"] == "create-total[ms]") & (df["bits"] == bits)]
    errors = _resolve_list(errors, sub, "errors")
    if errors is None:
        _skip(out.name,
              f"no create-total[ms] rows at bits={bits} with multiple capacities per errors")
        return
    sub = sub[sub["errors"].isin(errors)]
    if sub.empty: # or sub["capacity"].nunique() < 2:
        _skip(out.name,
              f"need create-total[ms] at bits={bits}, errors in {errors} across capacities "
              f"(try: for e in {' '.join(map(str, errors))}; do ./bench --sweep --errors $e; done)")
        return
    _plot_grouped(sub, x_col="capacity",
                  color_col="errors", color_label="errors",
                  xlabel="Capacity", ylabel="Creation time",
                  title=f"Creation time vs capacity (bits={bits}){title_suffix}",
                  out=out)


def plot_capacity_vs_decode_fixed(df: pd.DataFrame, *, bits: int,
                                  errors: list[int] | None, out: Path,
                                  title_suffix: str = "") -> None:
    """Plot #4: decode time vs capacity, at one or more fixed errors values.

    One line per (errors, implementation): color = errors, marker = impl."""
    sub = df[(df["metric"] == "recover[ms]") & (df["bits"] == bits)]
    errors = _resolve_list(errors, sub, "errors")
    if errors is None:
        _skip(out.name,
              f"no recover[ms] rows at bits={bits} with multiple capacities per errors")
        return
    sub = sub[sub["errors"].isin(errors)]
    if sub.empty or sub["capacity"].nunique() < 2:
        _skip(out.name,
              f"need recover[ms] at bits={bits}, errors in {errors} across capacities "
              f"(try: for e in {' '.join(map(str, errors))}; do ./bench --sweep --errors $e; done)")
        return
    _plot_grouped(sub, x_col="capacity",
                  color_col="errors", color_label="errors",
                  xlabel="Capacity", ylabel="Decode time",
                  title=f"Decode time vs capacity (bits={bits}){title_suffix}",
                  out=out)


# Load the benchmark output and compute various metrics, without plotting.
def compute_stats(df: pd.DataFrame) -> pd.DataFrame:
    generic_times = df[(df["metric"] == "recover[ms]") & (df["implementation"] == "GENERIC")]["value"]
    generic_times.index += 1
    clmul_times = df[(df["metric"] == "recover[ms]") & (df["implementation"] == "CLMUL")]["value"]
    clmul_speedup = generic_times.div(clmul_times)
    print(f"Average speedup from CLMUL: {clmul_speedup.mean():.3f}x")

    generic_create_op_time = df[(df["metric"] == "create[ns]") & (df["implementation"] == "GENERIC")]["value"]
    generic_create_op_time.index += 2
    clmul_create_op_time = df[(df["metric"] == "create[ns]") & (df["implementation"] == "CLMUL")]["value"]
    clmul_difference = generic_create_op_time.div(clmul_create_op_time)
    print(f"Average timing change from CLMUL for create[ns] (adding one element to sketch): {clmul_difference.mean():.3f}x")
    return clmul_speedup


PLOT_REGISTRY = {
    "diff": ("diff.png",
             lambda df, args, out: plot_diff_vs_decode(
                 df, capacities=args.fixed_capacity, bits=args.bits, out=out,
                 title_suffix=args.title_suffix)),
    "cap50": ("cap50.png",
              lambda df, args, out: plot_capacity_vs_decode_50pct(
                  df, bits=args.bits, out=out,
                  title_suffix=args.title_suffix)),
    "create": ("create.png",
               lambda df, args, out: plot_capacity_vs_create(
                   df, bits=args.bits, errors=args.fixed_errors, out=out,
                   title_suffix=args.title_suffix)),
    "cap-fixed": ("cap_fixed.png",
                  lambda df, args, out: plot_capacity_vs_decode_fixed(
                      df, bits=args.bits, errors=args.fixed_errors, out=out,
                      title_suffix=args.title_suffix)),
}


def _parse_int_list(s: str) -> list[int]:
    """Parse a comma-separated list of ints (e.g. '1024,2048,4096')."""
    return [int(v.strip()) for v in s.split(",") if v.strip()]


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot minisketch benchmark TSV output.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("inputs", nargs="+", type=Path,
                        help="one or more TSV files from src/bench")
    parser.add_argument("--output-dir", type=Path, default=Path("."),
                        help="directory for PNG output (default: .)")
    parser.add_argument("--bits", type=int, default=64,
                        help="field size to filter on (default: 64)")
    parser.add_argument("--fixed-capacity", type=_parse_int_list, default=None,
                        metavar="C1,C2,...",
                        help="comma-separated capacities for plot #1 "
                             "(default: auto-pick the largest one available)")
    parser.add_argument("--fixed-errors", type=_parse_int_list, default=None,
                        metavar="E1,E2,...",
                        help="comma-separated errors values for plots #3 and #4 "
                             "(default: auto-pick the largest one available)")
    parser.add_argument("--plots", nargs="+",
                        choices=list(PLOT_REGISTRY.keys()) + ["all","stats"],
                        default=["all"],
                        help="which plots to render (default: all)")
    parser.add_argument("--title-suffix", default="",
                        metavar="TEXT",
                        help="text appended to every plot title "
                             "(e.g. ', Cortex A72 @ 1.8 GHz')")
    args = parser.parse_args()

    df = load_tsv(args.inputs)
    print(f"loaded {len(df)} rows from {len(args.inputs)} file(s)")

    if "stats" in args.plots:
        df = compute_stats(df)
        print("computed stats")
        print(df)
        return

    selected = PLOT_REGISTRY.keys() if "all" in args.plots else args.plots
    for key in selected:
        filename, fn = PLOT_REGISTRY[key]
        fn(df, args, args.output_dir / filename)


if __name__ == "__main__":
    main()
