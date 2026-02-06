#!/usr/bin/env python3
"""Grid of cost-vs-speedup scatter plots, one subplot per script.

Usage:
    python plot_results.py <path/to/results.csv>

Produces grid_plot.png in the same directory as the input CSV.
"""

import sys
import os
import math

import pandas as pd
import matplotlib.pyplot as plt

MODE_STYLE = {
    "s3optdynamic": {"color": "#3498db", "marker": "o", "label": "s3optdynamic"},
    "single_shot": {"color": "#e67e22", "marker": "D", "label": "single_shot"},
    # Current runner mode names (post-refactor)
    "s3approxdynamic": {"color": "#3498db", "marker": "o", "label": "s3approxdynamic"},
    "s3approxsingleshot": {
        "color": "#e67e22",
        "marker": "D",
        "label": "s3approxsingleshot",
    },
}


def load(csv_path: str) -> pd.DataFrame:
    """Read CSV, skipping bare section-header lines (fewer fields than the header)."""
    df = pd.read_csv(csv_path, on_bad_lines="skip")
    # The header has 11 columns; any row that pandas couldn't fully parse is already
    # skipped.  Additionally drop rows where 'script' ended up NaN (safety net).
    df = df.dropna(subset=["script"])
    # Ensure numeric columns are actually numeric (section headers that slipped
    # through would make them object dtype).
    for col in ("wall_time_sec", "cost_usd"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    # Keep rows even when cost is unavailable (e.g. logs skipped); wall time is required.
    df = df.dropna(subset=["wall_time_sec"])
    return df


def main():
    if len(sys.argv) < 2:
        print("Usage: python plot_results.py <path/to/results.csv>", file=sys.stderr)
        sys.exit(1)

    csv_path = sys.argv[1]
    df = load(csv_path)

    # --------------- per-script baseline lookup ---------------
    # Group key is (script, input, width) so that if a CSV ever contains
    # multiple input sizes we still match the right baseline.
    baseline = (
        df[df["mode"] == "noopt"]
        .set_index(["script", "input", "width"])["wall_time_sec"]
    )
    baseline_cost = (
        df[df["mode"] == "noopt"]
        .set_index(["script", "input", "width"])["cost_usd"]
    )

    # Keep only the two optimised modes we want to plot.
    plot_df = df[df["mode"].isin(MODE_STYLE)].copy()
    if plot_df.empty:
        modes = ", ".join(sorted(df["mode"].dropna().astype(str).unique()))
        print(
            f"No plottable rows found for configured modes in {csv_path}. "
            f"Available modes: {modes}",
            file=sys.stderr,
        )
        sys.exit(1)

    # Attach the noopt wall_time for each row's group.
    plot_df["noopt_wall_time"] = plot_df.apply(
        lambda r: baseline.get((r["script"], r["input"], r["width"])), axis=1
    )
    plot_df["noopt_cost"] = plot_df.apply(
        lambda r: baseline_cost.get((r["script"], r["input"], r["width"])), axis=1
    )
    plot_df = plot_df.dropna(subset=["noopt_wall_time"])
    if plot_df.empty:
        print(
            f"No plottable rows remain after baseline matching in {csv_path}.",
            file=sys.stderr,
        )
        sys.exit(1)
    plot_df["speedup"] = plot_df["noopt_wall_time"] / plot_df["wall_time_sec"]

    # If cost is unavailable, fall back to wall-time scatter.
    if plot_df["cost_usd"].notna().any():
        x_col = "cost_usd"
        x_label = "Cost (USD)"
    else:
        x_col = "wall_time_sec"
        x_label = "Wall time (s)"

    # --------------- layout ---------------
    scripts = sorted(plot_df["script"].unique(), key=_script_sort_key)
    n = len(scripts)
    cols = math.ceil(math.sqrt(n))
    rows = math.ceil(n / cols)

    fig, axes = plt.subplots(rows, cols, figsize=(cols * 4.5, rows * 3.5))
    # Normalise axes to a flat list regardless of grid shape.
    if n == 1:
        axes = [axes]
    else:
        axes = list(axes.flat) if hasattr(axes, "flat") else axes

    handles = []  # collected once for the shared legend

    for idx, script in enumerate(scripts):
        ax = axes[idx]
        sub = plot_df[plot_df["script"] == script]

        # Baseline wall time and cost for title (take the first group's value).
        baseline_t = sub["noopt_wall_time"].iloc[0]
        baseline_cost_val = sub["noopt_cost"].iloc[0]
        if pd.notna(baseline_cost_val):
            subtitle = f"noopt: {baseline_t:.1f} s, ${baseline_cost_val:.4f}"
        else:
            subtitle = f"noopt: {baseline_t:.1f} s, cost N/A"
        ax.set_title(f"{script}  ({subtitle})")

        # Reference line ─ speedup == 1 means no improvement.
        ax.axhline(1.0, color="gray", linestyle="--", linewidth=0.9, alpha=0.7)

        for mode, style in MODE_STYLE.items():
            row = sub[sub["mode"] == mode]
            if row.empty:
                continue
            sc = ax.scatter(
                row[x_col],
                row["speedup"],
                c=style["color"],
                marker=style["marker"],
                s=100,
                edgecolors="black",
                linewidths=0.8,
                zorder=3,
                label=style["label"],
            )
            # Collect legend handles only on the first subplot.
            if idx == 0:
                handles.append(sc)

        ax.set_xlabel(x_label)
        ax.set_ylabel("Speedup vs noopt")
        ax.grid(True, alpha=0.3)

    # Hide any unused axes in the grid.
    for idx in range(n, len(axes)):
        axes[idx].set_visible(False)

    # Shared legend at the bottom centre.
    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=max(1, len(handles)),
        fontsize=11,
        frameon=True,
        bbox_to_anchor=(0.5, -0.02),
    )

    fig.suptitle(os.path.basename(csv_path), fontsize=13, fontweight="bold")
    fig.tight_layout(rect=[0, 0.04, 1, 0.96])

    # --------------- save ---------------
    out_path = os.path.join(os.path.dirname(os.path.abspath(csv_path)), "grid_plot.png")
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Saved: {out_path}")


def _script_sort_key(name: str):
    """Sort '6.sh' numerically by the leading integer, fall back to lexicographic."""
    try:
        return (0, int(name.split(".")[0]))
    except (ValueError, IndexError):
        return (1, name)


if __name__ == "__main__":
    main()
