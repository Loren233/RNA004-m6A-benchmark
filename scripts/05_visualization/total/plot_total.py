#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Plot a multi-tool, multi-depth ground-truth comparison grid."""

from __future__ import annotations

import argparse
import os
import re
import sys
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.ticker import LogFormatterMathtext, LogLocator
import numpy as np
import pandas as pd


DEFAULT_GROUPS = [
    "depth_gt40",
    "depth_31_40",
    "depth_21_30",
    "depth_11_20",
    "depth_1_10",
]

DEFAULT_GROUP_LABELS = [
    "depth > 40",
    "30 < depth ≤ 40",
    "20 < depth ≤ 30",
    "10 < depth ≤ 20",
    "1 ≤ depth ≤ 10",
]

CONFIG_COLUMNS = [
    "tool_name",
    "matched_table",
    "summary_table",
    "x_column",
    "y_column",
    "group_column",
]

PAPER_DENSITY_CMAP = LinearSegmentedColormap.from_list(
    "paper_density",
    [
        "#3b0f70",
        "#8c2981",
        "#de4968",
        "#f98e52",
        "#fee08b",
    ],
    N=256,
)


def log(message: str) -> None:
    print(message, file=sys.stderr)


def fail(message: str) -> None:
    raise ValueError(message)


def read_tsv(path: str) -> pd.DataFrame:
    if not Path(path).exists():
        raise FileNotFoundError(f"input file does not exist: {path}")
    return pd.read_csv(
        path,
        sep="\t",
        compression="infer",
        dtype=str,
        keep_default_na=False,
        na_values=[],
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a multi-tool depth-grid stoichiometry density figure."
    )
    parser.add_argument("--tool-config", required=True)
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--truth-name", required=True)
    parser.add_argument("--groups", nargs="+", default=DEFAULT_GROUPS)
    parser.add_argument("--group-labels", nargs="+", default=None)
    parser.add_argument("--bins", type=int, default=150)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def normalize_display_label(label: str) -> str:
    return (
        str(label)
        .replace("<=", "≤")
        .replace(">=", "≥")
        .replace(" <= ", " ≤ ")
        .replace(" >= ", " ≥ ")
    )


def resolve_group_labels(groups: list[str], labels: list[str] | None) -> list[str]:
    if labels is None:
        if groups == DEFAULT_GROUPS:
            return list(DEFAULT_GROUP_LABELS)
        return list(groups)
    if len(labels) != len(groups):
        fail("--groups and --group-labels must have the same length")
    return list(labels)


def validate_bins_and_dpi(bins: int, dpi: int) -> None:
    if bins <= 0:
        fail("--bins must be a positive integer")
    if dpi <= 0:
        fail("--dpi must be a positive integer")


def validate_config(df: pd.DataFrame) -> pd.DataFrame:
    for col in CONFIG_COLUMNS:
        if col not in df.columns:
            fail(f"tool config is missing required column: {col}")

    if len(df) == 0:
        fail("tool config must contain at least one tool row")

    out = df.copy()
    for col in CONFIG_COLUMNS:
        out[col] = out[col].astype(str).str.strip()
        if out[col].eq("").any():
            fail(f"tool config column {col} contains empty values")

    return out


def validate_matched(
    df: pd.DataFrame,
    x_col: str,
    y_col: str,
    group_col: str,
    source_name: str,
) -> pd.DataFrame:
    for col in [x_col, y_col, group_col]:
        if col not in df.columns:
            fail(
                f"{source_name} matched table is missing required column: {col}\n"
                f"existing columns: {list(df.columns)}"
            )

    out = df.copy()

    for col in [x_col, y_col]:
        out[col] = pd.to_numeric(out[col], errors="coerce")

    before = len(out)
    out = out.dropna(subset=[x_col, y_col, group_col]).copy()
    dropped = before - len(out)
    if dropped:
        log(f"drop_non_numeric_or_missing\t{source_name}\t{dropped}")

    for col in [x_col, y_col]:
        arr = out[col].to_numpy(dtype=float)
        if len(arr) and not np.isfinite(arr).all():
            fail(f"{source_name} matched table column {col} contains infinite values")
        out = out.loc[(out[col] >= 0.0) & (out[col] <= 1.0)].copy()

    out[group_col] = out[group_col].astype(str)
    if out[group_col].str.strip().eq("").any():
        fail(f"{source_name} matched table column {group_col} contains empty group labels")

    return out


def density_values(x: np.ndarray, y: np.ndarray, bins: int) -> np.ndarray:
    if len(x) == 0:
        return np.array([], dtype=float)

    counts, x_edges, y_edges = np.histogram2d(
        x,
        y,
        bins=bins,
        range=[[0.0, 1.0], [0.0, 1.0]],
    )
    x_bin = np.clip(np.searchsorted(x_edges, x, side="right") - 1, 0, bins - 1)
    y_bin = np.clip(np.searchsorted(y_edges, y, side="right") - 1, 0, bins - 1)

    return counts[x_bin, y_bin].astype(float)


def calc_stats(x: np.ndarray, y: np.ndarray) -> tuple[str, str, str]:
    sub = pd.DataFrame({"x": x, "y": y}).dropna()
    n = len(sub)

    if n < 2:
        return "nan", "nan", f"{n:,}"

    pearson_r = sub["x"].corr(sub["y"], method="pearson")
    spearman_rho = sub["x"].rank().corr(sub["y"].rank(), method="pearson")

    if pd.isna(pearson_r):
        pearson_r2 = "nan"
    else:
        pearson_r2 = f"{float(pearson_r) ** 2:.3f}"

    if pd.isna(spearman_rho):
        spearman_rho2 = "nan"
    else:
        spearman_rho2 = f"{float(spearman_rho) ** 2:.3f}"

    return pearson_r2, spearman_rho2, f"{n:,}"


def output_paths(prefix: str) -> tuple[Path, Path]:
    return Path(prefix + ".pdf"), Path(prefix + ".png")


def refuse_overwrite(paths: tuple[Path, Path], force: bool) -> None:
    if force:
        return

    existing = [str(path) for path in paths if path.exists()]
    if existing:
        fail("output exists; use --force to overwrite: " + ", ".join(existing))


def atomic_save(fig: plt.Figure, path: Path, dpi: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_name(f"{path.name}.tmp.{os.getpid()}{path.suffix}")
    fig.savefig(tmp_path, dpi=dpi, facecolor="white", bbox_inches="tight")
    os.replace(tmp_path, path)


def warn_if_deeprm_not_offset0(tool_name: str, matched_path: str) -> None:
    if tool_name.lower() != "deeprm":
        return

    path_text = str(matched_path)

    if re.search(r"offset[_-](p1|m1|\+1|-1)", path_text):
        log(
            "WARNING\tDeepRM matched table does not look like offset_0. "
            f"Current path: {path_text}"
        )

    if "offset_0" not in path_text and "offset0" not in path_text:
        log(
            "WARNING\tDeepRM path does not explicitly contain offset_0. "
            f"Please check wrapper DEEPRM_LABEL. Current path: {path_text}"
        )


def load_tools(config_path: str) -> list[dict[str, object]]:
    config = validate_config(read_tsv(config_path))
    tools: list[dict[str, object]] = []

    for row in config.to_dict("records"):
        tool_name = str(row["tool_name"])
        matched_path = str(row["matched_table"])

        warn_if_deeprm_not_offset0(tool_name, matched_path)

        matched = validate_matched(
            read_tsv(matched_path),
            row["x_column"],
            row["y_column"],
            row["group_column"],
            source_name=tool_name,
        )

        tools.append(
            {
                "config": row,
                "matched": matched,
            }
        )

        log(
            "loaded_tool\t"
            f"{tool_name}\tmatched={matched_path}\t"
            f"x={row['x_column']}\ty={row['y_column']}\tgroup={row['group_column']}\t"
            f"rows={len(matched)}"
        )

    return tools


def plot(args: argparse.Namespace) -> tuple[Path, Path]:
    groups = list(args.groups)
    labels = [
        normalize_display_label(label)
        for label in resolve_group_labels(groups, args.group_labels)
    ]

    validate_bins_and_dpi(args.bins, args.dpi)

    pdf_path, png_path = output_paths(args.output_prefix)
    refuse_overwrite((pdf_path, png_path), args.force)

    tools = load_tools(args.tool_config)
    log("tool_order\t" + ",".join(str(tool["config"]["tool_name"]) for tool in tools))

    row_count = len(groups)
    col_count = len(tools)

    panels: dict[tuple[int, int], tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    global_max = 1.0

    for col_index, tool in enumerate(tools):
        cfg = tool["config"]
        matched = tool["matched"]

        for row_index, group in enumerate(groups):
            sub = matched.loc[matched[cfg["group_column"]] == group]
            x = sub[cfg["x_column"]].to_numpy(dtype=float)
            y = sub[cfg["y_column"]].to_numpy(dtype=float)

            density = density_values(x, y, args.bins)
            if len(density):
                global_max = max(global_max, float(density.max()))

            panels[(row_index, col_index)] = (x, y, density)

    fig_width = max(2.70 * col_count, 7.0)
    fig_height = max(2.45 * row_count + 0.80, 7.5)

    fig = plt.figure(figsize=(fig_width, fig_height))

    height_ratios: list[float] = []
    for _ in range(row_count):
        height_ratios.extend([0.16, 1.00])
    height_ratios.extend([0.25, 0.16])

    gs = fig.add_gridspec(
        nrows=row_count * 2 + 2,
        ncols=col_count,
        height_ratios=height_ratios,
        hspace=0.20,
        wspace=0.28,
    )

    axes = np.empty((row_count, col_count), dtype=object)
    label_axes = []

    for row_index, label in enumerate(labels):
        label_ax = fig.add_subplot(gs[row_index * 2, :])
        label_ax.set_axis_off()
        label_ax.set_xlim(0.0, 1.0)
        label_ax.set_ylim(0.0, 1.0)
        label_axes.append((label_ax, label))

        for col_index in range(col_count):
            axes[row_index, col_index] = fig.add_subplot(
                gs[row_index * 2 + 1, col_index]
            )

    spacer_ax = fig.add_subplot(gs[-2, :])
    spacer_ax.set_axis_off()

    cbar_ax = fig.add_subplot(gs[-1, :])

    norm = LogNorm(vmin=1.0, vmax=max(1.0, global_max))
    cmap = PAPER_DENSITY_CMAP
    ticks = [0.0, 0.5, 1.0]

    for row_index, group in enumerate(groups):
        for col_index, tool in enumerate(tools):
            cfg = tool["config"]
            tool_name = str(cfg["tool_name"])
            ax = axes[row_index, col_index]

            x, y, density = panels[(row_index, col_index)]

            if len(x):
                order = np.argsort(density)
                ax.scatter(
                    x[order],
                    y[order],
                    c=density[order],
                    cmap=cmap,
                    norm=norm,
                    s=3.8,
                    linewidths=0,
                    rasterized=True,
                )
            else:
                ax.text(
                    0.5,
                    0.55,
                    "-",
                    transform=ax.transAxes,
                    ha="center",
                    va="center",
                    fontsize=16,
                    zorder=4,
                )

            pearson_r2, spearman_rho2, n_text = calc_stats(x, y)

            stats_text = (
                f"R² = {pearson_r2}\nρ² = {spearman_rho2}\nn = {n_text}"
                if row_index == 0 and col_index == 0
                else f"{pearson_r2}\n{spearman_rho2}\n{n_text}"
            )

            ax.text(
                0.04,
                0.94,
                stats_text,
                transform=ax.transAxes,
                ha="left",
                va="top",
                fontsize=8,
                bbox=dict(
                    boxstyle="round,pad=0.20",
                    facecolor="white",
                    edgecolor="none",
                    alpha=0.70,
                ),
                zorder=5,
            )

            ax.set_title(tool_name, fontsize=9, fontweight="bold", pad=2)

            ax.set_xlim(0.0, 1.0)
            ax.set_ylim(0.0, 1.0)
            ax.set_aspect("equal", adjustable="box")

            ax.set_xticks(ticks)
            ax.set_yticks(ticks)
            ax.set_xticklabels(["0.0", "0.5", "1.0"], fontsize=7)
            ax.set_yticklabels(["0.0", "0.5", "1.0"], fontsize=7)

            ax.tick_params(
                axis="both",
                which="both",
                labelsize=7,
                labelbottom=True,
                labelleft=True,
                pad=1,
            )

            if col_index == 0:
                ax.set_ylabel("model modification stoichiometry", fontsize=8)
            else:
                ax.set_ylabel("")

            if row_index == row_count - 1:
                ax.set_xlabel(
                    f"{args.truth_name} modification stoichiometry",
                    fontsize=8,
                )
            else:
                ax.set_xlabel("")

    for label_ax, label in label_axes:
        label_ax.axhline(
            0.42,
            xmin=0.0,
            xmax=1.0,
            color="0.25",
            linewidth=0.8,
            clip_on=False,
            zorder=1,
        )
        label_ax.text(
            0.5,
            0.47,
            label,
            ha="center",
            va="center",
            fontsize=10,
            fontweight="bold",
            bbox=dict(
                facecolor="white",
                edgecolor="none",
                alpha=0.95,
                pad=0.25,
            ),
            zorder=3,
        )

    scalar = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    scalar.set_array([])

    colorbar = fig.colorbar(
        scalar,
        cax=cbar_ax,
        orientation="horizontal",
    )
    colorbar.set_label("density", fontsize=8)
    colorbar.locator = LogLocator(base=10)
    colorbar.formatter = LogFormatterMathtext(base=10)
    colorbar.update_ticks()
    colorbar.ax.tick_params(labelsize=7)

    atomic_save(fig, pdf_path, args.dpi)
    atomic_save(fig, png_path, args.dpi)
    plt.close(fig)

    return pdf_path, png_path


def main() -> int:
    args = parse_args()

    try:
        log(f"start_time\t{datetime.now().isoformat(timespec='seconds')}")
        log(f"tool_config\t{args.tool_config}")
        log("groups\t" + ",".join(args.groups))

        pdf_path, png_path = plot(args)

        log(f"output_pdf\t{pdf_path}")
        log(f"output_png\t{png_path}")
        log(f"finish_time\t{datetime.now().isoformat(timespec='seconds')}")
        return 0

    except Exception as exc:
        log(f"ERROR: {exc}")
        return 1


if __name__ == "__main__":
    sys.exit(main())
