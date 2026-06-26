#!/usr/bin/env python3
"""Plot a multi-tool, multi-depth ground-truth comparison grid."""

from __future__ import annotations

import argparse
import os
import sys
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatterMathtext, LogLocator
import numpy as np
import pandas as pd


DEFAULT_GROUPS = ["depth_gt40", "depth_31_40", "depth_21_30", "depth_11_20", "depth_1_10"]
DEFAULT_GROUP_LABELS = ["depth > 40", "30 < depth <= 40", "20 < depth <= 30", "10 < depth <= 20", "1 <= depth <= 10"]
CONFIG_COLUMNS = ["tool_name", "matched_table", "summary_table", "x_column", "y_column", "group_column"]
SUMMARY_COLUMNS = ["group", "n", "pearson_r2", "spearman_rho2"]


def log(message: str) -> None:
    print(message, file=sys.stderr)


def fail(message: str) -> None:
    raise ValueError(message)


def read_tsv(path: str) -> pd.DataFrame:
    if not Path(path).exists():
        raise FileNotFoundError(f"input file does not exist: {path}")
    return pd.read_csv(path, sep="\t", compression="infer", dtype=str, keep_default_na=False, na_values=[])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create a multi-tool depth-grid stoichiometry density figure.")
    parser.add_argument("--tool-config", required=True)
    parser.add_argument("--output-prefix", required=True)
    parser.add_argument("--truth-name", required=True)
    parser.add_argument("--groups", nargs="+", default=DEFAULT_GROUPS)
    parser.add_argument("--group-labels", nargs="+", default=None)
    parser.add_argument("--bins", type=int, default=150)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def resolve_group_labels(groups: list[str], labels: list[str] | None) -> list[str]:
    if labels is None:
        if groups == DEFAULT_GROUPS:
            return list(DEFAULT_GROUP_LABELS)
        return list(groups)
    if len(labels) != len(groups):
        fail("--groups and --group-labels must have the same length")
    return labels


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


def validate_matched(df: pd.DataFrame, x_col: str, y_col: str, group_col: str) -> pd.DataFrame:
    for col in [x_col, y_col, group_col]:
        if col not in df.columns:
            fail(f"matched table is missing required column: {col}")
    out = df.copy()
    for col in [x_col, y_col]:
        values = pd.to_numeric(out[col], errors="coerce")
        if values.isna().any():
            fail(f"matched table column {col} contains missing or non-numeric values")
        arr = values.to_numpy(dtype=float)
        if not np.isfinite(arr).all():
            fail(f"matched table column {col} contains infinite values")
        if ((arr < 0.0) | (arr > 1.0)).any():
            fail(f"matched table column {col} contains values outside [0, 1]")
        out[col] = values
    group_values = out[group_col].astype(str)
    if group_values.str.strip().eq("").any():
        fail(f"matched table column {group_col} contains empty group labels")
    out[group_col] = group_values
    return out


def _is_na(value: object) -> bool:
    if pd.isna(value):
        return True
    return str(value).strip().upper() in {"", "NA", "NAN"}


def validate_summary(df: pd.DataFrame, groups: list[str]) -> pd.DataFrame:
    for col in SUMMARY_COLUMNS:
        if col not in df.columns:
            fail(f"summary table is missing required column: {col}")
    if df["group"].astype(str).duplicated().any():
        duplicates = sorted(df.loc[df["group"].astype(str).duplicated(), "group"].astype(str).unique())
        fail("summary table contains duplicate group rows: " + ", ".join(duplicates))
    missing = [group for group in groups if group not in set(df["group"].astype(str))]
    if missing:
        fail("summary table is missing requested group(s): " + ", ".join(missing))
    out = df.copy()
    out["group"] = out["group"].astype(str)
    for group in groups:
        row = out.loc[out["group"] == group].iloc[0]
        n_value = pd.to_numeric(pd.Series([row["n"]]), errors="coerce").iloc[0]
        if pd.isna(n_value) or not float(n_value).is_integer() or int(n_value) < 0:
            fail(f"summary table group {group} has invalid non-negative integer n")
        for col in ["pearson_r2", "spearman_rho2"]:
            value = row[col]
            if _is_na(value):
                continue
            numeric = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
            if pd.isna(numeric) or numeric < 0.0 or numeric > 1.0:
                fail(f"summary table group {group} column {col} must be within [0, 1] or NA")
    return out


def format_stat(value: object) -> str:
    if _is_na(value):
        return "NA"
    return f"{float(value):.3f}"


def format_n(value: object) -> str:
    return f"{int(float(value)):,}"


def stats_for(summary: pd.DataFrame, group: str) -> tuple[str, str, str]:
    row = summary.loc[summary["group"] == group].iloc[0]
    return format_stat(row["pearson_r2"]), format_stat(row["spearman_rho2"]), format_n(row["n"])


def density_values(x: np.ndarray, y: np.ndarray, bins: int) -> np.ndarray:
    if len(x) == 0:
        return np.array([], dtype=float)
    counts, x_edges, y_edges = np.histogram2d(x, y, bins=bins, range=[[0.0, 1.0], [0.0, 1.0]])
    x_bin = np.clip(np.searchsorted(x_edges, x, side="right") - 1, 0, bins - 1)
    y_bin = np.clip(np.searchsorted(y_edges, y, side="right") - 1, 0, bins - 1)
    return counts[x_bin, y_bin].astype(float)


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


def load_tools(config_path: str, groups: list[str]) -> list[dict[str, object]]:
    config = validate_config(read_tsv(config_path))
    tools = []
    for row in config.to_dict("records"):
        matched = validate_matched(read_tsv(row["matched_table"]), row["x_column"], row["y_column"], row["group_column"])
        summary = validate_summary(read_tsv(row["summary_table"]), groups)
        tools.append({"config": row, "matched": matched, "summary": summary})
    return tools


def plot(args: argparse.Namespace) -> tuple[Path, Path]:
    groups = list(args.groups)
    labels = resolve_group_labels(groups, args.group_labels)
    validate_bins_and_dpi(args.bins, args.dpi)
    pdf_path, png_path = output_paths(args.output_prefix)
    refuse_overwrite((pdf_path, png_path), args.force)
    tools = load_tools(args.tool_config, groups)
    log("tool_order\t" + ",".join(str(tool["config"]["tool_name"]) for tool in tools))

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

    row_count = len(groups)
    col_count = len(tools)
    fig_width = max(3.0 * col_count, 4.2)
    fig_height = max(2.25 * row_count, 3.2)
    fig, axes = plt.subplots(row_count, col_count, figsize=(fig_width, fig_height), squeeze=False)
    norm = LogNorm(vmin=1.0, vmax=global_max)
    cmap = "viridis"

    for row_index, (group, label) in enumerate(zip(groups, labels)):
        for col_index, tool in enumerate(tools):
            cfg = tool["config"]
            ax = axes[row_index, col_index]
            x, y, density = panels[(row_index, col_index)]
            if len(x):
                order = np.argsort(density)
                ax.scatter(x[order], y[order], c=density[order], cmap=cmap, norm=norm, s=4, linewidths=0, rasterized=True)
            else:
                ax.text(0.5, 0.55, "-", transform=ax.transAxes, ha="center", va="center", fontsize=16)
            pearson_r2, spearman_rho2, n_text = stats_for(tool["summary"], group)
            if len(x) == 0:
                pearson_r2, spearman_rho2, n_text = "NA", "NA", "0"
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
                bbox={"boxstyle": "round,pad=0.22", "facecolor": "white", "edgecolor": "0.75", "alpha": 0.9},
            )
            if row_index == 0:
                ax.set_title(str(cfg["tool_name"]), fontsize=11)
            ax.set_xlim(0.0, 1.0)
            ax.set_ylim(0.0, 1.0)
            ax.set_xticks([0.0, 0.5, 1.0])
            ax.set_yticks([0.0, 0.5, 1.0])
            if col_index == 0:
                ax.set_ylabel("model modification stoichiometry")
            else:
                ax.set_ylabel("")
                ax.tick_params(labelleft=False)
            if row_index == row_count - 1:
                ax.set_xlabel(f"{args.truth_name} modification stoichiometry")
            else:
                ax.set_xlabel("")
                ax.tick_params(labelbottom=False)

    fig.subplots_adjust(left=0.08, right=0.98, top=0.92, bottom=0.11, hspace=0.32, wspace=0.10)
    for row_index, label in enumerate(labels):
        row_axes = axes[row_index, :]
        left = row_axes[0].get_position().x0
        right = row_axes[-1].get_position().x1
        top = max(ax.get_position().y1 for ax in row_axes)
        y_text = min(top + 0.018, 0.985)
        fig.text((left + right) / 2.0, y_text, label, ha="center", va="bottom", fontsize=10)
        y_line = min(top + 0.014, 0.982)
        fig.add_artist(plt.Line2D([left, right], [y_line, y_line], transform=fig.transFigure, color="0.72", linewidth=0.8))

    scalar = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    scalar.set_array([])
    colorbar = fig.colorbar(scalar, ax=axes.ravel().tolist(), location="bottom", fraction=0.035, pad=0.055, aspect=40)
    colorbar.set_label("density")
    colorbar.locator = LogLocator(base=10)
    colorbar.formatter = LogFormatterMathtext(base=10)
    colorbar.update_ticks()

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
