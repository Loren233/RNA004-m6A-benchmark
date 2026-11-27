#!/usr/bin/env python3

"""
Depth-filtered cross-project Δβ consistency analysis.

Δβ definition:
    Δβ = β_HEK293T - β_HeLa

Cross-project comparison:
    x = Δβ_106788
    y = Δβ_80229

A site is retained at threshold T only when:
    min_depth_106788 > T
and
    min_depth_80229 > T

Since each project's min_depth is already:
    min(depth_HEK, depth_HeLa)

this means that all four measurements exceed T:
    106788 HEK
    106788 HeLa
    80229 HEK
    80229 HeLa
"""

import os as _os
from pathlib import Path as _Path
RNA004_STAT_ROOT = _os.environ.get(
    "RNA004_STAT_ROOT",
    str(_Path(__file__).resolve().parents[0])
)


from pathlib import Path
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr


ROOT = Path(f"{RNA004_STAT_ROOT}")

INPUT_DIR = (
    ROOT
    / "cross_project_delta_consistency"
    / "merged_tables"
)

OUTPUT_DIR = (
    ROOT
    / "cross_project_delta_consistency"
    / "depth_filtered"
)

TABLE_DIR = OUTPUT_DIR / "filtered_tables"
SINGLE_PLOT_DIR = OUTPUT_DIR / "single_plots"
PANEL_DIR = OUTPUT_DIR / "panels"

for directory in [
    OUTPUT_DIR,
    TABLE_DIR,
    SINGLE_PLOT_DIR,
    PANEL_DIR,
]:
    directory.mkdir(parents=True, exist_ok=True)


TOOLS = [
    "deeprm",
    "dorado",
    "m6anet",
    "singlemod",
]

DEPTH_THRESHOLDS = [10, 20, 40]

PLOT_SAMPLE_MAX = 200_000

REQUIRED_COLUMNS = [
    "chr",
    "pos_1based",
    "strand",
    "site_id",
    "delta_beta_106788",
    "delta_beta_80229",
    "min_depth_106788",
    "min_depth_80229",
]


def load_merged_table(tool):
    path = (
        INPUT_DIR
        / f"{tool}.cross_project.delta_beta.merged.tsv"
    )

    if not path.exists():
        raise FileNotFoundError(
            f"Input table not found: {path}"
        )

    df = pd.read_csv(
        path,
        sep="\t",
        low_memory=False,
    )

    missing = [
        column
        for column in REQUIRED_COLUMNS
        if column not in df.columns
    ]

    if missing:
        raise RuntimeError(
            f"{path} is missing columns: {missing}\n"
            f"Available columns: {list(df.columns)}"
        )

    numeric_columns = [
        "pos_1based",
        "delta_beta_106788",
        "delta_beta_80229",
        "min_depth_106788",
        "min_depth_80229",
    ]

    for column in numeric_columns:
        df[column] = pd.to_numeric(
            df[column],
            errors="coerce",
        )

    df = df.dropna(
        subset=[
            "delta_beta_106788",
            "delta_beta_80229",
            "min_depth_106788",
            "min_depth_80229",
        ]
    ).copy()

    df = df[
        df["delta_beta_106788"].between(-1, 1)
        & df["delta_beta_80229"].between(-1, 1)
        & (df["min_depth_106788"] > 0)
        & (df["min_depth_80229"] > 0)
    ].copy()

    duplicated = df.duplicated(
        ["chr", "pos_1based", "strand"],
        keep=False,
    )

    if duplicated.any():
        raise RuntimeError(
            f"{path} contains "
            f"{int(duplicated.sum()):,} duplicated genomic-key rows."
        )

    return df


def filter_by_depth(df, threshold):
    """
    Retain sites for which all four measurements have depth > threshold.

    min_depth_106788 = min(depth_HEK_106788, depth_HeLa_106788)
    min_depth_80229  = min(depth_HEK_80229, depth_HeLa_80229)
    """

    return df[
        (df["min_depth_106788"] > threshold)
        & (df["min_depth_80229"] > threshold)
    ].copy()


def safe_correlations(x, y):
    """
    Return Pearson and Spearman values safely.

    Correlation cannot be calculated when:
      - fewer than two observations remain;
      - all x values are identical;
      - all y values are identical.
    """

    if len(x) < 2:
        return np.nan, np.nan, np.nan, np.nan

    if np.all(x == x[0]) or np.all(y == y[0]):
        return np.nan, np.nan, np.nan, np.nan

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        pearson_r, pearson_p = pearsonr(x, y)
        spearman_rho, spearman_p = spearmanr(x, y)

    return (
        pearson_r,
        pearson_p,
        spearman_rho,
        spearman_p,
    )


def calculate_statistics(df, tool, threshold):
    x = df["delta_beta_106788"].to_numpy()
    y = df["delta_beta_80229"].to_numpy()

    n = len(df)

    (
        pearson_r,
        pearson_p,
        spearman_rho,
        spearman_p,
    ) = safe_correlations(x, y)

    if n == 0:
        mae = np.nan
        rmse = np.nan
        mean_difference = np.nan
        sign_concordance = np.nan
        nonzero_sign_concordance = np.nan
        same_direction_quadrants = np.nan
        both_zero_fraction = np.nan
    else:
        difference = y - x

        mae = np.mean(np.abs(difference))
        rmse = np.sqrt(np.mean(difference ** 2))
        mean_difference = np.mean(difference)

        sign_x = np.sign(x)
        sign_y = np.sign(y)

        sign_concordance = np.mean(
            sign_x == sign_y
        )

        both_nonzero = (
            (sign_x != 0)
            & (sign_y != 0)
        )

        if both_nonzero.any():
            nonzero_sign_concordance = np.mean(
                sign_x[both_nonzero]
                == sign_y[both_nonzero]
            )
        else:
            nonzero_sign_concordance = np.nan

        same_direction_quadrants = np.mean(
            ((x > 0) & (y > 0))
            | ((x < 0) & (y < 0))
        )

        both_zero_fraction = np.mean(
            (x == 0) & (y == 0)
        )

    return {
        "tool": tool,
        "depth_rule": f"all four depths > {threshold}",
        "depth_threshold": threshold,
        "n_sites": n,
        "pearson_r": pearson_r,
        "pearson_p": pearson_p,
        "spearman_rho": spearman_rho,
        "spearman_p": spearman_p,
        "mae": mae,
        "rmse": rmse,
        "mean_delta_80229_minus_106788": mean_difference,
        "sign_concordance_including_zero": sign_concordance,
        "sign_concordance_both_nonzero": nonzero_sign_concordance,
        "same_direction_quadrant_fraction": same_direction_quadrants,
        "both_zero_fraction": both_zero_fraction,
    }


def format_metric(value):
    if pd.isna(value):
        return "NA"
    return f"{value:.3f}"


def sample_for_plot(df):
    if len(df) <= PLOT_SAMPLE_MAX:
        return df

    return df.sample(
        n=PLOT_SAMPLE_MAX,
        random_state=42,
    )


def draw_scatter(ax, df, stats, tool, threshold):
    plot_df = sample_for_plot(df)

    if not plot_df.empty:
        ax.scatter(
            plot_df["delta_beta_106788"],
            plot_df["delta_beta_80229"],
            s=4,
            alpha=0.20,
            rasterized=True,
        )

    # Perfect cross-project agreement.
    ax.plot(
        [-1, 1],
        [-1, 1],
        linestyle="--",
        linewidth=1,
    )

    # Separate direction quadrants.
    ax.axhline(
        0,
        linestyle=":",
        linewidth=1,
    )

    ax.axvline(
        0,
        linestyle=":",
        linewidth=1,
    )

    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)

    ax.set_xlabel(
        "Δβ_106788 = β_HEK293T − β_HeLa"
    )

    ax.set_ylabel(
        "Δβ_80229 = β_HEK293T − β_HeLa"
    )

    ax.set_title(
        f"{tool} | all four depths > {threshold}"
    )

    text = (
        f"n = {stats['n_sites']:,}\n"
        f"Pearson r = {format_metric(stats['pearson_r'])}\n"
        f"Spearman ρ = {format_metric(stats['spearman_rho'])}\n"
        f"MAE = {format_metric(stats['mae'])}\n"
        f"RMSE = {format_metric(stats['rmse'])}\n"
        f"Sign agreement = "
        f"{format_metric(stats['sign_concordance_including_zero'])}\n"
        f"Nonzero sign agreement = "
        f"{format_metric(stats['sign_concordance_both_nonzero'])}"
    )

    ax.text(
        0.03,
        0.97,
        text,
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=9,
        bbox={
            "boxstyle": "round",
            "alpha": 0.15,
        },
    )

    ax.set_aspect(
        "equal",
        adjustable="box",
    )


def make_single_plot(
    df,
    stats,
    tool,
    threshold,
):
    fig, ax = plt.subplots(
        figsize=(7, 7)
    )

    draw_scatter(
        ax=ax,
        df=df,
        stats=stats,
        tool=tool,
        threshold=threshold,
    )

    fig.tight_layout()

    output = (
        SINGLE_PLOT_DIR
        / (
            f"{tool}."
            f"delta_beta_80229_vs_106788."
            f"depth_gt_{threshold}.png"
        )
    )

    fig.savefig(
        output,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)

    return output


def make_panel(
    filtered_data,
    statistics,
    threshold,
):
    fig, axes = plt.subplots(
        nrows=2,
        ncols=2,
        figsize=(14, 14),
    )

    axes = axes.flatten()

    for ax, tool in zip(axes, TOOLS):
        draw_scatter(
            ax=ax,
            df=filtered_data[tool],
            stats=statistics[tool],
            tool=tool,
            threshold=threshold,
        )

    fig.suptitle(
        (
            "Cross-project Δβ consistency "
            f"| all four depths > {threshold}"
        ),
        fontsize=16,
        y=0.995,
    )

    fig.tight_layout(
        rect=[0, 0, 1, 0.98]
    )

    output = (
        PANEL_DIR
        / (
            "cross_project_delta_consistency."
            f"depth_gt_{threshold}.panel.png"
        )
    )

    fig.savefig(
        output,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)

    return output


def main():
    all_tables = {}

    for tool in TOOLS:
        print(
            f"[INFO] Loading merged table: {tool}"
        )

        all_tables[tool] = load_merged_table(
            tool
        )

        print(
            f"[INFO] {tool}: "
            f"{len(all_tables[tool]):,} unfiltered sites"
        )

    summary_rows = []

    for threshold in DEPTH_THRESHOLDS:
        print()
        print(
            "============================================"
        )
        print(
            f"[INFO] Processing depth > {threshold}"
        )
        print(
            "============================================"
        )

        filtered_data = {}
        threshold_statistics = {}

        for tool in TOOLS:
            filtered = filter_by_depth(
                all_tables[tool],
                threshold,
            )

            filtered_data[tool] = filtered

            stats = calculate_statistics(
                filtered,
                tool,
                threshold,
            )

            threshold_statistics[tool] = stats
            summary_rows.append(stats)

            filtered_output = (
                TABLE_DIR
                / (
                    f"{tool}."
                    f"cross_project.delta_beta."
                    f"depth_gt_{threshold}.tsv"
                )
            )

            filtered.to_csv(
                filtered_output,
                sep="\t",
                index=False,
            )

            single_plot = make_single_plot(
                df=filtered,
                stats=stats,
                tool=tool,
                threshold=threshold,
            )

            print(
                f"[INFO] {tool}, depth > {threshold}: "
                f"n={len(filtered):,}, "
                f"Pearson={format_metric(stats['pearson_r'])}, "
                f"Spearman={format_metric(stats['spearman_rho'])}"
            )

            print(
                f"[INFO] Single plot: {single_plot}"
            )

        panel_output = make_panel(
            filtered_data=filtered_data,
            statistics=threshold_statistics,
            threshold=threshold,
        )

        print(
            f"[INFO] Panel: {panel_output}"
        )

    summary = pd.DataFrame(
        summary_rows
    )

    summary_output = (
        OUTPUT_DIR
        / "cross_project_delta_by_depth_summary.tsv"
    )

    summary.to_csv(
        summary_output,
        sep="\t",
        index=False,
    )

    print()
    print(
        f"[INFO] Summary written to: {summary_output}"
    )


if __name__ == "__main__":
    main()
