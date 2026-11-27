#!/usr/bin/env python3

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr


KEY = ["chr", "pos_1based", "strand"]

TOOLS = [
    "dorado",
    "m6anet",
    "singlemod",
]


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--glori", required=True)

    parser.add_argument("--dorado", required=True)
    parser.add_argument("--m6anet", required=True)
    parser.add_argument("--singlemod", required=True)

    parser.add_argument("--output-dir", required=True)

    parser.add_argument(
        "--sample-size",
        type=int,
        default=300000,
    )

    return parser.parse_args()


def normalize_strand(series):
    return (
        series.astype(str)
        .str.strip()
        .replace(
            {
                "1": "+",
                "1.0": "+",
                "plus": "+",
                "-1": "-",
                "-1.0": "-",
                "minus": "-",
            }
        )
    )


def load_glori(path):
    print("=" * 80)
    print("Loading GLORI")
    print(path)

    df = pd.read_csv(
        path,
        sep="\t",
        compression="infer",
        usecols=[
            "chr",
            "pos_1based",
            "strand",
            "glori_beta_hek",
            "glori_beta_hela",
            "glori_delta_beta",
            "glori_total_hek",
            "glori_total_hela",
            "glori_min_depth",
        ],
        low_memory=False,
    )

    df["chr"] = df["chr"].astype(str).str.strip()
    df["strand"] = normalize_strand(df["strand"])

    numeric_cols = [
        "pos_1based",
        "glori_beta_hek",
        "glori_beta_hela",
        "glori_delta_beta",
        "glori_total_hek",
        "glori_total_hela",
        "glori_min_depth",
    ]

    for col in numeric_cols:
        df[col] = pd.to_numeric(
            df[col],
            errors="coerce",
        )

    valid = (
        df["pos_1based"].notna()
        & df["strand"].isin(["+", "-"])
        & np.isfinite(df["glori_delta_beta"])
        & df["glori_beta_hek"].between(0, 1)
        & df["glori_beta_hela"].between(0, 1)
        & (df["glori_total_hek"] > 0)
        & (df["glori_total_hela"] > 0)
    )

    df = df.loc[valid].copy()

    df["pos_1based"] = (
        np.rint(df["pos_1based"])
        .astype(np.int64)
    )

    duplicated = df.duplicated(
        KEY,
        keep=False,
    )

    if duplicated.any():
        raise ValueError(
            "GLORI contains duplicated genomic keys: "
            f"{duplicated.sum():,} rows"
        )

    print(
        f"Valid GLORI sites with coverage in both cell lines: "
        f"{len(df):,}"
    )

    return df


def detect_column(columns, candidates, label):
    lookup = {
        str(c).lower(): c
        for c in columns
    }

    for candidate in candidates:
        if candidate.lower() in lookup:
            return lookup[candidate.lower()]

    raise ValueError(
        f"Cannot find {label} column.\n"
        f"Available: {columns}"
    )


def load_tool(path, tool):
    print()
    print("=" * 80)
    print(f"Loading {tool}")
    print(path)

    header = pd.read_csv(
        path,
        sep="\t",
        nrows=0,
        compression="infer",
    )

    columns = list(header.columns)

    chr_col = detect_column(
        columns,
        ["chr", "chrom", "chromosome"],
        "chromosome",
    )

    pos_col = detect_column(
        columns,
        [
            "pos_1based",
            "position_1based",
            "position",
            "pos",
        ],
        "position",
    )

    strand_col = detect_column(
        columns,
        ["strand"],
        "strand",
    )

    delta_col = detect_column(
        columns,
        [
            "delta_beta",
            "delta_beta_hek_minus_hela",
        ],
        "delta beta",
    )

    hek_depth_col = detect_column(
        columns,
        [
            "depth_hek",
            "hek_depth",
            "beta_hek_depth",
        ],
        "HEK depth",
    )

    hela_depth_col = detect_column(
        columns,
        [
            "depth_hela",
            "hela_depth",
            "beta_hela_depth",
        ],
        "HeLa depth",
    )

    usecols = [
        chr_col,
        pos_col,
        strand_col,
        delta_col,
        hek_depth_col,
        hela_depth_col,
    ]

    # Avoid duplicate read when naming happens to coincide.
    usecols = list(dict.fromkeys(usecols))

    df = pd.read_csv(
        path,
        sep="\t",
        usecols=usecols,
        compression="infer",
        low_memory=False,
    )

    df = df.rename(
        columns={
            chr_col: "chr",
            pos_col: "pos_1based",
            strand_col: "strand",
            delta_col: "delta_beta_tool",
            hek_depth_col: "tool_depth_hek",
            hela_depth_col: "tool_depth_hela",
        }
    )

    df["chr"] = df["chr"].astype(str).str.strip()
    df["strand"] = normalize_strand(df["strand"])

    for col in [
        "pos_1based",
        "delta_beta_tool",
        "tool_depth_hek",
        "tool_depth_hela",
    ]:
        df[col] = pd.to_numeric(
            df[col],
            errors="coerce",
        )

    valid = (
        df["pos_1based"].notna()
        & df["strand"].isin(["+", "-"])
        & np.isfinite(df["delta_beta_tool"])
        & (df["tool_depth_hek"] > 0)
        & (df["tool_depth_hela"] > 0)
    )

    df = df.loc[valid].copy()

    df["pos_1based"] = (
        np.rint(df["pos_1based"])
        .astype(np.int64)
    )

    df["tool_min_depth"] = np.minimum(
        df["tool_depth_hek"],
        df["tool_depth_hela"],
    )

    duplicated = df.duplicated(
        KEY,
        keep=False,
    )

    if duplicated.any():
        raise ValueError(
            f"{tool} contains duplicated genomic keys: "
            f"{duplicated.sum():,} rows"
        )

    print(f"{tool} delta sites: {len(df):,}")

    return df


def compute_metrics(df, tool, threshold_name):
    x = df["glori_delta_beta"].to_numpy(
        dtype=np.float64
    )

    y = df["delta_beta_tool"].to_numpy(
        dtype=np.float64
    )

    diff = y - x

    if len(df) >= 2:
        pearson = float(
            pearsonr(x, y).statistic
        )

        spearman = float(
            spearmanr(x, y).statistic
        )
    else:
        pearson = np.nan
        spearman = np.nan

    nonzero_both = (
        (x != 0)
        & (y != 0)
    )

    sign_concordance_all = float(
        np.mean(np.sign(x) == np.sign(y))
    )

    if nonzero_both.any():
        sign_concordance_nonzero = float(
            np.mean(
                np.sign(x[nonzero_both])
                == np.sign(y[nonzero_both])
            )
        )
    else:
        sign_concordance_nonzero = np.nan

    return {
        "tool": tool,
        "depth_filter": threshold_name,
        "n": int(len(df)),
        "pearson_r": pearson,
        "spearman_rho": spearman,
        "mae": float(
            np.mean(np.abs(diff))
        ),
        "rmse": float(
            np.sqrt(np.mean(diff ** 2))
        ),
        "mean_tool_minus_glori": float(
            np.mean(diff)
        ),
        "median_tool_minus_glori": float(
            np.median(diff)
        ),
        "sign_concordance_all": (
            sign_concordance_all
        ),
        "sign_concordance_nonzero_both": (
            sign_concordance_nonzero
        ),
        "glori_delta_mean": float(
            np.mean(x)
        ),
        "tool_delta_mean": float(
            np.mean(y)
        ),
    }


def plot_single(
    df,
    tool,
    threshold_name,
    output_path,
    sample_size,
    metrics,
):
    if len(df) > sample_size:
        plot_df = df.sample(
            n=sample_size,
            random_state=42,
        )
    else:
        plot_df = df

    x = plot_df[
        "glori_delta_beta"
    ].to_numpy()

    y = plot_df[
        "delta_beta_tool"
    ].to_numpy()

    fig, ax = plt.subplots(
        figsize=(7, 6)
    )

    hb = ax.hexbin(
        x,
        y,
        gridsize=160,
        bins="log",
        mincnt=1,
    )

    ax.plot(
        [-1, 1],
        [-1, 1],
        linestyle="--",
        linewidth=1,
    )

    ax.axhline(
        0,
        linestyle=":",
        linewidth=0.8,
    )

    ax.axvline(
        0,
        linestyle=":",
        linewidth=0.8,
    )

    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)

    ax.set_xlabel(
        "Δβ GLORI = β_HEK − β_HeLa"
    )

    ax.set_ylabel(
        f"Δβ {tool} PRJEB80229 = β_HEK − β_HeLa"
    )

    ax.set_title(
        f"{tool} PRJEB80229 vs GLORI\n"
        f"{threshold_name}"
    )

    text = (
        f"n = {metrics['n']:,}\n"
        f"Pearson r = {metrics['pearson_r']:.3f}\n"
        f"Spearman ρ = {metrics['spearman_rho']:.3f}\n"
        f"MAE = {metrics['mae']:.3f}\n"
        f"RMSE = {metrics['rmse']:.3f}\n"
        f"Sign concordance = "
        f"{metrics['sign_concordance_all']:.3f}"
    )

    ax.text(
        0.03,
        0.97,
        text,
        transform=ax.transAxes,
        va="top",
        ha="left",
        bbox={
            "boxstyle": "round",
            "alpha": 0.15,
        },
    )

    fig.colorbar(
        hb,
        ax=ax,
        label="Displayed-site count (log scale)",
    )

    fig.tight_layout()
    fig.savefig(
        output_path,
        dpi=300,
    )

    plt.close(fig)


def make_panel(
    plot_data,
    threshold_name,
    output_path,
    sample_size,
):
    fig, axes = plt.subplots(
        1,
        3,
        figsize=(18, 5.5),
    )

    for ax, tool in zip(
        axes,
        TOOLS,
    ):
        df, metrics = plot_data[tool]

        if len(df) > sample_size:
            plot_df = df.sample(
                sample_size,
                random_state=42,
            )
        else:
            plot_df = df

        hb = ax.hexbin(
            plot_df["glori_delta_beta"],
            plot_df["delta_beta_tool"],
            gridsize=130,
            bins="log",
            mincnt=1,
        )

        ax.plot(
            [-1, 1],
            [-1, 1],
            linestyle="--",
            linewidth=1,
        )

        ax.axhline(
            0,
            linestyle=":",
            linewidth=0.8,
        )

        ax.axvline(
            0,
            linestyle=":",
            linewidth=0.8,
        )

        ax.set_xlim(-1, 1)
        ax.set_ylim(-1, 1)

        ax.set_title(tool)

        ax.set_xlabel(
            "Δβ GLORI"
        )

        ax.set_ylabel(
            f"Δβ {tool}"
        )

        text = (
            f"n={metrics['n']:,}\n"
            f"Pearson={metrics['pearson_r']:.3f}\n"
            f"Spearman={metrics['spearman_rho']:.3f}\n"
            f"MAE={metrics['mae']:.3f}"
        )

        ax.text(
            0.03,
            0.97,
            text,
            transform=ax.transAxes,
            va="top",
            bbox={
                "boxstyle": "round",
                "alpha": 0.15,
            },
        )

    fig.suptitle(
        f"PRJEB80229 tool Δβ vs GLORI Δβ — {threshold_name}"
    )

    fig.tight_layout()

    fig.savefig(
        output_path,
        dpi=300,
    )

    plt.close(fig)


def main():
    args = parse_args()

    output_dir = Path(
        args.output_dir
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    glori = load_glori(
        args.glori
    )

    paths = {
        "dorado": args.dorado,
        "m6anet": args.m6anet,
        "singlemod": args.singlemod,
    }

    tool_tables = {
        tool: load_tool(
            paths[tool],
            tool,
        )
        for tool in TOOLS
    }

    thresholds = [
        ("all_common_sites", None),
        ("depth_ge_10", 10),
        ("depth_ge_20", 20),
        ("depth_ge_40", 40),
    ]

    all_metrics = []

    for threshold_name, threshold in thresholds:

        print()
        print("#" * 80)
        print(threshold_name)
        print("#" * 80)

        panel_data = {}

        for tool in TOOLS:

            merged = tool_tables[
                tool
            ].merge(
                glori,
                on=KEY,
                how="inner",
                validate="one_to_one",
                sort=False,
            )

            if threshold is not None:
                merged = merged.loc[
                    (
                        merged[
                            "tool_depth_hek"
                        ] >= threshold
                    )
                    & (
                        merged[
                            "tool_depth_hela"
                        ] >= threshold
                    )
                    & (
                        merged[
                            "glori_total_hek"
                        ] >= threshold
                    )
                    & (
                        merged[
                            "glori_total_hela"
                        ] >= threshold
                    )
                ].copy()

            metrics = compute_metrics(
                merged,
                tool,
                threshold_name,
            )

            all_metrics.append(
                metrics
            )

            print(
                f"{tool:10s} "
                f"n={metrics['n']:,} "
                f"Pearson={metrics['pearson_r']:.4f} "
                f"Spearman={metrics['spearman_rho']:.4f} "
                f"MAE={metrics['mae']:.4f} "
                f"RMSE={metrics['rmse']:.4f}"
            )

            matched_path = (
                output_dir
                / (
                    f"{tool}.80229_vs_GLORI."
                    f"{threshold_name}.matched.tsv.gz"
                )
            )

            merged[
                KEY
                + [
                    "delta_beta_tool",
                    "glori_delta_beta",
                    "tool_depth_hek",
                    "tool_depth_hela",
                    "tool_min_depth",
                    "glori_total_hek",
                    "glori_total_hela",
                    "glori_min_depth",
                ]
            ].to_csv(
                matched_path,
                sep="\t",
                index=False,
                compression={
                    "method": "gzip",
                    "compresslevel": 1,
                },
                float_format="%.8g",
            )

            plot_single(
                merged,
                tool,
                threshold_name,
                output_dir
                / (
                    f"{tool}.80229_vs_GLORI."
                    f"{threshold_name}.png"
                ),
                args.sample_size,
                metrics,
            )

            panel_data[tool] = (
                merged,
                metrics,
            )

        make_panel(
            panel_data,
            threshold_name,
            output_dir
            / (
                f"panel.80229_vs_GLORI."
                f"{threshold_name}.png"
            ),
            args.sample_size,
        )

    metrics_df = pd.DataFrame(
        all_metrics
    )

    metrics_df.to_csv(
        output_dir
        / "80229_tools_vs_GLORI.metrics.tsv",
        sep="\t",
        index=False,
        float_format="%.10g",
    )

    print()
    print("=" * 80)
    print("FINAL METRICS")
    print("=" * 80)

    print(
        metrics_df.to_string(
            index=False
        )
    )


if __name__ == "__main__":
    main()
