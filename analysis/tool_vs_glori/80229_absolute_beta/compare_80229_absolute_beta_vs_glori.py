#!/usr/bin/env python3

import argparse
import gc
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from matplotlib.colors import LinearSegmentedColormap
from scipy.stats import pearsonr, spearmanr


# ============================================================
# Global settings
# ============================================================

KEY = [
    "chr",
    "pos_1based",
    "strand",
]

TOOLS = [
    "dorado",
    "m6anet",
    "singlemod",
]

CELLS = [
    "HEK293T",
    "HeLa",
]


# ============================================================
# Custom paper-density colour map
# ============================================================

PAPER_DENSITY = LinearSegmentedColormap.from_list(
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


# ============================================================
# Arguments
# ============================================================

def parse_args():

    parser = argparse.ArgumentParser(
        description=(
            "Compare PRJEB80229 absolute m6A stoichiometry "
            "against pooled GLORI ground truth."
        )
    )

    parser.add_argument(
        "--glori",
        required=True,
    )

    parser.add_argument(
        "--dorado-hek",
        required=True,
    )

    parser.add_argument(
        "--dorado-hela",
        required=True,
    )

    parser.add_argument(
        "--m6anet-hek",
        required=True,
    )

    parser.add_argument(
        "--m6anet-hela",
        required=True,
    )

    parser.add_argument(
        "--singlemod-hek",
        required=True,
    )

    parser.add_argument(
        "--singlemod-hela",
        required=True,
    )

    parser.add_argument(
        "--output-dir",
        required=True,
    )

    parser.add_argument(
        "--chunk-size",
        type=int,
        default=500_000,
    )

    parser.add_argument(
        "--sample-size",
        type=int,
        default=300_000,
    )

    return parser.parse_args()


# ============================================================
# Utilities
# ============================================================

def detect_column(
    columns,
    candidates,
    label,
):

    lookup = {
        str(column).strip().lower(): column
        for column in columns
    }

    for candidate in candidates:

        key = candidate.lower()

        if key in lookup:
            return lookup[key]

    raise ValueError(
        f"Cannot detect {label} column.\n"
        f"Available columns: {columns}\n"
        f"Candidates: {candidates}"
    )


def normalize_strand(series):

    return (
        series
        .astype(str)
        .str.strip()
        .replace(
            {
                "1": "+",
                "1.0": "+",
                "plus": "+",
                "forward": "+",

                "-1": "-",
                "-1.0": "-",
                "minus": "-",
                "reverse": "-",
            }
        )
    )


# ============================================================
# Load GLORI pooled ground truth
# ============================================================

def load_glori(path):

    print("=" * 80)
    print("Loading GLORI")
    print(path)
    print("=" * 80)

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
            "glori_total_hek",
            "glori_total_hela",
        ],
        low_memory=False,
    )

    df["chr"] = (
        df["chr"]
        .astype(str)
        .str.strip()
    )

    df["strand"] = normalize_strand(
        df["strand"]
    )

    numeric_columns = [
        "pos_1based",
        "glori_beta_hek",
        "glori_beta_hela",
        "glori_total_hek",
        "glori_total_hela",
    ]

    for column in numeric_columns:

        df[column] = pd.to_numeric(
            df[column],
            errors="coerce",
        )

    valid_position = (
        df["pos_1based"].notna()
        & np.isfinite(df["pos_1based"])
    )

    df = df.loc[
        valid_position
        & df["strand"].isin(["+", "-"])
    ].copy()

    df["pos_1based"] = (
        np.rint(
            df["pos_1based"]
        )
        .astype(np.int64)
    )

    duplicated = df.duplicated(
        subset=KEY,
        keep=False,
    )

    if duplicated.any():

        raise ValueError(
            "GLORI contains duplicated genomic keys: "
            f"{int(duplicated.sum()):,} rows"
        )

    print(
        f"GLORI reference sites: "
        f"{len(df):,}"
    )

    return df


# ============================================================
# Determine tool columns
# ============================================================

def tool_columns(
    path,
    tool,
):

    header = pd.read_csv(
        path,
        sep="\t",
        nrows=0,
        compression="infer",
    )

    columns = list(
        header.columns
    )

    chr_col = detect_column(
        columns,
        [
            "chr",
            "chrom",
            "chromosome",
        ],
        f"{tool} chromosome",
    )

    pos_col = detect_column(
        columns,
        [
            "pos_1based",
            "position_1based",
            "genomic_position",
            "position",
            "pos",
        ],
        f"{tool} position",
    )

    strand_col = detect_column(
        columns,
        [
            "strand",
            "ref_strand",
        ],
        f"{tool} strand",
    )

    if tool == "dorado":

        beta_candidates = [
            "dorado_stoich",
            "stoich",
            "stoichiometry",
            "mod_ratio",
        ]

        depth_candidates = [
            "dorado_depth",
            "depth",
            "coverage",
            "count_all",
        ]

    elif tool == "m6anet":

        beta_candidates = [
            "m6anet_stoich",
            "mod_ratio",
            "stoich",
            "stoichiometry",
        ]

        depth_candidates = [
            "m6anet_depth",
            "n_reads",
            "depth",
            "coverage",
        ]

    elif tool == "singlemod":

        beta_candidates = [
            "singlemod_stoich",
            "stoich",
            "stoichiometry",
        ]

        depth_candidates = [
            "singlemod_depth",
            "depth",
            "coverage",
            "count_all",
        ]

    else:

        raise ValueError(
            f"Unsupported tool: {tool}"
        )

    beta_col = detect_column(
        columns,
        beta_candidates,
        f"{tool} beta",
    )

    depth_col = detect_column(
        columns,
        depth_candidates,
        f"{tool} depth",
    )

    print()
    print(
        f"{tool} input: {path}"
    )

    print(
        "Selected columns:"
    )

    print(
        f"  chr    = {chr_col}"
    )

    print(
        f"  pos    = {pos_col}"
    )

    print(
        f"  strand = {strand_col}"
    )

    print(
        f"  beta   = {beta_col}"
    )

    print(
        f"  depth  = {depth_col}"
    )

    return {
        "chr": chr_col,
        "pos": pos_col,
        "strand": strand_col,
        "beta": beta_col,
        "depth": depth_col,
    }


# ============================================================
# Match one tool / one cell line against GLORI
# ============================================================

def match_tool_to_glori(
    path,
    tool,
    cell,
    glori,
    chunk_size,
):

    columns = tool_columns(
        path,
        tool,
    )

    if cell == "HEK293T":

        gt_beta_column = (
            "glori_beta_hek"
        )

        gt_depth_column = (
            "glori_total_hek"
        )

    elif cell == "HeLa":

        gt_beta_column = (
            "glori_beta_hela"
        )

        gt_depth_column = (
            "glori_total_hela"
        )

    else:

        raise ValueError(
            f"Unknown cell line: {cell}"
        )

    # --------------------------------------------------------
    # Prepare GLORI table for this cell line
    # --------------------------------------------------------

    gt = glori[
        KEY
        + [
            gt_beta_column,
            gt_depth_column,
        ]
    ].copy()

    gt = gt.rename(
        columns={
            gt_beta_column:
                "glori_beta",

            gt_depth_column:
                "glori_depth",
        }
    )

    gt_valid = (
        gt["glori_beta"].notna()
        & np.isfinite(
            gt["glori_beta"]
        )
        & gt["glori_beta"].between(
            0,
            1,
        )
        & gt["glori_depth"].notna()
        & np.isfinite(
            gt["glori_depth"]
        )
        & (
            gt["glori_depth"]
            > 0
        )
    )

    gt = gt.loc[
        gt_valid
    ].copy()

    print()
    print("=" * 80)

    print(
        f"Matching {tool} {cell} "
        f"against GLORI {cell}"
    )

    print(
        f"Valid GLORI sites: "
        f"{len(gt):,}"
    )

    print("=" * 80)

    # --------------------------------------------------------
    # Read tool table in chunks
    # --------------------------------------------------------

    usecols = list(
        dict.fromkeys(
            [
                columns["chr"],
                columns["pos"],
                columns["strand"],
                columns["beta"],
                columns["depth"],
            ]
        )
    )

    matched_chunks = []

    total_input = 0
    total_valid = 0
    total_matched = 0

    for chunk in pd.read_csv(
        path,
        sep="\t",
        usecols=usecols,
        compression="infer",
        chunksize=chunk_size,
        low_memory=False,
    ):

        total_input += len(
            chunk
        )

        chunk = chunk.rename(
            columns={
                columns["chr"]:
                    "chr",

                columns["pos"]:
                    "pos_1based",

                columns["strand"]:
                    "strand",

                columns["beta"]:
                    "tool_beta",

                columns["depth"]:
                    "tool_depth",
            }
        )

        chunk["chr"] = (
            chunk["chr"]
            .astype(str)
            .str.strip()
        )

        chunk["strand"] = (
            normalize_strand(
                chunk["strand"]
            )
        )

        chunk["pos_1based"] = (
            pd.to_numeric(
                chunk["pos_1based"],
                errors="coerce",
            )
        )

        chunk["tool_beta"] = (
            pd.to_numeric(
                chunk["tool_beta"],
                errors="coerce",
            )
        )

        chunk["tool_depth"] = (
            pd.to_numeric(
                chunk["tool_depth"],
                errors="coerce",
            )
        )

        valid = (
            chunk["pos_1based"].notna()
            & np.isfinite(
                chunk["pos_1based"]
            )
            & chunk["strand"].isin(
                ["+", "-"]
            )
            & chunk["tool_beta"].notna()
            & np.isfinite(
                chunk["tool_beta"]
            )
            & chunk["tool_beta"].between(
                0,
                1,
            )
            & chunk["tool_depth"].notna()
            & np.isfinite(
                chunk["tool_depth"]
            )
            & (
                chunk["tool_depth"]
                > 0
            )
        )

        chunk = chunk.loc[
            valid
        ].copy()

        total_valid += len(
            chunk
        )

        if chunk.empty:
            continue

        chunk["pos_1based"] = (
            np.rint(
                chunk["pos_1based"]
            )
            .astype(np.int64)
        )

        merged = chunk.merge(
            gt,
            on=KEY,
            how="inner",
            validate="many_to_one",
            sort=False,
        )

        if not merged.empty:

            keep_columns = (
                KEY
                + [
                    "tool_beta",
                    "tool_depth",
                    "glori_beta",
                    "glori_depth",
                ]
            )

            matched_chunks.append(
                merged[
                    keep_columns
                ].copy()
            )

            total_matched += len(
                merged
            )

        print(
            f"  read={total_input:,}  "
            f"valid={total_valid:,}  "
            f"matched={total_matched:,}",
            flush=True,
        )

    if not matched_chunks:

        raise RuntimeError(
            f"No matched sites for "
            f"{tool} {cell}"
        )

    matched = pd.concat(
        matched_chunks,
        ignore_index=True,
    )

    del matched_chunks

    gc.collect()

    # --------------------------------------------------------
    # Final uniqueness check
    # --------------------------------------------------------

    duplicated = matched.duplicated(
        subset=KEY,
        keep=False,
    )

    if duplicated.any():

        examples = matched.loc[
            duplicated,
            KEY,
        ].head(10)

        raise ValueError(
            f"{tool} {cell}: "
            f"{int(duplicated.sum()):,} "
            "matched rows have duplicated genomic keys.\n"
            f"Examples:\n{examples}"
        )

    print()
    print(
        f"Final matched sites "
        f"{tool} {cell}: "
        f"{len(matched):,}"
    )

    return matched


# ============================================================
# Concordance correlation coefficient
# ============================================================

def concordance_correlation_coefficient(
    x,
    y,
):

    x = np.asarray(
        x,
        dtype=np.float64,
    )

    y = np.asarray(
        y,
        dtype=np.float64,
    )

    mean_x = np.mean(
        x
    )

    mean_y = np.mean(
        y
    )

    var_x = np.var(
        x
    )

    var_y = np.var(
        y
    )

    covariance = np.mean(
        (x - mean_x)
        * (y - mean_y)
    )

    denominator = (
        var_x
        + var_y
        + (
            mean_x
            - mean_y
        ) ** 2
    )

    if denominator == 0:
        return np.nan

    return (
        2
        * covariance
        / denominator
    )


# ============================================================
# Metrics
# ============================================================

def calculate_metrics(
    df,
    tool,
    cell,
    filter_name,
):

    if len(df) == 0:

        return {
            "tool": tool,
            "cell_line": cell,
            "depth_filter": filter_name,
            "n": 0,
            "pearson_r": np.nan,
            "spearman_rho": np.nan,
            "ccc": np.nan,
            "mae": np.nan,
            "rmse": np.nan,
            "mean_tool_minus_glori": np.nan,
            "median_tool_minus_glori": np.nan,
            "glori_mean_beta": np.nan,
            "tool_mean_beta": np.nan,
            "glori_zero_fraction": np.nan,
            "tool_zero_fraction": np.nan,
            "both_zero_fraction": np.nan,
            "glori_zero_tool_positive_fraction": np.nan,
            "tool_zero_glori_positive_fraction": np.nan,
        }

    x = df[
        "glori_beta"
    ].to_numpy(
        dtype=np.float64
    )

    y = df[
        "tool_beta"
    ].to_numpy(
        dtype=np.float64
    )

    difference = (
        y - x
    )

    if len(df) >= 2:

        pearson = float(
            pearsonr(
                x,
                y,
            ).statistic
        )

        spearman = float(
            spearmanr(
                x,
                y,
            ).statistic
        )

    else:

        pearson = np.nan
        spearman = np.nan

    return {
        "tool":
            tool,

        "cell_line":
            cell,

        "depth_filter":
            filter_name,

        "n":
            int(
                len(df)
            ),

        "pearson_r":
            pearson,

        "spearman_rho":
            spearman,

        "ccc":
            float(
                concordance_correlation_coefficient(
                    x,
                    y,
                )
            ),

        "mae":
            float(
                np.mean(
                    np.abs(
                        difference
                    )
                )
            ),

        "rmse":
            float(
                np.sqrt(
                    np.mean(
                        difference ** 2
                    )
                )
            ),

        "mean_tool_minus_glori":
            float(
                np.mean(
                    difference
                )
            ),

        "median_tool_minus_glori":
            float(
                np.median(
                    difference
                )
            ),

        "glori_mean_beta":
            float(
                np.mean(
                    x
                )
            ),

        "tool_mean_beta":
            float(
                np.mean(
                    y
                )
            ),

        "glori_zero_fraction":
            float(
                np.mean(
                    x == 0
                )
            ),

        "tool_zero_fraction":
            float(
                np.mean(
                    y == 0
                )
            ),

        "both_zero_fraction":
            float(
                np.mean(
                    (x == 0)
                    & (y == 0)
                )
            ),

        "glori_zero_tool_positive_fraction":
            float(
                np.mean(
                    (x == 0)
                    & (y > 0)
                )
            ),

        "tool_zero_glori_positive_fraction":
            float(
                np.mean(
                    (y == 0)
                    & (x > 0)
                )
            ),
    }


# ============================================================
# Statistics text
# ============================================================

def build_metrics_text(
    metrics,
):

    return (
        f"n = {metrics['n']:,}\n"
        f"Pearson r = "
        f"{metrics['pearson_r']:.3f}\n"
        f"Spearman ρ = "
        f"{metrics['spearman_rho']:.3f}\n"
        f"CCC = "
        f"{metrics['ccc']:.3f}\n"
        f"MAE = "
        f"{metrics['mae']:.3f}\n"
        f"RMSE = "
        f"{metrics['rmse']:.3f}"
    )


# ============================================================
# Single figure
# ============================================================

def plot_one(
    df,
    tool,
    cell,
    filter_name,
    metrics,
    output,
    sample_size,
):

    if len(df) == 0:
        return

    if len(df) > sample_size:

        plot_df = df.sample(
            n=sample_size,
            random_state=42,
        )

    else:

        plot_df = df

    fig, ax = plt.subplots(
        figsize=(7, 6)
    )

    # --------------------------------------------------------
    # Density layer
    # --------------------------------------------------------

    hb = ax.hexbin(
        plot_df["glori_beta"],
        plot_df["tool_beta"],
        gridsize=150,
        mincnt=1,
        bins="log",
        cmap=PAPER_DENSITY,
        zorder=1,
    )

    # --------------------------------------------------------
    # y = x
    # --------------------------------------------------------

    ax.plot(
        [0, 1],
        [0, 1],
        linestyle="--",
        linewidth=1.2,
        zorder=3,
    )

    ax.set_xlim(
        0,
        1,
    )

    ax.set_ylim(
        0,
        1,
    )

    ax.set_xlabel(
        f"GLORI β — {cell}"
    )

    ax.set_ylabel(
        f"{tool} PRJEB80229 β — {cell}"
    )

    ax.set_title(
        f"{tool} PRJEB80229 vs GLORI — {cell}\n"
        f"{filter_name}"
    )

    # --------------------------------------------------------
    # Statistics box
    #
    # White opaque background + high zorder makes the
    # box cover the points underneath.
    # --------------------------------------------------------

    text = build_metrics_text(
        metrics
    )

    ax.text(
        0.03,
        0.97,
        text,
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=10,
        zorder=20,
        bbox={
            "boxstyle":
                "round,pad=0.35",

            "facecolor":
                "white",

            "edgecolor":
                "#777777",

            "linewidth":
                0.8,

            "alpha":
                0.97,
        },
    )

    # --------------------------------------------------------
    # Colour bar
    # --------------------------------------------------------

    colorbar = fig.colorbar(
        hb,
        ax=ax,
    )

    colorbar.set_label(
        "Displayed-site count (log scale)"
    )

    fig.tight_layout()

    fig.savefig(
        output,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(
        fig
    )


# ============================================================
# 2 x 3 panel
# ============================================================

def make_panel(
    results,
    filter_name,
    output,
    sample_size,
):

    fig, axes = plt.subplots(
        2,
        3,
        figsize=(18, 11),
    )

    for row_index, cell in enumerate(
        CELLS
    ):

        for column_index, tool in enumerate(
            TOOLS
        ):

            ax = axes[
                row_index,
                column_index,
            ]

            df, metrics = results[
                (
                    tool,
                    cell,
                )
            ]

            if len(df) == 0:

                ax.axis(
                    "off"
                )

                continue

            if len(df) > sample_size:

                plot_df = df.sample(
                    n=sample_size,
                    random_state=42,
                )

            else:

                plot_df = df

            # ------------------------------------------------
            # Density
            # ------------------------------------------------

            ax.hexbin(
                plot_df["glori_beta"],
                plot_df["tool_beta"],
                gridsize=120,
                mincnt=1,
                bins="log",
                cmap=PAPER_DENSITY,
                zorder=1,
            )

            # ------------------------------------------------
            # y = x
            # ------------------------------------------------

            ax.plot(
                [0, 1],
                [0, 1],
                linestyle="--",
                linewidth=1.1,
                zorder=3,
            )

            ax.set_xlim(
                0,
                1,
            )

            ax.set_ylim(
                0,
                1,
            )

            ax.set_xlabel(
                f"GLORI β — {cell}"
            )

            ax.set_ylabel(
                f"{tool} β — {cell}"
            )

            ax.set_title(
                f"{tool} — {cell}"
            )

            text = (
                f"n = {metrics['n']:,}\n"
                f"Pearson r = "
                f"{metrics['pearson_r']:.3f}\n"
                f"Spearman ρ = "
                f"{metrics['spearman_rho']:.3f}\n"
                f"MAE = "
                f"{metrics['mae']:.3f}"
            )

            # ------------------------------------------------
            # Statistics box covers data points
            # ------------------------------------------------

            ax.text(
                0.03,
                0.97,
                text,
                transform=ax.transAxes,
                va="top",
                ha="left",
                fontsize=10,
                zorder=20,
                bbox={
                    "boxstyle":
                        "round,pad=0.35",

                    "facecolor":
                        "white",

                    "edgecolor":
                        "#777777",

                    "linewidth":
                        0.8,

                    "alpha":
                        0.97,
                },
            )

    fig.suptitle(
        "PRJEB80229 absolute β vs GLORI "
        f"— {filter_name}",
        fontsize=17,
        y=0.995,
    )

    fig.tight_layout(
        rect=[
            0,
            0,
            1,
            0.98,
        ]
    )

    fig.savefig(
        output,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(
        fig
    )


# ============================================================
# Main
# ============================================================

def main():

    args = parse_args()

    output_dir = Path(
        args.output_dir
    )

    output_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    if any(
        output_dir.iterdir()
    ):

        raise RuntimeError(
            "Output directory is not empty: "
            f"{output_dir}"
        )

    # --------------------------------------------------------
    # GLORI
    # --------------------------------------------------------

    glori = load_glori(
        args.glori
    )

    # --------------------------------------------------------
    # Input paths
    # --------------------------------------------------------

    input_paths = {
        (
            "dorado",
            "HEK293T",
        ):
            args.dorado_hek,

        (
            "dorado",
            "HeLa",
        ):
            args.dorado_hela,

        (
            "m6anet",
            "HEK293T",
        ):
            args.m6anet_hek,

        (
            "m6anet",
            "HeLa",
        ):
            args.m6anet_hela,

        (
            "singlemod",
            "HEK293T",
        ):
            args.singlemod_hek,

        (
            "singlemod",
            "HeLa",
        ):
            args.singlemod_hela,
    }

    # --------------------------------------------------------
    # Match every tool / cell line against GLORI
    # --------------------------------------------------------

    matched_tables = {}

    for tool in TOOLS:

        for cell in CELLS:

            path = input_paths[
                (
                    tool,
                    cell,
                )
            ]

            matched = match_tool_to_glori(
                path=path,
                tool=tool,
                cell=cell,
                glori=glori,
                chunk_size=args.chunk_size,
            )

            matched_tables[
                (
                    tool,
                    cell,
                )
            ] = matched

            matched_output = (
                output_dir
                / (
                    f"{tool}.{cell}."
                    "80229_vs_GLORI."
                    "matched.tsv.gz"
                )
            )

            matched.to_csv(
                matched_output,
                sep="\t",
                index=False,
                compression={
                    "method":
                        "gzip",

                    "compresslevel":
                        1,
                },
                float_format="%.8g",
            )

            print(
                f"Matched table written: "
                f"{matched_output}",
                flush=True,
            )

            gc.collect()

    # --------------------------------------------------------
    # Depth thresholds
    # --------------------------------------------------------

    thresholds = [
        (
            "all_common_sites",
            None,
        ),
        (
            "depth_ge_10",
            10,
        ),
        (
            "depth_ge_20",
            20,
        ),
        (
            "depth_ge_40",
            40,
        ),
    ]

    all_metrics = []

    # --------------------------------------------------------
    # Analyse each threshold
    # --------------------------------------------------------

    for (
        filter_name,
        threshold,
    ) in thresholds:

        print()
        print("#" * 80)
        print(
            f"Depth filter: {filter_name}"
        )
        print("#" * 80)

        panel_results = {}

        for tool in TOOLS:

            for cell in CELLS:

                matched = matched_tables[
                    (
                        tool,
                        cell,
                    )
                ]

                if threshold is None:

                    filtered = matched

                else:

                    filtered = matched.loc[
                        (
                            matched["tool_depth"]
                            >= threshold
                        )
                        & (
                            matched["glori_depth"]
                            >= threshold
                        )
                    ]

                metrics = calculate_metrics(
                    df=filtered,
                    tool=tool,
                    cell=cell,
                    filter_name=filter_name,
                )

                all_metrics.append(
                    metrics
                )

                panel_results[
                    (
                        tool,
                        cell,
                    )
                ] = (
                    filtered,
                    metrics,
                )

                print(
                    f"{tool:10s} "
                    f"{cell:8s} "
                    f"{filter_name:18s} "
                    f"n={metrics['n']:,} "
                    f"Pearson="
                    f"{metrics['pearson_r']:.4f} "
                    f"Spearman="
                    f"{metrics['spearman_rho']:.4f} "
                    f"CCC="
                    f"{metrics['ccc']:.4f} "
                    f"MAE="
                    f"{metrics['mae']:.4f} "
                    f"RMSE="
                    f"{metrics['rmse']:.4f}",
                    flush=True,
                )

                # --------------------------------------------
                # Individual plot
                # --------------------------------------------

                individual_output = (
                    output_dir
                    / (
                        f"{tool}."
                        f"{cell}."
                        f"{filter_name}.png"
                    )
                )

                plot_one(
                    df=filtered,
                    tool=tool,
                    cell=cell,
                    filter_name=filter_name,
                    metrics=metrics,
                    output=individual_output,
                    sample_size=args.sample_size,
                )

        # ----------------------------------------------------
        # Combined 2 x 3 panel
        # ----------------------------------------------------

        panel_output = (
            output_dir
            / (
                "panel.absolute_beta."
                f"{filter_name}.png"
            )
        )

        make_panel(
            results=panel_results,
            filter_name=filter_name,
            output=panel_output,
            sample_size=args.sample_size,
        )

    # --------------------------------------------------------
    # Metrics table
    # --------------------------------------------------------

    metrics_df = pd.DataFrame(
        all_metrics
    )

    metrics_output = (
        output_dir
        / (
            "80229_absolute_beta_vs_GLORI."
            "metrics.tsv"
        )
    )

    metrics_df.to_csv(
        metrics_output,
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

    print()
    print(
        f"Metrics table: "
        f"{metrics_output}"
    )

    print(
        f"Output directory: "
        f"{output_dir}"
    )


if __name__ == "__main__":
    main()
