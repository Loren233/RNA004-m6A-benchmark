
import os as _os
from pathlib import Path as _Path
RNA004_STAT_ROOT = _os.environ.get(
    "RNA004_STAT_ROOT",
    str(_Path(__file__).resolve().parents[3])
)

#!/usr/bin/env python3

from pathlib import Path
import warnings

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


# ============================================================
# 0. PATHS
# ============================================================

STAT802 = Path(
    f"{RNA004_STAT_ROOT}/80229"
)

WORK = Path(
    f"{RNA004_STAT_ROOT}/"
    "tool_vs_glori/80229_highconf_sensitivity"
)

RESULT_DIR = WORK / "results"
FIG_DIR = WORK / "figures"

FIG_DIR.mkdir(
    parents=True,
    exist_ok=True,
)

GLORI_FILE = (
    RESULT_DIR /
    "GLORI_HEK_HeLa.highconf.tsv.gz"
)

ABS_METRICS = (
    RESULT_DIR /
    "absolute_beta_strictness_metrics.tsv"
)

DELTA_METRICS = (
    RESULT_DIR /
    "delta_beta_strictness_metrics.tsv"
)


# ============================================================
# 1. TOOL FILES
# ============================================================

FILES = {

    "Dorado": {

        "HEK293T":
            STAT802 /
            "dorado_Hek293T.sitelevel.genomic.standardized.collapsed.tsv",

        "HeLa":
            STAT802 /
            "dorado_HeLa.sitelevel.genomic.standardized.recollapsed.tsv",
    },

    "m6Anet": {

        "HEK293T":
            STAT802 /
            "m6anet_Hek293T.sitelevel.genomic.standardized.collapsed.tsv",

        "HeLa":
            STAT802 /
            "m6anet_HeLa.sitelevel.genomic.standardized.collapsed.tsv",
    },

    "SingleMod": {

        "HEK293T":
            STAT802 /
            "singlemod_Hek293T.sitelevel.genomic.standardized.collapsed.tsv",

        "HeLa":
            STAT802 /
            "singlemod_HeLa.sitelevel.genomic.standardized.collapsed.tsv",
    },
}


# ============================================================
# 2. PARAMETERS
# ============================================================

KEY = [
    "chr",
    "pos_1based",
    "strand",
]

COVERAGE_THRESHOLD = 20

CHUNK_SIZE = 1_000_000

MAX_PLOT_POINTS = 600_000

RNG = np.random.default_rng(
    20260811
)


# ============================================================
# 3. S0-S7 DEFINITIONS
# ============================================================

STAGES = [

    "raw_common",
    "cov20",
    "p005",
    "p001",
    "p0001",
    "beta010",
    "beta020",
    "beta050",
]


STAGE_PREFIX = {

    "raw_common":
        "S0",

    "cov20":
        "S1",

    "p005":
        "S2",

    "p001":
        "S3",

    "p0001":
        "S4",

    "beta010":
        "S5",

    "beta020":
        "S6",

    "beta050":
        "S7",
}


STAGE_SHORT_LABEL = {

    "raw_common":
        "Raw",

    "cov20":
        "Coverage >20",

    "p005":
        "AdjP <0.05",

    "p001":
        "AdjP <0.01",

    "p0001":
        "AdjP <0.001",

    "beta010":
        "β >0.10",

    "beta020":
        "β >0.20",

    "beta050":
        "β >0.50",
}


STAGE_ABS_TITLE = {

    "raw_common":
        "All common sites",

    "cov20":
        "Tool coverage >20; GLORI coverage >20",

    "p005":
        "Coverage >20; both GLORI replicates AdjP <0.05",

    "p001":
        "Coverage >20; both GLORI replicates AdjP <0.01",

    "p0001":
        "Coverage >20; both GLORI replicates AdjP <0.001",

    "beta010":
        "S4 + GLORI β >0.10",

    "beta020":
        "S4 + GLORI β >0.20",

    "beta050":
        "S4 + GLORI β >0.50",
}


STAGE_DELTA_TITLE = {

    "raw_common":
        "All common HEK293T-HeLa sites",

    "cov20":
        "Coverage >20 in tool and GLORI for both cell lines",

    "p005":
        "S1 + GLORI AdjP <0.05 in both replicates of ≥1 cell line",

    "p001":
        "S1 + GLORI AdjP <0.01 in both replicates of ≥1 cell line",

    "p0001":
        "S1 + GLORI AdjP <0.001 in both replicates of ≥1 cell line",

    "beta010":
        "S4 + max(GLORI β_HEK, β_HeLa) >0.10",

    "beta020":
        "S4 + max(GLORI β_HEK, β_HeLa) >0.20",

    "beta050":
        "S4 + max(GLORI β_HEK, β_HeLa) >0.50",
}


# ============================================================
# 4. PLOT STYLE
# ============================================================

paper_density = [
    "#3b0f70",
    "#8c2981",
    "#de4968",
    "#f98e52",
    "#fee08b",
]

density_cmap = (
    LinearSegmentedColormap.from_list(
        "paper_density",
        paper_density,
    )
)


# ============================================================
# 5. FILE HELPERS
# ============================================================

def resolve_file(path):

    if path.exists():
        return path

    gz = Path(
        str(path) + ".gz"
    )

    if gz.exists():
        return gz

    raise FileNotFoundError(
        f"Cannot find:\n{path}\n"
        f"or:\n{gz}"
    )


def find_column(
    columns,
    candidates,
    label,
):

    for name in candidates:

        if name in columns:
            return name

    raise KeyError(
        f"Cannot identify {label}.\n"
        f"Candidates: {candidates}\n"
        f"Columns: {list(columns)}"
    )


def resolve_tool_columns(path):

    header = pd.read_csv(
        path,
        sep="\t",
        nrows=0,
    )

    columns = list(
        header.columns
    )

    for key in KEY:

        if key not in columns:

            raise KeyError(
                f"{key} missing from:\n"
                f"{path}"
            )


    depth_col = find_column(
        columns,
        [
            "depth",
            "dorado_depth",
            "m6anet_depth",
            "singlemod_depth",
            "coverage",
            "count_all",
            "n_reads",
        ],
        "depth",
    )


    stoich_col = find_column(
        columns,
        [
            "stoich",
            "dorado_stoich",
            "m6anet_stoich",
            "singlemod_stoich",
            "stoichiometry",
            "mod_ratio",
            "beta",
        ],
        "stoichiometry",
    )


    return (
        depth_col,
        stoich_col,
    )


# ============================================================
# 6. METRICS
# ============================================================

def concordance_ccc(
    x,
    y,
):

    x = np.asarray(
        x,
        dtype=float,
    )

    y = np.asarray(
        y,
        dtype=float,
    )


    if len(x) < 2:
        return np.nan


    vx = np.var(
        x,
        ddof=1,
    )

    vy = np.var(
        y,
        ddof=1,
    )

    mx = np.mean(x)
    my = np.mean(y)

    cov = np.cov(
        x,
        y,
        ddof=1,
    )[0, 1]


    denominator = (
        vx +
        vy +
        (mx - my) ** 2
    )


    if denominator == 0:
        return np.nan


    return (
        2 * cov /
        denominator
    )


def calculate_metrics(
    x,
    y,
):

    x = np.asarray(
        x,
        dtype=float,
    )

    y = np.asarray(
        y,
        dtype=float,
    )


    valid = (
        np.isfinite(x) &
        np.isfinite(y)
    )


    x = x[valid]
    y = y[valid]


    n = len(x)


    if n < 3:

        return {

            "n":
                n,

            "pearson":
                np.nan,

            "spearman":
                np.nan,

            "ccc":
                np.nan,

            "mae":
                np.nan,

            "rmse":
                np.nan,
        }


    with warnings.catch_warnings():

        warnings.simplefilter(
            "ignore"
        )

        pearson = pearsonr(
            x,
            y,
        ).statistic

        spearman = spearmanr(
            x,
            y,
        ).statistic


    error = (
        y - x
    )


    return {

        "n":
            n,

        "pearson":
            pearson,

        "spearman":
            spearman,

        "ccc":
            concordance_ccc(
                x,
                y,
            ),

        "mae":
            np.mean(
                np.abs(error)
            ),

        "rmse":
            np.sqrt(
                np.mean(
                    error ** 2
                )
            ),
    }


def sample_for_plot(df):

    if len(df) <= MAX_PLOT_POINTS:

        return df


    idx = RNG.choice(
        len(df),
        size=MAX_PLOT_POINTS,
        replace=False,
    )


    return df.iloc[
        idx
    ]


def format_stats(stat):

    return (
        f"n = {stat['n']:,}\n"
        f"Pearson = {stat['pearson']:.3f}\n"
        f"Spearman = {stat['spearman']:.3f}\n"
        f"MAE = {stat['mae']:.3f}\n"
        f"RMSE = {stat['rmse']:.3f}"
    )


# ============================================================
# 7. LOAD GLORI
# ============================================================

print(
    "=" * 80
)

print(
    "Loading GLORI..."
)

print(
    "=" * 80
)


gt = pd.read_csv(
    GLORI_FILE,
    sep="\t",
    low_memory=False,
)


required_gt = {

    "chr",
    "pos_1based",
    "strand",

    "glori_total_hek",
    "glori_beta_hek",

    "glori_total_hela",
    "glori_beta_hela",

    "glori_adj_p_hek_rep1",
    "glori_adj_p_hek_rep2",

    "glori_adj_p_hela_rep1",
    "glori_adj_p_hela_rep2",
}


missing = (
    required_gt -
    set(gt.columns)
)


if missing:

    raise KeyError(
        f"Missing GLORI columns:\n"
        f"{sorted(missing)}"
    )


gt["chr"] = (
    gt["chr"]
    .astype(str)
)

gt["strand"] = (
    gt["strand"]
    .astype(str)
)

gt["pos_1based"] = (
    pd.to_numeric(
        gt["pos_1based"],
        errors="raise",
    )
    .astype(np.int64)
)


if gt.duplicated(KEY).any():

    raise RuntimeError(
        "GLORI genomic keys are not unique."
    )


print(
    "GLORI rows:",
    f"{len(gt):,}"
)


# ============================================================
# 8. LOAD TOOL × CELL
# ============================================================

def load_tool_cell(
    tool,
    cell,
    path,
):

    path = resolve_file(
        path
    )


    (
        depth_col,
        stoich_col,
    ) = resolve_tool_columns(
        path
    )


    print()

    print(
        "=" * 80
    )

    print(
        f"{tool} -- {cell}"
    )

    print(
        path
    )

    print(
        f"depth  = {depth_col}"
    )

    print(
        f"stoich = {stoich_col}"
    )

    print(
        "=" * 80
    )


    if cell == "HEK293T":

        gt_cols = (
            KEY +
            [
                "glori_total_hek",
                "glori_beta_hek",

                "glori_adj_p_hek_rep1",
                "glori_adj_p_hek_rep2",
            ]
        )


    elif cell == "HeLa":

        gt_cols = (
            KEY +
            [
                "glori_total_hela",
                "glori_beta_hela",

                "glori_adj_p_hela_rep1",
                "glori_adj_p_hela_rep2",
            ]
        )


    else:

        raise ValueError(
            cell
        )


    gt_sub = (
        gt[
            gt_cols
        ]
    )


    usecols = (
        KEY +
        [
            depth_col,
            stoich_col,
        ]
    )


    pieces = []


    reader = pd.read_csv(
        path,
        sep="\t",
        usecols=usecols,
        chunksize=CHUNK_SIZE,
        low_memory=False,
    )


    for chunk_id, chunk in enumerate(
        reader,
        start=1,
    ):

        input_n = len(
            chunk
        )


        chunk = chunk.rename(
            columns={

                depth_col:
                    "tool_depth",

                stoich_col:
                    "tool_beta",
            }
        )


        chunk["chr"] = (
            chunk["chr"]
            .astype(str)
        )

        chunk["strand"] = (
            chunk["strand"]
            .astype(str)
        )


        chunk["pos_1based"] = (
            pd.to_numeric(
                chunk["pos_1based"],
                errors="coerce",
            )
        )


        chunk["tool_depth"] = (
            pd.to_numeric(
                chunk["tool_depth"],
                errors="coerce",
            )
        )


        chunk["tool_beta"] = (
            pd.to_numeric(
                chunk["tool_beta"],
                errors="coerce",
            )
        )


        chunk = chunk.dropna(
            subset=[
                "pos_1based",
                "tool_depth",
                "tool_beta",
            ]
        )


        chunk["pos_1based"] = (
            chunk["pos_1based"]
            .astype(np.int64)
        )


        merged = chunk.merge(
            gt_sub,
            on=KEY,
            how="inner",
            sort=False,
            validate="many_to_one",
        )


        if len(merged):

            pieces.append(
                merged
            )


        print(
            f"chunk {chunk_id:03d}: "
            f"input={input_n:,}, "
            f"matched={len(merged):,}",
            flush=True,
        )


    if not pieces:

        raise RuntimeError(
            f"No matched sites for "
            f"{tool} {cell}"
        )


    out = pd.concat(
        pieces,
        ignore_index=True,
    )


    dup = (
        out.duplicated(KEY)
        .sum()
    )


    if dup:

        raise RuntimeError(
            f"{tool} {cell}: "
            f"{dup:,} duplicated genomic keys."
        )


    print(
        f"TOTAL MATCHED: "
        f"{len(out):,}"
    )


    return out


# ============================================================
# 9. LOAD ALL TOOL DATA
# ============================================================

matched = {}


for tool, cell_files in FILES.items():

    matched[tool] = {}


    for cell, path in cell_files.items():

        matched[tool][cell] = (
            load_tool_cell(
                tool,
                cell,
                path,
            )
        )


# ============================================================
# 10. SIGNIFICANCE FUNCTIONS
# ============================================================

def absolute_sig(
    df,
    cell,
    cutoff,
):

    if cell == "HEK293T":

        p1 = (
            df[
                "glori_adj_p_hek_rep1"
            ]
        )

        p2 = (
            df[
                "glori_adj_p_hek_rep2"
            ]
        )


    elif cell == "HeLa":

        p1 = (
            df[
                "glori_adj_p_hela_rep1"
            ]
        )

        p2 = (
            df[
                "glori_adj_p_hela_rep2"
            ]
        )


    else:

        raise ValueError(
            cell
        )


    return (

        np.isfinite(p1)
        &
        np.isfinite(p2)
        &
        (p1 < cutoff)
        &
        (p2 < cutoff)
    )


def delta_sig(
    df,
    cutoff,
):

    hek_sig = (

        np.isfinite(
            df[
                "glori_adj_p_hek_rep1"
            ]
        )

        &
        np.isfinite(
            df[
                "glori_adj_p_hek_rep2"
            ]
        )

        &
        (
            df[
                "glori_adj_p_hek_rep1"
            ]
            <
            cutoff
        )

        &
        (
            df[
                "glori_adj_p_hek_rep2"
            ]
            <
            cutoff
        )
    )


    hela_sig = (

        np.isfinite(
            df[
                "glori_adj_p_hela_rep1"
            ]
        )

        &
        np.isfinite(
            df[
                "glori_adj_p_hela_rep2"
            ]
        )

        &
        (
            df[
                "glori_adj_p_hela_rep1"
            ]
            <
            cutoff
        )

        &
        (
            df[
                "glori_adj_p_hela_rep2"
            ]
            <
            cutoff
        )
    )


    return (
        hek_sig |
        hela_sig
    )


# ============================================================
# 11. ABSOLUTE MASKS
# ============================================================

def absolute_masks(
    df,
    cell,
):

    if cell == "HEK293T":

        gt_beta_col = (
            "glori_beta_hek"
        )

        gt_depth_col = (
            "glori_total_hek"
        )


    else:

        gt_beta_col = (
            "glori_beta_hela"
        )

        gt_depth_col = (
            "glori_total_hela"
        )


    raw = (

        np.isfinite(
            df["tool_beta"]
        )

        &
        np.isfinite(
            df[gt_beta_col]
        )
    )


    cov20 = (

        raw

        &
        (
            df["tool_depth"] >
            COVERAGE_THRESHOLD
        )

        &
        (
            df[gt_depth_col] >
            COVERAGE_THRESHOLD
        )
    )


    p005 = (

        cov20

        &
        absolute_sig(
            df,
            cell,
            0.05,
        )
    )


    p001 = (

        cov20

        &
        absolute_sig(
            df,
            cell,
            0.01,
        )
    )


    p0001 = (

        cov20

        &
        absolute_sig(
            df,
            cell,
            0.001,
        )
    )


    beta010 = (

        p0001

        &
        (
            df[gt_beta_col] >
            0.10
        )
    )


    beta020 = (

        p0001

        &
        (
            df[gt_beta_col] >
            0.20
        )
    )


    beta050 = (

        p0001

        &
        (
            df[gt_beta_col] >
            0.50
        )
    )


    masks = {

        "raw_common":
            raw,

        "cov20":
            cov20,

        "p005":
            p005,

        "p001":
            p001,

        "p0001":
            p0001,

        "beta010":
            beta010,

        "beta020":
            beta020,

        "beta050":
            beta050,
    }


    return (
        gt_beta_col,
        masks,
    )


# ============================================================
# 12. BUILD DELTA TABLES
# ============================================================

delta_tables = {}


for tool in FILES:

    hek = (
        matched[tool]["HEK293T"]
        [
            KEY +
            [
                "tool_beta",
                "tool_depth",

                "glori_beta_hek",
                "glori_total_hek",

                "glori_adj_p_hek_rep1",
                "glori_adj_p_hek_rep2",
            ]
        ]
        .copy()
    )


    hek = hek.rename(
        columns={

            "tool_beta":
                "tool_beta_hek",

            "tool_depth":
                "tool_depth_hek",
        }
    )


    hela = (
        matched[tool]["HeLa"]
        [
            KEY +
            [
                "tool_beta",
                "tool_depth",

                "glori_beta_hela",
                "glori_total_hela",

                "glori_adj_p_hela_rep1",
                "glori_adj_p_hela_rep2",
            ]
        ]
        .copy()
    )


    hela = hela.rename(
        columns={

            "tool_beta":
                "tool_beta_hela",

            "tool_depth":
                "tool_depth_hela",
        }
    )


    d = hek.merge(
        hela,
        on=KEY,
        how="inner",
        sort=False,
        validate="one_to_one",
    )


    d["tool_delta_beta"] = (

        d["tool_beta_hek"]
        -
        d["tool_beta_hela"]
    )


    d["glori_delta_beta"] = (

        d["glori_beta_hek"]
        -
        d["glori_beta_hela"]
    )


    d["glori_max_stoich"] = np.maximum(

        d["glori_beta_hek"],
        d["glori_beta_hela"],
    )


    delta_tables[tool] = d


# ============================================================
# 13. DELTA MASKS
# ============================================================

def delta_masks(d):

    raw = (

        np.isfinite(
            d["tool_delta_beta"]
        )

        &
        np.isfinite(
            d["glori_delta_beta"]
        )
    )


    cov20 = (

        raw

        &
        (
            d["tool_depth_hek"] >
            COVERAGE_THRESHOLD
        )

        &
        (
            d["tool_depth_hela"] >
            COVERAGE_THRESHOLD
        )

        &
        (
            d["glori_total_hek"] >
            COVERAGE_THRESHOLD
        )

        &
        (
            d["glori_total_hela"] >
            COVERAGE_THRESHOLD
        )
    )


    p005 = (

        cov20

        &
        delta_sig(
            d,
            0.05,
        )
    )


    p001 = (

        cov20

        &
        delta_sig(
            d,
            0.01,
        )
    )


    p0001 = (

        cov20

        &
        delta_sig(
            d,
            0.001,
        )
    )


    beta010 = (

        p0001

        &
        (
            d["glori_max_stoich"] >
            0.10
        )
    )


    beta020 = (

        p0001

        &
        (
            d["glori_max_stoich"] >
            0.20
        )
    )


    beta050 = (

        p0001

        &
        (
            d["glori_max_stoich"] >
            0.50
        )
    )


    return {

        "raw_common":
            raw,

        "cov20":
            cov20,

        "p005":
            p005,

        "p001":
            p001,

        "p0001":
            p0001,

        "beta010":
            beta010,

        "beta020":
            beta020,

        "beta050":
            beta050,
    }


# ============================================================
# 14. ABSOLUTE BETA S0-S7 SCATTERS
# ============================================================

def plot_absolute_stage(stage):

    prefix = (
        STAGE_PREFIX[stage]
    )


    fig, axes = plt.subplots(
        2,
        3,
        figsize=(15, 9),
        constrained_layout=True,
    )


    for row, cell in enumerate(
        [
            "HEK293T",
            "HeLa",
        ]
    ):

        for col, tool in enumerate(
            [
                "Dorado",
                "m6Anet",
                "SingleMod",
            ]
        ):

            ax = axes[
                row,
                col,
            ]


            df = (
                matched[tool][cell]
            )


            (
                gt_beta_col,
                masks,
            ) = absolute_masks(
                df,
                cell,
            )


            mask = (
                masks[stage]
            )


            sub = (
                df.loc[
                    mask,
                    [
                        gt_beta_col,
                        "tool_beta",
                    ]
                ]
                .dropna()
            )


            stat = calculate_metrics(

                sub[
                    gt_beta_col
                ],

                sub[
                    "tool_beta"
                ],
            )


            plot_sub = (
                sample_for_plot(
                    sub
                )
            )


            if len(plot_sub):

                ax.hexbin(

                    plot_sub[
                        gt_beta_col
                    ],

                    plot_sub[
                        "tool_beta"
                    ],

                    gridsize=75,
                    mincnt=1,
                    bins="log",
                    cmap=density_cmap,
                )


            ax.plot(
                [0, 1],
                [0, 1],
                "--",
                linewidth=1.0,
                alpha=0.8,
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
                "GLORI β"
            )

            ax.set_ylabel(
                f"{tool} β"
            )


            ax.set_title(
                f"{cell} — {tool}"
            )


            ax.text(

                0.03,
                0.97,

                format_stats(
                    stat
                ),

                transform=ax.transAxes,

                ha="left",
                va="top",

                bbox=dict(
                    boxstyle="round,pad=0.35",
                    facecolor="white",
                    edgecolor="0.6",
                    alpha=0.97,
                ),

                zorder=20,
            )


    fig.suptitle(
        f"{prefix} — Absolute β vs GLORI\n"
        f"{STAGE_ABS_TITLE[stage]}",
        fontsize=18,
    )


    output = (
        FIG_DIR /
        f"{prefix}_absolute_beta_{stage}.png"
    )


    fig.savefig(
        output,
        dpi=300,
    )

    plt.close(
        fig
    )


    print(
        "Written:",
        output
    )


# ============================================================
# 15. DELTA BETA S0-S7 SCATTERS
# ============================================================

def plot_delta_stage(stage):

    prefix = (
        STAGE_PREFIX[stage]
    )


    fig, axes = plt.subplots(
        1,
        3,
        figsize=(15, 5),
        constrained_layout=True,
    )


    for ax, tool in zip(
        axes,
        [
            "Dorado",
            "m6Anet",
            "SingleMod",
        ]
    ):

        d = (
            delta_tables[tool]
        )


        masks = (
            delta_masks(
                d
            )
        )


        mask = (
            masks[stage]
        )


        sub = (
            d.loc[
                mask,
                [
                    "glori_delta_beta",
                    "tool_delta_beta",
                ]
            ]
            .dropna()
        )


        stat = calculate_metrics(

            sub[
                "glori_delta_beta"
            ],

            sub[
                "tool_delta_beta"
            ],
        )


        plot_sub = (
            sample_for_plot(
                sub
            )
        )


        if len(plot_sub):

            ax.hexbin(

                plot_sub[
                    "glori_delta_beta"
                ],

                plot_sub[
                    "tool_delta_beta"
                ],

                gridsize=70,
                mincnt=1,
                bins="log",
                cmap=density_cmap,
            )


        ax.plot(
            [-1, 1],
            [-1, 1],
            "--",
            linewidth=1.0,
            alpha=0.8,
        )


        ax.axhline(
            0,
            linewidth=0.7,
            alpha=0.5,
        )

        ax.axvline(
            0,
            linewidth=0.7,
            alpha=0.5,
        )


        ax.set_xlim(
            -1,
            1,
        )

        ax.set_ylim(
            -1,
            1,
        )


        ax.set_xlabel(
            "GLORI Δβ (HEK − HeLa)"
        )

        ax.set_ylabel(
            f"{tool} Δβ"
        )

        ax.set_title(
            tool
        )


        ax.text(

            0.03,
            0.97,

            format_stats(
                stat
            ),

            transform=ax.transAxes,

            ha="left",
            va="top",

            bbox=dict(
                boxstyle="round,pad=0.35",
                facecolor="white",
                edgecolor="0.6",
                alpha=0.97,
            ),

            zorder=20,
        )


    fig.suptitle(
        f"{prefix} — Δβ vs GLORI\n"
        f"{STAGE_DELTA_TITLE[stage]}",
        fontsize=18,
    )


    output = (
        FIG_DIR /
        f"{prefix}_delta_beta_{stage}.png"
    )


    fig.savefig(
        output,
        dpi=300,
    )

    plt.close(
        fig
    )


    print(
        "Written:",
        output
    )


# ============================================================
# 16. GENERATE ALL S0-S7 SCATTERS
# ============================================================

print()
print(
    "=" * 80
)

print(
    "Generating S0-S7 scatter plots..."
)

print(
    "=" * 80
)


for stage in STAGES:

    plot_absolute_stage(
        stage
    )

    plot_delta_stage(
        stage
    )


# ============================================================
# 17. LOAD SUMMARY TABLES
# ============================================================

absolute_summary = pd.read_csv(
    ABS_METRICS,
    sep="\t",
)

delta_summary = pd.read_csv(
    DELTA_METRICS,
    sep="\t",
)


# ============================================================
# 18. SUMMARY PLOT FUNCTION
# ============================================================

def summary_plot(
    df,
    metric,
    title,
    ylabel,
    output,
    ylim=None,
    log_y=False,
):

    fig, ax = plt.subplots(
        figsize=(10, 6)
    )


    for tool in [
        "Dorado",
        "m6Anet",
        "SingleMod",
    ]:

        t = (
            df[
                df["tool"] ==
                tool
            ]
            .set_index(
                "stage"
            )
            .reindex(
                STAGES
            )
        )


        ax.plot(
            [
                f"{STAGE_PREFIX[s]}\n"
                f"{STAGE_SHORT_LABEL[s]}"
                for s in STAGES
            ],
            t[metric],
            marker="o",
            linewidth=2,
            label=tool,
        )


    if ylim is not None:

        ax.set_ylim(
            *ylim
        )


    if log_y:

        ax.set_yscale(
            "log"
        )


    ax.set_xlabel(
        "Increasing filter strictness"
    )

    ax.set_ylabel(
        ylabel
    )

    ax.set_title(
        title
    )

    ax.legend()

    fig.tight_layout()


    fig.savefig(
        output,
        dpi=300,
    )

    plt.close(
        fig
    )


    print(
        "Written:",
        output
    )


# ============================================================
# 19. ABSOLUTE BETA SUMMARY FIGURES
# ============================================================

for cell in [
    "HEK293T",
    "HeLa",
]:

    subset = (
        absolute_summary[
            absolute_summary[
                "cell"
            ]
            ==
            cell
        ]
    )


    summary_plot(

        subset,

        metric="pearson",

        title=(
            f"Absolute β — Pearson sensitivity ({cell})"
        ),

        ylabel=(
            "Pearson correlation"
        ),

        output=(
            FIG_DIR /
            f"SUMMARY_absolute_beta_pearson_{cell}.png"
        ),

        ylim=(-0.05, 1.0),
    )


    summary_plot(

        subset,

        metric="spearman",

        title=(
            f"Absolute β — Spearman sensitivity ({cell})"
        ),

        ylabel=(
            "Spearman correlation"
        ),

        output=(
            FIG_DIR /
            f"SUMMARY_absolute_beta_spearman_{cell}.png"
        ),

        ylim=(-0.05, 1.0),
    )


    summary_plot(

        subset,

        metric="mae",

        title=(
            f"Absolute β — MAE vs filter strictness ({cell})"
        ),

        ylabel=(
            "MAE"
        ),

        output=(
            FIG_DIR /
            f"SUMMARY_absolute_beta_mae_{cell}.png"
        ),
    )


    summary_plot(

        subset,

        metric="rmse",

        title=(
            f"Absolute β — RMSE vs filter strictness ({cell})"
        ),

        ylabel=(
            "RMSE"
        ),

        output=(
            FIG_DIR /
            f"SUMMARY_absolute_beta_rmse_{cell}.png"
        ),
    )


    summary_plot(

        subset,

        metric="n",

        title=(
            f"Absolute β — retained sites ({cell})"
        ),

        ylabel=(
            "Retained sites"
        ),

        output=(
            FIG_DIR /
            f"SUMMARY_absolute_beta_n_{cell}.png"
        ),

        log_y=True,
    )


# ============================================================
# 20. DELTA BETA SUMMARY FIGURES
# ============================================================

summary_plot(

    delta_summary,

    metric="pearson",

    title=(
        "Δβ — Pearson correlation vs filter strictness"
    ),

    ylabel=(
        "Pearson correlation"
    ),

    output=(
        FIG_DIR /
        "SUMMARY_delta_beta_pearson.png"
    ),

    ylim=(-0.1, 1.0),
)


summary_plot(

    delta_summary,

    metric="spearman",

    title=(
        "Δβ — Spearman correlation vs filter strictness"
    ),

    ylabel=(
        "Spearman correlation"
    ),

    output=(
        FIG_DIR /
        "SUMMARY_delta_beta_spearman.png"
    ),

    ylim=(-0.1, 1.0),
)


summary_plot(

    delta_summary,

    metric="mae",

    title=(
        "Δβ — MAE vs filter strictness"
    ),

    ylabel=(
        "MAE"
    ),

    output=(
        FIG_DIR /
        "SUMMARY_delta_beta_mae.png"
    ),
)


summary_plot(

    delta_summary,

    metric="rmse",

    title=(
        "Δβ — RMSE vs filter strictness"
    ),

    ylabel=(
        "RMSE"
    ),

    output=(
        FIG_DIR /
        "SUMMARY_delta_beta_rmse.png"
    ),
)


summary_plot(

    delta_summary,

    metric="n",

    title=(
        "Δβ — retained sites vs filter strictness"
    ),

    ylabel=(
        "Retained sites"
    ),

    output=(
        FIG_DIR /
        "SUMMARY_delta_beta_n.png"
    ),

    log_y=True,
)


if (
    "sign_concordance_nonzero"
    in delta_summary.columns
):

    summary_plot(

        delta_summary,

        metric="sign_concordance_nonzero",

        title=(
            "Δβ — directional concordance vs filter strictness"
        ),

        ylabel=(
            "Sign concordance"
        ),

        output=(
            FIG_DIR /
            "SUMMARY_delta_beta_sign_concordance.png"
        ),

        ylim=(0, 1),
    )


# ============================================================
# 21. DONE
# ============================================================

print()
print(
    "=" * 80
)

print(
    "DONE"
)

print(
    "=" * 80
)

print(
    "Figures written to:"
)

print(
    FIG_DIR
)
