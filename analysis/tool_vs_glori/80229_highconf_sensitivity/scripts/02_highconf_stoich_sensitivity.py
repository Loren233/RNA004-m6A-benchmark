
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

RESULT_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)

GLORI_FILE = (
    RESULT_DIR /
    "GLORI_HEK_HeLa.highconf.tsv.gz"
)


# ============================================================
# 1. INPUT FILES
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
# 2. ANALYSIS PARAMETERS
# ============================================================

KEY = [
    "chr",
    "pos_1based",
    "strand",
]


# Main coverage filter
MAIN_COVERAGE = 20

# Additional coverage sensitivity
STRICT_COVERAGE = 40


# Nested adjusted-p thresholds
P_THRESHOLDS = [
    0.05,
    0.01,
    0.001,
]


# After reaching p < 0.001,
# progressively increase GT stoichiometry threshold
STOICH_THRESHOLDS = [
    0.10,
    0.20,
    0.50,
]


CHUNK_SIZE = 1_000_000

MAX_PLOT_POINTS = 600_000

RNG = np.random.default_rng(
    20260811
)


# ============================================================
# 3. MAIN FILTER ORDER
# ============================================================

MAIN_STAGE_ORDER = [
    "raw_common",
    "cov20",
    "p005",
    "p001",
    "p0001",
    "beta010",
    "beta020",
    "beta050",
]


MAIN_STAGE_LABELS = [
    "Raw",
    "Cov >20",
    "p <0.05",
    "p <0.01",
    "p <0.001",
    "β >0.10",
    "β >0.20",
    "β >0.50",
]


FILTER_DESCRIPTIONS = {

    "raw_common":
        "All common sites",

    "cov20":
        "Tool and GLORI coverage >20",

    "p005":
        "Coverage >20 + both GLORI replicates AdjP <0.05",

    "p001":
        "Coverage >20 + both GLORI replicates AdjP <0.01",

    "p0001":
        "Coverage >20 + both GLORI replicates AdjP <0.001",

    "beta010":
        "p<0.001 + GLORI stoichiometry >0.10",

    "beta020":
        "p<0.001 + GLORI stoichiometry >0.20",

    "beta050":
        "p<0.001 + GLORI stoichiometry >0.50",
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

density_cmap = LinearSegmentedColormap.from_list(
    "paper_density",
    paper_density,
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
        f"Cannot find file:\n{path}\n"
        f"or:\n{gz}"
    )


def find_column(
    columns,
    candidates,
    label,
):

    for candidate in candidates:

        if candidate in columns:
            return candidate

    raise KeyError(
        f"Cannot identify {label} column.\n"
        f"Candidates: {candidates}\n"
        f"Available columns:\n"
        f"{list(columns)}"
    )


def resolve_tool_columns(path):

    header = pd.read_csv(
        path,
        sep="\t",
        nrows=0,
    )

    cols = list(
        header.columns
    )


    for required in KEY:

        if required not in cols:

            raise KeyError(
                f"Missing genomic key {required} in:\n"
                f"{path}\n"
                f"Columns:\n{cols}"
            )


    depth_col = find_column(
        cols,
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
        cols,
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


    return {
        "depth":
            depth_col,

        "stoich":
            stoich_col,
    }


# ============================================================
# 6. METRICS
# ============================================================

def concordance_ccc(x, y):

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


    mean_x = np.mean(x)
    mean_y = np.mean(y)

    var_x = np.var(
        x,
        ddof=1,
    )

    var_y = np.var(
        y,
        ddof=1,
    )

    cov_xy = np.cov(
        x,
        y,
        ddof=1,
    )[0, 1]


    denominator = (
        var_x +
        var_y +
        (mean_x - mean_y) ** 2
    )


    if denominator == 0:
        return np.nan


    return (
        2 * cov_xy /
        denominator
    )


def calculate_metrics(
    gt_values,
    tool_values,
):

    x = np.asarray(
        gt_values,
        dtype=float,
    )

    y = np.asarray(
        tool_values,
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

            "mean_bias_tool_minus_gt":
                np.nan,

            "gt_sd":
                np.nan,

            "tool_sd":
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

        "mean_bias_tool_minus_gt":
            np.mean(error),

        # Important for interpreting
        # possible range-restriction effects
        "gt_sd":
            np.std(
                x,
                ddof=1,
            ),

        "tool_sd":
            np.std(
                y,
                ddof=1,
            ),
    }


def calculate_sign_concordance(
    gt_delta,
    tool_delta,
):

    x = np.asarray(
        gt_delta,
        dtype=float,
    )

    y = np.asarray(
        tool_delta,
        dtype=float,
    )


    valid = (
        np.isfinite(x) &
        np.isfinite(y) &
        (x != 0) &
        (y != 0)
    )


    if valid.sum() == 0:
        return np.nan


    return np.mean(
        np.sign(
            x[valid]
        )
        ==
        np.sign(
            y[valid]
        )
    )


def sample_for_plot(df):

    if len(df) <= MAX_PLOT_POINTS:
        return df


    idx = RNG.choice(
        len(df),
        size=MAX_PLOT_POINTS,
        replace=False,
    )


    return df.iloc[idx]


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

print("=" * 80)
print("Loading GLORI")
print("=" * 80)


if not GLORI_FILE.exists():

    raise FileNotFoundError(
        f"GLORI file does not exist:\n"
        f"{GLORI_FILE}"
    )


gt = pd.read_csv(
    GLORI_FILE,
    sep="\t",
    low_memory=False,
)


# We now explicitly require replicate-level p-values.
required_gt = {

    "chr",
    "pos_1based",
    "strand",

    "glori_total_hek",
    "glori_beta_hek",

    "glori_total_hela",
    "glori_beta_hela",

    "glori_delta_beta",

    "glori_adj_p_hek_rep1",
    "glori_adj_p_hek_rep2",

    "glori_adj_p_hela_rep1",
    "glori_adj_p_hela_rep2",
}


missing_gt = (
    required_gt -
    set(gt.columns)
)


if missing_gt:

    raise KeyError(
        "Missing GLORI columns:\n"
        f"{sorted(missing_gt)}"
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


gt_dup = (
    gt.duplicated(KEY)
    .sum()
)


if gt_dup:

    raise RuntimeError(
        f"GLORI contains "
        f"{gt_dup:,} duplicate genomic keys"
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


    mapping = resolve_tool_columns(
        path
    )


    depth_col = (
        mapping["depth"]
    )

    stoich_col = (
        mapping["stoich"]
    )


    print()
    print("=" * 80)
    print(
        f"{tool} -- {cell}"
    )
    print(path)

    print(
        f"depth column  = {depth_col}"
    )

    print(
        f"stoich column = {stoich_col}"
    )

    print("=" * 80)


    if cell == "HEK293T":

        gt_cols = KEY + [

            "glori_total_hek",
            "glori_beta_hek",

            "glori_adj_p_hek_rep1",
            "glori_adj_p_hek_rep2",
        ]


    elif cell == "HeLa":

        gt_cols = KEY + [

            "glori_total_hela",
            "glori_beta_hela",

            "glori_adj_p_hela_rep1",
            "glori_adj_p_hela_rep2",
        ]


    else:

        raise ValueError(
            cell
        )


    gt_sub = (
        gt[gt_cols]
    )


    usecols = (
        KEY +
        [
            depth_col,
            stoich_col,
        ]
    )


    matched_parts = []


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


        chunk["pos_1based"] = pd.to_numeric(
            chunk["pos_1based"],
            errors="coerce",
        )


        chunk["tool_depth"] = pd.to_numeric(
            chunk["tool_depth"],
            errors="coerce",
        )


        chunk["tool_beta"] = pd.to_numeric(
            chunk["tool_beta"],
            errors="coerce",
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

            matched_parts.append(
                merged
            )


        print(
            f"chunk {chunk_id:03d}: "
            f"input={input_n:,}, "
            f"valid={len(chunk):,}, "
            f"matched={len(merged):,}",
            flush=True,
        )


    if not matched_parts:

        raise RuntimeError(
            f"No GLORI matches for "
            f"{tool} {cell}"
        )


    out = pd.concat(
        matched_parts,
        ignore_index=True,
    )


    duplicate_count = (
        out.duplicated(KEY)
        .sum()
    )


    if duplicate_count:

        raise RuntimeError(
            f"{tool} {cell}: "
            f"{duplicate_count:,} duplicated genomic keys"
        )


    print(
        f"TOTAL MATCHED {tool} {cell}: "
        f"{len(out):,}"
    )


    return out


# ============================================================
# 9. LOAD ALL TOOLS
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
# 10. GLORI SIGNIFICANCE HELPERS
# ============================================================

def absolute_glori_sig(
    df,
    cell,
    cutoff,
):
    """
    STRICT RULE:
    BOTH GLORI replicates must have
    precomputed adjusted p-value < cutoff.
    """

    if cell == "HEK293T":

        p1 = (
            df["glori_adj_p_hek_rep1"]
        )

        p2 = (
            df["glori_adj_p_hek_rep2"]
        )


    elif cell == "HeLa":

        p1 = (
            df["glori_adj_p_hela_rep1"]
        )

        p2 = (
            df["glori_adj_p_hela_rep2"]
        )


    else:

        raise ValueError(
            cell
        )


    return (
        np.isfinite(p1) &
        np.isfinite(p2) &
        (p1 < cutoff) &
        (p2 < cutoff)
    )


def delta_glori_sig(
    df,
    cutoff,
):
    """
    For differential analysis:

    At least ONE cell line must have
    significant m6A evidence.

    But within that cell line,
    BOTH GLORI replicates must pass.

    This preserves strong gain/loss sites.
    """

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
# 11. ABSOLUTE BETA NESTED MASKS
# ============================================================

def get_absolute_columns(
    cell,
):

    if cell == "HEK293T":

        return (
            "glori_beta_hek",
            "glori_total_hek",
        )


    elif cell == "HeLa":

        return (
            "glori_beta_hela",
            "glori_total_hela",
        )


    raise ValueError(
        cell
    )


def absolute_main_masks(
    df,
    cell,
):

    (
        gt_beta_col,
        gt_depth_col,
    ) = get_absolute_columns(
        cell
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
            df["tool_depth"]
            >
            MAIN_COVERAGE
        )

        &
        (
            df[gt_depth_col]
            >
            MAIN_COVERAGE
        )
    )


    p005 = (
        cov20
        &
        absolute_glori_sig(
            df,
            cell,
            0.05,
        )
    )


    p001 = (
        cov20
        &
        absolute_glori_sig(
            df,
            cell,
            0.01,
        )
    )


    p0001 = (
        cov20
        &
        absolute_glori_sig(
            df,
            cell,
            0.001,
        )
    )


    beta010 = (
        p0001
        &
        (
            df[gt_beta_col]
            >
            0.10
        )
    )


    beta020 = (
        p0001
        &
        (
            df[gt_beta_col]
            >
            0.20
        )
    )


    beta050 = (
        p0001
        &
        (
            df[gt_beta_col]
            >
            0.50
        )
    )


    return (
        gt_beta_col,
        {

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
    )


# ============================================================
# 12. ABSOLUTE BETA ANALYSIS
# ============================================================

absolute_rows = []


for tool in FILES:

    for cell in [
        "HEK293T",
        "HeLa",
    ]:

        df = (
            matched[tool][cell]
        )


        (
            gt_beta_col,
            masks,
        ) = absolute_main_masks(
            df,
            cell,
        )


        for stage in MAIN_STAGE_ORDER:

            mask = (
                masks[stage]
            )


            stat = calculate_metrics(
                df.loc[
                    mask,
                    gt_beta_col,
                ],

                df.loc[
                    mask,
                    "tool_beta",
                ],
            )


            absolute_rows.append(
                {

                    "tool":
                        tool,

                    "cell":
                        cell,

                    "stage":
                        stage,

                    "filter_description":
                        FILTER_DESCRIPTIONS[
                            stage
                        ],

                    **stat,
                }
            )


absolute_summary = pd.DataFrame(
    absolute_rows
)


absolute_output = (
    RESULT_DIR /
    "absolute_beta_strictness_metrics.tsv"
)


absolute_summary.to_csv(
    absolute_output,
    sep="\t",
    index=False,
)


# ============================================================
# 13. BUILD HEK-HeLa DELTA TABLES
# ============================================================

delta_tables = {}


for tool in FILES:

    hek = (
        matched[tool]["HEK293T"]
        [
            KEY + [

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
            KEY + [

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


    # For differential filtering:
    # retain site if either cell line has high GT stoichiometry.
    d["glori_max_stoich"] = np.maximum(

        d["glori_beta_hek"],
        d["glori_beta_hela"],
    )


    delta_tables[tool] = d


    print(
        f"{tool} paired HEK-HeLa sites: "
        f"{len(d):,}"
    )


# ============================================================
# 14. DELTA NESTED MASKS
# ============================================================

def delta_main_masks(d):

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
            d["tool_depth_hek"]
            >
            MAIN_COVERAGE
        )

        &
        (
            d["tool_depth_hela"]
            >
            MAIN_COVERAGE
        )

        &
        (
            d["glori_total_hek"]
            >
            MAIN_COVERAGE
        )

        &
        (
            d["glori_total_hela"]
            >
            MAIN_COVERAGE
        )
    )


    p005 = (
        cov20
        &
        delta_glori_sig(
            d,
            0.05,
        )
    )


    p001 = (
        cov20
        &
        delta_glori_sig(
            d,
            0.01,
        )
    )


    p0001 = (
        cov20
        &
        delta_glori_sig(
            d,
            0.001,
        )
    )


    beta010 = (
        p0001
        &
        (
            d["glori_max_stoich"]
            >
            0.10
        )
    )


    beta020 = (
        p0001
        &
        (
            d["glori_max_stoich"]
            >
            0.20
        )
    )


    beta050 = (
        p0001
        &
        (
            d["glori_max_stoich"]
            >
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
# 15. DELTA ANALYSIS
# ============================================================

delta_rows = []


for tool, d in delta_tables.items():

    masks = delta_main_masks(
        d
    )


    for stage in MAIN_STAGE_ORDER:

        mask = (
            masks[stage]
        )


        stat = calculate_metrics(

            d.loc[
                mask,
                "glori_delta_beta",
            ],

            d.loc[
                mask,
                "tool_delta_beta",
            ],
        )


        sign_concordance = (
            calculate_sign_concordance(

                d.loc[
                    mask,
                    "glori_delta_beta",
                ],

                d.loc[
                    mask,
                    "tool_delta_beta",
                ],
            )
        )


        delta_rows.append(
            {

                "tool":
                    tool,

                "stage":
                    stage,

                "filter_description":
                    FILTER_DESCRIPTIONS[
                        stage
                    ],

                **stat,

                "sign_concordance_nonzero":
                    sign_concordance,
            }
        )


delta_summary = pd.DataFrame(
    delta_rows
)


delta_output = (
    RESULT_DIR /
    "delta_beta_strictness_metrics.tsv"
)


delta_summary.to_csv(
    delta_output,
    sep="\t",
    index=False,
)


# ============================================================
# 16. EXTRA COVERAGE >40 SENSITIVITY
# ============================================================

coverage40_rows = []


# ------------------------------------------------------------
# Absolute beta
# ------------------------------------------------------------

for tool in FILES:

    for cell in [
        "HEK293T",
        "HeLa",
    ]:

        df = (
            matched[tool][cell]
        )


        (
            gt_beta_col,
            gt_depth_col,
        ) = get_absolute_columns(
            cell
        )


        base = (

            np.isfinite(
                df["tool_beta"]
            )
            &
            np.isfinite(
                df[gt_beta_col]
            )
        )


        for coverage in [
            20,
            40,
        ]:

            cov = (

                base

                &
                (
                    df["tool_depth"]
                    >
                    coverage
                )

                &
                (
                    df[gt_depth_col]
                    >
                    coverage
                )
            )


            strict_sig = (

                cov

                &
                absolute_glori_sig(
                    df,
                    cell,
                    0.001,
                )
            )


            strict_beta = (

                strict_sig

                &
                (
                    df[gt_beta_col]
                    >
                    0.50
                )
            )


            for stage_name, mask in [

                (
                    "p0001",
                    strict_sig,
                ),

                (
                    "p0001_beta050",
                    strict_beta,
                ),
            ]:

                stat = calculate_metrics(

                    df.loc[
                        mask,
                        gt_beta_col,
                    ],

                    df.loc[
                        mask,
                        "tool_beta",
                    ],
                )


                coverage40_rows.append(
                    {

                        "analysis":
                            "absolute_beta",

                        "tool":
                            tool,

                        "cell":
                            cell,

                        "coverage_threshold":
                            coverage,

                        "stage":
                            stage_name,

                        **stat,
                    }
                )


# ------------------------------------------------------------
# Delta beta
# ------------------------------------------------------------

for tool, d in delta_tables.items():

    base = (

        np.isfinite(
            d["tool_delta_beta"]
        )
        &
        np.isfinite(
            d["glori_delta_beta"]
        )
    )


    for coverage in [
        20,
        40,
    ]:

        cov = (

            base

            &
            (
                d["tool_depth_hek"]
                >
                coverage
            )

            &
            (
                d["tool_depth_hela"]
                >
                coverage
            )

            &
            (
                d["glori_total_hek"]
                >
                coverage
            )

            &
            (
                d["glori_total_hela"]
                >
                coverage
            )
        )


        strict_sig = (

            cov

            &
            delta_glori_sig(
                d,
                0.001,
            )
        )


        strict_beta = (

            strict_sig

            &
            (
                d["glori_max_stoich"]
                >
                0.50
            )
        )


        for stage_name, mask in [

            (
                "p0001",
                strict_sig,
            ),

            (
                "p0001_beta050",
                strict_beta,
            ),
        ]:

            stat = calculate_metrics(

                d.loc[
                    mask,
                    "glori_delta_beta",
                ],

                d.loc[
                    mask,
                    "tool_delta_beta",
                ],
            )


            sign_concordance = (
                calculate_sign_concordance(

                    d.loc[
                        mask,
                        "glori_delta_beta",
                    ],

                    d.loc[
                        mask,
                        "tool_delta_beta",
                    ],
                )
            )


            coverage40_rows.append(
                {

                    "analysis":
                        "delta_beta",

                    "tool":
                        tool,

                    "cell":
                        "HEK_minus_HeLa",

                    "coverage_threshold":
                        coverage,

                    "stage":
                        stage_name,

                    **stat,

                    "sign_concordance_nonzero":
                        sign_concordance,
                }
            )


coverage40_summary = pd.DataFrame(
    coverage40_rows
)


coverage40_output = (
    RESULT_DIR /
    "coverage20_vs_40_sensitivity.tsv"
)


coverage40_summary.to_csv(
    coverage40_output,
    sep="\t",
    index=False,
)


# ============================================================
# 17. PRINT MAIN RESULTS
# ============================================================

print()
print("=" * 110)
print("ABSOLUTE BETA — NESTED FILTERING")
print("=" * 110)


print(
    absolute_summary[
        [
            "tool",
            "cell",
            "stage",
            "n",
            "pearson",
            "spearman",
            "ccc",
            "mae",
            "rmse",
            "gt_sd",
        ]
    ]
    .to_string(
        index=False
    )
)


print()
print("=" * 110)
print("DELTA BETA — NESTED FILTERING")
print("=" * 110)


print(
    delta_summary[
        [
            "tool",
            "stage",
            "n",
            "pearson",
            "spearman",
            "ccc",
            "mae",
            "rmse",
            "gt_sd",
            "sign_concordance_nonzero",
        ]
    ]
    .to_string(
        index=False
    )
)


# ============================================================
# 18. GENERIC SENSITIVITY CURVE
# ============================================================

def plot_metric_curve(
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
                df["tool"] == tool
            ]
            .set_index(
                "stage"
            )
            .reindex(
                MAIN_STAGE_ORDER
            )
        )


        ax.plot(
            MAIN_STAGE_LABELS,
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


    ax.set_ylabel(
        ylabel
    )

    ax.set_title(
        title
    )

    ax.tick_params(
        axis="x",
        rotation=30,
    )

    ax.legend()

    fig.tight_layout()


    fig.savefig(
        output,
        dpi=300,
    )

    plt.close(fig)


# ============================================================
# 19. ABSOLUTE BETA SUMMARY CURVES
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


    plot_metric_curve(

        subset,

        metric="pearson",

        title=(
            f"Absolute β sensitivity — {cell}"
        ),

        ylabel=(
            "Pearson correlation"
        ),

        output=(
            FIG_DIR /
            f"absolute_beta_strictness_pearson_{cell}.png"
        ),

        ylim=(-0.05, 1.0),
    )


    plot_metric_curve(

        subset,

        metric="mae",

        title=(
            f"Absolute β error sensitivity — {cell}"
        ),

        ylabel=(
            "MAE"
        ),

        output=(
            FIG_DIR /
            f"absolute_beta_strictness_mae_{cell}.png"
        ),
    )


    plot_metric_curve(

        subset,

        metric="n",

        title=(
            f"Absolute β retained sites — {cell}"
        ),

        ylabel=(
            "Number of retained sites"
        ),

        output=(
            FIG_DIR /
            f"absolute_beta_strictness_n_{cell}.png"
        ),

        log_y=True,
    )


# ============================================================
# 20. DELTA BETA SUMMARY CURVES
# ============================================================

plot_metric_curve(

    delta_summary,

    metric="pearson",

    title=(
        "Δβ sensitivity to increasing filter strictness"
    ),

    ylabel=(
        "Pearson correlation"
    ),

    output=(
        FIG_DIR /
        "delta_beta_strictness_pearson.png"
    ),

    ylim=(-0.1, 1.0),
)


plot_metric_curve(

    delta_summary,

    metric="mae",

    title=(
        "Δβ error vs increasing filter strictness"
    ),

    ylabel=(
        "MAE"
    ),

    output=(
        FIG_DIR /
        "delta_beta_strictness_mae.png"
    ),
)


plot_metric_curve(

    delta_summary,

    metric="n",

    title=(
        "Δβ retained sites vs increasing filter strictness"
    ),

    ylabel=(
        "Number of retained sites"
    ),

    output=(
        FIG_DIR /
        "delta_beta_strictness_n.png"
    ),

    log_y=True,
)


plot_metric_curve(

    delta_summary,

    metric="sign_concordance_nonzero",

    title=(
        "Δβ directional concordance vs filter strictness"
    ),

    ylabel=(
        "Sign concordance"
    ),

    output=(
        FIG_DIR /
        "delta_beta_strictness_sign_concordance.png"
    ),

    ylim=(0, 1),
)


# ============================================================
# 21. SCATTER PLOT HELPERS
# ============================================================

def make_absolute_scatter(
    stage,
):

    fig, axes = plt.subplots(
        2,
        3,
        figsize=(14, 9),
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
            ) = absolute_main_masks(
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


            plot_sub = sample_for_plot(
                sub
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
                linewidth=1,
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
        "Absolute β vs GLORI\n"
        +
        FILTER_DESCRIPTIONS[
            stage
        ]
    )


    fig.savefig(
        FIG_DIR /
        f"absolute_beta_{stage}.png",
        dpi=300,
    )

    plt.close(fig)


def make_delta_scatter(
    stage,
):

    fig, axes = plt.subplots(
        1,
        3,
        figsize=(14, 4.8),
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


        masks = delta_main_masks(
            d
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


        plot_sub = sample_for_plot(
            sub
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
            linewidth=1,
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
        "Δβ vs GLORI\n"
        +
        FILTER_DESCRIPTIONS[
            stage
        ]
    )


    fig.savefig(
        FIG_DIR /
        f"delta_beta_{stage}.png",
        dpi=300,
    )

    plt.close(fig)


# ============================================================
# 22. KEY SCATTER PLOTS
# ============================================================

# Don't create eight giant panels.
# Keep only biologically/statistically meaningful checkpoints.

SCATTER_STAGES = [
    "cov20",
    "p005",
    "p0001",
    "beta050",
]


for stage in SCATTER_STAGES:

    make_absolute_scatter(
        stage
    )

    make_delta_scatter(
        stage
    )


# ============================================================
# 23. DONE
# ============================================================

print()
print("=" * 110)
print("DONE")
print("=" * 110)

print()
print(
    "Absolute beta metrics:"
)
print(
    absolute_output
)

print()
print(
    "Delta beta metrics:"
)
print(
    delta_output
)

print()
print(
    "Coverage 20 vs 40 sensitivity:"
)
print(
    coverage40_output
)

print()
print(
    "Figures:"
)
print(
    FIG_DIR
)
