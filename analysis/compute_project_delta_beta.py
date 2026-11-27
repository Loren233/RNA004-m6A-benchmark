#!/usr/bin/env python3

"""
Calculate within-project cell-line methylation changes:

    delta_beta = beta_HEK293T - beta_HeLa

The script processes PRJEB106788 and PRJEB80229 separately.

For each project and tool:
  1. Read the collapsed HeLa table.
  2. Read the collapsed HEK293T table in chunks.
  3. Match sites by chr + pos_1based + strand.
  4. Calculate delta_beta.
  5. Preserve HEK/HeLa depth and modified counts for later
     depth filtering and cross-project consensus analysis.

Only sites observed in both cell lines are retained.
"""

import os as _os
from pathlib import Path as _Path
RNA004_STAT_ROOT = _os.environ.get(
    "RNA004_STAT_ROOT",
    str(_Path(__file__).resolve().parents[0])
)


import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


# ------------------------------------------------------------------
# Basic configuration
# ------------------------------------------------------------------

KEY_COLUMNS = ["chr", "pos_1based", "strand"]

ROOT = Path(f"{RNA004_STAT_ROOT}")

PROJECTS = {
    "106788": {
        "input_dir": ROOT / "106788" / "collapsed",
        "output_dir": ROOT / "106788" / "delta_beta",
        "files": {
            "deeprm": {
                "hek": "deeprm_Hek293T.sitelevel.genome.collapsed.tsv",
                "hela": "deeprm_HeLa.sitelevel.genome.collapsed.tsv",
            },
            "dorado": {
                "hek": "dorado_Hek293T.sitelevel.genome.collapsed.tsv",
                "hela": "dorado_HeLa.sitelevel.genome.collapsed.tsv",
            },
            "m6anet": {
                "hek": "m6anet_Hek293t.sitelevel_proba.genome.collapsed.tsv",
                "hela": "m6anet_HeLa.sitelevel_proba.genome.collapsed.tsv",
            },
            "singlemod": {
                "hek": "singlemod_Hek293t.sitelevel.genome.collapsed.tsv",
                "hela": "singlemod_HeLa.sitelevel.genome.collapsed.tsv",
            },
        },
    },

    "80229": {
        "input_dir": ROOT / "80229",
        "output_dir": ROOT / "80229" / "delta_beta",
        "files": {
            "deeprm": {
                "hek": "deeprm_Hek293T.sitelevel.genomic.standardized.collapsed.tsv",
                "hela": "deeprm_HeLa.sitelevel.genomic.standardized.collapsed.tsv",
            },
            "dorado": {
                "hek": "dorado_Hek293T.sitelevel.genomic.standardized.collapsed.tsv",
                "hela": "dorado_HeLa.sitelevel.genomic.standardized.recollapsed.tsv",
            },
            "m6anet": {
                "hek": "m6anet_Hek293T.sitelevel.genomic.standardized.collapsed.tsv",
                "hela": "m6anet_HeLa.sitelevel.genomic.standardized.collapsed.tsv",
            },
            "singlemod": {
                "hek": "singlemod_Hek293T.sitelevel.genomic.standardized.collapsed.tsv",
                "hela": "singlemod_HeLa.sitelevel.genomic.standardized.collapsed.tsv",
            },
        },
    },
}


# Candidate columns are listed in priority order.
COLUMN_CANDIDATES = {
    "deeprm": {
        "beta": [
            "deeprm_stoich",
            "stoich",
            "m6A_fraction",
            "mod_ratio",
        ],
        "depth": [
            "deeprm_depth",
            "depth",
            "total",
            "Total",
        ],
        "modified_count": [
            "deeprm_mod_count",
            "deeprm_estimated_m6A_count",
            "modified_count",
            "m6A_count",
        ],
    },

    "dorado": {
        "beta": [
            "dorado_stoich",
            "stoich",
            "m6A_fraction",
            "mod_ratio",
        ],
        "depth": [
            "dorado_depth",
            "depth",
            "total",
            "Total",
        ],
        "modified_count": [
            "dorado_m6A_count",
            "modified_count",
            "m6A_count",
        ],
    },

    "m6anet": {
        "beta": [
            "m6anet_stoich",
            "stoich",
            "mod_ratio",
            "m6A_fraction",
        ],
        "depth": [
            "n_reads",
            "depth",
            "m6anet_depth",
            "total",
            "Total",
        ],
        "modified_count": [
            "m6anet_modified_count",
            "modified_count",
            "m6A_count",
        ],
    },

    "singlemod": {
        "beta": [
            "singlemod_stoich",
            "stoich",
            "m6A_fraction",
            "mod_ratio",
        ],
        "depth": [
            "singlemod_depth",
            "depth",
            "total",
            "Total",
        ],
        "modified_count": [
            "singlemod_m6A_count",
            "modified_count",
            "m6A_count",
        ],
    },
}


OUTPUT_COLUMNS = [
    "project",
    "tool",
    "chr",
    "pos_1based",
    "start_0based",
    "end_0based_exclusive",
    "strand",
    "site_id",
    "beta_hek",
    "beta_hela",
    "delta_beta",
    "depth_hek",
    "depth_hela",
    "min_depth",
    "modified_count_hek",
    "modified_count_hela",
    "direction_raw",
]


# ------------------------------------------------------------------
# Command-line arguments
# ------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Calculate beta_HEK293T - beta_HeLa within "
            "PRJEB106788 and PRJEB80229."
        )
    )

    parser.add_argument(
        "--project",
        choices=["106788", "80229", "all"],
        default="all",
        help="Project to process. Default: all",
    )

    parser.add_argument(
        "--chunk-size",
        type=int,
        default=500_000,
        help="Number of HEK rows read per chunk. Default: 500000",
    )

    return parser.parse_args()


# ------------------------------------------------------------------
# Column detection
# ------------------------------------------------------------------

def select_column(columns, candidates, field, path, required=True):
    for candidate in candidates:
        if candidate in columns:
            return candidate

    if required:
        raise RuntimeError(
            f"Cannot identify {field} column in {path}\n"
            f"Available columns: {list(columns)}\n"
            f"Tried: {candidates}"
        )

    return None


def inspect_table(path, tool):
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")

    header = pd.read_csv(path, sep="\t", nrows=0)
    columns = list(header.columns)

    missing_keys = [column for column in KEY_COLUMNS if column not in columns]

    if missing_keys:
        raise RuntimeError(
            f"Missing genomic key columns in {path}: {missing_keys}"
        )

    candidates = COLUMN_CANDIDATES[tool]

    beta_column = select_column(
        columns,
        candidates["beta"],
        "beta/stoichiometry",
        path,
        required=True,
    )

    depth_column = select_column(
        columns,
        candidates["depth"],
        "depth",
        path,
        required=True,
    )

    modified_count_column = select_column(
        columns,
        candidates["modified_count"],
        "modified count",
        path,
        required=False,
    )

    return {
        "beta": beta_column,
        "depth": depth_column,
        "modified_count": modified_count_column,
    }


# ------------------------------------------------------------------
# Data cleaning
# ------------------------------------------------------------------

def normalize_strand(series):
    normalized = (
        series.astype("string")
        .str.strip()
        .replace(
            {
                "1": "+",
                "1.0": "+",
                "+1": "+",
                "+1.0": "+",
                "plus": "+",
                "forward": "+",

                "-1": "-",
                "-1.0": "-",
                "minus": "-",
                "reverse": "-",
            }
        )
    )

    invalid = ~normalized.isin(["+", "-"])

    if invalid.any():
        examples = normalized[invalid].drop_duplicates().head(10).tolist()
        raise RuntimeError(
            f"Unexpected strand values: {examples}"
        )

    return normalized


def prepare_frame(df, schema, suffix):
    """
    Convert tool-specific column names to common names.

    Example:
        deeprm_stoich -> beta_hek
        deeprm_depth  -> depth_hek
    """

    rename_map = {
        schema["beta"]: f"beta_{suffix}",
        schema["depth"]: f"depth_{suffix}",
    }

    if schema["modified_count"] is not None:
        rename_map[schema["modified_count"]] = (
            f"modified_count_{suffix}"
        )

    df = df.rename(columns=rename_map)

    df["chr"] = df["chr"].astype("string").str.strip()

    df["strand"] = normalize_strand(df["strand"])

    df["pos_1based"] = pd.to_numeric(
        df["pos_1based"],
        errors="coerce",
    )

    df[f"beta_{suffix}"] = pd.to_numeric(
        df[f"beta_{suffix}"],
        errors="coerce",
    )

    df[f"depth_{suffix}"] = pd.to_numeric(
        df[f"depth_{suffix}"],
        errors="coerce",
    )

    if schema["modified_count"] is not None:
        df[f"modified_count_{suffix}"] = pd.to_numeric(
            df[f"modified_count_{suffix}"],
            errors="coerce",
        )
    else:
        df[f"modified_count_{suffix}"] = np.nan

    required = [
        "chr",
        "pos_1based",
        "strand",
        f"beta_{suffix}",
        f"depth_{suffix}",
    ]

    df = df.dropna(subset=required).copy()

    df["pos_1based"] = df["pos_1based"].astype("int64")

    # Only retain sites with positive depth.
    df = df[df[f"depth_{suffix}"] > 0].copy()

    invalid_beta = (
        (df[f"beta_{suffix}"] < -1e-9)
        | (df[f"beta_{suffix}"] > 1 + 1e-9)
    )

    if invalid_beta.any():
        number = int(invalid_beta.sum())
        minimum = df.loc[invalid_beta, f"beta_{suffix}"].min()
        maximum = df.loc[invalid_beta, f"beta_{suffix}"].max()

        raise RuntimeError(
            f"{number} invalid beta values in {suffix}: "
            f"range {minimum} to {maximum}"
        )

    return df


def read_selected_table(path, schema):
    usecols = KEY_COLUMNS + [
        schema["beta"],
        schema["depth"],
    ]

    if schema["modified_count"] is not None:
        usecols.append(schema["modified_count"])

    # Remove repeated names while preserving order.
    usecols = list(dict.fromkeys(usecols))

    return pd.read_csv(
        path,
        sep="\t",
        usecols=usecols,
        dtype={
            "chr": "string",
            "strand": "string",
        },
        low_memory=False,
    )


# ------------------------------------------------------------------
# Main project/tool processing
# ------------------------------------------------------------------

def process_pair(
    project,
    tool,
    hek_path,
    hela_path,
    output_path,
    chunk_size,
):
    print(
        f"\n[{project} | {tool}]",
        file=sys.stderr,
    )
    print(f"HEK:  {hek_path}", file=sys.stderr)
    print(f"HeLa: {hela_path}", file=sys.stderr)

    hek_schema = inspect_table(hek_path, tool)
    hela_schema = inspect_table(hela_path, tool)

    print(
        "HEK columns: "
        f"beta={hek_schema['beta']}, "
        f"depth={hek_schema['depth']}, "
        f"modified_count={hek_schema['modified_count']}",
        file=sys.stderr,
    )

    print(
        "HeLa columns: "
        f"beta={hela_schema['beta']}, "
        f"depth={hela_schema['depth']}, "
        f"modified_count={hela_schema['modified_count']}",
        file=sys.stderr,
    )

    # --------------------------------------------------------------
    # Load the smaller HeLa table into memory.
    # --------------------------------------------------------------

    hela = read_selected_table(hela_path, hela_schema)
    hela_input_rows = len(hela)

    hela = prepare_frame(
        hela,
        hela_schema,
        suffix="hela",
    )

    duplicated_hela = hela.duplicated(
        KEY_COLUMNS,
        keep=False,
    )

    if duplicated_hela.any():
        number = int(duplicated_hela.sum())

        raise RuntimeError(
            f"{hela_path} still contains {number} duplicated "
            "genomic-key rows after collapse."
        )

    # --------------------------------------------------------------
    # Configure chunked HEK reading.
    # --------------------------------------------------------------

    hek_usecols = KEY_COLUMNS + [
        hek_schema["beta"],
        hek_schema["depth"],
    ]

    if hek_schema["modified_count"] is not None:
        hek_usecols.append(hek_schema["modified_count"])

    hek_usecols = list(dict.fromkeys(hek_usecols))

    output_path.parent.mkdir(parents=True, exist_ok=True)

    if output_path.exists():
        output_path.unlink()

    wrote_header = False

    hek_input_rows = 0
    hek_valid_rows = 0
    matched_sites = 0

    delta_sum = 0.0
    absolute_delta_sum = 0.0

    positive_delta = 0
    negative_delta = 0
    zero_delta = 0

    depth_ge_10 = 0
    depth_ge_20 = 0
    depth_ge_40 = 0

    reader = pd.read_csv(
        hek_path,
        sep="\t",
        usecols=hek_usecols,
        dtype={
            "chr": "string",
            "strand": "string",
        },
        chunksize=chunk_size,
        low_memory=False,
    )

    for chunk_number, chunk in enumerate(reader, start=1):
        hek_input_rows += len(chunk)

        chunk = prepare_frame(
            chunk,
            hek_schema,
            suffix="hek",
        )

        hek_valid_rows += len(chunk)

        duplicated_chunk = chunk.duplicated(
            KEY_COLUMNS,
            keep=False,
        )

        if duplicated_chunk.any():
            number = int(duplicated_chunk.sum())

            raise RuntimeError(
                f"{hek_path} contains at least {number} duplicated "
                f"keys in HEK chunk {chunk_number}."
            )

        paired = chunk.merge(
            hela,
            on=KEY_COLUMNS,
            how="inner",
            validate="one_to_one",
        )

        if paired.empty:
            if chunk_number % 10 == 0:
                print(
                    f"  chunks={chunk_number}, "
                    f"HEK rows read={hek_input_rows:,}, "
                    f"matched={matched_sites:,}",
                    file=sys.stderr,
                )
            continue

        paired["delta_beta"] = (
            paired["beta_hek"] - paired["beta_hela"]
        )

        paired["min_depth"] = np.minimum(
            paired["depth_hek"],
            paired["depth_hela"],
        )

        paired["start_0based"] = (
            paired["pos_1based"] - 1
        )

        paired["end_0based_exclusive"] = (
            paired["pos_1based"]
        )

        paired["site_id"] = (
            paired["chr"].astype(str)
            + "|"
            + paired["pos_1based"].astype(str)
            + "|"
            + paired["strand"].astype(str)
        )

        paired["project"] = project
        paired["tool"] = tool

        zero_mask = np.isclose(
            paired["delta_beta"],
            0.0,
            atol=1e-12,
        )

        paired["direction_raw"] = np.select(
            [
                paired["delta_beta"] > 0,
                paired["delta_beta"] < 0,
            ],
            [
                "HEK_higher",
                "HeLa_higher",
            ],
            default="equal",
        )

        number_matched = len(paired)
        matched_sites += number_matched

        delta_sum += paired["delta_beta"].sum()
        absolute_delta_sum += paired["delta_beta"].abs().sum()

        positive_delta += int(
            (paired["delta_beta"] > 0).sum()
        )

        negative_delta += int(
            (paired["delta_beta"] < 0).sum()
        )

        zero_delta += int(zero_mask.sum())

        depth_ge_10 += int(
            (paired["min_depth"] >= 10).sum()
        )

        depth_ge_20 += int(
            (paired["min_depth"] >= 20).sum()
        )

        depth_ge_40 += int(
            (paired["min_depth"] >= 40).sum()
        )

        paired = paired[OUTPUT_COLUMNS]

        paired.to_csv(
            output_path,
            sep="\t",
            index=False,
            mode="a",
            header=not wrote_header,
        )

        wrote_header = True

        if chunk_number % 10 == 0:
            print(
                f"  chunks={chunk_number}, "
                f"HEK rows read={hek_input_rows:,}, "
                f"matched={matched_sites:,}",
                file=sys.stderr,
            )

    # Write an empty table with a header when no sites matched.
    if not wrote_header:
        pd.DataFrame(
            columns=OUTPUT_COLUMNS
        ).to_csv(
            output_path,
            sep="\t",
            index=False,
        )

    mean_delta = (
        delta_sum / matched_sites
        if matched_sites > 0
        else np.nan
    )

    mean_absolute_delta = (
        absolute_delta_sum / matched_sites
        if matched_sites > 0
        else np.nan
    )

    summary = {
        "project": project,
        "tool": tool,
        "hek_input_rows": hek_input_rows,
        "hek_valid_rows": hek_valid_rows,
        "hela_input_rows": hela_input_rows,
        "hela_valid_rows": len(hela),
        "matched_sites": matched_sites,
        "delta_beta_gt_0_HEK_higher": positive_delta,
        "delta_beta_lt_0_HeLa_higher": negative_delta,
        "delta_beta_equal_0": zero_delta,
        "matched_min_depth_ge_10": depth_ge_10,
        "matched_min_depth_ge_20": depth_ge_20,
        "matched_min_depth_ge_40": depth_ge_40,
        "mean_delta_beta": mean_delta,
        "mean_absolute_delta_beta": mean_absolute_delta,
        "output_file": str(output_path),
    }

    print(
        f"Completed {project} {tool}: "
        f"{matched_sites:,} matched sites",
        file=sys.stderr,
    )

    return summary


# ------------------------------------------------------------------
# Program entry
# ------------------------------------------------------------------

def main():
    args = parse_args()

    if args.chunk_size <= 0:
        raise ValueError("--chunk-size must be greater than zero")

    if args.project == "all":
        selected_projects = ["106788", "80229"]
    else:
        selected_projects = [args.project]

    all_summaries = []

    for project in selected_projects:
        project_config = PROJECTS[project]

        input_dir = project_config["input_dir"]
        output_dir = project_config["output_dir"]

        output_dir.mkdir(
            parents=True,
            exist_ok=True,
        )

        project_summaries = []

        for tool, filenames in project_config["files"].items():
            hek_path = input_dir / filenames["hek"]
            hela_path = input_dir / filenames["hela"]

            output_path = (
                output_dir
                / f"{tool}.HEK_minus_HeLa.delta_beta.tsv"
            )

            summary = process_pair(
                project=project,
                tool=tool,
                hek_path=hek_path,
                hela_path=hela_path,
                output_path=output_path,
                chunk_size=args.chunk_size,
            )

            project_summaries.append(summary)
            all_summaries.append(summary)

        project_summary_path = (
            output_dir / "delta_beta_summary.tsv"
        )

        pd.DataFrame(project_summaries).to_csv(
            project_summary_path,
            sep="\t",
            index=False,
        )

        print(
            f"Project summary: {project_summary_path}",
            file=sys.stderr,
        )

    global_summary_path = (
        ROOT / "delta_beta_all_projects_summary.tsv"
    )

    pd.DataFrame(all_summaries).to_csv(
        global_summary_path,
        sep="\t",
        index=False,
    )

    print(
        f"\nGlobal summary: {global_summary_path}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
