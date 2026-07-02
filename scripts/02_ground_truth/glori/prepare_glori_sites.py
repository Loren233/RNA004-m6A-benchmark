#!/usr/bin/env python3
"""Prepare final HEK293T GLORI ground-truth sites for benchmarking."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List

import pandas as pd


EXPECTED_COLUMNS = [
    "chr",
    "pos",
    "strand",
    "glori_depth_rep1",
    "glori_stoich_rep1",
    "glori_depth_rep2",
    "glori_stoich_rep2",
    "glori_depth_mean",
    "glori_stoich_mean",
]

OUTPUT_COLUMNS = [
    "chr",
    "pos_1based",
    "strand",
    "depth",
    "stoich",
    "glori_depth_rep1",
    "glori_stoich_rep1",
    "glori_depth_rep2",
    "glori_stoich_rep2",
    "glori_depth_mean",
    "glori_stoich_mean",
    "source",
    "cell_line",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Final GLORI replicate-intersection TSV")
    parser.add_argument("--output", required=True, help="Prepared output TSV")
    parser.add_argument("--force", action="store_true", help="Overwrite an existing output")
    return parser.parse_args()


def require_columns(df: pd.DataFrame) -> None:
    missing = [column for column in EXPECTED_COLUMNS if column not in df.columns]
    if missing:
        raise SystemExit("ERROR: missing required input column(s): %s" % ", ".join(missing))


def numeric_column(df: pd.DataFrame, column: str) -> pd.Series:
    return pd.to_numeric(df[column], errors="coerce")


def validate(df: pd.DataFrame) -> None:
    errors = []  # type: List[str]
    if df["chr"].isna().any() or df["chr"].astype(str).str.strip().eq("").any():
        errors.append("chr contains missing or empty values")

    pos = numeric_column(df, "pos_1based")
    if pos.isna().any() or (pos % 1 != 0).any() or (pos <= 0).any():
        errors.append("pos_1based must contain positive integers")

    if (~df["strand"].isin(["+", "-"])).any():
        errors.append("strand must be either + or -")

    depth = numeric_column(df, "depth")
    if depth.isna().any() or (depth < 0).any():
        errors.append("depth must be numeric and non-negative")

    stoich = numeric_column(df, "stoich")
    if stoich.isna().any() or (stoich < 0).any() or (stoich > 1).any():
        errors.append("stoich must be numeric and between 0 and 1")

    duplicate_rows = int(df.duplicated(["chr", "pos_1based", "strand"], keep=False).sum())
    if duplicate_rows:
        errors.append("duplicate genomic keys found: %d rows" % duplicate_rows)

    if errors:
        raise SystemExit("ERROR: invalid GLORI ground-truth table:\n- " + "\n- ".join(errors))


def main() -> int:
    args = parse_args()
    input_path = Path(args.input)
    output_path = Path(args.output)

    if output_path.exists() and not args.force:
        raise SystemExit("ERROR: output exists; use --force to overwrite: %s" % output_path)

    df = pd.read_csv(input_path, sep="\t")
    input_rows = len(df)
    require_columns(df)

    df = df.rename(
        columns={
            "pos": "pos_1based",
            "glori_depth_mean": "depth",
            "glori_stoich_mean": "stoich",
        }
    )
    df["glori_depth_mean"] = df["depth"]
    df["glori_stoich_mean"] = df["stoich"]
    df["source"] = "GLORI"
    df["cell_line"] = "HEK293T"

    validate(df)

    df["pos_1based"] = pd.to_numeric(df["pos_1based"], errors="raise").astype("int64")
    duplicate_rows = int(df.duplicated(["chr", "pos_1based", "strand"], keep=False).sum())
    depth = pd.to_numeric(df["depth"], errors="raise")
    stoich = pd.to_numeric(df["stoich"], errors="raise")

    output = df[OUTPUT_COLUMNS].sort_values(["chr", "pos_1based", "strand"])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(output_path, sep="\t", index=False)

    print("input_rows\t%d" % input_rows)
    print("output_rows\t%d" % len(output))
    print("duplicate_key_count\t%d" % duplicate_rows)
    print("depth_range\t%s\t%s" % (depth.min(), depth.max()))
    print("stoichiometry_range\t%s\t%s" % (stoich.min(), stoich.max()))
    print("output_path\t%s" % output_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
