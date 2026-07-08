#!/usr/bin/env python3
"""Match standardized prediction sites to a ground-truth table by genomic offset.

This script is ground-truth agnostic. It can be used for GLORI, eTAM-seq,
m6A-SAC-seq, m6ACE-seq, miCLIP-derived site tables, or any standardized
site-level truth table with chr + position + strand + stoichiometry columns.

Convention:
  prediction_pos_1based = truth_pos_1based + offset

Outputs:
  <output-prefix>_offset_<label>.tsv
  <output-prefix>_offset_summary.tsv
"""

from __future__ import annotations

import argparse
import os
import sys
from datetime import datetime
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


def log(message: str) -> None:
    print(message, file=sys.stderr, flush=True)


def unique_join(values: Iterable[object]) -> str:
    out = []
    for value in values:
        if pd.isna(value):
            continue
        text = str(value).strip()
        if text == "":
            continue
        out.append(text)
    return ",".join(sorted(set(out)))


def read_table(path: str, sep: str) -> pd.DataFrame:
    if not Path(path).exists():
        raise FileNotFoundError(f"input file does not exist: {path}")
    return pd.read_csv(path, sep=sep, compression="infer", dtype=str, keep_default_na=False, na_values=[])


def require_columns(df: pd.DataFrame, columns: list[str], label: str) -> None:
    missing = [col for col in columns if col and col not in df.columns]
    if missing:
        raise SystemExit(f"ERROR: {label} missing required columns: {missing}\nExisting columns: {list(df.columns)}")


def parse_offsets(offsets: list[int], labels: list[str] | None) -> list[tuple[int, str]]:
    if labels is None:
        out = []
        for offset in offsets:
            if offset < 0:
                label = f"m{abs(offset)}"
            elif offset > 0:
                label = f"p{offset}"
            else:
                label = "0"
            out.append((offset, label))
        return out
    if len(offsets) != len(labels):
        raise SystemExit("ERROR: --offsets and --offset-labels must have the same length")
    return list(zip(offsets, labels))


def safe_numeric(series: pd.Series, label: str) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    bad = values.isna() & series.astype(str).str.strip().ne("")
    if bad.any():
        log(f"WARNING\tnon_numeric_values\t{label}\t{int(bad.sum())}")
    return values


def corr_stats(df: pd.DataFrame, x_col: str, y_col: str) -> dict[str, object]:
    sub = df[[x_col, y_col]].copy()
    sub[x_col] = pd.to_numeric(sub[x_col], errors="coerce")
    sub[y_col] = pd.to_numeric(sub[y_col], errors="coerce")
    sub = sub.dropna()
    n = len(sub)

    if n < 2:
        return {
            "n": n,
            "pearson_r": np.nan,
            "pearson_r2": np.nan,
            "spearman_rho": np.nan,
            "spearman_rho2": np.nan,
        }

    pearson_r = sub[x_col].corr(sub[y_col], method="pearson")
    spearman_rho = sub[x_col].rank().corr(sub[y_col].rank(), method="pearson")

    return {
        "n": n,
        "pearson_r": pearson_r,
        "pearson_r2": pearson_r ** 2 if pd.notna(pearson_r) else np.nan,
        "spearman_rho": spearman_rho,
        "spearman_rho2": spearman_rho ** 2 if pd.notna(spearman_rho) else np.nan,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("--prediction-table", required=True)
    parser.add_argument("--truth-table", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--output-prefix", required=True)

    parser.add_argument("--prediction-name", default="prediction")
    parser.add_argument("--truth-name", default="truth")

    parser.add_argument("--prediction-chrom-column", default="chr")
    parser.add_argument("--prediction-position-column", default="pos_1based")
    parser.add_argument("--prediction-strand-column", default="strand")
    parser.add_argument("--prediction-depth-column", default="depth")
    parser.add_argument("--prediction-modified-count-column", default="modified_count")
    parser.add_argument("--prediction-stoich-column", default="stoich")
    parser.add_argument("--prediction-motif-column", default="motif")

    parser.add_argument("--truth-chrom-column", default="chr")
    parser.add_argument("--truth-position-column", default="pos_1based")
    parser.add_argument("--truth-strand-column", default="strand")
    parser.add_argument("--truth-depth-column", default="")
    parser.add_argument("--truth-stoich-column", required=True)
    parser.add_argument("--truth-motif-column", default="motif")

    parser.add_argument(
        "--truth-extra-columns",
        default="",
        help="Comma-separated extra truth columns to carry into matched output."
    )
    parser.add_argument(
        "--prediction-extra-columns",
        default="",
        help="Comma-separated extra prediction columns to carry into matched output."
    )

    parser.add_argument("--offsets", nargs="+", type=int, default=[-1, 0, 1])
    parser.add_argument("--offset-labels", nargs="+", default=None)
    parser.add_argument("--separator", default="\t")
    parser.add_argument("--force", action="store_true")

    return parser.parse_args()


def main() -> int:
    args = parse_args()

    started = datetime.now().isoformat(timespec="seconds")
    log(f"start_time\t{started}")
    log(f"prediction_table\t{args.prediction_table}")
    log(f"truth_table\t{args.truth_table}")
    log(f"output_dir\t{args.output_dir}")
    log(f"output_prefix\t{args.output_prefix}")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    offsets = parse_offsets(args.offsets, args.offset_labels)

    pred = read_table(args.prediction_table, args.separator)
    truth = read_table(args.truth_table, args.separator)

    require_columns(
        pred,
        [
            args.prediction_chrom_column,
            args.prediction_position_column,
            args.prediction_strand_column,
            args.prediction_depth_column,
            args.prediction_stoich_column,
        ],
        "prediction table",
    )
    require_columns(
        truth,
        [
            args.truth_chrom_column,
            args.truth_position_column,
            args.truth_strand_column,
            args.truth_stoich_column,
        ] + ([args.truth_depth_column] if args.truth_depth_column else []),
        "truth table",
    )

    has_pred_modified = args.prediction_modified_count_column in pred.columns
    has_pred_motif = args.prediction_motif_column in pred.columns
    has_truth_depth = bool(args.truth_depth_column) and args.truth_depth_column in truth.columns
    has_truth_motif = bool(args.truth_motif_column) and args.truth_motif_column in truth.columns

    pred_extra = [c.strip() for c in args.prediction_extra_columns.split(",") if c.strip() and c.strip() in pred.columns]
    truth_extra = [c.strip() for c in args.truth_extra_columns.split(",") if c.strip() and c.strip() in truth.columns]

    pred2 = pd.DataFrame(
        {
            "chr": pred[args.prediction_chrom_column].astype(str),
            "prediction_pos_1based": safe_numeric(pred[args.prediction_position_column], "prediction position"),
            "strand": pred[args.prediction_strand_column].astype(str),
            "prediction_depth": safe_numeric(pred[args.prediction_depth_column], "prediction depth"),
            "prediction_stoich": safe_numeric(pred[args.prediction_stoich_column], "prediction stoich"),
        }
    )

    if has_pred_modified:
        pred2["prediction_modified_count"] = safe_numeric(pred[args.prediction_modified_count_column], "prediction modified count")
    else:
        pred2["prediction_modified_count"] = pred2["prediction_depth"] * pred2["prediction_stoich"]

    if has_pred_motif:
        pred2["prediction_motif"] = pred[args.prediction_motif_column].astype(str)
    else:
        pred2["prediction_motif"] = "NA"

    for col in pred_extra:
        pred2[f"prediction_{col}"] = pred[col].astype(str)

    truth2 = pd.DataFrame(
        {
            "chr": truth[args.truth_chrom_column].astype(str),
            "truth_pos_1based": safe_numeric(truth[args.truth_position_column], "truth position"),
            "strand": truth[args.truth_strand_column].astype(str),
            "truth_stoich": safe_numeric(truth[args.truth_stoich_column], "truth stoich"),
        }
    )

    if has_truth_depth:
        truth2["truth_depth"] = safe_numeric(truth[args.truth_depth_column], "truth depth")
    else:
        truth2["truth_depth"] = np.nan

    if has_truth_motif:
        truth2["truth_motif"] = truth[args.truth_motif_column].astype(str)
    else:
        truth2["truth_motif"] = "NA"

    for col in truth_extra:
        truth2[f"truth_{col}"] = truth[col].astype(str)

    pred2 = pred2.dropna(
        subset=["chr", "prediction_pos_1based", "strand", "prediction_depth", "prediction_stoich"]
    ).copy()
    truth2 = truth2.dropna(
        subset=["chr", "truth_pos_1based", "strand", "truth_stoich"]
    ).copy()

    pred2["prediction_pos_1based"] = pred2["prediction_pos_1based"].astype(int)
    truth2["truth_pos_1based"] = truth2["truth_pos_1based"].astype(int)

    pred_agg = {
        "prediction_depth": ("prediction_depth", "sum"),
        "prediction_modified_count": ("prediction_modified_count", "sum"),
        "prediction_stoich": ("prediction_stoich", "mean"),
        "prediction_motif": ("prediction_motif", unique_join),
    }
    for col in pred_extra:
        pred_agg[f"prediction_{col}"] = (f"prediction_{col}", unique_join)

    pred2 = (
        pred2.groupby(["chr", "prediction_pos_1based", "strand"], as_index=False)
        .agg(**pred_agg)
    )
    mask = pred2["prediction_depth"] > 0
    pred2.loc[mask, "prediction_stoich"] = (
        pred2.loc[mask, "prediction_modified_count"] / pred2.loc[mask, "prediction_depth"]
    )

    truth_agg = {
        "truth_depth": ("truth_depth", "mean"),
        "truth_stoich": ("truth_stoich", "mean"),
        "truth_motif": ("truth_motif", unique_join),
    }
    for col in truth_extra:
        truth_agg[f"truth_{col}"] = (f"truth_{col}", unique_join)

    truth2 = (
        truth2.groupby(["chr", "truth_pos_1based", "strand"], as_index=False)
        .agg(**truth_agg)
    )

    summary_rows = []

    for offset, label in offsets:
        out_path = output_dir / f"{args.output_prefix}_offset_{label}.tsv"
        if out_path.exists() and not args.force:
            raise SystemExit(f"ERROR: output exists; use --force to overwrite: {out_path}")

        pred_shift = pred2.copy()
        pred_shift["truth_pos_1based"] = pred_shift["prediction_pos_1based"] - offset

        matched = truth2.merge(
            pred_shift,
            on=["chr", "truth_pos_1based", "strand"],
            how="inner",
        )

        matched["offset"] = offset
        matched["prediction_name"] = args.prediction_name
        matched["truth_name"] = args.truth_name

        preferred = [
            "chr",
            "truth_pos_1based",
            "prediction_pos_1based",
            "strand",
            "offset",
            "truth_name",
            "prediction_name",
            "truth_depth",
            "truth_stoich",
            "prediction_depth",
            "prediction_modified_count",
            "prediction_stoich",
            "truth_motif",
            "prediction_motif",
        ]
        remaining = [c for c in matched.columns if c not in preferred]
        matched = matched[preferred + remaining]

        matched = matched.sort_values(["chr", "truth_pos_1based", "strand"], kind="mergesort")
        matched.to_csv(out_path, sep=args.separator, index=False)

        stats = corr_stats(matched, "truth_stoich", "prediction_stoich")
        summary_rows.append(
            {
                "offset": offset,
                "offset_label": label,
                "overlap_n": stats["n"],
                "pearson_r": stats["pearson_r"],
                "pearson_r2": stats["pearson_r2"],
                "spearman_rho": stats["spearman_rho"],
                "spearman_rho2": stats["spearman_rho2"],
                "matched_table": str(out_path),
            }
        )

        log(f"offset\t{offset}\tlabel\t{label}\toverlap\t{stats['n']}\toutput\t{out_path}")

    summary = pd.DataFrame(summary_rows)
    summary_path = output_dir / f"{args.output_prefix}_offset_summary.tsv"
    if summary_path.exists() and not args.force:
        raise SystemExit(f"ERROR: output exists; use --force to overwrite: {summary_path}")
    summary.to_csv(summary_path, sep=args.separator, index=False)

    log(f"summary\t{summary_path}")
    log(f"finish_time\t{datetime.now().isoformat(timespec='seconds')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
