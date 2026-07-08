#!/usr/bin/env python3
"""Calculate Pearson R² and Spearman rho² for matched site-level tables."""

from __future__ import annotations

import argparse
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd


def log(message: str) -> None:
    print(message, file=sys.stderr, flush=True)


def read_tsv(path: str) -> pd.DataFrame:
    if not Path(path).exists():
        raise FileNotFoundError(f"input file does not exist: {path}")
    return pd.read_csv(path, sep="\t", compression="infer", dtype=str, keep_default_na=False, na_values=[])


def corr_stats(df: pd.DataFrame, x_col: str, y_col: str) -> dict[str, object]:
    sub = df[[x_col, y_col]].copy()
    sub[x_col] = pd.to_numeric(sub[x_col], errors="coerce")
    sub[y_col] = pd.to_numeric(sub[y_col], errors="coerce")
    sub = sub.dropna()
    sub = sub[np.isfinite(sub[x_col]) & np.isfinite(sub[y_col])].copy()

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
    parser.add_argument("--input", required=True, help="Matched TSV table")
    parser.add_argument("--output", required=True, help="Output summary TSV")
    parser.add_argument("--x-column", default="truth_stoich")
    parser.add_argument("--y-column", default="prediction_stoich")
    parser.add_argument("--group-column", default="", help="Optional group column, e.g. depth_group")
    parser.add_argument("--groups", nargs="+", default=None, help="Optional group order")
    parser.add_argument("--allow-empty-groups", action="store_true")
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    output = Path(args.output)
    if output.exists() and not args.force:
        raise SystemExit(f"ERROR: output exists; use --force to overwrite: {output}")
    output.parent.mkdir(parents=True, exist_ok=True)

    log(f"start_time\t{datetime.now().isoformat(timespec='seconds')}")
    log(f"input\t{args.input}")
    log(f"output\t{args.output}")
    log(f"x_column\t{args.x_column}")
    log(f"y_column\t{args.y_column}")
    log(f"group_column\t{args.group_column}")

    df = read_tsv(args.input)

    for col in [args.x_column, args.y_column]:
        if col not in df.columns:
            raise SystemExit(f"ERROR: missing required column: {col}\nExisting columns: {list(df.columns)}")

    rows = []

    if args.group_column:
        if args.group_column not in df.columns:
            raise SystemExit(f"ERROR: missing group column: {args.group_column}")

        if args.groups:
            groups = list(args.groups)
        else:
            groups = list(dict.fromkeys(df[args.group_column].astype(str).tolist()))

        for group in groups:
            sub = df.loc[df[args.group_column].astype(str) == str(group)].copy()
            if len(sub) == 0 and not args.allow_empty_groups:
                raise SystemExit(f"ERROR: group has no rows: {group}")
            stats = corr_stats(sub, args.x_column, args.y_column)
            rows.append(
                {
                    "group": group,
                    "n": stats["n"],
                    "pearson_r": stats["pearson_r"],
                    "pearson_r2": stats["pearson_r2"],
                    "spearman_rho": stats["spearman_rho"],
                    "spearman_rho2": stats["spearman_rho2"],
                }
            )
    else:
        stats = corr_stats(df, args.x_column, args.y_column)
        rows.append(
            {
                "group": "all",
                "n": stats["n"],
                "pearson_r": stats["pearson_r"],
                "pearson_r2": stats["pearson_r2"],
                "spearman_rho": stats["spearman_rho"],
                "spearman_rho2": stats["spearman_rho2"],
            }
        )

    out = pd.DataFrame(rows)
    out.to_csv(output, sep="\t", index=False)

    log(out.to_string(index=False))
    log(f"finish_time\t{datetime.now().isoformat(timespec='seconds')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
