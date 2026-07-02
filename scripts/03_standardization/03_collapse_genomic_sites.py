#!/usr/bin/env python3

"""
Collapse duplicated genomic sites after standardization.

When to use:
  Use this script AFTER tool-specific results have been standardized to genomic
  coordinates, and BEFORE benchmarking against ground truth datasets.

Why:
  Transcript-level tools may output multiple transcript/isoform rows that map to
  the same genomic position. If these duplicate genomic sites are not collapsed,
  overlap and correlation analyses can be biased.

Typical use cases:
  - SingleMod transcript-level output converted to genome coordinates
  - m6Anet transcript-level output converted to genome coordinates
  - Any Dorado / DeepRM / other site-level table that contains duplicated
    chr + position + strand rows after standardization

Input requirement:
  The input table should already contain genomic coordinate columns, such as:
    chr, pos_1based, strand

Collapse rule:
  Group by:
    chr + pos_1based + strand

  Then:
    depth: sum
    modified_count: sum
    stoich: modified_count / depth if both are available
    motif: unique motifs joined by comma
    source_rows: number of input rows collapsed into the genomic site
"""

import argparse
import sys
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Collapse duplicated genomic sites in standardized site-level tables."
    )

    parser.add_argument("--input", required=True, help="Input standardized TSV file")
    parser.add_argument("--output", required=True, help="Output collapsed TSV file")

    parser.add_argument("--chrom-column", default="chr", help="Chromosome column name")
    parser.add_argument("--position-column", default="pos_1based", help="1-based genomic position column name")
    parser.add_argument("--strand-column", default="strand", help="Strand column name")

    parser.add_argument("--depth-column", default="depth", help="Depth / total reads column name")
    parser.add_argument("--modified-count-column", default="modified_count", help="Modified read count column name")
    parser.add_argument("--stoich-column", default="stoich", help="Modification stoichiometry column name")

    parser.add_argument("--motif-column", default="motif", help="Motif column name")
    parser.add_argument("--site-id-column", default="site_id", help="Output site ID column name")

    parser.add_argument(
        "--separator",
        default="\t",
        help="Input/output separator. Default: tab"
    )

    parser.add_argument(
        "--keep-columns",
        default="",
        help=(
            "Optional comma-separated columns to keep by unique joining. "
            "Example: transcript_id,tool,input_site_id"
        )
    )

    parser.add_argument(
        "--drop-zero-depth",
        action="store_true",
        help="Drop collapsed sites with depth <= 0"
    )

    return parser.parse_args()


def require_columns(df, columns):
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise SystemExit(f"ERROR: missing required columns: {missing}")


def unique_join(series):
    values = []
    for x in series:
        if pd.isna(x):
            continue
        s = str(x)
        if s == "":
            continue
        values.append(s)
    return ",".join(sorted(set(values)))


def main():
    args = parse_args()

    df = pd.read_csv(args.input, sep=args.separator)

    key_cols = [args.chrom_column, args.position_column, args.strand_column]
    require_columns(df, key_cols)

    has_depth = args.depth_column in df.columns
    has_modified = args.modified_count_column in df.columns
    has_stoich = args.stoich_column in df.columns
    has_motif = args.motif_column in df.columns

    if not has_depth and not has_stoich:
        raise SystemExit(
            "ERROR: input must contain at least a depth column or a stoich column."
        )

    if has_depth:
        df[args.depth_column] = pd.to_numeric(df[args.depth_column], errors="coerce").fillna(0)

    if has_modified:
        df[args.modified_count_column] = pd.to_numeric(df[args.modified_count_column], errors="coerce").fillna(0)

    if has_stoich:
        df[args.stoich_column] = pd.to_numeric(df[args.stoich_column], errors="coerce").fillna(0)

    if not has_modified and has_depth and has_stoich:
        args.modified_count_column = "__inferred_modified_count"
        df[args.modified_count_column] = df[args.depth_column] * df[args.stoich_column]
        has_modified = True

    keep_columns = [c.strip() for c in args.keep_columns.split(",") if c.strip()]
    keep_columns = [c for c in keep_columns if c in df.columns and c not in key_cols]

    agg = {}

    if has_depth:
        agg[args.depth_column] = "sum"

    if has_modified:
        agg[args.modified_count_column] = "sum"

    if has_motif:
        agg[args.motif_column] = unique_join

    for col in keep_columns:
        agg[col] = unique_join

    collapsed = (
        df.groupby(key_cols, as_index=False)
          .agg(**{col: pd.NamedAgg(column=col, aggfunc=func) for col, func in agg.items()})
    )

    source_rows = (
        df.groupby(key_cols, as_index=False)
          .size()
          .rename(columns={"size": "source_rows"})
    )

    collapsed = collapsed.merge(source_rows, on=key_cols, how="left")

    if has_depth and has_modified:
        collapsed[args.stoich_column] = collapsed.apply(
            lambda r: r[args.modified_count_column] / r[args.depth_column]
            if r[args.depth_column] > 0 else 0,
            axis=1
        )
    elif has_stoich:
        mean_stoich = (
            df.groupby(key_cols, as_index=False)[args.stoich_column]
              .mean()
        )
        collapsed = collapsed.merge(mean_stoich, on=key_cols, how="left")

    if args.drop_zero_depth and has_depth:
        collapsed = collapsed[collapsed[args.depth_column] > 0].copy()

    collapsed["start_0based"] = pd.to_numeric(
        collapsed[args.position_column],
        errors="coerce"
    ).astype("Int64") - 1

    collapsed["end_0based_exclusive"] = pd.to_numeric(
        collapsed[args.position_column],
        errors="coerce"
    ).astype("Int64")

    if has_motif:
        motif_for_id = collapsed[args.motif_column].astype(str)
    else:
        motif_for_id = "NA"

    collapsed[args.site_id_column] = (
        collapsed[args.chrom_column].astype(str) + "|" +
        collapsed[args.position_column].astype(str) + "|" +
        collapsed[args.strand_column].astype(str) + "|" +
        motif_for_id
    )

    preferred = [
        args.chrom_column,
        args.position_column,
        "start_0based",
        "end_0based_exclusive",
        args.strand_column,
        args.site_id_column,
    ]

    if has_modified:
        preferred.append(args.modified_count_column)

    if has_depth:
        preferred.append(args.depth_column)

    if args.stoich_column in collapsed.columns:
        preferred.append(args.stoich_column)

    if has_motif:
        preferred.append(args.motif_column)

    preferred.append("source_rows")

    for col in keep_columns:
        if col not in preferred:
            preferred.append(col)

    remaining = [c for c in collapsed.columns if c not in preferred]
    collapsed = collapsed[preferred + remaining]

    collapsed = collapsed.sort_values(
        by=[args.chrom_column, args.position_column, args.strand_column],
        kind="mergesort"
    )

    input_rows = len(df)
    output_rows = len(collapsed)
    removed_rows = input_rows - output_rows

    collapsed.to_csv(args.output, sep=args.separator, index=False)

    print(f"Input rows: {input_rows}", file=sys.stderr)
    print(f"Collapsed rows: {output_rows}", file=sys.stderr)
    print(f"Removed duplicate rows: {removed_rows}", file=sys.stderr)
    print(f"Output: {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
