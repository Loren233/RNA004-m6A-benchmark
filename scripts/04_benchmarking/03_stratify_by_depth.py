#!/usr/bin/env python3
"""Assign matched sites to non-overlapping prediction-depth groups."""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import sys
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Sequence, TextIO, Tuple


def log(message: str) -> None:
    print(message, file=sys.stderr, flush=True)


def now() -> str:
    return datetime.now().isoformat(timespec="seconds")


def open_text(path: str) -> TextIO:
    if path.endswith(".gz"):
        return gzip.open(path, "rt", newline="")
    return open(path, "r", newline="")


def atomic_path(path: Path) -> Path:
    return Path(str(path) + ".tmp.%d" % os.getpid())


def summary_path(output: Path) -> Path:
    return output.with_name(output.stem + "_depth_summary.tsv")


def labels_and_rules(edges: Sequence[int]) -> List[Tuple[str, str, int, Optional[int]]]:
    if not edges:
        return [("depth_gt0", "depth >= 1", 1, None)]
    sorted_edges = sorted(edges)
    if sorted_edges[0] < 1:
        raise SystemExit("ERROR: bin edges must be positive integers")
    groups = []  # type: List[Tuple[str, str, int, Optional[int]]]
    start = 1
    for edge in sorted_edges:
        if edge < start:
            raise SystemExit("ERROR: bin edges must be strictly increasing")
        label = "depth_%d_%d" % (start, edge)
        rule = "%d <= depth <= %d" % (start, edge)
        groups.append((label, rule, start, edge))
        start = edge + 1
    groups.append(("depth_gt%d" % sorted_edges[-1], "depth > %d" % sorted_edges[-1], start, None))
    return groups


def assign_group(depth: int, groups: Sequence[Tuple[str, str, int, Optional[int]]]) -> str:
    for label, _rule, start, end in groups:
        if depth >= start and (end is None or depth <= end):
            return label
    raise SystemExit("ERROR: failed to assign depth %d to a group" % depth)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--depth-column", default="prediction_depth")
    parser.add_argument("--group-column", default="depth_group")
    parser.add_argument("--bin-edges", nargs="+", type=int, default=[10, 20, 30, 40])
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)

    output = Path(args.output)
    summary = summary_path(output)
    if not args.force:
        existing = [str(path) for path in (output, summary) if path.exists()]
        if existing:
            raise SystemExit("ERROR: output exists; use --force to overwrite: %s" % ", ".join(existing))
    output.parent.mkdir(parents=True, exist_ok=True)
    groups = labels_and_rules(args.bin_edges)
    stats = {
        label: {"count": 0, "min": None, "max": None, "rule": rule}
        for label, rule, _start, _end in groups
    }  # type: Dict[str, Dict[str, object]]
    total = 0
    tmp_output = atomic_path(output)
    tmp_summary = atomic_path(summary)

    log("start_time\t%s" % now())
    log("input_path\t%s" % args.input)
    log("output_path\t%s" % output)
    log("depth_column\t%s" % args.depth_column)
    log("bin_edges\t%s" % ",".join(str(edge) for edge in sorted(args.bin_edges)))

    try:
        with open_text(args.input) as in_handle, open(tmp_output, "w", newline="") as out_handle:
            reader = csv.DictReader(in_handle, delimiter="\t")
            if reader.fieldnames is None:
                raise SystemExit("ERROR: input table is missing a header")
            if args.depth_column not in reader.fieldnames:
                raise SystemExit("ERROR: missing depth column: %s" % args.depth_column)
            if args.group_column in reader.fieldnames:
                raise SystemExit("ERROR: group column already exists: %s" % args.group_column)
            fieldnames = list(reader.fieldnames) + [args.group_column]
            writer = csv.DictWriter(out_handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for line_number, row in enumerate(reader, start=2):
                try:
                    depth = int(row[args.depth_column])
                except ValueError:
                    raise SystemExit("ERROR: %s line %d column %s: depth must be an integer" % (args.input, line_number, args.depth_column))
                if depth <= 0:
                    raise SystemExit("ERROR: %s line %d column %s: depth must be positive" % (args.input, line_number, args.depth_column))
                label = assign_group(depth, groups)
                row[args.group_column] = label
                writer.writerow(row)
                total += 1
                item = stats[label]
                item["count"] = int(item["count"]) + 1
                item["min"] = depth if item["min"] is None else min(int(item["min"]), depth)
                item["max"] = depth if item["max"] is None else max(int(item["max"]), depth)

        with open(tmp_summary, "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["depth_group", "depth_rule", "site_count", "fraction_of_all_sites", "minimum_depth", "maximum_depth"])
            check_total = 0
            for label, rule, _start, _end in groups:
                item = stats[label]
                count = int(item["count"])
                check_total += count
                writer.writerow([
                    label,
                    rule,
                    str(count),
                    format(float(count) / total if total else 0.0, ".10g"),
                    "" if item["min"] is None else str(item["min"]),
                    "" if item["max"] is None else str(item["max"]),
                ])
        if check_total != total:
            raise SystemExit("ERROR: depth-group counts do not sum to total input rows")
        os.replace(tmp_output, output)
        os.replace(tmp_summary, summary)
        log("total_input_site_count\t%d" % total)
        for label, _rule, _start, _end in groups:
            log("%s\t%d" % (label, int(stats[label]["count"])))
        log("finish_time\t%s" % now())
        return 0
    finally:
        for path in (tmp_output, tmp_summary):
            try:
                if path.exists():
                    path.unlink()
            except OSError:
                pass


if __name__ == "__main__":
    raise SystemExit(main())
