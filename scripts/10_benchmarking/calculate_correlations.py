#!/usr/bin/env python3
"""Calculate stoichiometry correlations for matched benchmark sites."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import os
import sys
from collections import defaultdict
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


def parse_stoich(path: str, line_number: int, column: str, value: str) -> float:
    try:
        parsed = float(value)
    except ValueError:
        raise SystemExit("ERROR: %s line %d column %s: must be numeric" % (path, line_number, column))
    if not math.isfinite(parsed) or parsed < 0.0 or parsed > 1.0:
        raise SystemExit("ERROR: %s line %d column %s: must be finite and between 0 and 1" % (path, line_number, column))
    return parsed


def average_ranks(values: Sequence[float]) -> List[float]:
    indexed = sorted((value, index) for index, value in enumerate(values))
    ranks = [0.0] * len(values)
    i = 0
    while i < len(indexed):
        j = i + 1
        while j < len(indexed) and indexed[j][0] == indexed[i][0]:
            j += 1
        avg_rank = (i + 1 + j) / 2.0
        for k in range(i, j):
            ranks[indexed[k][1]] = avg_rank
        i = j
    return ranks


def pearson(x: Sequence[float], y: Sequence[float]) -> Optional[float]:
    n = len(x)
    if n < 2:
        return None
    mean_x = sum(x) / n
    mean_y = sum(y) / n
    dx = [value - mean_x for value in x]
    dy = [value - mean_y for value in y]
    ss_x = sum(value * value for value in dx)
    ss_y = sum(value * value for value in dy)
    if ss_x == 0.0 or ss_y == 0.0:
        return None
    return sum(a * b for a, b in zip(dx, dy)) / math.sqrt(ss_x * ss_y)


def spearman(x: Sequence[float], y: Sequence[float]) -> Optional[float]:
    if len(x) < 2:
        return None
    return pearson(average_ranks(x), average_ranks(y))


def format_value(value: Optional[float]) -> str:
    if value is None:
        return "NA"
    return format(value, ".10g")


def group_sort_key(group: str, rows: Sequence[Tuple[float, float, Optional[int]]]) -> Tuple[int, str]:
    depths = [depth for _truth, _pred, depth in rows if depth is not None]
    if depths:
        return (min(depths), group)
    if group.startswith("depth_") and "_" in group[6:]:
        try:
            return (int(group.split("_")[1]), group)
        except Exception:
            pass
    if group.startswith("depth_gt"):
        try:
            return (int(group.replace("depth_gt", "")) + 1, group)
        except Exception:
            pass
    return (10 ** 12, group)


def depth_rule(group: str) -> str:
    if group == "all":
        return "all matched sites"
    if group.startswith("depth_gt"):
        return "depth > %s" % group.replace("depth_gt", "")
    parts = group.split("_")
    if len(parts) == 3 and parts[0] == "depth":
        return "%s <= depth <= %s" % (parts[1], parts[2])
    return group


def summarize(group: str, rows: Sequence[Tuple[float, float, Optional[int]]]) -> List[str]:
    n = len(rows)
    truth = [row[0] for row in rows]
    pred = [row[1] for row in rows]
    pr = pearson(truth, pred)
    sr = spearman(truth, pred)
    if n < 2:
        log("warning\t%s\tfewer than two observations; correlations set to NA" % group)
    elif pr is None:
        log("warning\t%s\tconstant vector for Pearson; correlations set to NA" % group)
    elif sr is None:
        log("warning\t%s\tconstant rank vector for Spearman; correlations set to NA" % group)
    return [
        group,
        depth_rule(group),
        str(n),
        format_value(pr),
        format_value(None if pr is None else pr * pr),
        format_value(sr),
        format_value(None if sr is None else sr * sr),
        format_value(sum(truth) / n if n else None),
        format_value(sum(pred) / n if n else None),
        format_value(min(truth) if n else None),
        format_value(max(truth) if n else None),
        format_value(min(pred) if n else None),
        format_value(max(pred) if n else None),
    ]


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--truth-column", default="truth_stoich")
    parser.add_argument("--prediction-column", default="prediction_stoich")
    parser.add_argument("--group-column", default="depth_group")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)

    output = Path(args.output)
    if output.exists() and not args.force:
        raise SystemExit("ERROR: output exists; use --force to overwrite: %s" % output)
    output.parent.mkdir(parents=True, exist_ok=True)
    tmp = atomic_path(output)

    groups = defaultdict(list)  # type: Dict[str, List[Tuple[float, float, Optional[int]]]]
    all_rows = []  # type: List[Tuple[float, float, Optional[int]]]
    log("start_time\t%s" % now())
    log("input_path\t%s" % args.input)
    with open_text(args.input) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise SystemExit("ERROR: input table is missing a header")
        for col in (args.truth_column, args.prediction_column):
            if col not in reader.fieldnames:
                raise SystemExit("ERROR: missing required column: %s" % col)
        has_group = args.group_column in reader.fieldnames
        has_depth = "prediction_depth" in reader.fieldnames
        for line_number, row in enumerate(reader, start=2):
            truth = parse_stoich(args.input, line_number, args.truth_column, row[args.truth_column])
            pred = parse_stoich(args.input, line_number, args.prediction_column, row[args.prediction_column])
            depth = None
            if has_depth and row["prediction_depth"] != "":
                try:
                    depth = int(row["prediction_depth"])
                except ValueError:
                    depth = None
            item = (truth, pred, depth)
            all_rows.append(item)
            if has_group:
                group = row[args.group_column]
                if not group:
                    raise SystemExit("ERROR: %s line %d column %s: group label is empty" % (args.input, line_number, args.group_column))
                groups[group].append(item)
    try:
        with open(tmp, "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow([
                "group", "depth_rule", "n", "pearson_r", "pearson_r2", "spearman_rho",
                "spearman_rho2", "truth_stoich_mean", "prediction_stoich_mean",
                "truth_stoich_min", "truth_stoich_max", "prediction_stoich_min",
                "prediction_stoich_max",
            ])
            rows = [summarize("all", all_rows)]
            for group in sorted(groups, key=lambda g: group_sort_key(g, groups[g])):
                rows.append(summarize(group, groups[group]))
            writer.writerows(rows)
        os.replace(tmp, output)
        log("total_matched_site_count\t%d" % len(all_rows))
        for row in rows:
            log("group\t%s\tn\t%s\tpearson_r\t%s\tspearman_rho\t%s" % (row[0], row[2], row[3], row[5]))
        log("output_path\t%s" % output)
        log("finish_time\t%s" % now())
        return 0
    finally:
        try:
            if tmp.exists():
                tmp.unlink()
        except OSError:
            pass


if __name__ == "__main__":
    raise SystemExit(main())
