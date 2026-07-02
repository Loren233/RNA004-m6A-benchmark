#!/usr/bin/env python3
"""Match prediction sites to ground-truth sites under genomic offsets."""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import sqlite3
import sys
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, TextIO, Tuple


MatchKey = Tuple[str, int, str]


def log(message: str) -> None:
    print(message, file=sys.stderr, flush=True)


def now() -> str:
    return datetime.now().isoformat(timespec="seconds")


def open_text(path: str) -> TextIO:
    if path.endswith(".gz"):
        return gzip.open(path, "rt", newline="")
    return open(path, "r", newline="")


def offset_label(offset: int) -> str:
    if offset < 0:
        return "m%d" % abs(offset)
    if offset > 0:
        return "p%d" % offset
    return "0"


def output_paths(prefix: str, offsets: Sequence[int]) -> Dict[int, Path]:
    return {offset: Path("%s_offset_%s.tsv" % (prefix, offset_label(offset))) for offset in offsets}


def summary_path(prefix: str) -> Path:
    return Path("%s_offset_summary.tsv" % prefix)


def atomic_path(path: Path) -> Path:
    return Path(str(path) + ".tmp.%d" % os.getpid())


def require_no_overwrite(paths: Iterable[Path], force: bool) -> None:
    if force:
        return
    existing = [str(path) for path in paths if path.exists()]
    if existing:
        raise SystemExit("ERROR: output exists; use --force to overwrite: %s" % ", ".join(existing))


def reader_with_lines(path: str):
    handle = open_text(path)
    try:
        header_line = handle.readline()
        if not header_line:
            raise SystemExit("ERROR: empty input file: %s" % path)
        header = header_line.rstrip("\n\r").split("\t")
        for line_number, line in enumerate(handle, start=2):
            stripped = line.rstrip("\n\r")
            if not stripped:
                yield header, line_number, None
                continue
            fields = stripped.split("\t")
            if len(fields) != len(header):
                raise SystemExit(
                    "ERROR: %s line %d: expected %d fields, found %d"
                    % (path, line_number, len(header), len(fields))
                )
            yield header, line_number, dict(zip(header, fields))
    finally:
        handle.close()


def require_columns(path: str, header: Sequence[str], columns: Iterable[str]) -> None:
    missing = [column for column in columns if column not in header]
    if missing:
        raise SystemExit("ERROR: %s missing required column(s): %s" % (path, ", ".join(missing)))


def parse_int(path: str, line_number: int, column: str, value: str, positive: bool, nonnegative: bool = False) -> int:
    try:
        parsed = int(value)
    except ValueError:
        raise SystemExit("ERROR: %s line %d column %s: must be an integer" % (path, line_number, column))
    if positive and parsed <= 0:
        raise SystemExit("ERROR: %s line %d column %s: must be positive" % (path, line_number, column))
    if nonnegative and parsed < 0:
        raise SystemExit("ERROR: %s line %d column %s: must be non-negative" % (path, line_number, column))
    return parsed


def parse_float(path: str, line_number: int, column: str, value: str, minimum: float, maximum: float) -> float:
    try:
        parsed = float(value)
    except ValueError:
        raise SystemExit("ERROR: %s line %d column %s: must be numeric" % (path, line_number, column))
    if parsed != parsed or parsed < minimum or parsed > maximum:
        raise SystemExit(
            "ERROR: %s line %d column %s: must be between %s and %s"
            % (path, line_number, column, minimum, maximum)
        )
    return parsed


def load_truth(args: argparse.Namespace) -> Tuple[Dict[MatchKey, Dict[str, str]], int, int]:
    truth = {}  # type: Dict[MatchKey, Dict[str, str]]
    blank_lines = 0
    total = 0
    first = True
    for header, line_number, row in reader_with_lines(args.truth):
        if first:
            require_columns(
                args.truth,
                header,
                [
                    args.truth_chromosome_column,
                    args.truth_position_column,
                    args.truth_strand_column,
                    args.truth_depth_column,
                    args.truth_stoich_column,
                ],
            )
            first = False
        if row is None:
            blank_lines += 1
            continue
        chrom = row[args.truth_chromosome_column]
        if not chrom:
            raise SystemExit("ERROR: %s line %d column %s: must be non-empty" % (args.truth, line_number, args.truth_chromosome_column))
        pos = parse_int(args.truth, line_number, args.truth_position_column, row[args.truth_position_column], True)
        strand = row[args.truth_strand_column]
        if strand not in ("+", "-"):
            raise SystemExit("ERROR: %s line %d column %s: must be + or -" % (args.truth, line_number, args.truth_strand_column))
        parse_float(args.truth, line_number, args.truth_depth_column, row[args.truth_depth_column], 0.0, float("inf"))
        parse_float(args.truth, line_number, args.truth_stoich_column, row[args.truth_stoich_column], 0.0, 1.0)
        key = (chrom, pos, strand)
        if key in truth:
            raise SystemExit("ERROR: %s line %d: duplicate truth key %s:%d:%s" % (args.truth, line_number, chrom, pos, strand))
        truth[key] = row
        total += 1
    return truth, total, blank_lines


def validate_prediction(row: Dict[str, str], line_number: int, args: argparse.Namespace) -> Tuple[MatchKey, Dict[str, str]]:
    chrom = row[args.prediction_chromosome_column]
    if not chrom:
        raise SystemExit("ERROR: %s line %d column %s: must be non-empty" % (args.prediction, line_number, args.prediction_chromosome_column))
    pos = parse_int(args.prediction, line_number, args.prediction_position_column, row[args.prediction_position_column], True)
    strand = row[args.prediction_strand_column]
    if strand not in ("+", "-"):
        raise SystemExit("ERROR: %s line %d column %s: must be + or -" % (args.prediction, line_number, args.prediction_strand_column))
    depth = parse_int(args.prediction, line_number, args.prediction_depth_column, row[args.prediction_depth_column], True)
    parse_float(args.prediction, line_number, args.prediction_stoich_column, row[args.prediction_stoich_column], 0.0, 1.0)
    if args.prediction_modified_count_column in row:
        modified = parse_int(args.prediction, line_number, args.prediction_modified_count_column, row[args.prediction_modified_count_column], False, True)
        if modified > depth:
            raise SystemExit("ERROR: %s line %d column %s: cannot exceed depth" % (args.prediction, line_number, args.prediction_modified_count_column))
    if args.prediction_transcript_count_column in row:
        transcript_count = parse_int(args.prediction, line_number, args.prediction_transcript_count_column, row[args.prediction_transcript_count_column], True)
    else:
        transcript_count = None
    if args.prediction_record_count_column in row:
        record_count = parse_int(args.prediction, line_number, args.prediction_record_count_column, row[args.prediction_record_count_column], True)
        if transcript_count is not None and record_count < transcript_count:
            raise SystemExit("ERROR: %s line %d column %s: cannot be smaller than transcript_count" % (args.prediction, line_number, args.prediction_record_count_column))
    return (chrom, pos, strand), row


def write_header(writer: csv.writer) -> None:
    writer.writerow([
        "chr", "truth_pos_1based", "prediction_pos_1based", "strand", "offset",
        "truth_depth", "truth_stoich", "prediction_depth", "prediction_modified_count",
        "prediction_stoich", "transcript_count", "transcript_record_count",
        "prediction_source", "truth_source",
    ])


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prediction", required=True)
    parser.add_argument("--truth", required=True)
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--offsets", nargs="+", type=int, default=[-1, 0, 1])
    parser.add_argument("--prediction-chromosome-column", default="chr")
    parser.add_argument("--prediction-position-column", default="pos_1based")
    parser.add_argument("--prediction-strand-column", default="strand")
    parser.add_argument("--prediction-depth-column", default="depth")
    parser.add_argument("--prediction-stoich-column", default="stoich")
    parser.add_argument("--prediction-modified-count-column", default="modified_count")
    parser.add_argument("--prediction-transcript-count-column", default="transcript_count")
    parser.add_argument("--prediction-record-count-column", default="transcript_record_count")
    parser.add_argument("--truth-chromosome-column", default="chr")
    parser.add_argument("--truth-position-column", default="pos_1based")
    parser.add_argument("--truth-strand-column", default="strand")
    parser.add_argument("--truth-depth-column", default="depth")
    parser.add_argument("--truth-stoich-column", default="stoich")
    parser.add_argument("--prediction-source", default="prediction")
    parser.add_argument("--truth-source", default="truth")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)

    offsets = sorted(set(args.offsets))
    paths = output_paths(args.out_prefix, offsets)
    sum_path = summary_path(args.out_prefix)
    require_no_overwrite(list(paths.values()) + [sum_path], args.force)
    for path in list(paths.values()) + [sum_path]:
        path.parent.mkdir(parents=True, exist_ok=True)

    log("start_time\t%s" % now())
    log("prediction_path\t%s" % args.prediction)
    log("truth_path\t%s" % args.truth)
    log("offsets\t%s" % ",".join(str(o) for o in offsets))

    temp_paths = {offset: atomic_path(path) for offset, path in paths.items()}
    temp_summary = atomic_path(sum_path)
    handles = []
    db_path = str(Path(tempfile.gettempdir()) / ("prediction_keys_%d.sqlite" % os.getpid()))
    try:
        truth, truth_total, truth_blank = load_truth(args)
        log("loaded_truth_site_count\t%d" % truth_total)
        writers = {}
        for offset in offsets:
            handle = open(temp_paths[offset], "w", newline="")
            handles.append(handle)
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            write_header(writer)
            writers[offset] = writer

        conn = sqlite3.connect(db_path)
        conn.execute("CREATE TABLE keys (chr TEXT, pos INTEGER, strand TEXT, UNIQUE(chr, pos, strand))")
        prediction_total = 0
        prediction_blank = 0
        overlap_counts = {offset: 0 for offset in offsets}
        first = True
        for header, line_number, row in reader_with_lines(args.prediction):
            if first:
                require_columns(
                    args.prediction,
                    header,
                    [
                        args.prediction_chromosome_column,
                        args.prediction_position_column,
                        args.prediction_strand_column,
                        args.prediction_depth_column,
                        args.prediction_stoich_column,
                    ],
                )
                first = False
            if row is None:
                prediction_blank += 1
                continue
            key, row = validate_prediction(row, line_number, args)
            try:
                conn.execute("INSERT INTO keys VALUES (?, ?, ?)", key)
            except sqlite3.IntegrityError:
                raise SystemExit("ERROR: %s line %d: duplicate prediction key %s:%d:%s" % (args.prediction, line_number, key[0], key[1], key[2]))
            prediction_total += 1
            if prediction_total % 5000000 == 0:
                log("prediction_rows_processed\t%d" % prediction_total)
            chrom, pred_pos, strand = key
            for offset in offsets:
                truth_key = (chrom, pred_pos - offset, strand)
                truth_row = truth.get(truth_key)
                if truth_row is None:
                    continue
                writers[offset].writerow([
                    chrom,
                    str(truth_key[1]),
                    str(pred_pos),
                    strand,
                    str(offset),
                    truth_row[args.truth_depth_column],
                    truth_row[args.truth_stoich_column],
                    row[args.prediction_depth_column],
                    row.get(args.prediction_modified_count_column, ""),
                    row[args.prediction_stoich_column],
                    row.get(args.prediction_transcript_count_column, ""),
                    row.get(args.prediction_record_count_column, ""),
                    args.prediction_source,
                    args.truth_source,
                ])
                overlap_counts[offset] += 1
        conn.close()
        for handle in handles:
            handle.close()

        with open(temp_summary, "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow([
                "offset", "definition", "overlap_count", "truth_total_sites",
                "overlap_fraction_of_truth", "prediction_total_sites",
                "overlap_fraction_of_prediction", "matched_file",
            ])
            for offset in offsets:
                count = overlap_counts[offset]
                writer.writerow([
                    str(offset),
                    "prediction_pos=truth_pos+offset",
                    str(count),
                    str(truth_total),
                    format(float(count) / truth_total if truth_total else 0.0, ".10g"),
                    str(prediction_total),
                    format(float(count) / prediction_total if prediction_total else 0.0, ".10g"),
                    str(paths[offset]),
                ])

        for offset, path in paths.items():
            os.replace(temp_paths[offset], path)
        os.replace(temp_summary, sum_path)
        log("total_prediction_site_count\t%d" % prediction_total)
        for offset in offsets:
            log("overlap_count_offset_%s\t%d" % (offset_label(offset), overlap_counts[offset]))
        log("prediction_blank_lines\t%d" % prediction_blank)
        log("truth_blank_lines\t%d" % truth_blank)
        log("summary_path\t%s" % sum_path)
        log("finish_time\t%s" % now())
        return 0
    finally:
        for handle in handles:
            try:
                handle.close()
            except Exception:
                pass
        for path in list(temp_paths.values()) + [temp_summary, Path(db_path)]:
            try:
                if Path(path).exists():
                    Path(path).unlink()
            except OSError:
                pass


if __name__ == "__main__":
    raise SystemExit(main())
