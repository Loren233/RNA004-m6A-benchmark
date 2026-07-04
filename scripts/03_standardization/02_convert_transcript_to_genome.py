#!/usr/bin/env python3
"""
Convert transcript-coordinate site tables to genomic coordinates.

This script only performs coordinate conversion.

It does not standardize tool-specific columns. Tool-specific column standardization
should be handled before this step by:

    01_standardize_site_tables.py

Expected workflow:

    tool-specific output
        -> 01_standardize_site_tables.py
        -> 02_convert_transcript_to_genome.py
        -> 03_collapse_genomic_sites.py

Input requirement:
    The input table must contain transcript-level coordinates, usually:

        transcript_id
        transcript_pos_1based

Output:
    The output table contains genomic coordinate columns:

        chr
        pos_1based
        start_0based
        end_0based_exclusive
        strand
        transcript_id
        transcript_position_0based

    All non-coordinate input columns are preserved whenever possible.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import sys
from bisect import bisect_right
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, TextIO, Tuple


Exon = Tuple[int, int]


class ConversionError(Exception):
    pass


class TranscriptModel:
    def __init__(self, chrom: str, strand: str, exons: List[Exon]) -> None:
        self.chrom = chrom
        self.strand = strand
        self.exons = exons
        self.cumulative_ends: List[int] = []

        total = 0
        for exon_start, exon_end in exons:
            total += exon_end - exon_start + 1
            self.cumulative_ends.append(total)

    @property
    def length(self) -> int:
        if not self.cumulative_ends:
            return 0
        return self.cumulative_ends[-1]

    def map_position(self, transcript_position_0based: int) -> Optional[int]:
        if transcript_position_0based < 0 or transcript_position_0based >= self.length:
            return None

        exon_index = bisect_right(self.cumulative_ends, transcript_position_0based)
        previous_end = 0 if exon_index == 0 else self.cumulative_ends[exon_index - 1]
        offset = transcript_position_0based - previous_end
        exon_start, exon_end = self.exons[exon_index]

        if self.strand == "+":
            return exon_start + offset

        if self.strand == "-":
            return exon_end - offset

        raise ConversionError("unsupported transcript strand: %s" % self.strand)


def log(message: str) -> None:
    print(message, file=sys.stderr, flush=True)


def open_text(path: str) -> TextIO:
    if path.endswith(".gz"):
        return gzip.open(path, "rt", newline="")
    return open(path, "r", newline="")


def delimiter_from_option(value: str) -> str:
    if value in ("\\t", "tab", "tsv"):
        return "\t"
    if value in (",", "comma", "csv"):
        return ","
    if len(value) == 1:
        return value
    raise ConversionError("unsupported separator: %s" % value)


def parse_gtf_attributes(attributes: str) -> Dict[str, str]:
    parsed: Dict[str, str] = {}

    for field in attributes.strip().rstrip(";").split(";"):
        field = field.strip()
        if not field:
            continue

        parts = field.split(None, 1)
        if len(parts) != 2:
            continue

        key, value = parts
        parsed[key] = value.strip().strip('"')

    return parsed


def load_transcripts(gtf_path: str) -> Dict[str, TranscriptModel]:
    raw: Dict[str, List[Tuple[str, str, int, int]]] = defaultdict(list)

    with open_text(gtf_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue

            fields = line.rstrip("\n\r").split("\t")
            if len(fields) < 9:
                continue

            if fields[2] != "exon":
                continue

            chrom = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]

            attrs = parse_gtf_attributes(fields[8])
            transcript_id = attrs.get("transcript_id")

            if not transcript_id:
                continue

            raw[transcript_id].append((chrom, strand, start, end))

    models: Dict[str, TranscriptModel] = {}

    for transcript_id, records in raw.items():
        chroms = {record[0] for record in records}
        strands = {record[1] for record in records}

        if len(chroms) != 1 or len(strands) != 1:
            raise ConversionError("inconsistent exon records for transcript: %s" % transcript_id)

        chrom = records[0][0]
        strand = records[0][1]

        if strand == "+":
            records.sort(key=lambda record: record[2])
        elif strand == "-":
            records.sort(key=lambda record: record[2], reverse=True)
        else:
            raise ConversionError("unsupported strand for transcript %s: %s" % (transcript_id, strand))

        exons = [(record[2], record[3]) for record in records]
        models[transcript_id] = TranscriptModel(chrom, strand, exons)

    if not models:
        raise ConversionError("no transcript models loaded from GTF: %s" % gtf_path)

    return models


def require_columns(fieldnames: Optional[Sequence[str]], columns: Iterable[str]) -> None:
    if fieldnames is None:
        raise ConversionError("input table is missing a header row")

    missing = [column for column in columns if column not in fieldnames]
    if missing:
        raise ConversionError("missing required input column(s): %s" % ", ".join(missing))


def build_output_header(
    input_header: Sequence[str],
    transcript_column: str,
    position_column: str,
) -> List[str]:
    standard = [
        "chr",
        "pos_1based",
        "start_0based",
        "end_0based_exclusive",
        "strand",
        "transcript_id",
        "transcript_position_0based",
    ]

    header = list(standard)
    consumed = {transcript_column, position_column}

    for column in input_header:
        if column in consumed:
            continue

        output_column = column

        if output_column in header:
            output_column = "input_" + output_column

        if output_column not in header:
            header.append(output_column)

    return header


def parse_transcript_position(value: str, line_number: int, column: str) -> int:
    try:
        parsed = int(value)
    except ValueError:
        raise ConversionError(
            "line %d: non-integer transcript coordinate in column %s: %s"
            % (line_number, column, value)
        )

    return parsed


def convert(args: argparse.Namespace) -> int:
    input_path = Path(args.input)
    output_path = Path(args.output)

    if not input_path.exists() or input_path.stat().st_size == 0:
        raise ConversionError("input file missing or empty: %s" % input_path)

    if output_path.exists() and not args.force:
        raise ConversionError("output exists; use --force to overwrite: %s" % output_path)

    output_path.parent.mkdir(parents=True, exist_ok=True)

    delimiter = delimiter_from_option(args.input_separator)
    models = load_transcripts(args.gtf)

    total_sites = 0
    converted_sites = 0
    missing_transcripts: Dict[str, int] = defaultdict(int)
    out_of_range: Dict[str, int] = defaultdict(int)

    output_tmp = str(output_path) + ".tmp.%d" % os.getpid()

    log("input_path\t%s" % input_path)
    log("gtf_path\t%s" % args.gtf)
    log("output_path\t%s" % output_path)
    log("transcript_column\t%s" % args.transcript_column)
    log("position_column\t%s" % args.position_column)
    log("position_base\t%s" % args.position_base)

    success = False

    try:
        with open_text(str(input_path)) as in_handle:
            reader = csv.DictReader(in_handle, delimiter=delimiter)

            require_columns(
                reader.fieldnames,
                [args.transcript_column, args.position_column],
            )

            input_header = reader.fieldnames or []
            output_header = build_output_header(
                input_header=input_header,
                transcript_column=args.transcript_column,
                position_column=args.position_column,
            )

            with open(output_tmp, "w", newline="") as out_handle:
                writer = csv.DictWriter(
                    out_handle,
                    fieldnames=output_header,
                    delimiter="\t",
                    lineterminator="\n",
                    extrasaction="ignore",
                )
                writer.writeheader()

                for line_number, row in enumerate(reader, start=2):
                    total_sites += 1

                    transcript_id = row[args.transcript_column].strip()
                    if not transcript_id:
                        raise ConversionError("line %d: transcript ID is empty" % line_number)

                    model = models.get(transcript_id)
                    if model is None:
                        missing_transcripts[transcript_id] += 1
                        continue

                    input_position = parse_transcript_position(
                        row[args.position_column],
                        line_number,
                        args.position_column,
                    )

                    if args.position_base == 0:
                        transcript_position_0based = input_position
                    else:
                        transcript_position_0based = input_position - 1

                    genomic_position_1based = model.map_position(transcript_position_0based)

                    if genomic_position_1based is None:
                        out_of_range[transcript_id] += 1
                        continue

                    out_row: Dict[str, str] = {
                        "chr": model.chrom,
                        "pos_1based": str(genomic_position_1based),
                        "start_0based": str(genomic_position_1based - 1),
                        "end_0based_exclusive": str(genomic_position_1based),
                        "strand": model.strand,
                        "transcript_id": transcript_id,
                        "transcript_position_0based": str(transcript_position_0based),
                    }

                    consumed = {args.transcript_column, args.position_column}

                    for column in input_header:
                        if column in consumed:
                            continue

                        output_column = column
                        if output_column in out_row:
                            output_column = "input_" + output_column

                        out_row[output_column] = row.get(column, "")

                    writer.writerow(out_row)
                    converted_sites += 1

                    if converted_sites % 1000000 == 0:
                        log("converted_sites\t%d" % converted_sites)

        os.replace(output_tmp, output_path)
        success = True

        log("input_sites\t%d" % total_sites)
        log("converted_sites\t%d" % converted_sites)
        log("missing_transcript_ids\t%d" % len(missing_transcripts))
        log("out_of_range_transcript_ids\t%d" % len(out_of_range))

        for transcript_id, count in sorted(missing_transcripts.items()):
            log("missing_transcript_id\t%s\t%d" % (transcript_id, count))

        for transcript_id, count in sorted(out_of_range.items()):
            log("out_of_range_transcript_id\t%s\t%d" % (transcript_id, count))

        log("output_file_path\t%s" % output_path)

        return 0

    finally:
        if not success and os.path.exists(output_tmp):
            try:
                os.remove(output_tmp)
            except OSError:
                pass


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("--input", required=True, help="Input transcript-level TSV/CSV")
    parser.add_argument("--gtf", required=True, help="GTF annotation file")
    parser.add_argument("--output", required=True, help="Output genomic-coordinate TSV")

    parser.add_argument(
        "--transcript-column",
        required=True,
        help="Column containing transcript IDs",
    )
    parser.add_argument(
        "--position-column",
        required=True,
        help="Column containing transcript coordinate",
    )
    parser.add_argument(
        "--position-base",
        required=True,
        type=int,
        choices=(0, 1),
        help="Base of input transcript coordinate: 0 or 1",
    )
    parser.add_argument(
        "--input-separator",
        default="\\t",
        help="Input separator: tab, tsv, comma, csv, or a single character. Output is always TSV.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overwrite output if it already exists",
    )

    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        return convert(args)
    except ConversionError as exc:
        print("ERROR: %s" % exc, file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
