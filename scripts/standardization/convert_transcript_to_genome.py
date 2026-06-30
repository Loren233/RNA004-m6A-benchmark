#!/usr/bin/env python3
"""Convert transcript-coordinate site tables to genomic coordinates."""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from bisect import bisect_right
from collections import defaultdict
from typing import Dict, Iterable, List, Optional, Sequence, TextIO, Tuple


Exon = Tuple[int, int]


class TranscriptModel(object):
    def __init__(self, chrom: str, strand: str, exons: List[Exon]) -> None:
        self.chrom = chrom
        self.strand = strand
        self.exons = exons
        self.cumulative_ends = []  # type: List[int]
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
            genomic_position_1based = exon_start + offset
        else:
            genomic_position_1based = exon_end - offset
        return genomic_position_1based


def open_text(path: str) -> TextIO:
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r", newline="")


def parse_gtf_attributes(attributes: str) -> Dict[str, str]:
    parsed = {}  # type: Dict[str, str]
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
    raw = defaultdict(list)  # type: Dict[str, List[Tuple[str, str, int, int]]]
    with open_text(gtf_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
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

    models = {}  # type: Dict[str, TranscriptModel]
    for transcript_id, records in raw.items():
        chroms = {record[0] for record in records}
        strands = {record[1] for record in records}
        if len(chroms) != 1 or len(strands) != 1:
            raise ValueError("Inconsistent GTF exon records for transcript: %s" % transcript_id)
        chrom = records[0][0]
        strand = records[0][1]
        if strand == "+":
            records.sort(key=lambda record: record[2])
        elif strand == "-":
            records.sort(key=lambda record: record[2], reverse=True)
        else:
            raise ValueError("Unsupported strand for transcript %s: %s" % (transcript_id, strand))
        exons = [(record[2], record[3]) for record in records]
        models[transcript_id] = TranscriptModel(chrom, strand, exons)
    return models


def delimiter_from_option(value: str) -> str:
    if value in ("\\t", "tab", "tsv"):
        return "\t"
    if value in (",", "comma", "csv"):
        return ","
    if len(value) == 1:
        return value
    raise ValueError("Unsupported input separator: %s" % value)


def require_columns(fieldnames: Optional[Sequence[str]], columns: Iterable[Optional[str]]) -> None:
    if fieldnames is None:
        raise SystemExit("ERROR: input table is missing a header row")
    missing = [column for column in columns if column and column not in fieldnames]
    if missing:
        raise SystemExit("ERROR: missing required input column(s): %s" % ", ".join(missing))


def normalized_columns(args: argparse.Namespace, row: Dict[str, str]) -> Dict[str, str]:
    normalized = {}  # type: Dict[str, str]
    tool = (args.tool_name or "").lower()

    if tool == "dorado":
        if args.modified_count_column:
            normalized["dorado_m6A_count"] = row.get(args.modified_count_column, "")
        if args.depth_column:
            normalized["dorado_depth"] = row.get(args.depth_column, "")
        if args.percent_column:
            normalized["dorado_m6A_percent"] = row.get(args.percent_column, "")
        if args.stoich_column:
            normalized["dorado_stoich"] = row.get(args.stoich_column, "")
    elif tool == "m6anet":
        if args.depth_column:
            normalized["m6anet_depth"] = row.get(args.depth_column, "")
        if "probability_modified" in row:
            normalized["m6anet_site_probability"] = row.get("probability_modified", "")
        if args.stoich_column:
            normalized["m6anet_stoich"] = row.get(args.stoich_column, "")
        if "kmer" in row:
            normalized["kmer"] = row.get("kmer", "")
    else:
        if args.modified_count_column:
            normalized["modified_count"] = row.get(args.modified_count_column, "")
        if args.depth_column:
            normalized["depth"] = row.get(args.depth_column, "")
        if args.percent_column:
            normalized["percent"] = row.get(args.percent_column, "")
        if args.stoich_column:
            normalized["stoich"] = row.get(args.stoich_column, "")

    return normalized


def make_output_header(
    input_header: Sequence[str], args: argparse.Namespace, first_normalized: Dict[str, str]
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
    for column in first_normalized.keys():
        if column not in header:
            header.append(column)

    consumed = {
        args.transcript_column,
        args.position_column,
        args.depth_column,
        args.stoich_column,
        args.modified_count_column,
        args.percent_column,
    }
    for column in input_header:
        if column in consumed:
            continue
        output_column = column
        if output_column in header:
            output_column = "input_" + output_column
        if output_column not in header:
            header.append(output_column)
    return header


def convert(args: argparse.Namespace) -> int:
    delimiter = delimiter_from_option(args.input_separator)
    models = load_transcripts(args.gtf)

    total_sites = 0
    converted_sites = 0
    missing_transcripts = defaultdict(int)  # type: Dict[str, int]
    out_of_range = defaultdict(int)  # type: Dict[str, int]

    with open_text(args.input) as in_handle:
        reader = csv.DictReader(in_handle, delimiter=delimiter)
        required = [
            args.transcript_column,
            args.position_column,
            args.depth_column,
            args.stoich_column,
            args.modified_count_column,
            args.percent_column,
        ]
        require_columns(reader.fieldnames, required)
        input_header = reader.fieldnames or []
        rows = list(reader)

    output_rows = []  # type: List[Dict[str, str]]
    output_header = None  # type: Optional[List[str]]

    for row in rows:
        total_sites += 1
        transcript_id = row[args.transcript_column]
        model = models.get(transcript_id)
        if model is None:
            missing_transcripts[transcript_id] += 1
            continue

        try:
            input_position = int(row[args.position_column])
        except ValueError:
            raise SystemExit(
                "ERROR: non-integer transcript coordinate in column %s: %s"
                % (args.position_column, row[args.position_column])
            )
        transcript_position_0based = input_position if args.position_base == 0 else input_position - 1
        genomic_position_1based = model.map_position(transcript_position_0based)
        if genomic_position_1based is None:
            out_of_range[transcript_id] += 1
            continue

        normalized = normalized_columns(args, row)
        if output_header is None:
            output_header = make_output_header(input_header, args, normalized)

        out_row = {}  # type: Dict[str, str]
        out_row["chr"] = model.chrom
        out_row["pos_1based"] = str(genomic_position_1based)
        out_row["start_0based"] = str(genomic_position_1based - 1)
        out_row["end_0based_exclusive"] = str(genomic_position_1based)
        out_row["strand"] = model.strand
        out_row["transcript_id"] = transcript_id
        out_row["transcript_position_0based"] = str(transcript_position_0based)
        out_row.update(normalized)

        consumed = {
            args.transcript_column,
            args.position_column,
            args.depth_column,
            args.stoich_column,
            args.modified_count_column,
            args.percent_column,
        }
        for column in input_header:
            if column in consumed:
                continue
            output_column = column if column not in out_row else "input_" + column
            out_row[output_column] = row.get(column, "")

        output_rows.append(out_row)
        converted_sites += 1

    if output_header is None:
        output_header = make_output_header(input_header, args, {})

    with open(args.output, "w", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, fieldnames=output_header, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for out_row in output_rows:
            writer.writerow(out_row)

    print("input_sites\t%d" % total_sites, file=sys.stderr)
    print("converted_sites\t%d" % converted_sites, file=sys.stderr)
    print("missing_transcript_ids\t%d" % len(missing_transcripts), file=sys.stderr)
    for transcript_id, count in sorted(missing_transcripts.items()):
        print("missing_transcript_id\t%s\t%d" % (transcript_id, count), file=sys.stderr)
    print("out_of_range_transcript_ids\t%d" % len(out_of_range), file=sys.stderr)
    for transcript_id, count in sorted(out_of_range.items()):
        print("out_of_range_transcript_id\t%s\t%d" % (transcript_id, count), file=sys.stderr)

    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--transcript-column", required=True)
    parser.add_argument("--position-column", required=True)
    parser.add_argument("--position-base", required=True, type=int, choices=(0, 1))
    parser.add_argument("--input-separator", default="\\t")
    parser.add_argument("--depth-column")
    parser.add_argument("--stoich-column")
    parser.add_argument("--modified-count-column")
    parser.add_argument("--percent-column")
    parser.add_argument("--tool-name")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return convert(args)


if __name__ == "__main__":
    raise SystemExit(main())
