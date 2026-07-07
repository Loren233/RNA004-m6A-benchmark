#!/usr/bin/env python3
"""Standardize tool-specific site tables into shared genomic site tables."""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import shutil
import signal
import subprocess
import sys
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, TextIO, Tuple


DORADO_REQUIRED_COLUMNS = [
    "chr",
    "pos_1based",
    "strand",
    "transcript_id",
    "dorado_m6A_count",
    "dorado_depth",
]

DORADO_OUTPUT_COLUMNS = [
    "chr",
    "pos_1based",
    "strand",
    "depth",
    "modified_count",
    "stoich",
    "transcript_count",
    "transcript_record_count",
    "source",
]

M6ANET_REQUIRED_COLUMNS = [
    "transcript_id",
    "transcript_position",
    "n_reads",
    "probability_modified",
    "kmer",
    "mod_ratio",
]

M6ANET_OUTPUT_COLUMNS = [
    "ref",
    "transcript_id",
    "pos_1based",
    "transcript_pos_1based",
    "start_0based",
    "end",
    "site_id",
    "m6anet_stoich",
    "depth",
    "modified_count",
    "probability_modified",
    "motif",
    "source",
]


DEEPRM_REQUIRED_COLUMNS = [
    "transcript_id",
    "pos_1based",
    "depth",
    "mod_count",
    "stoich",
    "modified_base",
]

DEEPRM_OUTPUT_COLUMNS = [
    "ref",
    "transcript_id",
    "pos_1based",
    "transcript_pos_1based",
    "start_0based",
    "end",
    "site_id",
    "deeprm_stoich",
    "stoich",
    "depth",
    "modified_count",
    "probability_modified",
    "modified_base",
    "source",
]

class StandardizationError(Exception):
    pass


def log(message: str) -> None:
    print(message, file=sys.stderr, flush=True)


def now() -> str:
    return datetime.now().isoformat(timespec="seconds")


def open_text(path: str) -> TextIO:
    if path.endswith(".gz"):
        return gzip.open(path, "rt", newline="")
    return open(path, "r", newline="")


def read_header(path: str) -> List[str]:
    with open_text(path) as handle:
        header_line = handle.readline()
    if not header_line:
        raise StandardizationError("input file is empty: %s" % path)
    return header_line.rstrip("\n\r").split("\t")


def require_columns(header: Sequence[str], required: Iterable[str]) -> Dict[str, int]:
    indexes = {column: i for i, column in enumerate(header)}
    missing = [column for column in required if column not in indexes]
    if missing:
        raise StandardizationError("missing required input column(s): %s" % ", ".join(missing))
    return indexes


def parse_positive_int(value: str, line_number: int, column: str) -> int:
    try:
        parsed = int(value)
    except ValueError:
        raise StandardizationError("line %d: %s must be an integer" % (line_number, column))
    if parsed <= 0:
        raise StandardizationError("line %d: %s must be positive" % (line_number, column))
    return parsed


def parse_nonnegative_int(value: str, line_number: int, column: str) -> int:
    try:
        parsed = int(value)
    except ValueError:
        raise StandardizationError("line %d: %s must be an integer" % (line_number, column))
    if parsed < 0:
        raise StandardizationError("line %d: %s must be non-negative" % (line_number, column))
    return parsed


def validate_dorado_fields(fields: Sequence[str], indexes: Dict[str, int], line_number: int) -> None:
    chrom = fields[indexes["chr"]]
    pos = fields[indexes["pos_1based"]]
    strand = fields[indexes["strand"]]
    transcript_id = fields[indexes["transcript_id"]]
    modified = fields[indexes["dorado_m6A_count"]]
    depth = fields[indexes["dorado_depth"]]

    if not chrom:
        raise StandardizationError("line %d: chr must be non-empty" % line_number)
    parse_positive_int(pos, line_number, "pos_1based")
    if strand not in ("+", "-"):
        raise StandardizationError("line %d: strand must be + or -" % line_number)
    if not transcript_id:
        raise StandardizationError("line %d: transcript_id must be non-empty" % line_number)
    modified_count = parse_nonnegative_int(modified, line_number, "dorado_m6A_count")
    depth_count = parse_positive_int(depth, line_number, "dorado_depth")
    if modified_count > depth_count:
        raise StandardizationError(
            "line %d: dorado_m6A_count must be <= dorado_depth" % line_number
        )


def validate_dorado_input(path: str, header: Sequence[str], indexes: Dict[str, int]) -> Tuple[int, int]:
    expected_fields = len(header)
    record_count = 0
    blank_lines = 0
    with open_text(path) as handle:
        handle.readline()
        for line_number, line in enumerate(handle, start=2):
            stripped = line.rstrip("\n\r")
            if not stripped:
                blank_lines += 1
                continue
            fields = stripped.split("\t")
            if len(fields) != expected_fields:
                raise StandardizationError(
                    "line %d: expected %d fields, found %d"
                    % (line_number, expected_fields, len(fields))
                )
            validate_dorado_fields(fields, indexes, line_number)
            record_count += 1
            if record_count % 10000000 == 0:
                log("parsed_input_records\t%d" % record_count)
    return record_count, blank_lines


def emit_body(input_path: str) -> int:
    with open_text(input_path) as handle:
        handle.readline()
        for line in handle:
            if line.strip():
                sys.stdout.write(line)
    return 0


def find_gnu_sort() -> str:
    env_sort = os.environ.get("GNU_SORT")
    candidates = []  # type: List[str]
    if env_sort:
        candidates.append(env_sort)
    for name in ("gsort", "sort"):
        path = shutil.which(name)
        if path:
            candidates.append(path)

    # Common bundled Git-for-Windows coreutils location used by this workspace.
    home = Path.home()
    candidates.append(
        str(
            home
            / ".cache"
            / "codex-runtimes"
            / "codex-primary-runtime"
            / "dependencies"
            / "native"
            / "git"
            / "usr"
            / "bin"
            / "sort.exe"
        )
    )

    seen = set()
    for candidate in candidates:
        if candidate in seen:
            continue
        seen.add(candidate)
        if not Path(candidate).exists() and os.path.sep in candidate:
            continue
        try:
            result = subprocess.run(
                [candidate, "--version"], text=True, capture_output=True, check=False
            )
        except OSError:
            continue
        version_text = (result.stdout or "") + (result.stderr or "")
        if result.returncode == 0 and ("GNU coreutils" in version_text or "GNU sort" in version_text):
            return candidate
    raise StandardizationError(
        "GNU sort is required but was not found. Install coreutils sort or set GNU_SORT."
    )


def build_sort_command(
    sort_bin: str,
    tmp_dir: str,
    sort_threads: int,
    sort_memory: str,
    indexes: Dict[str, int],
) -> List[str]:
    chr_col = indexes["chr"] + 1
    pos_col = indexes["pos_1based"] + 1
    strand_col = indexes["strand"] + 1
    transcript_col = indexes["transcript_id"] + 1
    return [
        sort_bin,
        "-t",
        "\t",
        "-k%d,%dV" % (chr_col, chr_col),
        "-k%d,%dn" % (pos_col, pos_col),
        "-k%d,%d" % (strand_col, strand_col),
        "-k%d,%d" % (transcript_col, transcript_col),
        "--parallel",
        str(sort_threads),
        "-S",
        sort_memory,
        "-T",
        tmp_dir,
    ]


def default_tmp_dir(output: str, tmp_dir: Optional[str]) -> str:
    if tmp_dir:
        return tmp_dir
    output_parent = str(Path(output).resolve().parent)
    if output_parent and Path(output_parent).exists():
        return output_parent
    return tempfile.gettempdir()


def emit_site(writer: csv.writer, state: Dict[str, object], stats: Dict[str, object]) -> None:
    if not state:
        return
    chrom = str(state["chr"])
    pos = int(state["pos"])
    strand = str(state["strand"])
    depth = int(state["depth"])
    modified_count = int(state["modified_count"])
    transcripts = state["transcripts"]  # type: ignore[assignment]
    record_count = int(state["record_count"])

    transcript_count = len(transcripts)  # type: ignore[arg-type]
    if depth <= 0:
        raise StandardizationError("aggregated depth must be > 0 for %s:%s:%s" % (chrom, pos, strand))
    if modified_count < 0 or modified_count > depth:
        raise StandardizationError("aggregated modified_count out of range for %s:%s:%s" % (chrom, pos, strand))
    stoich = float(modified_count) / float(depth)
    if stoich < 0 or stoich > 1:
        raise StandardizationError("aggregated stoich out of range for %s:%s:%s" % (chrom, pos, strand))
    if transcript_count < 1:
        raise StandardizationError("transcript_count must be >= 1 for %s:%s:%s" % (chrom, pos, strand))
    if record_count < transcript_count:
        raise StandardizationError(
            "transcript_record_count must be >= transcript_count for %s:%s:%s" % (chrom, pos, strand)
        )

    last_key = stats.get("last_output_key")
    current_key = (chrom, pos, strand)
    if last_key == current_key:
        raise StandardizationError("duplicate output genomic key: %s:%s:%s" % current_key)
    stats["last_output_key"] = current_key

    writer.writerow(
        [
            chrom,
            str(pos),
            strand,
            str(depth),
            str(modified_count),
            format(stoich, ".10g"),
            str(transcript_count),
            str(record_count),
            "Dorado",
        ]
    )

    stats["output_sites"] = int(stats["output_sites"]) + 1
    stats["collapsed_records"] = int(stats["collapsed_records"]) + record_count - 1
    if transcript_count > 1:
        stats["multi_transcript_sites"] = int(stats["multi_transcript_sites"]) + 1
    stats["max_transcript_count"] = max(int(stats["max_transcript_count"]), transcript_count)
    stats["min_depth"] = depth if stats["min_depth"] is None else min(int(stats["min_depth"]), depth)
    stats["max_depth"] = depth if stats["max_depth"] is None else max(int(stats["max_depth"]), depth)
    stats["min_stoich"] = stoich if stats["min_stoich"] is None else min(float(stats["min_stoich"]), stoich)
    stats["max_stoich"] = stoich if stats["max_stoich"] is None else max(float(stats["max_stoich"]), stoich)


def aggregate_sorted_dorado(
    sorted_stream: TextIO,
    header: Sequence[str],
    indexes: Dict[str, int],
    output_tmp: str,
) -> Dict[str, object]:
    stats = {
        "input_records_from_sort": 0,
        "output_sites": 0,
        "collapsed_records": 0,
        "multi_transcript_sites": 0,
        "max_transcript_count": 0,
        "min_depth": None,
        "max_depth": None,
        "min_stoich": None,
        "max_stoich": None,
        "last_output_key": None,
    }  # type: Dict[str, object]

    state = {}  # type: Dict[str, object]
    expected_fields = len(header)
    with open(output_tmp, "w", newline="") as out_handle:
        writer = csv.writer(out_handle, delimiter="\t", lineterminator="\n")
        writer.writerow(DORADO_OUTPUT_COLUMNS)
        for sorted_line_number, line in enumerate(sorted_stream, start=1):
            stripped = line.rstrip("\n\r")
            if not stripped:
                continue
            fields = stripped.split("\t")
            if len(fields) != expected_fields:
                raise StandardizationError(
                    "sorted stream line %d: expected %d fields, found %d"
                    % (sorted_line_number, expected_fields, len(fields))
                )
            stats["input_records_from_sort"] = int(stats["input_records_from_sort"]) + 1
            chrom = fields[indexes["chr"]]
            pos = int(fields[indexes["pos_1based"]])
            strand = fields[indexes["strand"]]
            transcript_id = fields[indexes["transcript_id"]]
            modified_count = int(fields[indexes["dorado_m6A_count"]])
            depth = int(fields[indexes["dorado_depth"]])
            key = (chrom, pos, strand)

            if state and state["key"] != key:
                emit_site(writer, state, stats)
                state = {}

            if not state:
                state = {
                    "key": key,
                    "chr": chrom,
                    "pos": pos,
                    "strand": strand,
                    "depth": 0,
                    "modified_count": 0,
                    "transcripts": set(),
                    "record_count": 0,
                }

            state["depth"] = int(state["depth"]) + depth
            state["modified_count"] = int(state["modified_count"]) + modified_count
            transcripts = state["transcripts"]  # type: ignore[assignment]
            transcripts.add(transcript_id)  # type: ignore[attr-defined]
            state["record_count"] = int(state["record_count"]) + 1

        emit_site(writer, state, stats)
    return stats


def run_dorado(args: argparse.Namespace) -> int:
    input_path = str(Path(args.input))
    output_path = Path(args.output)
    if output_path.exists() and not args.force:
        raise StandardizationError("output exists; use --force to overwrite: %s" % output_path)

    tmp_dir = default_tmp_dir(args.output, args.tmp_dir)
    Path(tmp_dir).mkdir(parents=True, exist_ok=True)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_tmp = str(output_path) + ".tmp.%d" % os.getpid()

    sort_bin = find_gnu_sort()
    header = read_header(input_path)
    indexes = require_columns(header, DORADO_REQUIRED_COLUMNS)
    sort_cmd = build_sort_command(sort_bin, tmp_dir, args.sort_threads, args.sort_memory, indexes)

    log("start_time\t%s" % now())
    log("input_path\t%s" % input_path)
    log("output_path\t%s" % output_path)
    log("temporary_directory\t%s" % tmp_dir)
    log("sort_command\t%s" % " ".join(sort_cmd))

    producer = None  # type: Optional[subprocess.Popen[str]]
    sorter = None  # type: Optional[subprocess.Popen[str]]
    success = False
    try:
        input_records, blank_lines = validate_dorado_input(input_path, header, indexes)
        if blank_lines:
            log("blank_lines_ignored\t%d" % blank_lines)

        producer_cmd = [sys.executable, __file__, "_emit-body", input_path]
        env = os.environ.copy()
        env["LC_ALL"] = "C"
        producer = subprocess.Popen(producer_cmd, stdout=subprocess.PIPE, text=True)
        if producer.stdout is None:
            raise StandardizationError("failed to open producer stdout")
        sorter = subprocess.Popen(
            sort_cmd,
            stdin=producer.stdout,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env,
        )
        producer.stdout.close()
        if sorter.stdout is None:
            raise StandardizationError("failed to open sort stdout")

        stats = aggregate_sorted_dorado(sorter.stdout, header, indexes, output_tmp)
        sort_stderr = sorter.stderr.read() if sorter.stderr is not None else ""
        sort_rc = sorter.wait()
        producer_rc = producer.wait()
        if producer_rc != 0:
            raise StandardizationError("input body emitter failed with exit code %d" % producer_rc)
        if sort_rc != 0:
            raise StandardizationError(
                "GNU sort failed with exit code %d%s"
                % (sort_rc, (": " + sort_stderr.strip()) if sort_stderr.strip() else "")
            )
        if int(stats["input_records_from_sort"]) != input_records:
            raise StandardizationError(
                "sorted record count mismatch: validated %d, sorted %d"
                % (input_records, int(stats["input_records_from_sort"]))
            )

        os.replace(output_tmp, output_path)
        success = True

        log("input_record_count\t%d" % input_records)
        log("output_genomic_site_count\t%d" % int(stats["output_sites"]))
        log("collapsed_record_count\t%d" % int(stats["collapsed_records"]))
        log("multi_transcript_site_count\t%d" % int(stats["multi_transcript_sites"]))
        log("maximum_transcript_count\t%d" % int(stats["max_transcript_count"]))
        log("depth_range\t%s\t%s" % (stats["min_depth"], stats["max_depth"]))
        log("stoichiometry_range\t%s\t%s" % (stats["min_stoich"], stats["max_stoich"]))
        log("output_file_path\t%s" % output_path)
        log("finish_time\t%s" % now())
        return 0
    finally:
        if not success and os.path.exists(output_tmp):
            try:
                os.remove(output_tmp)
            except OSError:
                pass
        for proc in (sorter, producer):
            if proc is not None and proc.poll() is None:
                proc.terminate()

def parse_probability(value: str, line_number: int, column: str) -> float:
    try:
        parsed = float(value)
    except ValueError:
        raise StandardizationError("line %d: %s must be numeric" % (line_number, column))
    if parsed < 0.0 or parsed > 1.0:
        raise StandardizationError("line %d: %s must be between 0 and 1" % (line_number, column))
    return parsed


def run_m6anet(args: argparse.Namespace) -> int:
    input_path = Path(args.input)
    output_path = Path(args.output)

    if output_path.exists() and not args.force:
        raise StandardizationError("output exists; use --force to overwrite: %s" % output_path)

    if not input_path.exists() or input_path.stat().st_size == 0:
        raise StandardizationError("input file missing or empty: %s" % input_path)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_tmp = str(output_path) + ".tmp.%d" % os.getpid()

    input_records = 0
    output_records = 0
    skipped_zero_depth = 0

    log("start_time\t%s" % now())
    log("tool\tm6Anet")
    log("input_path\t%s" % input_path)
    log("output_path\t%s" % output_path)
    log("position_base\t%s" % args.position_base)

    success = False

    try:
        with open_text(str(input_path)) as in_handle:
            reader = csv.DictReader(in_handle)

            if reader.fieldnames is None:
                raise StandardizationError("input has no header: %s" % input_path)

            missing = [column for column in M6ANET_REQUIRED_COLUMNS if column not in reader.fieldnames]
            if missing:
                raise StandardizationError(
                    "missing required input column(s): %s" % ", ".join(missing)
                )

            with open(output_tmp, "w", newline="") as out_handle:
                writer = csv.DictWriter(
                    out_handle,
                    fieldnames=M6ANET_OUTPUT_COLUMNS,
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()

                for line_number, row in enumerate(reader, start=2):
                    input_records += 1

                    transcript_id = row["transcript_id"].strip()
                    if not transcript_id:
                        raise StandardizationError(
                            "line %d: transcript_id must be non-empty" % line_number
                        )

                    transcript_position = parse_nonnegative_int(
                        row["transcript_position"], line_number, "transcript_position"
                    )
                    depth = parse_positive_int(row["n_reads"], line_number, "n_reads")
                    probability_modified = parse_probability(
                        row["probability_modified"], line_number, "probability_modified"
                    )
                    stoich = parse_probability(row["mod_ratio"], line_number, "mod_ratio")
                    motif = row["kmer"].strip()

                    if depth <= 0:
                        skipped_zero_depth += 1
                        continue

                    if args.position_base == "0":
                        pos_1based = transcript_position + 1
                    else:
                        pos_1based = transcript_position

                    if pos_1based <= 0:
                        raise StandardizationError(
                            "line %d: converted pos_1based must be positive" % line_number
                        )

                    start_0based = pos_1based - 1
                    end = pos_1based

                    modified_count = int(round(float(depth) * stoich))
                    if modified_count < 0:
                        modified_count = 0
                    if modified_count > depth:
                        modified_count = depth

                    site_id = "%s|%s|%s" % (transcript_id, pos_1based, motif)

                    writer.writerow(
                        {
                            "ref": transcript_id,
                            "transcript_id": transcript_id,
                            "pos_1based": str(pos_1based),
                            "transcript_pos_1based": str(pos_1based),
                            "start_0based": str(start_0based),
                            "end": str(end),
                            "site_id": site_id,
                            "m6anet_stoich": format(stoich, ".10g"),
                            "depth": str(depth),
                            "modified_count": str(modified_count),
                            "probability_modified": format(probability_modified, ".10g"),
                            "motif": motif,
                            "source": "m6Anet",
                        }
                    )

                    output_records += 1

                    if output_records % 1000000 == 0:
                        log("written_records\t%d" % output_records)

        os.replace(output_tmp, output_path)
        success = True

        log("input_record_count\t%d" % input_records)
        log("output_record_count\t%d" % output_records)
        log("skipped_zero_depth\t%d" % skipped_zero_depth)
        log("output_file_path\t%s" % output_path)
        log("finish_time\t%s" % now())

        return 0

    finally:
        if not success and os.path.exists(output_tmp):
            try:
                os.remove(output_tmp)
            except OSError:
                pass


def parse_nonnegative_count_like(value: str, line_number: int, column: str) -> int:
    try:
        parsed = float(value)
    except ValueError:
        raise StandardizationError("line %d: %s must be numeric" % (line_number, column))
    if parsed < 0:
        raise StandardizationError("line %d: %s must be non-negative" % (line_number, column))
    return int(round(parsed))


def run_deeprm(args: argparse.Namespace) -> int:
    input_path = Path(args.input)
    output_path = Path(args.output)

    if output_path.exists() and not args.force:
        raise StandardizationError("output exists; use --force to overwrite: %s" % output_path)

    if not input_path.exists() or input_path.stat().st_size == 0:
        raise StandardizationError("input file missing or empty: %s" % input_path)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_tmp = str(output_path) + ".tmp.%d" % os.getpid()

    input_records = 0
    output_records = 0
    skipped_zero_depth = 0

    log("start_time\t%s" % now())
    log("tool\tDeepRM")
    log("input_path\t%s" % input_path)
    log("output_path\t%s" % output_path)

    success = False

    try:
        with open_text(str(input_path)) as in_handle:
            reader = csv.DictReader(in_handle, delimiter="\t")

            if reader.fieldnames is None:
                raise StandardizationError("input has no header: %s" % input_path)

            missing = [column for column in DEEPRM_REQUIRED_COLUMNS if column not in reader.fieldnames]
            if missing:
                raise StandardizationError(
                    "missing required input column(s): %s" % ", ".join(missing)
                )

            with open(output_tmp, "w", newline="") as out_handle:
                writer = csv.DictWriter(
                    out_handle,
                    fieldnames=DEEPRM_OUTPUT_COLUMNS,
                    delimiter="\t",
                    lineterminator="\n",
                )
                writer.writeheader()

                for line_number, row in enumerate(reader, start=2):
                    input_records += 1

                    transcript_id = row["transcript_id"].strip()
                    if not transcript_id:
                        raise StandardizationError(
                            "line %d: transcript_id must be non-empty" % line_number
                        )

                    pos_1based = parse_positive_int(
                        row["pos_1based"], line_number, "pos_1based"
                    )
                    depth = parse_positive_int(row["depth"], line_number, "depth")
                    modified_count = parse_nonnegative_count_like(
                        row["mod_count"], line_number, "mod_count"
                    )
                    stoich = parse_probability(row["stoich"], line_number, "stoich")
                    modified_base = row["modified_base"].strip()

                    if depth <= 0:
                        skipped_zero_depth += 1
                        continue

                    if modified_count > depth:
                        raise StandardizationError(
                            "line %d: mod_count must be <= depth" % line_number
                        )

                    start_0based = pos_1based - 1
                    end = pos_1based
                    site_id = "%s|%s|%s" % (transcript_id, pos_1based, modified_base)

                    writer.writerow(
                        {
                            "ref": transcript_id,
                            "transcript_id": transcript_id,
                            "pos_1based": str(pos_1based),
                            "transcript_pos_1based": str(pos_1based),
                            "start_0based": str(start_0based),
                            "end": str(end),
                            "site_id": site_id,
                            "deeprm_stoich": format(stoich, ".10g"),
                            "stoich": format(stoich, ".10g"),
                            "depth": str(depth),
                            "modified_count": str(modified_count),
                            "probability_modified": format(stoich, ".10g"),
                            "modified_base": modified_base,
                            "source": "DeepRM",
                        }
                    )

                    output_records += 1

                    if output_records % 10000000 == 0:
                        log("written_records\t%d" % output_records)

        os.replace(output_tmp, output_path)
        success = True

        log("input_record_count\t%d" % input_records)
        log("output_record_count\t%d" % output_records)
        log("skipped_zero_depth\t%d" % skipped_zero_depth)
        log("output_file_path\t%s" % output_path)
        log("finish_time\t%s" % now())

        return 0

    finally:
        if not success and os.path.exists(output_tmp):
            try:
                os.remove(output_tmp)
            except OSError:
                pass

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command")

    dorado = subparsers.add_parser("dorado", help="standardize Dorado genomic site table")
    dorado.add_argument("--input", required=True)
    dorado.add_argument("--output", required=True)
    dorado.add_argument("--tmp-dir")
    dorado.add_argument("--sort-threads", type=int, default=1)
    dorado.add_argument("--sort-memory", default="2G")
    dorado.add_argument("--force", action="store_true")
    dorado.set_defaults(func=run_dorado)
	
    m6anet = subparsers.add_parser("m6anet", help="standardize m6Anet transcript-level site table")
    m6anet.add_argument("--input", required=True)
    m6anet.add_argument("--output", required=True)
    m6anet.add_argument(
        "--position-base",
        choices=["0", "1"],
        default="0",
        help="base of m6Anet transcript_position; default assumes 0-based input",
    )
    m6anet.add_argument("--force", action="store_true")
    m6anet.set_defaults(func=run_m6anet)


    deeprm = subparsers.add_parser("deeprm", help="standardize DeepRM transcript-level site table")
    deeprm.add_argument("--input", required=True)
    deeprm.add_argument("--output", required=True)
    deeprm.add_argument("--force", action="store_true")
    deeprm.set_defaults(func=run_deeprm)

    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    if argv is None:
        argv = sys.argv[1:]
    if argv and argv[0] == "_emit-body":
        if len(argv) != 2:
            print("ERROR: internal _emit-body expects exactly one input path", file=sys.stderr)
            return 2
        return emit_body(argv[1])

    parser = build_parser()
    args = parser.parse_args(argv)
    if not hasattr(args, "func"):
        parser.print_help(sys.stderr)
        return 2
    if getattr(args, "sort_threads", 1) < 1:
        raise StandardizationError("--sort-threads must be >= 1")
    try:
        return args.func(args)
    except KeyboardInterrupt:
        raise
    except StandardizationError as exc:
        print("ERROR: %s" % exc, file=sys.stderr)
        return 1


if __name__ == "__main__":
    signal.signal(signal.SIGINT, signal.default_int_handler)
    raise SystemExit(main())
