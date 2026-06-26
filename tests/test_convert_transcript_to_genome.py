#!/usr/bin/env python3
"""Lightweight synthetic checks for transcript-to-genome conversion."""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CONVERTER = ROOT / "scripts" / "09_standardization" / "convert_transcript_to_genome.py"
FIXTURE_DIR = ROOT / "tests" / "fixtures" / "transcript_to_genome"


def main() -> int:
    with tempfile.TemporaryDirectory() as tmpdir:
        output = Path(tmpdir) / "converted.tsv"
        cmd = [
            sys.executable,
            str(CONVERTER),
            "--input",
            str(FIXTURE_DIR / "sites.tsv"),
            "--gtf",
            str(FIXTURE_DIR / "synthetic.gtf"),
            "--output",
            str(output),
            "--transcript-column",
            "transcript_id",
            "--position-column",
            "transcript_position_0based",
            "--position-base",
            "0",
            "--depth-column",
            "depth",
            "--stoich-column",
            "stoich",
            "--tool-name",
            "synthetic",
        ]
        result = subprocess.run(cmd, check=True, text=True, capture_output=True)
        stderr = result.stderr
        assert "input_sites\t6" in stderr
        assert "converted_sites\t4" in stderr
        assert "missing_transcript_id\tTXMISSING.1\t1" in stderr
        assert "out_of_range_transcript_id\tTXPLUS.1\t1" in stderr

        with output.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))

    by_label = {row["label"]: row for row in rows}
    assert by_label["plus_first_base"]["pos_1based"] == "100"
    assert by_label["plus_exon_junction"]["pos_1based"] == "200"
    assert by_label["minus_first_base"]["pos_1based"] == "209"
    assert by_label["minus_exon_junction"]["pos_1based"] == "109"
    print("synthetic transcript-to-genome validation ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
