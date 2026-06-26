#!/usr/bin/env python3
"""Tests for Dorado site-table standardization."""

from __future__ import annotations

import csv
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "09_standardization" / "standardize_site_tables.py"
HEADER = [
    "chr",
    "pos_1based",
    "start_0based",
    "end_0based_exclusive",
    "strand",
    "transcript_id",
    "transcript_position_0based",
    "dorado_m6A_count",
    "dorado_depth",
    "dorado_m6A_percent",
    "dorado_stoich",
]


def find_test_gnu_sort() -> str:
    candidates = []
    if os.environ.get("GNU_SORT"):
        candidates.append(os.environ["GNU_SORT"])
    for name in ("gsort", "sort"):
        path = shutil.which(name)
        if path:
            candidates.append(path)
    candidates.append(
        str(
            Path.home()
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
    for candidate in candidates:
        if not Path(candidate).exists() and os.path.sep in candidate:
            continue
        try:
            result = subprocess.run([candidate, "--version"], text=True, capture_output=True)
        except OSError:
            continue
        text = result.stdout + result.stderr
        if result.returncode == 0 and ("GNU coreutils" in text or "GNU sort" in text):
            return candidate
    raise unittest.SkipTest("GNU sort not available")


class DoradoStandardizationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.gnu_sort = find_test_gnu_sort()

    def run_standardize(self, rows, extra_args=None, expect_success=True):
        extra_args = extra_args or []
        tmpdir = tempfile.TemporaryDirectory()
        self.addCleanup(tmpdir.cleanup)
        tmp = Path(tmpdir.name)
        input_path = tmp / "input.tsv"
        output_path = tmp / "output.tsv"
        with input_path.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(HEADER)
            writer.writerows(rows)
        env = os.environ.copy()
        env["GNU_SORT"] = self.gnu_sort
        cmd = [
            sys.executable,
            str(SCRIPT),
            "dorado",
            "--input",
            str(input_path),
            "--output",
            str(output_path),
            "--tmp-dir",
            str(tmp),
        ] + extra_args
        result = subprocess.run(cmd, text=True, capture_output=True, env=env)
        if expect_success:
            self.assertEqual(result.returncode, 0, msg=result.stderr)
            with output_path.open(newline="") as handle:
                output_rows = list(csv.DictReader(handle, delimiter="\t"))
            return result, output_path, output_rows
        self.assertNotEqual(result.returncode, 0)
        return result, output_path, []

    def test_aggregates_unsorted_transcript_records(self):
        rows = [
            ["chr1", "200", "199", "200", "+", "ENST2", "5", "2", "5", "40", "0.4"],
            ["chr1", "100", "99", "100", "+", "ENST1", "1", "1", "4", "25", "0.25"],
            ["chr1", "100", "99", "100", "+", "ENST2", "2", "3", "6", "50", "0.5"],
        ]
        result, _, output = self.run_standardize(rows)
        self.assertIn("input_record_count\t3", result.stderr)
        self.assertEqual(len(output), 2)
        first = output[0]
        self.assertEqual((first["chr"], first["pos_1based"], first["strand"]), ("chr1", "100", "+"))
        self.assertEqual(first["depth"], "10")
        self.assertEqual(first["modified_count"], "4")
        self.assertEqual(first["stoich"], "0.4")
        self.assertEqual(first["transcript_count"], "2")
        self.assertEqual(first["transcript_record_count"], "2")

    def test_zero_modified_count_site(self):
        rows = [["chr1", "300", "299", "300", "-", "ENST3", "9", "0", "7", "0", "0"]]
        _, _, output = self.run_standardize(rows)
        self.assertEqual(output[0]["modified_count"], "0")
        self.assertEqual(output[0]["stoich"], "0")

    def test_duplicate_transcript_records_count_as_records(self):
        rows = [
            ["chr1", "100", "99", "100", "+", "ENST1", "1", "1", "4", "25", "0.25"],
            ["chr1", "100", "99", "100", "+", "ENST1", "1", "2", "5", "40", "0.4"],
        ]
        _, _, output = self.run_standardize(rows)
        self.assertEqual(output[0]["depth"], "9")
        self.assertEqual(output[0]["modified_count"], "3")
        self.assertEqual(output[0]["transcript_count"], "1")
        self.assertEqual(output[0]["transcript_record_count"], "2")

    def test_invalid_modified_count_gt_depth(self):
        rows = [["chr1", "100", "99", "100", "+", "ENST1", "1", "5", "4", "125", "1.25"]]
        result, output_path, _ = self.run_standardize(rows, expect_success=False)
        self.assertIn("dorado_m6A_count must be <= dorado_depth", result.stderr)
        self.assertFalse(output_path.exists())

    def test_invalid_strand(self):
        rows = [["chr1", "100", "99", "100", ".", "ENST1", "1", "1", "4", "25", "0.25"]]
        result, _, _ = self.run_standardize(rows, expect_success=False)
        self.assertIn("strand must be + or -", result.stderr)

    def test_refuses_overwrite_without_force(self):
        rows = [["chr1", "100", "99", "100", "+", "ENST1", "1", "1", "4", "25", "0.25"]]
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            input_path = tmp / "input.tsv"
            output_path = tmp / "output.tsv"
            output_path.write_text("existing\n", encoding="utf-8")
            with input_path.open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                writer.writerow(HEADER)
                writer.writerows(rows)
            env = os.environ.copy()
            env["GNU_SORT"] = self.gnu_sort
            result = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "dorado",
                    "--input",
                    str(input_path),
                    "--output",
                    str(output_path),
                    "--tmp-dir",
                    str(tmp),
                ],
                text=True,
                capture_output=True,
                env=env,
            )
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("output exists", result.stderr)
            self.assertEqual(output_path.read_text(encoding="utf-8"), "existing\n")


if __name__ == "__main__":
    unittest.main()
