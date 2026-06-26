#!/usr/bin/env python3
"""Focused tests for benchmarking scripts."""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MATCH = ROOT / "scripts" / "10_benchmarking" / "match_sites_by_offset.py"
STRATIFY = ROOT / "scripts" / "10_benchmarking" / "stratify_by_depth.py"
CORR = ROOT / "scripts" / "10_benchmarking" / "calculate_correlations.py"


PRED_HEADER = ["chr", "pos_1based", "strand", "depth", "modified_count", "stoich", "transcript_count", "transcript_record_count", "source"]
TRUTH_HEADER = ["chr", "pos", "strand", "glori_depth_mean", "glori_stoich_mean"]


def write_tsv(path, header, rows):
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def read_tsv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class BenchmarkingTests(unittest.TestCase):
    def run_match(self, tmp, pred_rows, truth_rows, extra=None, expect_success=True):
        pred = tmp / "pred.tsv"
        truth = tmp / "truth.tsv"
        prefix = tmp / "matched"
        write_tsv(pred, PRED_HEADER, pred_rows)
        write_tsv(truth, TRUTH_HEADER, truth_rows)
        cmd = [
            sys.executable, str(MATCH),
            "--prediction", str(pred),
            "--truth", str(truth),
            "--out-prefix", str(prefix),
            "--truth-position-column", "pos",
            "--truth-depth-column", "glori_depth_mean",
            "--truth-stoich-column", "glori_stoich_mean",
            "--prediction-source", "Dorado",
            "--truth-source", "GLORI",
        ] + (extra or [])
        result = subprocess.run(cmd, text=True, capture_output=True)
        if expect_success:
            self.assertEqual(result.returncode, 0, result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0)
        return result, prefix

    def test_offsets_and_mismatches_and_zero_stoich(self):
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            pred_rows = [
                ["chr1", "100", "+", "10", "0", "0", "1", "1", "Dorado"],
                ["chr1", "101", "+", "11", "3", "0.27", "1", "1", "Dorado"],
                ["chr1", "99", "+", "12", "4", "0.33", "1", "1", "Dorado"],
                ["chr1", "100", "-", "13", "5", "0.38", "1", "1", "Dorado"],
                ["chr2", "100", "+", "14", "6", "0.42", "1", "1", "Dorado"],
            ]
            truth_rows = [["chr1", "100", "+", "20", "0.5"]]
            _result, prefix = self.run_match(tmp, pred_rows, truth_rows)
            off0 = read_tsv(Path(str(prefix) + "_offset_0.tsv"))
            offp1 = read_tsv(Path(str(prefix) + "_offset_p1.tsv"))
            offm1 = read_tsv(Path(str(prefix) + "_offset_m1.tsv"))
            self.assertEqual(len(off0), 1)
            self.assertEqual(off0[0]["prediction_stoich"], "0")
            self.assertEqual(len(offp1), 1)
            self.assertEqual(offp1[0]["prediction_pos_1based"], "101")
            self.assertEqual(len(offm1), 1)
            self.assertEqual(offm1[0]["prediction_pos_1based"], "99")

    def test_duplicate_truth_and_prediction_fail(self):
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            result, _ = self.run_match(
                tmp,
                [["chr1", "100", "+", "10", "1", "0.1", "1", "1", "Dorado"]],
                [["chr1", "100", "+", "20", "0.5"], ["chr1", "100", "+", "21", "0.6"]],
                expect_success=False,
            )
            self.assertIn("duplicate truth key", result.stderr)
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            result, _ = self.run_match(
                tmp,
                [
                    ["chr1", "100", "+", "10", "1", "0.1", "1", "1", "Dorado"],
                    ["chr1", "100", "+", "11", "2", "0.2", "1", "1", "Dorado"],
                ],
                [["chr1", "100", "+", "20", "0.5"]],
                expect_success=False,
            )
            self.assertIn("duplicate prediction key", result.stderr)

    def test_malformed_depth_fails(self):
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            result, _ = self.run_match(
                tmp,
                [["chr1", "100", "+", "bad", "1", "0.1", "1", "1", "Dorado"]],
                [["chr1", "100", "+", "20", "0.5"]],
                expect_success=False,
            )
            self.assertIn("depth", result.stderr)

    def test_depth_bins_and_summary(self):
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            inp = tmp / "matched.tsv"
            out = tmp / "stratified.tsv"
            header = ["chr", "truth_pos_1based", "prediction_pos_1based", "strand", "offset", "truth_depth", "truth_stoich", "prediction_depth", "prediction_modified_count", "prediction_stoich", "transcript_count", "transcript_record_count", "prediction_source", "truth_source"]
            depths = [1, 10, 11, 20, 21, 30, 31, 40, 41]
            rows = [["chr1", str(100+d0), str(100+d0), "+", "0", "10", "0.1", str(d0), "1", "0.1", "1", "1", "Dorado", "GLORI"] for d0 in depths]
            write_tsv(inp, header, rows)
            result = subprocess.run([sys.executable, str(STRATIFY), "--input", str(inp), "--output", str(out)], text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            groups = [row["depth_group"] for row in read_tsv(out)]
            self.assertEqual(groups, ["depth_1_10", "depth_1_10", "depth_11_20", "depth_11_20", "depth_21_30", "depth_21_30", "depth_31_40", "depth_31_40", "depth_gt40"])
            summary = read_tsv(tmp / "stratified_depth_summary.tsv")
            self.assertEqual(sum(int(row["site_count"]) for row in summary), len(depths))

    def test_correlations_known_values_and_constant(self):
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            inp = tmp / "stratified.tsv"
            out = tmp / "corr.tsv"
            header = ["truth_stoich", "prediction_stoich", "depth_group", "prediction_depth"]
            rows = [["0.1", "0.1", "depth_1_10", "1"], ["0.2", "0.2", "depth_1_10", "2"], ["0.2", "0.3", "depth_11_20", "11"], ["0.4", "0.4", "depth_11_20", "12"]]
            write_tsv(inp, header, rows)
            result = subprocess.run([sys.executable, str(CORR), "--input", str(inp), "--output", str(out)], text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            data = {row["group"]: row for row in read_tsv(out)}
            self.assertAlmostEqual(float(data["all"]["pearson_r"]), 0.9233805169, places=9)
            self.assertAlmostEqual(float(data["all"]["spearman_rho"]), 0.9486832981, places=9)
            inp2 = tmp / "constant.tsv"
            out2 = tmp / "constant_corr.tsv"
            write_tsv(inp2, header, [["0.1", "0.5", "depth_1_10", "1"], ["0.2", "0.5", "depth_1_10", "2"]])
            result2 = subprocess.run([sys.executable, str(CORR), "--input", str(inp2), "--output", str(out2)], text=True, capture_output=True)
            self.assertEqual(result2.returncode, 0, result2.stderr)
            self.assertEqual(read_tsv(out2)[0]["pearson_r"], "NA")

    def test_overwrite_refusal(self):
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            inp = tmp / "matched.tsv"
            out = tmp / "stratified.tsv"
            out.write_text("existing\n", encoding="utf-8")
            write_tsv(inp, ["prediction_depth"], [["1"]])
            result = subprocess.run([sys.executable, str(STRATIFY), "--input", str(inp), "--output", str(out)], text=True, capture_output=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("output exists", result.stderr)


if __name__ == "__main__":
    unittest.main()
