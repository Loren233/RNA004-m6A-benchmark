#!/usr/bin/env python3
"""Focused tests for visualization scripts."""

from __future__ import annotations

import csv
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SINGLE = ROOT / "scripts" / "11_visualization" / "single" / "plot_single.py"
TOTAL = ROOT / "scripts" / "11_visualization" / "total" / "plot_total.py"


GROUPS = ["depth_gt40", "depth_31_40", "depth_21_30", "depth_11_20", "depth_1_10"]


def write_tsv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def matched_rows(groups: list[str] = GROUPS) -> list[list[object]]:
    rows = []
    for group_index, group in enumerate(groups):
        for point_index in range(6):
            x = 0.05 + 0.12 * point_index
            y = min(1.0, x + 0.02 * group_index)
            rows.append([f"chr{group_index + 1}", 100 + point_index, x, y, group])
    return rows


def summary_rows(groups: list[str] = GROUPS, n: int = 6) -> list[list[object]]:
    return [[group, n, 0.8, 0.7] for group in groups]


class VisualizationTests(unittest.TestCase):
    def run_cmd(self, cmd: list[str], tmp: Path, expect_success: bool = True) -> subprocess.CompletedProcess[str]:
        env = os.environ.copy()
        env["MPLCONFIGDIR"] = str(tmp / "mplconfig")
        result = subprocess.run(cmd, text=True, capture_output=True, env=env)
        if expect_success:
            self.assertEqual(result.returncode, 0, result.stderr)
        else:
            self.assertNotEqual(result.returncode, 0)
        return result

    def make_single_inputs(self, tmp: Path, groups: list[str] = GROUPS) -> tuple[Path, Path]:
        matched = tmp / "matched.tsv"
        summary = tmp / "summary.tsv"
        write_tsv(matched, ["chr", "pos", "truth_stoich", "prediction_stoich", "depth_group"], matched_rows(groups))
        write_tsv(summary, ["group", "n", "pearson_r2", "spearman_rho2"], summary_rows(groups))
        return matched, summary

    def test_plot_single_creates_pdf_and_png(self):
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            matched, summary = self.make_single_inputs(tmp)
            prefix = tmp / "single"
            self.run_cmd([
                sys.executable, str(SINGLE),
                "--matched-table", str(matched),
                "--summary-table", str(summary),
                "--output-prefix", str(prefix),
                "--tool-name", "Dorado",
                "--truth-name", "GLORI",
            ], tmp)
            self.assertTrue(Path(str(prefix) + ".pdf").exists())
            self.assertTrue(Path(str(prefix) + ".png").exists())

    def test_single_accepts_custom_group_list(self):
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            groups = ["depth_hi", "depth_lo"]
            matched, summary = self.make_single_inputs(tmp, groups)
            prefix = tmp / "custom_single"
            self.run_cmd([
                sys.executable, str(SINGLE),
                "--matched-table", str(matched),
                "--summary-table", str(summary),
                "--output-prefix", str(prefix),
                "--tool-name", "Dorado",
                "--truth-name", "GLORI",
                "--groups", "depth_hi", "depth_lo",
                "--group-labels", "high", "low",
            ], tmp)
            self.assertTrue(Path(str(prefix) + ".pdf").exists())
            self.assertTrue(Path(str(prefix) + ".png").exists())

    def make_tool_inputs(self, tmp: Path, name: str, empty_group: str | None = None) -> tuple[Path, Path]:
        matched = tmp / f"{name}_matched.tsv"
        summary = tmp / f"{name}_summary.tsv"
        rows = [row for row in matched_rows() if row[-1] != empty_group]
        write_tsv(matched, ["chr", "pos", "truth_stoich", "prediction_stoich", "depth_group"], rows)
        summary_data = []
        for group in GROUPS:
            if group == empty_group:
                summary_data.append([group, 0, "NA", "NA"])
            else:
                summary_data.append([group, 6, 0.8, 0.7])
        write_tsv(summary, ["group", "n", "pearson_r2", "spearman_rho2"], summary_data)
        return matched, summary

    def make_config(self, tmp: Path, tools: list[str], empty_group: str | None = None) -> Path:
        config = tmp / "tools.tsv"
        rows = []
        for tool in tools:
            matched, summary = self.make_tool_inputs(tmp, tool, empty_group=empty_group if tool == tools[0] else None)
            rows.append([tool, matched, summary, "truth_stoich", "prediction_stoich", "depth_group"])
        write_tsv(config, ["tool_name", "matched_table", "summary_table", "x_column", "y_column", "group_column"], rows)
        return config

    def test_plot_total_creates_pdf_and_png_for_two_tools(self):
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            config = self.make_config(tmp, ["ToolA", "ToolB"])
            prefix = tmp / "total"
            self.run_cmd([
                sys.executable, str(TOTAL),
                "--tool-config", str(config),
                "--output-prefix", str(prefix),
                "--truth-name", "GLORI",
            ], tmp)
            self.assertTrue(Path(str(prefix) + ".pdf").exists())
            self.assertTrue(Path(str(prefix) + ".png").exists())

    def test_total_tool_column_order_follows_config(self):
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            config = self.make_config(tmp, ["ToolB", "ToolA"])
            prefix = tmp / "ordered_total"
            result = self.run_cmd([
                sys.executable, str(TOTAL),
                "--tool-config", str(config),
                "--output-prefix", str(prefix),
                "--truth-name", "GLORI",
            ], tmp)
            self.assertIn("tool_order\tToolB,ToolA", result.stderr)

    def test_missing_requested_group_fails_clearly(self):
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            matched = tmp / "matched.tsv"
            summary = tmp / "summary.tsv"
            write_tsv(matched, ["truth_stoich", "prediction_stoich", "depth_group"], matched_rows(["depth_present"]))
            write_tsv(summary, ["group", "n", "pearson_r2", "spearman_rho2"], summary_rows(["depth_present"]))
            result = self.run_cmd([
                sys.executable, str(SINGLE),
                "--matched-table", str(matched),
                "--summary-table", str(summary),
                "--output-prefix", str(tmp / "missing"),
                "--tool-name", "Dorado",
                "--truth-name", "GLORI",
                "--groups", "depth_absent",
            ], tmp, expect_success=False)
            self.assertIn("requested group", result.stderr)

    def test_empty_group_in_total_creates_empty_panel(self):
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            config = self.make_config(tmp, ["ToolA", "ToolB"], empty_group="depth_gt40")
            prefix = tmp / "empty_total"
            self.run_cmd([
                sys.executable, str(TOTAL),
                "--tool-config", str(config),
                "--output-prefix", str(prefix),
                "--truth-name", "GLORI",
            ], tmp)
            self.assertTrue(Path(str(prefix) + ".pdf").exists())
            self.assertTrue(Path(str(prefix) + ".png").exists())

    def test_missing_required_matched_columns_fail(self):
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            matched = tmp / "matched.tsv"
            summary = tmp / "summary.tsv"
            write_tsv(matched, ["truth_stoich", "depth_group"], [[0.1, "depth_gt40"]])
            write_tsv(summary, ["group", "n", "pearson_r2", "spearman_rho2"], summary_rows())
            result = self.run_cmd([
                sys.executable, str(SINGLE),
                "--matched-table", str(matched),
                "--summary-table", str(summary),
                "--output-prefix", str(tmp / "bad_columns"),
                "--tool-name", "Dorado",
                "--truth-name", "GLORI",
            ], tmp, expect_success=False)
            self.assertIn("prediction_stoich", result.stderr)

    def test_missing_required_summary_columns_fail(self):
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            matched, _summary = self.make_single_inputs(tmp)
            bad_summary = tmp / "bad_summary.tsv"
            write_tsv(bad_summary, ["group", "n", "pearson_r2"], [[group, 6, 0.8] for group in GROUPS])
            result = self.run_cmd([
                sys.executable, str(SINGLE),
                "--matched-table", str(matched),
                "--summary-table", str(bad_summary),
                "--output-prefix", str(tmp / "bad_summary"),
                "--tool-name", "Dorado",
                "--truth-name", "GLORI",
            ], tmp, expect_success=False)
            self.assertIn("spearman_rho2", result.stderr)

    def test_existing_outputs_refused_without_force(self):
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            matched, summary = self.make_single_inputs(tmp)
            prefix = tmp / "exists"
            Path(str(prefix) + ".pdf").write_text("existing\n", encoding="utf-8")
            result = self.run_cmd([
                sys.executable, str(SINGLE),
                "--matched-table", str(matched),
                "--summary-table", str(summary),
                "--output-prefix", str(prefix),
                "--tool-name", "Dorado",
                "--truth-name", "GLORI",
            ], tmp, expect_success=False)
            self.assertIn("output exists", result.stderr)


if __name__ == "__main__":
    unittest.main()
