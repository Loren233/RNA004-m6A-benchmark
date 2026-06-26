#!/usr/bin/env python3
"""Repository privacy and portability checks."""

from __future__ import annotations

import re
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SKIP_DIRS = {
    ".git",
    "__pycache__",
    ".pytest_cache",
    "archive",
    "results",
}
SKIP_SUFFIXES = {
    ".bam",
    ".bai",
    ".sam",
    ".fastq",
    ".fq",
    ".pod5",
    ".fast5",
    ".blow5",
    ".slow5",
    ".docx",
    ".pdf",
    ".png",
    ".jpg",
    ".jpeg",
    ".gif",
    ".zip",
    ".gz",
    ".pt",
    ".npy",
    ".npz",
}
LARGE_TEXT_LIMIT = 1_000_000


def forbidden_literals() -> list[str]:
    return [
        "/" + "gpfs" + "/",
        "zhike" + "fu24",
        "PRJEB" + "80229_dorado",
        "HEK293T_" + "GLORI_rep1_rep2_intersection_mean.tsv",
        "HeLa_" + "eTAM_FTOm_FTOp_rep1_rep2_rep3_intersection_mean.tsv",
    ]


def forbidden_regexes() -> list[re.Pattern[str]]:
    unix_home = "/" + "home" + "/"
    mac_home = "/" + "Users" + "/"
    return [
        re.compile(re.escape(unix_home) + r"[^/\s]+/"),
        re.compile(re.escape(mac_home) + r"[^/\s]+/"),
        re.compile(r"C:" + re.escape("\\") + r"Users" + re.escape("\\") + r"[^" + re.escape("\\") + r"\s]+" + re.escape("\\")),
    ]


def should_skip(path: Path) -> bool:
    relative = path.relative_to(ROOT)
    if any(part in SKIP_DIRS for part in relative.parts):
        return True
    if path.suffix.lower() in SKIP_SUFFIXES:
        return True
    if path.suffix.lower() in {".tsv", ".csv", ".bed"} and path.stat().st_size > LARGE_TEXT_LIMIT:
        return True
    return False


def read_text(path: Path) -> str | None:
    data = path.read_bytes()
    if b"\x00" in data:
        return None
    try:
        return data.decode("utf-8")
    except UnicodeDecodeError:
        try:
            return data.decode("utf-8-sig")
        except UnicodeDecodeError:
            return None


class RepositoryHygieneTests(unittest.TestCase):
    def test_no_private_or_machine_specific_paths_in_text(self):
        failures: list[str] = []
        literals = forbidden_literals()
        regexes = forbidden_regexes()
        for path in ROOT.rglob("*"):
            if not path.is_file() or should_skip(path):
                continue
            text = read_text(path)
            if text is None:
                continue
            for line_number, line in enumerate(text.splitlines(), start=1):
                for literal in literals:
                    if literal in line:
                        failures.append(f"{path.relative_to(ROOT)}:{line_number}: contains forbidden literal")
                for pattern in regexes:
                    if pattern.search(line):
                        failures.append(f"{path.relative_to(ROOT)}:{line_number}: contains personal absolute path")
        if failures:
            self.fail("Forbidden private or machine-specific repository text:\n" + "\n".join(failures))


if __name__ == "__main__":
    unittest.main()
