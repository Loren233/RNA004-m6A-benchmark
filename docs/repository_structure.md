# Repository Structure

This repository separates metadata, scripts, environment definitions,
documentation, and selected generated outputs.

## Local configuration

No committed `configs/` directory is present in the current repository. Shell
wrappers read local paths and tool settings from environment variables. Do not
commit machine-specific paths, accounts, credentials, or HPC job identifiers.

## metadata/

Place verified metadata tables here. Current templates cover nanopore datasets,
ground-truth datasets, and software versions. Leave fields blank until the
values are known and verified.

## envs/

Conda environment definitions for the documented tool routes live here.

## scripts/

Command-line scripts and shell wrappers live here. The current top-level script
groups are:

- `scripts/01_deeprm/`
- `scripts/01_dorado/`
- `scripts/01_m6anet/`
- `scripts/01_singlemod/`
- `scripts/02_ground_truth/`
- `scripts/03_standardization/`
- `scripts/04_benchmarking/`
- `scripts/05_visualization/`
- `scripts/utilities/`

## results/

Place generated outputs here during local or HPC runs. Large result files are
ignored by Git. Small final summary tables or example outputs may be committed
later only if they are intentionally curated and documented.
