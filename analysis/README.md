# Downstream Statistical Analysis

This directory contains downstream statistical analyses
of RNA004-based m6A quantification.

## Analysis scope

- Absolute modification stoichiometry (beta)
- Differential modification stoichiometry (delta beta)
- Cross-project consistency
- Read-depth-dependent performance
- Orthogonal comparison against GLORI
- High-confidence filtering and sensitivity analysis
- DeepRM-specific sensitivity analysis

## Contents

Python and R scripts, Slurm job scripts, summary statistics,
and selected figures are provided.

Large intermediate files are excluded from this repository.
HPC-specific paths and input files may require configuration.

## Configuration

Run commands from the repository root.

The analysis root can be configured through the
RNA004_STAT_ROOT environment variable:

    export RNA004_STAT_ROOT="$(pwd)/analysis"

Python, R and Slurm scripts use this directory
to locate analysis inputs and outputs.

Slurm accounts must be specified according to
the local HPC configuration.

Large intermediate datasets are not included
and must be generated or supplied before running
the downstream analyses.
