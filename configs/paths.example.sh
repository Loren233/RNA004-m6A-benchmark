#!/usr/bin/env bash

# Example local path configuration for RNA004-m6A-benchmark.
# Copy this file to configs/paths.sh and edit the copy.
# Do not commit personal, cluster-specific, or credential-bearing paths.

export PROJECT_ROOT="/path/to/RNA004-m6A-benchmark"
export RAW_DATA_ROOT="/path/to/raw_nanopore_data"
export GROUND_TRUTH_ROOT="/path/to/ground_truth_tables"
export REFERENCE_ROOT="/path/to/reference_files"
export RESULTS_ROOT="/path/to/output_results"
export SCRATCH_ROOT="/path/to/scratch_space"

# Planned reference inputs.
export GENOME_FASTA="/path/to/reference/genome.fa"
export TRANSCRIPTOME_FASTA="/path/to/reference/transcriptome.fa"
export GENE_ANNOTATION_GTF="/path/to/reference/annotation.gtf"

# Planned executable locations when not available on PATH.
export DORADO_BIN="/path/to/dorado"
export MINIMAP2_BIN="/path/to/minimap2"
export SAMTOOLS_BIN="/path/to/samtools"
