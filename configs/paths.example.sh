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

# Planned sample inputs.
export SAMPLE_ID="sample_id"
export CELL_LINE="HEK293T_or_HeLa"
export POD5="/path/to/input.pod5_or_pod5_directory"
export BAM="/path/to/aligned.sorted.bam"
export FASTQ="/path/to/reads.fastq"
export BLOW5="/path/to/reads.blow5"
export GROUND_TRUTH="/path/to/ground_truth.tsv"
export OUTDIR="/path/to/output_directory"
export THREADS="8"

# Planned reference inputs.
export GENOME_FASTA="/path/to/reference/genome.fa"
export TRANSCRIPTOME_FASTA="/path/to/reference/transcriptome.fa"
export TRANSCRIPT_FASTA="/path/to/reference/transcriptome.fa"
export GENE_ANNOTATION_GTF="/path/to/reference/annotation.gtf"
export GTF="/path/to/reference/annotation.gtf"
export KMER_MODEL="/path/to/rna004.nucleotide.5mer.model"

# Planned executable locations when not available on PATH.
export DORADO="/path/to/dorado"
export MINIMAP2="/path/to/minimap2"
export SAMTOOLS="/path/to/samtools"
export MODKIT="/path/to/modkit"
export F5C="/path/to/f5c"
export BLUECRAB="/path/to/blue-crab"
export M6ANET="/path/to/m6anet"
export DEEPRM="/path/to/deeprm"
export SINGLEMOD_REPO="/path/to/SingleMod"
export BEDTOOLS="/path/to/bedtools"
export PYTHON="python"

# Planned model inputs.
export DORADO_MODEL="/path/to/dorado/rna004/model"
export DEEPRM_MODEL_FILE="/path/to/deeprm_weights.pt"
