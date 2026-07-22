# Pipeline Overview

The workflow was adapted from `docs/source/pipeline.docx` and synchronized with
the current script layout under `scripts/`.

## Benchmark Scope

Tools:

- Dorado
- DeepRM
- m6Anet
- SingleMod

Cell-line comparisons:

- HEK293/HEK293T: primary ground truth is GLORI.
- HeLa: primary ground truth is eTAM-seq.

## Shared Stages

1. Set the local environment variables required by the selected shell wrappers.
2. Confirm references and coordinate conventions.
3. Run basecalling/alignment or register existing compatible BAM/FASTQ files.
4. Run tool-specific m6A detection workflows.
5. Convert transcript-level outputs to genomic coordinates where needed.
6. Standardize site-level tables. Dorado standardization streams the converted
   genomic table through GNU sort and aggregates transcript-derived records by
   genomic site.
7. Match predictions to cell-line-specific ground truth.
8. Stratify by depth and calculate correlations.
9. Generate figures.

No basecalling, eventalign, model inference, Slurm job, software installation,
or large file conversion should be run during repository organization.
