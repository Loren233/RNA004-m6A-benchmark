# RNA004-m6A-benchmark

Benchmarking nanopore direct RNA sequencing tools for m6A detection using RNA004 datasets.

## Project Overview

This project will compare site-level m6A calls from RNA004 direct RNA sequencing
against orthogonal ground-truth assays:

- HEK293T nanopore datasets benchmarked against GLORI.
- HeLa nanopore datasets benchmarked against eTAM-seq.

The benchmark will evaluate multiple m6A detection approaches after converting
their outputs into a shared site-level representation and matching sites by
explicit coordinate-offset rules.

## Planned Workflow

1. Dorado basecalling from raw nanopore signal files.
2. Minimap2 splice-aware alignment of RNA reads.
3. BAM quality control.
4. DeepRM m6A detection.
5. m6Anet m6A detection.
6. Dorado modification calling.
7. SingleMod m6A detection.
8. GLORI preparation for HEK293T ground truth.
9. eTAM-seq preparation for HeLa ground truth.
10. Site-table standardization across tools and assays.
11. Coordinate-offset matching between nanopore calls and ground truth.
12. Depth-stratified correlation analysis.
13. Visualization and summary reporting.

## Repository Structure

```text
.
├── AGENTS.md
├── README.md
├── .gitignore
├── configs/
├── docs/
├── environments/
├── metadata/
├── results/
├── scripts/
└── slurm/
```

`configs/` contains example configuration files and local path templates.
Personal or cluster-specific paths belong in untracked local copies.

`metadata/` contains table templates for nanopore datasets, ground-truth
datasets, and software versions. These templates intentionally do not contain
accession numbers, completed sample metadata, or software versions yet.

`environments/` will contain reproducible software environment definitions, such
as Conda, Mamba, Apptainer, Docker, or module documentation. Environment files
are planned work.

`scripts/` will contain future analysis entry points for basecalling,
alignment, quality control, tool-specific output conversion, coordinate
matching, correlation analysis, and visualization. No runnable analysis scripts
are included yet.

`slurm/` will contain future cluster job templates for each pipeline stage. No
personal HPC paths should be committed.

`docs/` contains method and contributor documentation.

`results/` is reserved for derived benchmark outputs. Large outputs are ignored
by default and should be regenerated from metadata, configs, scripts, and raw
data stored outside the repository.

## Supported Cell Lines

Planned benchmark support:

- HEK293T, using GLORI as the ground-truth m6A reference.
- HeLa, using eTAM-seq as the ground-truth m6A reference.

Additional cell lines may be added later if matching nanopore and orthogonal
ground-truth datasets are available.

## Planned Tools

- Dorado for basecalling and modification calling.
- Minimap2 for splice-aware RNA alignment.
- DeepRM for m6A detection.
- m6Anet for m6A detection.
- SingleMod for m6A detection.
- Custom standardization, coordinate matching, depth-stratified correlation,
  and visualization code to be added under `scripts/`.

Software versions are not declared yet. Future versions should be recorded in
`metadata/software_versions.tsv` and, where possible, locked in files under
`environments/`.

## Data Availability

Raw POD5, FAST5, FASTQ, BAM, SAM, BLOW5, eventalign outputs, and large
intermediate files are not stored in this repository. Future users should place
or link local data outside version control and document sample-level paths using
untracked local configuration files derived from `configs/paths.example.sh`.

This scaffold does not invent accession numbers, data locations, completed
results, or software versions.

## Current Status

This repository currently contains project organization, templates, and
documentation only. Pipeline implementation, validation data preparation,
analysis execution, and figures are planned work.
