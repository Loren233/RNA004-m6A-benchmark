# HeLa eTAM-seq Workflow

This is a planned workflow document for benchmarking HeLa RNA004 nanopore m6A
calls against eTAM-seq.

Planned steps:

1. Register verified HeLa nanopore samples in `metadata/nanopore_datasets.tsv`.
2. Register verified eTAM-seq ground-truth tables in
   `metadata/ground_truth_datasets.tsv`.
3. Configure placeholder-free local paths in an untracked `configs/paths.sh`.
4. Run basecalling, alignment, BAM quality control, and m6A detection tools.
5. Prepare eTAM-seq sites with
   `scripts/08_ground_truth/etam_seq/prepare_etam_sites.py`.
6. Standardize site tables.
7. Match nanopore and eTAM-seq sites using explicit coordinate-offset rules.
8. Calculate depth-stratified correlations and generate figures.

No completed HeLa eTAM-seq results are included yet.
