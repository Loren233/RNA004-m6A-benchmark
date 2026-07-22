# HeLa eTAM-seq Workflow

This is a planned workflow document for benchmarking HeLa RNA004 nanopore m6A
calls against eTAM-seq.

Planned steps:

1. Register verified HeLa nanopore samples in `metadata/nanopore_datasets.tsv`.
2. Use the verified eTAM-seq ground-truth table, represented generically as
   `metadata/ground_truth_datasets/ground_truth.tsv`.
3. Set the local environment variables required by the selected shell wrappers.
4. Run basecalling, alignment, BAM quality control, and m6A detection tools.
5. Prepare eTAM-seq sites with
   `scripts/02_ground_truth/etam_seq/prepare_etam_sites.py`.
6. Standardize site tables.
7. Match nanopore and eTAM-seq sites using explicit coordinate-offset rules.
8. Calculate depth-stratified correlations and generate figures.

The selected completed HeLa eTAM-seq benchmark figure is listed in
`results/README.md`.
