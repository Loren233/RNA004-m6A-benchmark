# HEK293T GLORI Workflow

This is a planned workflow document for benchmarking HEK293T RNA004 nanopore
m6A calls against GLORI.

Planned steps:

1. Register verified HEK293T nanopore samples in `metadata/nanopore_datasets.tsv`.
2. Use the verified GLORI ground-truth table, represented generically as
   `metadata/ground_truth_datasets/ground_truth.tsv`.
3. Set the local environment variables required by the selected shell wrappers.
4. Run basecalling, alignment, BAM quality control, and m6A detection tools.
5. Prepare GLORI sites with `scripts/02_ground_truth/glori/prepare_glori_sites.py`.
6. Standardize site tables.
7. Match nanopore and GLORI sites using explicit coordinate-offset rules.
8. Calculate depth-stratified correlations and generate figures.

Selected completed HEK293T GLORI benchmark figures are listed in
`results/README.md`.
