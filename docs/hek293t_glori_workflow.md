# HEK293T GLORI Workflow

This is a planned workflow document for benchmarking HEK293T RNA004 nanopore
m6A calls against GLORI.

Planned steps:

1. Register verified HEK293T nanopore samples in `metadata/nanopore_datasets.tsv`.
2. Use the verified GLORI ground-truth table in
   `metadata/ground_truth_datasets/HEK293T_GLORI_rep1_rep2_intersection_mean.tsv`.
3. Configure placeholder-free local paths in an untracked `configs/paths.sh`.
4. Run basecalling, alignment, BAM quality control, and m6A detection tools.
5. Prepare GLORI sites with `scripts/08_ground_truth/glori/prepare_glori_sites.py`.
6. Standardize site tables.
7. Match nanopore and GLORI sites using explicit coordinate-offset rules.
8. Calculate depth-stratified correlations and generate figures.

No completed HEK293T GLORI results are included yet.
