# HEK293T GLORI Workflow

HEK293T RNA004 predictions from the selected tools are standardized, converted
to genomic coordinates when required, and matched to
`metadata/ground_truth_datasets/HEK293T_GLORI_rep1_rep2_intersection_mean.tsv`.
The ground-truth table is validated with
`scripts/02_ground_truth/glori/prepare_glori_sites.py`; site matching, depth
stratification, correlation calculation, and visualization use
`scripts/04_benchmarking/` and `scripts/05_visualization/`.

The table is derived from GSE210563 samples GSM6432590 and GSM6432591. Citation
and sample provenance are recorded in `docs/citations.md`. Completed HEK293T
figures are listed in `results/README.md`.
