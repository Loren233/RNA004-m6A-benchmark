# HeLa eTAM-seq Workflow

HeLa RNA004 predictions from the selected tools are standardized, converted to
genomic coordinates when required, and matched to
`metadata/ground_truth_datasets/HeLa_eTAM_FTOm_FTOp_rep1_rep2_rep3_intersection_mean.tsv`.
The ground-truth table is validated with
`scripts/02_ground_truth/etam_seq/prepare_etam_sites.py`; site matching,
depth stratification, correlation calculation, and visualization use
`scripts/04_benchmarking/` and `scripts/05_visualization/`.

The table is derived from GSE211303 HeLa poly(A) WT FTO-/FTO+ rep1-rep3
results. Citation and sample provenance are recorded in `docs/citations.md`.
The completed HeLa benchmark figure is listed in `results/README.md`.
