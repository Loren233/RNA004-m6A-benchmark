# Benchmarking

The benchmarking scripts are standalone Python programs intended to run on an
HPC server even when the server cannot access this repository directly. Copy the
three scripts in this directory to the server if needed.

Workflow:

```text
standardized prediction table + genomic ground truth
  -> match_sites_by_offset.py
  -> choose the globally supported offset
  -> stratify_by_depth.py
  -> calculate_correlations.py
```

The formal Dorado versus GLORI benchmark uses offset 0 when it has the clearly
largest overlap. Do not merge the -1, 0, and +1 matched tables, and do not
choose offsets separately per site.

## Example 1: Match Offsets

These paths are examples and are not hardcoded.

```bash
python scripts/10_benchmarking/match_sites_by_offset.py \
    --prediction /gpfs/work/bio/zhikefu24/surf/PRJEB80229/dorado/PRJEB80229_dorado_m6A.primary_standardized_genomic.tsv \
    --truth /gpfs/work/bio/zhikefu24/deeprm/ground_truth/glori/HEK293T_GLORI_rep1_rep2_intersection_mean.tsv \
    --out-prefix /gpfs/work/bio/zhikefu24/surf/PRJEB80229/dorado/PRJEB80229_dorado_vs_GLORI \
    --offsets -1 0 1 \
    --prediction-position-column pos_1based \
    --truth-position-column pos \
    --prediction-depth-column depth \
    --prediction-stoich-column stoich \
    --truth-depth-column glori_depth_mean \
    --truth-stoich-column glori_stoich_mean \
    --prediction-source Dorado \
    --truth-source GLORI
```

Offset convention:

```text
prediction_pos = truth_pos + offset
```

## Example 2: Stratify Offset Zero

```bash
python scripts/10_benchmarking/stratify_by_depth.py \
    --input /gpfs/work/bio/zhikefu24/surf/PRJEB80229/dorado/PRJEB80229_dorado_vs_GLORI_offset_0.tsv \
    --output /gpfs/work/bio/zhikefu24/surf/PRJEB80229/dorado/PRJEB80229_dorado_vs_GLORI_offset_0_stratified.tsv \
    --depth-column prediction_depth \
    --bin-edges 10 20 30 40
```

Default groups are `depth_1_10`, `depth_11_20`, `depth_21_30`,
`depth_31_40`, and `depth_gt40`.

## Example 3: Calculate Correlations

```bash
python scripts/10_benchmarking/calculate_correlations.py \
    --input /gpfs/work/bio/zhikefu24/surf/PRJEB80229/dorado/PRJEB80229_dorado_vs_GLORI_offset_0_stratified.tsv \
    --output /gpfs/work/bio/zhikefu24/surf/PRJEB80229/dorado/PRJEB80229_dorado_vs_GLORI_correlation_summary.tsv \
    --truth-column truth_stoich \
    --prediction-column prediction_stoich \
    --group-column depth_group
```

Correlations use `x = truth_stoich` and `y = prediction_stoich`. Pearson
`r2` and Spearman `rho2` are reported separately.
