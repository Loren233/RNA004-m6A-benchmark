# Benchmarking

The benchmarking scripts are standalone Python programs. They may be copied to a
compute server and run independently of the repository checkout.

Workflow:

```text
standardized prediction table + genomic ground truth
  -> match_sites_by_offset.py
  -> choose one globally supported offset
  -> stratify_by_depth.py
  -> calculate_correlations.py
```

Do not merge matched tables from different offsets, and do not choose a
different offset independently for each site.

## Example 1: Match Offsets

```bash
PREDICTION=/path/to/prediction_standardized.tsv
TRUTH=/path/to/ground_truth.tsv
OUT_PREFIX=/path/to/results/tool_vs_truth

python scripts/10_benchmarking/match_sites_by_offset.py \
    --prediction "$PREDICTION" \
    --truth "$TRUTH" \
    --out-prefix "$OUT_PREFIX" \
    --offsets -1 0 1 \
    --prediction-position-column pos_1based \
    --truth-position-column pos \
    --prediction-depth-column depth \
    --prediction-stoich-column stoich \
    --truth-depth-column truth_depth \
    --truth-stoich-column truth_stoich \
    --prediction-source prediction_tool \
    --truth-source ground_truth
```

Offset convention:

```text
prediction_pos = truth_pos + offset
```

The truth column arguments must be adjusted to the actual ground-truth schema.

## Example 2: Stratify the Selected Offset

```bash
MATCHED=/path/to/results/tool_vs_truth_offset_0.tsv
STRATIFIED=/path/to/results/tool_vs_truth_offset_0_stratified.tsv

python scripts/10_benchmarking/stratify_by_depth.py \
    --input "$MATCHED" \
    --output "$STRATIFIED" \
    --depth-column prediction_depth \
    --bin-edges 10 20 30 40
```

Default groups:

```text
depth_1_10
depth_11_20
depth_21_30
depth_31_40
depth_gt40
```

## Example 3: Calculate Correlations

```bash
STRATIFIED=/path/to/results/tool_vs_truth_offset_0_stratified.tsv
SUMMARY=/path/to/results/tool_vs_truth_correlation_summary.tsv

python scripts/10_benchmarking/calculate_correlations.py \
    --input "$STRATIFIED" \
    --output "$SUMMARY" \
    --truth-column truth_stoich \
    --prediction-column prediction_stoich \
    --group-column depth_group
```

Correlations use ground-truth stoichiometry as x and prediction stoichiometry
as y. Pearson r-squared and Spearman rho-squared must remain separate metrics.
