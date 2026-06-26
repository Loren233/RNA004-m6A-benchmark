# Visualization Layer

This directory contains the planned visualization layer for the RNA004 m6A
benchmark. It intentionally has only two figure types:

- `single`: one prediction tool versus one ground-truth dataset, with
  non-overlapping depth ranges stacked vertically.
- `total`: multiple prediction tools versus one ground-truth dataset, with
  non-overlapping depth ranges as rows and tools as columns.

There is no separate depth-performance line plot in this visualization design.
The cumulative `depth_ge*` groups overlap and are not part of the main single or
total figures.

The HPC server cannot pull directly from the GitHub repository, so upload the
standalone scripts manually before execution.

## Single-Tool Figure

`single/plot_single.py` creates one vertically stacked density figure for a
single prediction tool compared with one ground-truth dataset. Example:

```bash
MATCHED=/path/to/results/tool_vs_truth_offset_0_stratified.tsv
SUMMARY=/path/to/results/tool_vs_truth_correlation_summary.tsv
OUT_PREFIX=/path/to/results/figures/tool_vs_truth_single

python scripts/11_visualization/single/plot_single.py \
    --matched-table "$MATCHED" \
    --summary-table "$SUMMARY" \
    --output-prefix "$OUT_PREFIX" \
    --tool-name prediction_tool \
    --truth-name ground_truth
```

The default top-to-bottom depth groups are:

```text
depth_gt40   depth > 40
depth_31_40  30 < depth <= 40
depth_21_30  20 < depth <= 30
depth_11_20  10 < depth <= 20
depth_1_10   1 <= depth <= 10
```

Use `--groups` and `--group-labels` to supply another non-overlapping depth
scheme, such as a bottom row of `depth_5_10`.

## Total Figure

`total/plot_total.py` creates the final multi-tool grid only after the
standardized and benchmarked outputs for all desired tools are available.
Columns follow the row order in a TSV tool config:

```text
tool_name	matched_table	summary_table	x_column	y_column	group_column
DeepRM	/path/deeprm.tsv	/path/deeprm_summary.tsv	truth_stoich	prediction_stoich	depth_group
Dorado	/path/dorado.tsv	/path/dorado_summary.tsv	truth_stoich	prediction_stoich	depth_group
```

Example command:

```bash
python scripts/11_visualization/total/plot_total.py \
    --tool-config /path/to/tool_config.tsv \
    --output-prefix /path/to/figures/all_tools_vs_GLORI_total \
    --truth-name GLORI
```

Both scripts write `OUTPUT_PREFIX.pdf` and `OUTPUT_PREFIX.png`, refuse to
overwrite existing outputs unless `--force` is supplied, and use one shared
log-density scale per complete figure.
