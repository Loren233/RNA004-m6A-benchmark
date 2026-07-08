# Benchmarking scripts

This directory contains generic benchmarking scripts for comparing standardized
site-level prediction tables against standardized ground-truth tables.

The scripts are ground-truth agnostic. They can be used with GLORI, eTAM-seq,
m6A-SAC-seq, m6ACE-seq, miCLIP-derived site tables, or other site-level truth
sets, as long as the required columns are provided through command-line options.

## Scripts

### 01_match_sites_by_offset.py

Match prediction sites to ground-truth sites by genomic coordinate and strand.

Offset convention:

    prediction_pos_1based = truth_pos_1based + offset

Typical offsets:

    -1 0 1

Required prediction input columns are usually:

    chr
    pos_1based
    strand
    depth
    stoich
    modified_count
    motif

Required ground-truth input columns are usually:

    chr
    pos_1based
    strand
    truth stoichiometry column
    truth depth column

For eTAM-seq HeLa, use:

    truth stoich column: etam_stoich_mean
    truth depth column:  etam_depth_mean

The output matched tables contain generic downstream columns:

    truth_stoich
    prediction_stoich
    prediction_depth

These generic columns are used by downstream correlation, stratification, and
visualization scripts.

Example eTAM-seq usage:

    python scripts/04_benchmarking/01_match_sites_by_offset.py \
      --prediction-table /path/to/dorado_standardized_genomic_collapsed.tsv \
      --truth-table /path/to/HeLa_eTAM_FTOm_FTOp_rep1_rep2_rep3_intersection_mean.tsv \
      --output-dir /path/to/benchmark_vs_etam \
      --output-prefix PRJEB80229_HeLa_dorado_vs_eTAM \
      --prediction-name Dorado \
      --truth-name eTAM-seq \
      --prediction-stoich-column stoich \
      --prediction-depth-column depth \
      --prediction-modified-count-column modified_count \
      --truth-stoich-column etam_stoich_mean \
      --truth-depth-column etam_depth_mean \
      --truth-extra-columns motif,motif_type,source,etam_percent_mean,etam_FDR_max,etam_accessibility_mean \
      --offsets -1 0 1 \
      --offset-labels m1 0 p1 \
      --force

### 02_calculate_correlations.py

Calculate Pearson R-squared and Spearman rho-squared for one matched table.

It can calculate either:

- one overall correlation for the whole table
- grouped correlations using a column such as depth_group

Example grouped usage:

    python scripts/04_benchmarking/02_calculate_correlations.py \
      --input matched_depth_groups.tsv \
      --output matched_depth_correlation_summary.tsv \
      --x-column truth_stoich \
      --y-column prediction_stoich \
      --group-column depth_group \
      --groups depth_gt40 depth_31_40 depth_21_30 depth_11_20 depth_1_10 \
      --allow-empty-groups \
      --force

### 03_stratify_by_depth.py

Assign matched sites into non-overlapping prediction-depth groups.

Default bin edges:

    10 20 30 40

This creates:

    depth_1_10
    depth_11_20
    depth_21_30
    depth_31_40
    depth_gt40

Example usage:

    python scripts/04_benchmarking/03_stratify_by_depth.py \
      --input matched.tsv \
      --output matched_depth_groups.tsv \
      --depth-column prediction_depth \
      --group-column depth_group \
      --bin-edges 10 20 30 40 \
      --force

## Recommended workflow

    standardized prediction table
      ↓
    01_match_sites_by_offset.py
      ↓
    choose best offset
      ↓
    03_stratify_by_depth.py
      ↓
    02_calculate_correlations.py
      ↓
    plot_single.py or plot_total.py
