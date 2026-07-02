# 04_benchmarking

This folder contains scripts for benchmarking standardized RNA modification prediction results against ground truth datasets such as GLORI or eTAM.

The input files should already be standardized into genomic site-level coordinates before using these scripts.

The recommended order is:

1. `01_match_sites_by_offset.py`
2. `02_calculate_correlations.py`
3. `03_stratify_by_depth.py`

## Required input format

Prediction files should use genomic coordinates and include columns such as:

- `chr`
- `pos_1based`
- `strand`
- `stoich`
- `depth`
- `modified_count`
- `motif`

Ground truth files should also use genomic coordinates.

Do not use transcript-level coordinates directly for benchmarking against GLORI or eTAM. If the prediction file contains transcript IDs such as `ENST00000000233.10`, first run the scripts in `03_standardization`.

## 01_match_sites_by_offset.py

This script matches prediction sites to ground truth sites by genomic position.

It is used to find overlapping sites between a prediction table and a ground truth table.

This step can also test small coordinate offsets, such as:

- direct match
- prediction position plus 1
- prediction position minus 1

This is useful because different tools or references may report the modified base position with a one-base shift.

The output should be a joined table containing both prediction values and ground truth values for overlapping sites.

## 02_calculate_correlations.py

This script calculates overall correlation between predicted modification stoichiometry and ground truth stoichiometry.

Typical metrics include:

- number of matched sites
- Pearson correlation
- Pearson R squared
- Spearman correlation
- Spearman rho squared

This script should be run after site matching.

## 03_stratify_by_depth.py

This script calculates benchmarking statistics within sequencing-depth bins.

This is useful because low-depth sites usually have noisier stoichiometry estimates.

Example depth bins include:

- `depth > 40`
- `30 < depth <= 40`
- `20 < depth <= 30`
- `10 < depth <= 20`
- `5 < depth <= 10`

This step helps evaluate whether a tool performs better at higher-confidence, higher-depth sites.

## Recommended workflow

For each tool:

standardized genomic prediction  
→ `01_match_sites_by_offset.py`  
→ `02_calculate_correlations.py`  
→ `03_stratify_by_depth.py`  
→ `05_visualization`

## Notes

Benchmarking should be performed only after the prediction files and ground truth files are in compatible genomic coordinate systems.

If a prediction file was generated from transcript-level coordinates, first run:

1. `03_standardization/02_convert_transcript_to_genome.py`
2. `03_standardization/03_collapse_genomic_sites.py`, if duplicated genomic sites are present

The collapsed genomic result should be preferred for benchmarking when multiple transcript isoforms map to the same genomic site.
