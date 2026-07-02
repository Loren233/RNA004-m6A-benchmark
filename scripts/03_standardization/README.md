# 03_standardization

This folder contains scripts for converting tool-specific RNA modification outputs into a common genomic site-level format before benchmarking.

The recommended order is:

1. `01_standardize_site_tables.py`
2. `02_convert_transcript_to_genome.py`
3. `03_collapse_genomic_sites.py`

## 01_standardize_site_tables.py

This script converts raw site-level outputs from different tools into a shared column format.

Use this step when a tool output has tool-specific column names or formats.

Typical standardized columns include:

- `chr` or `ref`
- `pos_1based`
- `start_0based`
- `end`
- `strand`
- `modified_count`
- `depth`
- `stoich`
- `motif`
- `tool`

This step can be used for outputs from Dorado, DeepRM, m6Anet, SingleMod, or other m6A prediction tools.

## 02_convert_transcript_to_genome.py

This script converts transcript-level coordinates to genomic coordinates using a GTF annotation.

Use this step when the standardized output is still based on transcript IDs, for example:

- transcript ID: `ENST00000000233.10`
- transcript position: `122`

After conversion, the output should contain genomic coordinates, for example:

- chromosome: `chr7`
- genomic position: `127588532`

This step is needed before comparing transcript-level predictions with genomic ground truth datasets such as GLORI or eTAM.

This step is not needed if the tool output is already in genomic coordinates.

## 03_collapse_genomic_sites.py

This script collapses duplicated genomic sites after coordinate conversion.

This is needed because multiple transcript isoforms can map to the same genomic site. If duplicate genomic sites are not collapsed, the same genomic position may be counted multiple times during benchmarking.

Rows are grouped by:

- `chr`
- `pos_1based`
- `strand`

The collapse rule is:

- `depth`: sum
- `modified_count`: sum
- `stoich`: `modified_count / depth`
- `motif`: unique motifs joined by comma
- `source_rows`: number of input rows collapsed into the genomic site

Use this step after `02_convert_transcript_to_genome.py` if duplicated genomic sites are present.

If the standardized genomic output has no duplicated genomic sites, this step is optional.

## Recommended workflow

For transcript-level tools:

raw tool output  
→ `01_standardize_site_tables.py`  
→ `02_convert_transcript_to_genome.py`  
→ `03_collapse_genomic_sites.py`  
→ `04_benchmarking`

For tools already in genomic coordinates:

raw tool output  
→ `01_standardize_site_tables.py`  
→ `03_collapse_genomic_sites.py`, only if duplicated genomic sites exist  
→ `04_benchmarking`

## Notes

Do not benchmark transcript-level coordinates directly against genomic ground truth datasets.

For example, a SingleMod site reported as `ENST00000000233.10:122` must first be converted to a genomic coordinate before it can be compared with GLORI or eTAM.

If multiple standardized rows map to the same genomic site, use `03_collapse_genomic_sites.py` before overlap and correlation analysis.
