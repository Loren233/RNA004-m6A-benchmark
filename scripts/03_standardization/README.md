# 03_standardization

This folder contains scripts for converting tool-specific RNA modification outputs into a common genomic site-level format before benchmarking.

The general order is:

1. `01_standardize_site_tables.py`
2. `02_convert_transcript_to_genome.py`
3. `03_collapse_genomic_sites.py`

However, different tools enter this workflow at different steps because their intermediate outputs are not in the same format.

## Tool-specific entry points

| Tool | Entry point | Reason |
|---|---|---|
| Dorado / Modkit | `01_standardize_site_tables.py` | Dorado/Modkit outputs are already in genomic coordinates, but their column names and count fields are tool-specific. |
| m6Anet | `01_standardize_site_tables.py` | m6Anet `data.site_proba.csv` files must first be merged and standardized into the shared site-level format. |
| SingleMod | `02_convert_transcript_to_genome.py` | SingleMod skips `01_standardize_site_tables.py` because `singlemod_03_collect_sitelevel.sh` already creates `SingleMod_sitelevel_m6A_Total.tsv`. |

## Dorado / Modkit

Dorado/Modkit outputs are already in genomic coordinates. Therefore, transcript-to-genome conversion is usually not needed.

Recommended workflow:

1. `01_standardize_site_tables.py`
2. `03_collapse_genomic_sites.py`, only if duplicated genomic sites exist
3. `04_benchmarking`

## m6Anet

m6Anet inference produces chunk-level `data.site_proba.csv` files. These files are first merged into one transcript-level raw table, such as:

`m6anet_site_proba.transcript_raw.csv`

This merged table still has m6Anet-specific columns, so it must be standardized before transcript-to-genome conversion.

Recommended workflow:

1. Merge chunk-level `data.site_proba.csv` files
2. `01_standardize_site_tables.py`
3. `02_convert_transcript_to_genome.py`
4. `03_collapse_genomic_sites.py`
5. `04_benchmarking`

## SingleMod

SingleMod skips `01_standardize_site_tables.py`.

The SingleMod collection step already converts chunk-level `mr.bed` files into a transcript-level site table:

`SingleMod_sitelevel_m6A_Total.tsv`

Recommended workflow:

1. `singlemod_03_collect_sitelevel.sh`
2. `02_convert_transcript_to_genome.py`
3. `03_collapse_genomic_sites.py`
4. `04_benchmarking`

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
- `source`

This step is required for Dorado/Modkit and m6Anet.

This step can be skipped if an earlier tool-specific collection step already produces a compatible transcript-level or genomic site-level table. SingleMod is the current example.

## 02_convert_transcript_to_genome.py

This script converts transcript-level coordinates to genomic coordinates using a GTF annotation.

Use this step when the standardized output is still based on transcript IDs.

Example transcript-level coordinate:

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

## Notes

Do not benchmark transcript-level coordinates directly against genomic ground truth datasets.

For example, a SingleMod or m6Anet site reported on a transcript such as `ENST00000000233.10:122` must first be converted to a genomic coordinate before it can be compared with GLORI or eTAM.

If multiple standardized rows map to the same genomic site, use `03_collapse_genomic_sites.py` before overlap and correlation analysis.
