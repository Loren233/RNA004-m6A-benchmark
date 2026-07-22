# Citations

This file records citations for software, methods, reference resources, and
ground-truth datasets used or currently documented in this RNA004 m6A benchmark.
Unresolved items are marked `TO VERIFY` rather than guessed.

## Software

### Dorado

- Citation: Oxford Nanopore Technologies. Dorado: Oxford Nanopore's basecaller.
- Repository: https://github.com/nanoporetech/dorado
- Local version: `0.7.3+6e6c45c`
- RNA004 model: `rna004_130bps_sup@v5.0.0`
- Used for: basecalling, RNA004 modified-base calling, alignment when
  `--reference` is provided, and MM/ML tag generation.
- Note: Modkit is cited separately; it is not treated as part of Dorado.

### Modkit

- Citation: Oxford Nanopore Technologies. Modkit: a bioinformatics tool for
  working with modified bases.
- Repository: https://github.com/nanoporetech/modkit
- Local version: `0.6.1`
- Used for: aggregation of modified-base BAM calls into site-level bedMethyl
  output.

### minimap2

- Citation: Li, H. Minimap2: pairwise alignment for nucleotide sequences.
  Bioinformatics 34, 3094-3100 (2018).
- DOI: https://doi.org/10.1093/bioinformatics/bty191
- Additional official citation: Li, H. New strategies to improve minimap2
  alignment accuracy. Bioinformatics 37, 4572-4574 (2021).
- DOI: https://doi.org/10.1093/bioinformatics/btab705
- Repository: https://github.com/lh3/minimap2
- Local version: `2.30-r1287`
- Used for: transcriptome alignment in `m6anet_03_align_transcriptome.sh`.

### SAMtools

- Citation: Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J.,
  Homer, N., Marth, G., Abecasis, G., Durbin, R. and 1000 Genome Project Data
  Processing Subgroup. The Sequence Alignment/Map format and SAMtools.
  Bioinformatics 25, 2078-2079 (2009).
- DOI: https://doi.org/10.1093/bioinformatics/btp352
- Current recommended SAMtools/BCFtools citation: Danecek, P. et al. Twelve
  years of SAMtools and BCFtools. GigaScience 10, giab008 (2021).
- DOI: https://doi.org/10.1093/gigascience/giab008
- Official documentation: https://www.htslib.org/doc/
- Local version: `1.16.1`
- Used for: BAM sorting, indexing, filtering, and FASTQ extraction.

### f5c

- Citation: Gamaarachchi, H., Lam, C. W., Jayatilaka, G., Samarakoon, H.,
  Simpson, J. T., Smith, M. A. and Parameswaran, S. GPU accelerated adaptive
  banded event alignment for rapid comparative nanopore signal analysis.
  BMC Bioinformatics 21, 343 (2020).
- DOI: https://doi.org/10.1186/s12859-020-03697-x
- Repository: https://github.com/hasindu2008/f5c
- Local version: `1.6`
- Used for: RNA004 event alignment in m6Anet and SingleMod workflows.

### blue-crab

- Citation: blue-crab: lossless nanopore POD5 to SLOW5/BLOW5 file conversion.
- Repository: https://github.com/Psy-Fer/blue-crab
- Local release directory: `blue-crab-v0.4.0`
- Version note: version output unavailable locally because the launcher requires
  `python3.9`.
- Used for: POD5 to BLOW5 conversion.
- Note: cite blue-crab as software and cite SLOW5/BLOW5 separately as the file
  format.

## Modified-Base Methods

### DeepRM

- Local version: `1.0.9`
- Model: `deeprm_weights.pt`
- Model SHA256:
  `6c2f6444a1b5e29fe653441a3592830a1b29e8ede41f2e868346398c26779522`
- Used for: DeepRM m6A detection.
- Citation: TO VERIFY: original DeepRM m6A method publication.
- Repository: TO VERIFY: official DeepRM repository.

### m6Anet

- Citation: Hendra, C., Pratanwanich, P. N., Wan, Y. K., Goh, W. S. S.,
  Thiery, A. and Göke, J. Detection of m6A from direct RNA sequencing using a
  multiple instance learning framework. Nature Methods 19, 1590-1598 (2022).
- DOI: https://doi.org/10.1038/s41592-022-01666-1
- Repository: https://github.com/GoekeLab/m6anet
- Local version: `2.1.0`
- Pretrained model: `HEK293T_RNA004`
- Used for: site-level m6A prediction from direct RNA sequencing eventalign
  features.
- Note: in HeLa, this model is a cross-cell-line use of the HEK293T RNA004
  model, not a HeLa-specific model.

### SingleMod

- Citation: Xie, Y.-Y., Zhong, Z.-D., Chen, H.-X. et al. Single-molecule direct
  RNA sequencing reveals the shaping of epitranscriptome across multiple
  species. Nature Communications 16, 5119 (2025).
- DOI: https://doi.org/10.1038/s41467-025-60447-4
- Repository: https://github.com/xieyy46/SingleMod-v1
- Local commit: `bd7df5043fa1669d44e504b807fbcc37ed07daaf`
- Local git state: `bd7df50-dirty`
- RNA004 models: 36 motif-specific model files.
- Local modification: `SingleMod/organize_from_eventalign.py` includes
  `if can_motif not in idx: continue`.
- Used for: single-molecule and site-level m6A prediction from RNA004 direct RNA
  data.

## Signal Formats

### SLOW5/BLOW5

- Citation: Gamaarachchi, H., Samarakoon, H., Jenner, S. P., Ferguson, J. M.,
  Amos, T. G., Hammond, J. M., Saadat, H., Smith, M. A., Parameswaran, S. and
  Deveson, I. W. Fast nanopore sequencing data analysis with SLOW5.
  Nature Biotechnology 40, 1026-1029 (2022).
- DOI: https://doi.org/10.1038/s41587-021-01147-4
- Format documentation: https://hasindu2008.github.io/slow5specs
- Used for: BLOW5 signal files in m6Anet and SingleMod workflows.
- Note: SLOW5/BLOW5 is cited as a signal file format, not as a software version.

## Reference Annotation

### GENCODE v49

- Release page: https://www.gencodegenes.org/human/release_49.html
- GENCODE release: v49
- Ensembl release: 115
- Assembly annotation: GRCh38 / GRCh38.p14 on the GENCODE release page
- GTF date: 2025-07-08
- Files:
  - `gencode.v49.annotation.gtf.gz`
    MD5 `0ef4a024ea2d35b1b88c12447b0b70b9`
  - `gencode.v49.transcripts.fa`
    MD5 `7e3f327f3da54f5beef110a2e4a5df59`
  - `gencode.v49.transcripts.clean.fa`
    MD5 `000e70bbdf1c42fe42f4f9cce1f74944`
    with 533740 transcript sequences
- Local derived resource: `gencode.v49.transcripts.clean.fa`, with headers
  cleaned to retain only transcript IDs.
- Project publication: Mudge, J. M. et al. GENCODE 2025: reference gene
  annotation for human and mouse. Nucleic Acids Research 53, D966-D975 (2025).
- DOI: https://doi.org/10.1093/nar/gkae1078
- Note: RefSeq is not cited because it is not used in the current workflow.
- Note: this project does not claim that a full
  `GRCh38.primary_assembly.genome.fa` was used.

## Ground-Truth Methods and Datasets

### GLORI HEK293T Ground Truth

- Method citation: Liu, C., Sun, H., Yi, Y., Shen, W., Li, K. et al. Absolute
  quantification of single-base m6A methylation in the mammalian transcriptome
  using GLORI. Nature Biotechnology 41, 355-366 (2023).
- DOI: https://doi.org/10.1038/s41587-022-01487-9
- GEO series: GSE210563, "Absolute quantification of single-base m6A
  methylation in the mammalian transcriptome".
- BioProject: PRJNA866163
- GEO samples used:
  - GSM6432590, HEK293T_GLORI+_rep1, SRA SRX17361939,
    BioSample SAMN30156017.
  - GSM6432591, HEK293T_GLORI+_rep2, SRA SRX17361940,
    BioSample SAMN30156016.
- Source processed files:
  - `GSM6432590_HEK293T_GLORI_plus_rep1_totalm6A_FDR.csv.gz`
  - `GSM6432591_HEK293T_GLORI_plus_rep2_totalm6A_FDR.csv.gz`
- Intermediate site-level files:
  - `GSM6432590_HEK293T_GLORI_plus_rep1_site_level.tsv`
  - `GSM6432591_HEK293T_GLORI_plus_rep2_site_level.tsv`
- Final benchmark table: local replicate-intersection ground-truth table.
- Final benchmark file properties:
  - Total lines: 176650
  - Data rows: 176649
  - MD5: `f738a525c029499a376b6b0dd4cd3042`
- Interpretation: intersection of two HEK293T GLORI replicates with mean depth
  and mean stoichiometry.

### HeLa / eTAM-seq

The repository contains a HeLa eTAM-seq replicate-intersection ground-truth table
under `metadata/ground_truth_datasets/` and selected HeLa benchmark figures
under `results/`.

- Citation: TO VERIFY official eTAM-seq publication.
- Dataset accession: TO VERIFY source accession for the eTAM-seq data used to
  construct the local ground-truth table.

## Unresolved Citations

- TO VERIFY: original DeepRM m6A method publication and official repository.
- TO VERIFY: official eTAM-seq publication and source dataset accession.
- TO VERIFY: whether a separate Dorado model documentation citation is needed
  for `rna004_130bps_sup@v5.0.0` beyond the official Dorado repository/docs.
