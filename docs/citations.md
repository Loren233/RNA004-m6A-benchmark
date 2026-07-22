# Citations

## Detection software and methods

### Dorado

- Software: Oxford Nanopore Technologies, Dorado.
- Repository and documentation: https://github.com/nanoporetech/dorado
- Project version: `0.7.3+6e6c45c`.
- Workflow metadata: basecalling model `rna004_130bps_sup@v5.0.0`; modified-base
  calling is invoked with `--modified-bases m6A`.
- Citation policy: Dorado is cited as software using the official repository
  and documentation. The exact RNA004 basecalling model and modified-base
  option are recorded as workflow metadata; no separate model publication is
  cited for this benchmark.

### Modkit

- Software: Oxford Nanopore Technologies, Modkit.
- Repository: https://github.com/nanoporetech/modkit
- Project version: `0.6.1`.

### DeepRM

- Kang, G., Hwang, H., Jeon, H. et al. Comprehensive discovery of m6A sites in
  the human transcriptome at single-molecule resolution. *Nature
  Communications* 17, 664 (2026).
- DOI: https://doi.org/10.1038/s41467-025-67417-w
- Publication: https://www.nature.com/articles/s41467-025-67417-w
- Repository: https://github.com/vadanamu/DeepRM
- Documentation: https://deeprm.readthedocs.io/
- Project environment version: `1.0.9`.
- Model: `deeprm_weights.pt`; SHA256
  `6c2f6444a1b5e29fe653441a3592830a1b29e8ede41f2e868346398c26779522`.

### m6Anet

- Hendra, C. et al. Detection of m6A from direct RNA sequencing using a
  multiple instance learning framework. *Nature Methods* 19, 1590-1598 (2022).
- DOI: https://doi.org/10.1038/s41592-022-01666-1
- Repository: https://github.com/GoekeLab/m6anet
- Project version: `2.1.0`; pretrained model `HEK293T_RNA004`.

### SingleMod

- Xie, Y.-Y. et al. Single-molecule direct RNA sequencing reveals the shaping
  of epitranscriptome across multiple species. *Nature Communications* 16,
  5119 (2025).
- DOI: https://doi.org/10.1038/s41467-025-60447-4
- Repository: https://github.com/xieyy46/SingleMod-v1
- Project commit: `bd7df5043fa1669d44e504b807fbcc37ed07daaf` (`bd7df50-dirty`).

## Supporting software and formats

- minimap2: Li, H. *Bioinformatics* 34, 3094-3100 (2018), DOI
  https://doi.org/10.1093/bioinformatics/bty191; project version `2.30-r1287`.
- SAMtools: Danecek, P. et al. *GigaScience* 10, giab008 (2021), DOI
  https://doi.org/10.1093/gigascience/giab008; project version `1.16.1`.
- f5c: Gamaarachchi, H. et al. *BMC Bioinformatics* 21, 343 (2020), DOI
  https://doi.org/10.1186/s12859-020-03697-x; project version `1.6`.
- blue-crab: https://github.com/Psy-Fer/blue-crab; project release directory
  `blue-crab-v0.4.0`.
- SLOW5/BLOW5: Gamaarachchi, H. et al. *Nature Biotechnology* 40, 1026-1029
  (2022), DOI https://doi.org/10.1038/s41587-021-01147-4.

## Reference annotation

- Mudge, J. M. et al. GENCODE 2025: reference gene annotation for human and
  mouse. *Nucleic Acids Research* 53, D966-D975 (2025).
- DOI: https://doi.org/10.1093/nar/gkae1078
- Release: GENCODE v49, Ensembl 115, GRCh38.p14.
- Release page: https://www.gencodegenes.org/human/release_49.html

## Ground-truth methods and datasets

### GLORI, HEK293T

- Liu, C. et al. Absolute quantification of single-base m6A methylation in the
  mammalian transcriptome using GLORI. *Nature Biotechnology* 41, 355-366
  (2023).
- DOI: https://doi.org/10.1038/s41587-022-01487-9
- GEO series: GSE210563; BioProject: PRJNA866163.
- Samples: GSM6432590 (`HEK293T_GLORI+_rep1`) and GSM6432591
  (`HEK293T_GLORI+_rep2`).
- Local benchmark table: replicate intersection with mean depth and mean
  stoichiometry.

### eTAM-seq, HeLa

- Xiao, Y.-L., Liu, S., Ge, R. et al. Transcriptome-wide profiling and
  quantification of N6-methyladenosine by enzyme-assisted adenosine
  deamination. *Nature Biotechnology* 41, 993-1003 (2023).
- DOI: https://doi.org/10.1038/s41587-022-01587-6
- Publication: https://www.nature.com/articles/s41587-022-01587-6
- GEO SuperSeries: GSE201064.
- eTAM-seq SubSeries: GSE211303; BioProject: PRJNA869740.
- HeLa poly(A) WT samples used: GSM6463599-GSM6463601 (FTO-, rep1-rep3) and
  GSM6463602-GSM6463604 (FTO+, rep1-rep3). GSM6463605-GSM6463607 are the
  corresponding IVT samples in the SubSeries and are not represented as
  biological FTO-/FTO+ replicates in the local benchmark table.
- Local benchmark table: a project-derived genomic replicate-intersection table
  with mean stoichiometry, not an unchanged GEO download.
