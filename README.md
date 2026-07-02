# RNA004-m6A-benchmark

Reproducible benchmarking scaffold for RNA004 nanopore direct RNA sequencing
m6A detection. The project compares four detection routes:

- Dorado
- DeepRM
- m6Anet
- SingleMod

Supported cell-line benchmarks:

- HEK293/HEK293T RNA004 direct RNA data, primarily compared with GLORI.
- HeLa RNA004 direct RNA data, primarily compared with eTAM-seq.

Do not label GLORI as the HeLa ground truth, and do not label eTAM-seq as the
HEK293T ground truth.

## Data Flow by Tool

### Dorado

Citation: Oxford Nanopore Technologies. Dorado: Oxford Nanopore's basecaller. GitHub/software documentation. https://github.com/nanoporetech/dorado

Citation: Oxford Nanopore Technologies. Modkit: a bioinformatics tool for working with modified bases from Oxford Nanopore, including conversion of modBAM to bedMethyl. GitHub/software documentation. https://github.com/nanoporetech/modkit

```text
POD5
  -> Dorado basecaller with m6A modified-base calling and alignment
sorted modBAM with MM/ML tags
  -> samtools index
BAM + BAI
  -> modkit pileup
site-level bedMethyl
  -> standardize and filter columns
site-level TSV
  -> transcript-to-genome conversion only if transcript FASTA was used
GRCh38 genomic site-level TSV
  -> GLORI comparison for HEK293T or eTAM-seq comparison for HeLa
```

Dorado read-level modified-base information is stored in MM/ML tags. Modkit is
used for site-level aggregation. Each run must document whether `REFERENCE` was
a transcript FASTA or a GRCh38 genome FASTA.

### DeepRM

Citation: Kang, G. et al. Comprehensive discovery of m6A sites in the human transcriptome at single-molecule resolution. Nature Communications 17, 664 (2026). https://doi.org/10.1038/s41467-025-67417-w

```text
POD5 + compatible aligned BAM
  -> deeprm call prep
prepared DeepRM input
  -> deeprm call run
read-level prediction
  -> site-level aggregation
transcript-level or genomic site table
  -> transcript-to-genome conversion only when required
GRCh38 genomic site-level table
  -> GLORI comparison for HEK293T or eTAM-seq comparison for HeLa
```

The DeepRM model path should be configured with `DEEPRM_MODEL_FILE`. Do not
commit personal installation paths.

### m6Anet

Citation: Hendra, C. et al. Detection of m6A from direct RNA sequencing using a multiple instance learning framework. Nature Methods 19, 1590–1598 (2022). https://doi.org/10.1038/s41592-022-01666-1

```text
POD5
  -> Dorado basecalling/alignment -> sorted BAM -> samtools fastq -> FASTQ
  -> blue-crab p2s -> BLOW5 -> f5c index

sorted BAM + FASTQ + BLOW5 + transcript FASTA + RNA004 5-mer model
  -> f5c eventalign
eventalign.tsv
  -> m6anet dataprep
prepared input
  -> m6anet inference with --pretrained_model HEK293T_RNA004
data.site_proba.csv
  -> transcript-to-genome conversion
GRCh38 genomic site-level TSV
  -> GLORI comparison for HEK293T or eTAM-seq comparison for HeLa
```

This project uses the BLOW5/f5c RNA004 signal route, not a FAST5-based m6Anet
workflow.

For HEK293T, `HEK293T_RNA004` is documented as the matched RNA004 pretrained
model. For HeLa, the same model may be run only as a cross-cell-line
generalization experiment; it is not a HeLa-specific m6Anet model.

### SingleMod

Citation: Xie, Y.-Y. et al. Single-molecule direct RNA sequencing reveals the shaping of epitranscriptome across multiple species. Nature Communications 16, 5119 (2025). https://doi.org/10.1038/s41467-025-60447-4

```text
POD5
  -> blue-crab p2s
BLOW5

genome-aligned BAM
  -> samtools fastq -> FASTQ
  -> split by chromosome -> chromosome BAM -> bedtools bamtobed

chromosome BAM + BED + FASTQ + BLOW5 + GRCh38 FASTA
  -> f5c eventalign
chromosome eventalign
  -> organize_from_eventalign.py
temporary motif features
  -> merge_motif_npy.py
merged features
  -> 36 RNA004 motif models
read-level predictions
  -> probability > 0.5 and site aggregation
chromosome mr.bed
  -> merge chromosomes
genomic site-level TSV
  -> GLORI comparison for HEK293T or eTAM-seq comparison for HeLa
```

SingleMod output is already genomic when it uses chromosome-level genome-aligned
BAM and GRCh38 FASTA. Do not add a GTF transcript-to-genome conversion step for
that route. Site-level stoichiometry is calculated as `m6A_count / Total`.

## Reference and Coordinate Conventions

Citation: GENCODE Consortium. GENCODE Human Release 49 (GRCh38.p14; corresponding to Ensembl 115). https://www.gencodegenes.org/human/release_49.html

Citation: Frankish, A. et al. GENCODE: reference annotation for the human and mouse genomes in 2023. Nucleic Acids Research 51(D1), D942–D949 (2023). https://doi.org/10.1093/nar/gkac1071

See `docs/references_and_coordinates.md`. In short:

- GENCODE v49 transcript FASTA headers should be cleaned to retain transcript
  IDs.
- Transcript FASTA names should match GTF `transcript_id`.
- BED start is 0-based.
- GTF genomic coordinates are 1-based inclusive.
- m6Anet `transcript_position` is treated as 0-based by the current planned
  conversion.
- Site tables should use explicit names such as `start_0based`, `pos_1based`,
  and `transcript_position_0based`.

## Repository Layout

```text
RNA004-m6A-benchmark/
├── README.md
├── docs/
│   ├── pipeline_overview.md
│   ├── references_and_coordinates.md
│   ├── ground_truths.md
│   ├── dorado.md
│   ├── deeprm.md
│   ├── m6anet.md
│   ├── singlemod.md
│   └── source/pipeline.docx
├── configs/
│   ├── example.env
│   ├── hek293t.example.env
│   ├── hela.example.env
│   ├── paths.example.sh
│   ├── hek293t.yaml
│   └── hela.yaml
├── scripts/
│   ├── dorado/
│   ├── deeprm/
│   ├── m6anet/
│   ├── singlemod/
│   ├── comparison/
│   └── utilities/
├── metadata/
├── slurm/
├── results/
└── archive/
```

The older numbered scaffold directories under `scripts/` are retained as empty
placeholders from the initial project layout. The tool-oriented directories are
the preferred location for workflow templates derived from `pipeline.docx`.

## Script Execution Order

### HEK293T + GLORI

Citation: the current HEK293T ground truth is derived from GLORI samples
GSM6432590 and GSM6432591 in GSE210563; cite the GLORI method paper and GEO
records listed in `docs/citations.md`.

1. Configure local paths from `configs/example.env` and
   `configs/hek293t.example.env`.
2. Run Dorado basecalling/alignment or provide compatible existing BAM/FASTQ
   inputs.
3. Run each selected tool: Dorado, DeepRM, m6Anet, and/or SingleMod.
4. Convert transcript-level outputs to GRCh38 genomic coordinates where needed.
5. Compare genomic predictions to GLORI.
6. Run depth-stratified correlation and visualization steps.

### HeLa + eTAM-seq

1. Configure local paths from `configs/example.env` and
   `configs/hela.example.env`.
2. Run Dorado basecalling/alignment or provide compatible existing BAM/FASTQ
   inputs.
3. Run selected tools. For m6Anet, document use of `HEK293T_RNA004` as a
   cross-cell-line generalization experiment.
4. Convert transcript-level outputs to GRCh38 genomic coordinates where needed.
5. Compare genomic predictions to eTAM-seq.
6. Run depth-stratified correlation and visualization steps.

## Software and Environments

Expected tools include Dorado, Samtools, Modkit, f5c, blue-crab, DeepRM,
m6Anet, SingleMod, bedtools, Python, and shell utilities. Verified versions
should be recorded in `metadata/software_versions.tsv`; do not invent versions.

## Data Policy

Large data files and generated outputs must not be committed to Git. This
includes POD5, FAST5, FASTQ, BAM, BAI, SAM, BLOW5/SLOW5, eventalign files,
model weights, large TSV/CSV/BED outputs, logs, and generated results.

Current status: scaffold and documentation only. Scientific pipeline stages
must be run manually after local paths, references, models, and ground-truth
formats are verified.
