# SingleMod Workflow

Planned practical flow:

```text
POD5 -> blue-crab p2s -> BLOW5
genome-aligned BAM -> samtools fastq -> FASTQ
genome-aligned BAM -> split by chromosome -> chromosome BAM -> BED
chromosome BAM + BED + FASTQ + BLOW5 + GRCh38 FASTA
  -> f5c eventalign
  -> organize_from_eventalign.py
  -> merge_motif_npy.py
  -> 36 RNA004 motif models
  -> read-level predictions
  -> probability > 0.5 and site aggregation
  -> chromosome mr.bed
  -> merged genomic site-level TSV
  -> GLORI or eTAM-seq comparison
```

The 36 RNA004 motifs from `pipeline.docx` must be preserved when implementing
the chromosome-array stage.

SingleMod output is already genomic when this route uses chromosome-level
genome-aligned BAM and GRCh38 FASTA. Do not add a GTF transcript-to-genome
conversion step for that route.

Site-level stoichiometry:

```text
m6A_count / Total
```
