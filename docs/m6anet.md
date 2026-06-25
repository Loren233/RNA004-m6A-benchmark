# m6Anet Workflow

This project uses the BLOW5/f5c route for RNA004 signal handling.

```text
POD5
  -> Dorado basecalling/alignment -> sorted BAM
  -> samtools fastq -> FASTQ
  -> blue-crab p2s -> BLOW5
  -> f5c index

sorted BAM + FASTQ + BLOW5 + transcript FASTA + RNA004 5-mer model
  -> f5c eventalign
  -> m6anet dataprep
  -> m6anet inference
  -> transcript-to-genome conversion
  -> GLORI or eTAM-seq comparison
```

The preserved eventalign route is:

```bash
f5c eventalign \
  --rna \
  --slow5 "$BLOW5" \
  --kmer-model "$KMER_MODEL" \
  -b "$BAM" \
  -r "$FASTQ" \
  -g "$TRANSCRIPT_FASTA" \
  -o "$EVENTALIGN" \
  --summary "$SUMMARY" \
  --scale-events \
  --signal-index
```

Inference uses:

```bash
--pretrained_model HEK293T_RNA004
```

For HEK293T, this is the matched RNA004 pretrained model. For HeLa, it is a
cross-cell-line generalization experiment and must not be described as a
HeLa-specific model.
