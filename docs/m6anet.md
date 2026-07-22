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

After inference, `scripts/01_m6anet/m6anet_06_merge_site_proba.sh` merges
chunk-level `data.site_proba.csv` files. The merged table is standardized with
`scripts/03_standardization/01_standardize_site_tables.py m6anet`, then
converted with `scripts/03_standardization/02_convert_transcript_to_genome.py`.
The standardizer treats `transcript_position` as 0-based by default and
preserves m6Anet fields such as `probability_modified`, `mod_ratio`, `kmer`,
and read depth.
