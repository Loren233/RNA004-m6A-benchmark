# Dorado Workflow

Planned source flow:

```text
POD5 -> Dorado basecaller with m6A calling and alignment
     -> sorted modBAM with MM/ML tags
     -> samtools index
     -> modkit pileup
     -> site-level bedMethyl
     -> standardized site-level TSV
     -> transcript-to-genome conversion only if transcript reference was used
     -> GLORI or eTAM-seq comparison
```

The command template preserves these options from `pipeline.docx`:

- `--device cuda:0`
- `--reference "$REFERENCE"`
- `--modified-bases m6A`
- `--emit-moves`
- `--estimate-poly-a`
- `--min-qscore 0`

If `REFERENCE` is a transcript FASTA, downstream site tables may need
transcript-to-genome conversion. If `REFERENCE` is a GRCh38 genome FASTA, do not
convert again.

When Dorado is aligned to the cleaned GENCODE transcript FASTA, the site-level
table column named `chr` contains transcript IDs such as `ENST00000832841.1`.
Use `scripts/dorado/dorado_04_convert_tx_to_genome.sh` to call the generic
converter before genomic standardization. The wrapper maps `chr` to
`transcript_id` and uses `start_0based` as the 0-based transcript coordinate.
