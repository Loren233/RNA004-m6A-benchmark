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
