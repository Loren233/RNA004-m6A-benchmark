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
No Dorado-specific transcript-to-genome wrapper is present in the current
repository. If Dorado was aligned to transcript coordinates, use the generic
converter at `scripts/03_standardization/02_convert_transcript_to_genome.py`
after confirming the transcript ID and coordinate columns for that run.

## Site-Table Standardization

Production Dorado standardization is implemented in
`scripts/03_standardization/01_standardize_site_tables.py dorado`.

Workflow:

```text
Dorado primary-only genomic site records
  -> external genomic-coordinate sort
  -> genomic-site count aggregation
  -> shared standardized site table
```

The standardizer aggregates transcript-derived records by
`chr + pos_1based + strand`. For each genomic site it sums
`dorado_m6A_count` and `dorado_depth`, then calculates
`stoich = modified_count / depth`. It does not average `dorado_stoich`.

Usage example:

```bash
INPUT=/path/to/dorado_sitelevel_genomic.tsv
OUTPUT=/path/to/prediction_standardized.tsv
TMPDIR=/path/to/tmp_standardize

python scripts/03_standardization/01_standardize_site_tables.py dorado \
    --input "$INPUT" \
    --output "$OUTPUT" \
    --tmp-dir "$TMPDIR" \
    --sort-threads 4 \
    --sort-memory 4G
```

These are documentation example paths only; they are not hardcoded in the
standardization script.
