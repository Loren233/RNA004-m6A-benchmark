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

## Site-Table Standardization

Production Dorado standardization is implemented in
`scripts/09_standardization/standardize_site_tables.py dorado`.

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

Real usage example:

```bash
python scripts/09_standardization/standardize_site_tables.py dorado \
    --input /gpfs/work/bio/zhikefu24/surf/PRJEB80229/dorado/PRJEB80229_dorado_m6A.primary_sitelevel_genomic.tsv \
    --output /gpfs/work/bio/zhikefu24/surf/PRJEB80229/dorado/PRJEB80229_dorado_m6A.primary_standardized_genomic.tsv \
    --tmp-dir /gpfs/work/bio/zhikefu24/surf/PRJEB80229/dorado/tmp_standardize \
    --sort-threads 4 \
    --sort-memory 4G
```

These are documentation example paths only; they are not hardcoded in the
standardization script.
