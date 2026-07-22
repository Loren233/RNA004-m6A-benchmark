# DeepRM Workflow

The DeepRM route uses three wrappers under `scripts/01_deeprm/`:

```text
POD5 + compatible aligned BAM
  -> deeprm_01_call_prep.sh
  -> deeprm_02_call_run.sh
  -> deeprm_03_make_sitelevel_table.sh
  -> DeepRM transcript-level site table
```

The site-level wrapper writes `transcript_id`, `pos_1based`, `strand`, `depth`,
`mod_count`, `stoich`, `percent_modified`, `modified_base`, `score`, and the
remaining bedMethyl count fields. Standardize it with:

```bash
python scripts/03_standardization/01_standardize_site_tables.py deeprm \
    --input /path/to/deeprm_sitelevel.tsv \
    --output /path/to/deeprm_standardized_transcript.tsv
```

The standardized output contains `transcript_id` and
`transcript_pos_1based`. DeepRM has no separate transcript-to-genome wrapper;
use the generic converter and its column-mapping options:

```bash
python scripts/03_standardization/02_convert_transcript_to_genome.py \
    --input /path/to/deeprm_standardized_transcript.tsv \
    --gtf /path/to/gencode.v49.annotation.gtf \
    --output /path/to/deeprm_standardized_genomic.tsv \
    --transcript-column transcript_id \
    --position-column transcript_pos_1based \
    --position-base 1
```

The converter derives genomic `strand` from the GTF and preserves the
non-coordinate DeepRM fields, including stoichiometry and depth; those fields
do not require separate CLI mappings. `DEEPRM_MODEL_FILE` supplies the model
path to the calling wrappers.
