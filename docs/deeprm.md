# DeepRM Workflow

Planned source flow:

```text
POD5 + compatible aligned BAM
  -> deeprm call prep
  -> deeprm call run
  -> read-level predictions
  -> site-level aggregation
  -> transcript-to-genome conversion when required
  -> GLORI or eTAM-seq comparison
```

The model path should be supplied with `DEEPRM_MODEL_FILE`. The source document
included a personal installation path; this repository replaces it with a
portable variable and does not commit that path.

DeepRM inference is a GPU task according to the current command template.

The current DeepRM route contains three shell wrappers under
`scripts/01_deeprm/`:

- `deeprm_01_call_prep.sh`
- `deeprm_02_call_run.sh`
- `deeprm_03_make_sitelevel_table.sh`

`deeprm_03_make_sitelevel_table.sh` writes a transcript-level site table with
these columns:

```text
transcript_id
pos_1based
strand
depth
mod_count
stoich
percent_modified
modified_base
score
bed_start_0based
bed_end
canonical_count
other_mod_count
delete_count
fail_count
diff_count
nocall_count
```

The current standardizer supports DeepRM through:

```bash
python scripts/03_standardization/01_standardize_site_tables.py deeprm \
    --input /path/to/deeprm_sitelevel.tsv \
    --output /path/to/deeprm_standardized_transcript.tsv
```

If the standardized DeepRM table remains transcript-coordinate based, convert it
with the generic converter after confirming the GTF path:

```bash
python scripts/03_standardization/02_convert_transcript_to_genome.py \
    --input /path/to/deeprm_standardized_transcript.tsv \
    --gtf /path/to/gencode.v49.annotation.gtf \
    --output /path/to/deeprm_standardized_genomic.tsv \
    --transcript-column transcript_id \
    --position-column transcript_pos_1based \
    --position-base 1
```

No DeepRM-specific transcript-to-genome wrapper is present in the current
repository.
