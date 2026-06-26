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

Transcript-to-genome conversion for DeepRM is not enabled yet because the active
DeepRM site-level transcript table schema still needs verification. The
placeholder `scripts/deeprm/deeprm_04_convert_tx_to_genome.sh` records the
required unresolved fields: transcript ID column, transcript coordinate column,
coordinate base, depth column, and stoichiometry column.
