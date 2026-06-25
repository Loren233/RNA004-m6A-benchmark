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
