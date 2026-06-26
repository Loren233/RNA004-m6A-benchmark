# Ground Truths

Ground-truth processing is cell-line specific.

## HEK293/HEK293T

Primary benchmark: GLORI.

HEK293/HEK293T nanopore predictions should be compared primarily with GLORI
after both prediction and ground-truth tables are represented in compatible
GRCh38 genomic coordinates.

Current ground-truth input pattern:

```text
metadata/ground_truth_datasets/ground_truth.tsv
```

Processing order:

```text
Uploaded final HEK293T GLORI ground-truth table
  -> scripts/08_ground_truth/glori/prepare_glori_sites.py
validated prepared genomic site table
  -> standardization / benchmarking
```

The input table is already a genomic replicate-intersection ground truth and
should not undergo transcript-to-genome conversion.

## HeLa

Primary benchmark: eTAM-seq.

HeLa nanopore predictions should be compared primarily with eTAM-seq. Do not
label GLORI as the HeLa ground truth.

Current ground-truth input pattern:

```text
metadata/ground_truth_datasets/ground_truth.tsv
```

Processing order:

```text
Uploaded final HeLa eTAM-seq ground-truth table
  -> scripts/08_ground_truth/etam_seq/prepare_etam_sites.py
validated prepared genomic site table
  -> standardization / benchmarking
```

The input table is already a genomic replicate-intersection ground truth and
should not undergo transcript-to-genome conversion.

## Comparison Interface

Use separate preprocessing scripts or parameterized logic that explicitly names
the ground-truth type:

```bash
compare_ground_truth.sh \
  --prediction prediction_genomic.tsv \
  --ground-truth ground_truth.tsv \
  --ground-truth-type GLORI \
  --cell-line HEK293T
```

```bash
compare_ground_truth.sh \
  --prediction prediction_genomic.tsv \
  --ground-truth ground_truth.tsv \
  --ground-truth-type eTAM \
  --cell-line HeLa
```

Do not assume GLORI and eTAM-seq column names are identical. Inspect and
document their preprocessing separately before unifying comparison code.
