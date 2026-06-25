#!/usr/bin/env bash
set -euo pipefail

cat >&2 <<'USAGE'
Planned comparison interface:

compare_ground_truth.sh \
  --prediction prediction_genomic.tsv \
  --ground-truth ground_truth.tsv \
  --ground-truth-type GLORI_or_eTAM \
  --cell-line HEK293T_or_HeLa

Do not assume GLORI and eTAM column names are identical.
USAGE

exit 2
