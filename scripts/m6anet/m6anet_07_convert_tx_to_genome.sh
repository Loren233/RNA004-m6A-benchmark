#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/../utilities/common.sh"

PY="${PY:-${PYTHON:-python}}"
CONVERTER="${CONVERTER:-$(dirname "$0")/../09_standardization/convert_transcript_to_genome.py}"

require_var SITE_CSV
require_var GTF
require_var OUT_TSV
require_file "$SITE_CSV"
require_file "$GTF"

"$PY" "$CONVERTER" \
    --input "$SITE_CSV" \
    --gtf "$GTF" \
    --output "$OUT_TSV" \
    --input-separator "," \
    --transcript-column transcript_id \
    --position-column transcript_position \
    --position-base 0 \
    --depth-column n_reads \
    --stoich-column mod_ratio \
    --tool-name m6anet
