#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/../utilities/common.sh"

PY="${PY:-${PYTHON:-python}}"
CONVERTER="${CONVERTER:-$(dirname "$0")/../09_standardization/convert_transcript_to_genome.py}"

require_var SITE_TSV
require_var GTF
require_var OUT_TSV
require_file "$SITE_TSV"
require_file "$GTF"

"$PY" "$CONVERTER" \
    --input "$SITE_TSV" \
    --gtf "$GTF" \
    --output "$OUT_TSV" \
    --transcript-column chr \
    --position-column start_0based \
    --position-base 0 \
    --depth-column Total \
    --stoich-column m6A_fraction \
    --modified-count-column m6A_count \
    --percent-column m6A_percent \
    --tool-name dorado
