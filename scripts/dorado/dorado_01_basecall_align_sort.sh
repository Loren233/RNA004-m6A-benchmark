#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/../utilities/common.sh"

require_var DORADO
require_var SAMTOOLS
require_var DORADO_MODEL
require_var POD5
require_var REFERENCE
require_var OUTBAM
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-8}}"

"$DORADO" basecaller \
    --device "${DORADO_DEVICE:-cuda:0}" \
    --reference "$REFERENCE" \
    --modified-bases m6A \
    --emit-moves \
    --estimate-poly-a \
    --min-qscore 0 \
    "$DORADO_MODEL" \
    "$POD5" \
    | "$SAMTOOLS" view \
        -u \
        -F 2304 \
        - \
    | "$SAMTOOLS" sort \
        -@ "$THREADS" \
        -o "$OUTBAM" \
        -

"$SAMTOOLS" index -@ "$THREADS" "$OUTBAM"
