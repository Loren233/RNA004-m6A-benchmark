#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/../utilities/common.sh"

require_var DEEPRM
require_var DEEPRM_MODEL_FILE
require_var PREP_DIR
require_var BAM
require_var CALL_DIR
require_file "$BAM"
require_file "$DEEPRM_MODEL_FILE"

"$DEEPRM" call run \
    --model "$DEEPRM_MODEL_FILE" \
    --input "$PREP_DIR" \
    --bam "$BAM" \
    --output "$CALL_DIR" \
    --gpu "${DEEPRM_GPU:-1}" \
    --gpu-pool "${DEEPRM_GPU_POOL:-0}" \
    --batch "${DEEPRM_BATCH:-5000}"
