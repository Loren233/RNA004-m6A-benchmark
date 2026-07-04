#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/../utilities/common.sh"

require_var M6ANET
require_var INPUT_DIR
require_var OUTDIR

"$M6ANET" inference \
    --input_dir "$INPUT_DIR" \
    --out_dir "$OUTDIR" \
    --pretrained_model "${M6ANET_PRETRAINED_MODEL:-HEK293T_RNA004}" \
    --n_processes "${M6ANET_INFERENCE_PROCESSES:-8}" \
    --batch_size "${M6ANET_BATCH_SIZE:-64}"
