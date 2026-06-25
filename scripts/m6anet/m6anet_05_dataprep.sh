#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/../utilities/common.sh"

require_var M6ANET
require_var EVENTALIGN
require_var INPUT_DIR
require_file "$EVENTALIGN"

"$M6ANET" dataprep \
    --eventalign "$EVENTALIGN" \
    --out_dir "$INPUT_DIR" \
    --n_processes "${M6ANET_DATAPREP_PROCESSES:-4}"
