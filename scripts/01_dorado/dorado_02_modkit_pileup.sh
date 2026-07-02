#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/../utilities/common.sh"

require_var MODKIT
require_var BAM
require_var OUTBED
require_var REFERENCE
require_file "$BAM"

"$MODKIT" pileup \
    "$BAM" \
    "$OUTBED" \
    --ref "$REFERENCE" \
    --log-filepath "${LOG:-${OUTBED}.log}"
