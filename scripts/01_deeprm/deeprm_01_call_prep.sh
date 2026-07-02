#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/../utilities/common.sh"

require_var DEEPRM
require_var POD5
require_var BAM
require_var PREP_DIR
require_file "$BAM"

"$DEEPRM" call prep \
    -p "$POD5" \
    -b "$BAM" \
    -o "$PREP_DIR"
