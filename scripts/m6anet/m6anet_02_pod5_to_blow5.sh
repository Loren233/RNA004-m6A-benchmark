#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/../utilities/common.sh"

require_var BLUECRAB
require_var POD5
require_var BLOW5
THREADS="${THREADS:-8}"

"$BLUECRAB" p2s -t "$THREADS" "$POD5" -o "$BLOW5"
