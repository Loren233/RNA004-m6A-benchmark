#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/../utilities/common.sh"

require_var F5C
require_var BLOW5
require_var FASTQ
THREADS="${THREADS:-8}"
require_file "$BLOW5"
require_file "$FASTQ"

"$F5C" index -t "$THREADS" --slow5 "$BLOW5" "$FASTQ"
