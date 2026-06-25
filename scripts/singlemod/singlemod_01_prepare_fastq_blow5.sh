#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/../utilities/common.sh"

require_var SAMTOOLS
require_var BLUECRAB
require_var F5C
require_var BAM
require_var POD5
require_var FASTQ
require_var BLOW5
THREADS="${THREADS:-8}"

"$SAMTOOLS" fastq -@ "$THREADS" -n -F 2304 "$BAM" > "$FASTQ"
"$BLUECRAB" p2s -t "$THREADS" "$POD5" -o "$BLOW5"
"$F5C" index -t "$THREADS" --slow5 "$BLOW5" "$FASTQ"
