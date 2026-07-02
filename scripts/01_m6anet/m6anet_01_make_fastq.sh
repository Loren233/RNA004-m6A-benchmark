#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/../utilities/common.sh"

require_var SAMTOOLS
require_var BAM
require_var FASTQ
THREADS="${THREADS:-8}"
require_file "$BAM"

"$SAMTOOLS" fastq -@ "$THREADS" -n -F 2304 "$BAM" > "$FASTQ"
