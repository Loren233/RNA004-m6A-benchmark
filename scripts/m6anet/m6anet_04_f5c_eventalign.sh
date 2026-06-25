#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/../utilities/common.sh"

require_var F5C
require_var BLOW5
require_var KMER_MODEL
require_var BAM
require_var FASTQ
require_var TRANSCRIPT_FASTA
require_var EVENTALIGN
require_var SUMMARY
THREADS="${THREADS:-16}"

"$F5C" eventalign \
    --rna \
    --slow5 "$BLOW5" \
    --kmer-model "$KMER_MODEL" \
    -b "$BAM" \
    -r "$FASTQ" \
    -g "$TRANSCRIPT_FASTA" \
    -o "$EVENTALIGN" \
    --summary "$SUMMARY" \
    --scale-events \
    --signal-index \
    -t "$THREADS"
