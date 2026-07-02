#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/../utilities/common.sh"

require_var BED
require_var TSV
require_file "$BED"

awk 'BEGIN{
OFS="\t";
print "chr","pos_1based","start_0based","end","strand","m6A_count","Total","m6A_percent","m6A_fraction"
}
{
print $1,$2+1,$2,$3,$6,$12,$10,$11,$11/100
}' "$BED" > "$TSV"
