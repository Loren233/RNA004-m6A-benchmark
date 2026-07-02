#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/../utilities/common.sh"

require_var SAMTOOLS
require_var BAM
require_var WORK

N_CHUNKS="${N_CHUNKS:-24}"

OUT_DIR="${WORK}/ref_chunks"

mkdir -p "$OUT_DIR"

echo "========================================"
echo "m6Anet transcript chunk generation"
echo "BAM=$BAM"
echo "WORK=$WORK"
echo "OUT_DIR=$OUT_DIR"
echo "N_CHUNKS=$N_CHUNKS"
echo "Start: $(date)"
echo "========================================"

[[ -x "$SAMTOOLS" ]] || { echo "ERROR: samtools not executable: $SAMTOOLS" >&2; exit 1; }
[[ -s "$BAM" ]] || { echo "ERROR: BAM missing or empty: $BAM" >&2; exit 1; }

if [[ ! -s "${BAM}.bai" && ! -s "${BAM%.bam}.bai" ]]; then
    echo "ERROR: BAM index missing for $BAM" >&2
    exit 1
fi

rm -rf "$OUT_DIR"
mkdir -p "$OUT_DIR"

echo "[1/3] Extract transcript/reference records with mapped reads..."

"$SAMTOOLS" idxstats "$BAM" \
    | awk 'BEGIN{OFS="\t"} $1!="*" && $3>0 {print $1,0,$2,$3}' \
    | sort -k4,4nr \
    > "${OUT_DIR}/all_transcripts_with_counts.bed"

[[ -s "${OUT_DIR}/all_transcripts_with_counts.bed" ]] || {
    echo "ERROR: no mapped transcript/reference records found from BAM idxstats." >&2
    exit 1
}

echo "[2/3] Split transcript records into ${N_CHUNKS} chunks..."

awk -v n="$N_CHUNKS" -v out="$OUT_DIR" '
BEGIN{
    for(i=1;i<=n;i++){
        files[i]=sprintf("%s/chunk_%02d.bed", out, i)
    }
}
{
    idx=((NR-1)%n)+1
    print $1"\t"$2"\t"$3 > files[idx]
}
' "${OUT_DIR}/all_transcripts_with_counts.bed"

echo "[3/3] Summary..."

wc -l "${OUT_DIR}"/chunk_*.bed
echo
echo "First lines of chunk_01.bed:"
head "${OUT_DIR}/chunk_01.bed"

echo
echo "Output directory:"
du -sh "$OUT_DIR"
ls -lh "$OUT_DIR" | head

echo "Finished: $(date)"
