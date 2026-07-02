#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/../utilities/common.sh"

require_var WORK

N_CHUNKS="${N_CHUNKS:-24}"
ALLOW_MISSING="${ALLOW_MISSING:-0}"

SITE="${WORK}/sitelevel"
STATUS="${WORK}/status"
LOG="${WORK}/logs"

mkdir -p "$SITE" "$STATUS" "$LOG"

OUTBED="${SITE}/SingleMod_m6A_mr.all.bed"
OUTTSV="${SITE}/SingleMod_sitelevel_m6A_Total.tsv"
LIST="${SITE}/SingleMod_chunk_mr_bed_list.txt"
MISSING="${SITE}/SingleMod_missing_chunks.txt"
STATS="${SITE}/SingleMod_chunk_site_counts.tsv"

echo "========================================"
echo "SingleMod collect site-level outputs"
echo "WORK=${WORK}"
echo "SITE=${SITE}"
echo "N_CHUNKS=${N_CHUNKS}"
echo "ALLOW_MISSING=${ALLOW_MISSING}"
echo "Start: $(date)"
echo "========================================"

rm -f "$OUTBED" "$OUTTSV" "$LIST" "$MISSING" "$STATS"

echo -e "chunk\tmr_bed\tline_count" > "$STATS"

N_OK=0
N_MISSING=0

for i in $(seq -w 1 "$N_CHUNKS"); do
    chunk="chunk_${i}"
    bed="${SITE}/${chunk}/${chunk}.mr.bed"

    if [[ -s "$bed" ]]; then
        n=$(wc -l < "$bed")
        echo "$bed" >> "$LIST"
        echo -e "${chunk}\t${bed}\t${n}" >> "$STATS"
        echo "OK ${chunk}: lines=${n}"
        N_OK=$((N_OK + 1))
    else
        echo "$chunk" >> "$MISSING"
        echo -e "${chunk}\tMISSING\t0" >> "$STATS"
        echo "MISSING ${chunk}"
        N_MISSING=$((N_MISSING + 1))
    fi
done

echo "========================================"
echo "Completed chunk files: ${N_OK}"
echo "Missing chunk files: ${N_MISSING}"
echo "========================================"

if [[ "$N_MISSING" -gt 0 && "$ALLOW_MISSING" -eq 0 ]]; then
    echo "ERROR: some chunk mr.bed files are missing:"
    cat "$MISSING"
    echo "Do not collect until all chunks are completed, or set ALLOW_MISSING=1."
    exit 1
fi

if [[ ! -s "$LIST" ]]; then
    echo "ERROR: no chunk mr.bed files found." >&2
    exit 1
fi

echo "Merging chunk-level mr.bed files..."

cat $(cat "$LIST") \
    | sort -k1,1 -k2,2n \
    > "$OUTBED"

if [[ ! -s "$OUTBED" ]]; then
    echo "ERROR: merged BED is empty: $OUTBED" >&2
    exit 1
fi

echo "Creating TSV with header..."

awk 'BEGIN{
    OFS=FS="\t";
    print "ref","pos_1based","start_0based","end","site_id","singlemod_stoich","strand","Total","m6A_count","motif"
}
NF >= 9 {
    print $1,$3,$2,$3,$4,$5,$6,$7,$8,$9
}' "$OUTBED" > "$OUTTSV"

if [[ ! -s "$OUTTSV" ]]; then
    echo "ERROR: output TSV is empty: $OUTTSV" >&2
    exit 1
fi

echo "========================================"
echo "Final outputs:"
ls -lh "$OUTBED" "$OUTTSV" "$LIST" "$STATS"

if [[ -s "$MISSING" ]]; then
    ls -lh "$MISSING"
fi

echo
echo "Line counts:"
wc -l "$OUTBED" "$OUTTSV"

echo
echo "Preview TSV:"
head "$OUTTSV"

echo
echo "Chunk site counts:"
cat "$STATS"

echo
echo "Job finished: $(date)"
echo "========================================"
