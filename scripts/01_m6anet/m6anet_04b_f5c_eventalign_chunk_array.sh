#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/../utilities/common.sh"

require_var SAMTOOLS
require_var F5C
require_var BLOW5
require_var KMER_MODEL
require_var BAM
require_var FASTQ
require_var TRANSCRIPT_FASTA
require_var WORK

THREADS="${THREADS:-5}"
N_CHUNKS="${N_CHUNKS:-24}"
CLEAN_CHUNK_BAM="${CLEAN_CHUNK_BAM:-1}"

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    CHUNK_N="${SLURM_ARRAY_TASK_ID}"
else
    require_var CHUNK_N
fi

if (( CHUNK_N < 1 || CHUNK_N > N_CHUNKS )); then
    echo "ERROR: CHUNK_N=$CHUNK_N outside 1-${N_CHUNKS}" >&2
    exit 1
fi

CHUNK_ID=$(printf "chunk_%02d" "$CHUNK_N")

REF_CHUNK="${WORK}/ref_chunks/${CHUNK_ID}.bed"
SPLIT_DIR="${WORK}/split_bam/${CHUNK_ID}"
EVENT_DIR="${WORK}/eventalign/${CHUNK_ID}"
STATUS_DIR="${WORK}/status"

mkdir -p "$SPLIT_DIR" "$EVENT_DIR" "$STATUS_DIR" "$WORK/logs"

CHUNK_BAM="${SPLIT_DIR}/${CHUNK_ID}.bam"
EVENTALIGN="${EVENT_DIR}/${CHUNK_ID}.eventalign.tsv"
SUMMARY="${EVENT_DIR}/${CHUNK_ID}.eventalign.summary.tsv"

echo "========================================"
echo "m6Anet f5c eventalign chunk array"
echo "CHUNK_ID=$CHUNK_ID"
echo "THREADS=$THREADS"
echo "WORK=$WORK"
echo "REF_CHUNK=$REF_CHUNK"
echo "CHUNK_BAM=$CHUNK_BAM"
echo "EVENTALIGN=$EVENTALIGN"
echo "SUMMARY=$SUMMARY"
echo "Start: $(date)"
echo "========================================"

for f in "$BLOW5" "${BLOW5}.idx" "$FASTQ" "$BAM" "$TRANSCRIPT_FASTA" "$KMER_MODEL" "$REF_CHUNK"; do
    [[ -s "$f" ]] || { echo "ERROR: missing or empty input: $f" >&2; exit 1; }
done

if [[ ! -s "${BAM}.bai" && ! -s "${BAM%.bam}.bai" ]]; then
    echo "ERROR: BAM index missing for $BAM" >&2
    exit 1
fi

if [[ -s "$EVENTALIGN" && -s "$SUMMARY" ]]; then
    echo "[SKIP] Existing eventalign and summary found for $CHUNK_ID"
    touch "$STATUS_DIR/${CHUNK_ID}.eventalign.ok"
    exit 0
fi

echo "[1/3] Create chunk BAM from transcript BED..."

rm -f "${CHUNK_BAM}.tmp" "$CHUNK_BAM" "${CHUNK_BAM}.bai"

"$SAMTOOLS" view -@ "$THREADS" -bh -L "$REF_CHUNK" "$BAM" \
    | "$SAMTOOLS" sort -@ "$THREADS" -o "${CHUNK_BAM}.tmp"

mv "${CHUNK_BAM}.tmp" "$CHUNK_BAM"
"$SAMTOOLS" index -@ "$THREADS" "$CHUNK_BAM"

N_BAM=$("$SAMTOOLS" view -c "$CHUNK_BAM")
echo "BAM records for $CHUNK_ID: $N_BAM"

if [[ "$N_BAM" -eq 0 ]]; then
    echo "No BAM records for $CHUNK_ID. Mark as empty."
    touch "$STATUS_DIR/${CHUNK_ID}.empty.ok"
    rm -f "$CHUNK_BAM" "${CHUNK_BAM}.bai"
    rmdir "$SPLIT_DIR" 2>/dev/null || true
    exit 0
fi

echo "[2/3] f5c eventalign for $CHUNK_ID..."

rm -f "${EVENTALIGN}.tmp" "${SUMMARY}.tmp"

export BAM="$CHUNK_BAM"
export EVENTALIGN="${EVENTALIGN}.tmp"
export SUMMARY="${SUMMARY}.tmp"
export THREADS="$THREADS"

bash "$(dirname "$0")/m6anet_04_f5c_eventalign.sh"

[[ -s "${EVENTALIGN}.tmp" ]] || {
    echo "ERROR: eventalign output is empty for $CHUNK_ID" >&2
    exit 1
}

mv "${EVENTALIGN}.tmp" "$EVENTALIGN"
mv "${SUMMARY}.tmp" "$SUMMARY"

echo "[3/3] Clean temporary chunk BAM..."

if [[ "$CLEAN_CHUNK_BAM" == "1" ]]; then
    rm -f "$CHUNK_BAM" "${CHUNK_BAM}.bai" "${CHUNK_BAM}.csi"
    rmdir "$SPLIT_DIR" 2>/dev/null || true
fi

touch "$STATUS_DIR/${CHUNK_ID}.eventalign.ok"

echo "Finished $CHUNK_ID: $(date)"
echo "Output:"
ls -lh "$EVENTALIGN" "$SUMMARY"
