#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/../utilities/common.sh"

require_var SAMTOOLS
require_var F5C
require_var M6ANET
require_var BLOW5
require_var FASTQ
require_var BAM
require_var TRANSCRIPT_FASTA
require_var KMER_MODEL
require_var WORK

THREADS="${THREADS:-5}"
N_CHUNKS="${N_CHUNKS:-24}"
M6ANET_DATAPREP_PROCESSES="${M6ANET_DATAPREP_PROCESSES:-$THREADS}"
CLEAN_EVENTALIGN_AFTER_DATAPREP="${CLEAN_EVENTALIGN_AFTER_DATAPREP:-1}"

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
DATAPREP_DIR="${WORK}/dataprep/${CHUNK_ID}"
STATUS_DIR="${WORK}/status"

mkdir -p "$SPLIT_DIR" "$EVENT_DIR" "$DATAPREP_DIR" "$STATUS_DIR" "$WORK/logs"

CHUNK_BAM="${SPLIT_DIR}/${CHUNK_ID}.bam"

EVENTALIGN_FINAL="${EVENT_DIR}/${CHUNK_ID}.eventalign.tsv"
SUMMARY_FINAL="${EVENT_DIR}/${CHUNK_ID}.eventalign.summary.tsv"

EVENTALIGN_TMP="${EVENTALIGN_FINAL}.tmp"
SUMMARY_TMP="${SUMMARY_FINAL}.tmp"

echo "========================================"
echo "m6Anet eventalign + dataprep chunk"
echo "CHUNK_ID=$CHUNK_ID"
echo "THREADS=$THREADS"
echo "WORK=$WORK"
echo "REF_CHUNK=$REF_CHUNK"
echo "CHUNK_BAM=$CHUNK_BAM"
echo "EVENTALIGN_FINAL=$EVENTALIGN_FINAL"
echo "SUMMARY_FINAL=$SUMMARY_FINAL"
echo "DATAPREP_DIR=$DATAPREP_DIR"
echo "Start: $(date)"
echo "========================================"

[[ -x "$SAMTOOLS" ]] || { echo "ERROR: samtools not executable: $SAMTOOLS" >&2; exit 1; }
[[ -x "$F5C" ]] || { echo "ERROR: f5c not executable: $F5C" >&2; exit 1; }
[[ -x "$M6ANET" ]] || { echo "ERROR: m6anet not executable: $M6ANET" >&2; exit 1; }

for f in \
    "$REF_CHUNK" \
    "$BLOW5" \
    "${BLOW5}.idx" \
    "$FASTQ" \
    "$BAM" \
    "$TRANSCRIPT_FASTA" \
    "$KMER_MODEL"
do
    [[ -s "$f" ]] || { echo "ERROR: missing or empty input: $f" >&2; exit 1; }
done

if [[ ! -s "${BAM}.bai" && ! -s "${BAM%.bam}.bai" ]]; then
    echo "ERROR: BAM index missing for $BAM" >&2
    exit 1
fi

if [[ -s "$STATUS_DIR/${CHUNK_ID}.dataprep.ok" ]] \
    && [[ -d "$DATAPREP_DIR" ]] \
    && find "$DATAPREP_DIR" -type f -size +0c | grep -q .
then
    echo "[SKIP] Existing successful dataprep output found for $CHUNK_ID"
    exit 0
fi

cleanup_chunk_bam() {
    rm -rf "$SPLIT_DIR"
}
trap cleanup_chunk_bam EXIT

echo "[1/4] Create temporary chunk BAM..."

rm -rf "$SPLIT_DIR"
mkdir -p "$SPLIT_DIR"

"$SAMTOOLS" view -@ "$THREADS" -bh -L "$REF_CHUNK" "$BAM" \
    | "$SAMTOOLS" sort -@ "$THREADS" -o "$CHUNK_BAM"

"$SAMTOOLS" index -@ "$THREADS" "$CHUNK_BAM"

N_BAM=$("$SAMTOOLS" view -c "$CHUNK_BAM")
echo "BAM records for $CHUNK_ID: $N_BAM"

if [[ "$N_BAM" -eq 0 ]]; then
    echo "No BAM records for $CHUNK_ID. Mark as empty."
    touch "$STATUS_DIR/${CHUNK_ID}.empty.ok"
    exit 0
fi

echo "[2/4] f5c eventalign..."

rm -f "$EVENTALIGN_TMP" "$SUMMARY_TMP"

"$F5C" eventalign \
    --rna \
    --slow5 "$BLOW5" \
    --kmer-model "$KMER_MODEL" \
    -b "$CHUNK_BAM" \
    -r "$FASTQ" \
    -g "$TRANSCRIPT_FASTA" \
    -o "$EVENTALIGN_TMP" \
    --summary "$SUMMARY_TMP" \
    --scale-events \
    --signal-index \
    -t "$THREADS"

[[ -s "$EVENTALIGN_TMP" ]] || {
    echo "ERROR: eventalign output is empty for $CHUNK_ID" >&2
    exit 1
}

mv "$EVENTALIGN_TMP" "$EVENTALIGN_FINAL"
mv "$SUMMARY_TMP" "$SUMMARY_FINAL"

echo "[3/4] Remove temporary chunk BAM..."

rm -rf "$SPLIT_DIR"
trap - EXIT

echo "[4/4] m6anet dataprep..."

rm -rf "${DATAPREP_DIR}.tmp"
mkdir -p "${DATAPREP_DIR}.tmp"

"$M6ANET" dataprep \
    --eventalign "$EVENTALIGN_FINAL" \
    --out_dir "${DATAPREP_DIR}.tmp" \
    --n_processes "$M6ANET_DATAPREP_PROCESSES"

if ! find "${DATAPREP_DIR}.tmp" -type f -size +0c | grep -q .; then
    echo "ERROR: dataprep output appears empty for $CHUNK_ID" >&2
    exit 1
fi

rm -rf "$DATAPREP_DIR"
mv "${DATAPREP_DIR}.tmp" "$DATAPREP_DIR"

if [[ "$CLEAN_EVENTALIGN_AFTER_DATAPREP" == "1" ]]; then
    echo "[CLEANUP] dataprep succeeded, removing large eventalign TSV..."
    rm -f "$EVENTALIGN_FINAL"
else
    echo "[KEEP] CLEAN_EVENTALIGN_AFTER_DATAPREP=$CLEAN_EVENTALIGN_AFTER_DATAPREP, keeping eventalign TSV."
fi

touch "$STATUS_DIR/${CHUNK_ID}.dataprep.ok"

echo "Finished $CHUNK_ID: $(date)"
echo "Kept files:"
du -sh "$DATAPREP_DIR" "$EVENT_DIR" 2>/dev/null || true
ls -lh "$EVENT_DIR" || true
