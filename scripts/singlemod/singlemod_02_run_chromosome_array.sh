#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/../utilities/common.sh"

require_var PYTHON
require_var SAMTOOLS
require_var BEDTOOLS
require_var F5C
require_var SINGLEMOD_REPO
require_var FASTQ
require_var BLOW5
require_var BAM
require_var REF
require_var WORK

THREADS="${THREADS:-8}"
S="${S:-500000}"

CHRS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    CHR="${CHRS[$((SLURM_ARRAY_TASK_ID - 1))]}"
else
    require_var CHR
fi

PREFIX="full_${CHR}"

SPLIT="${WORK}/split_bam/${CHR}"
EVENT="${WORK}/eventalign/${CHR}"
TMP="${WORK}/tmp_features/${CHR}"
FEAT="${WORK}/features/${CHR}"
PRED="${WORK}/prediction/${CHR}"
SITE="${WORK}/sitelevel/${CHR}"
STATUS="${WORK}/status"

mkdir -p "$SPLIT" "$EVENT" "$TMP" "$FEAT" "$PRED" "$SITE" "$STATUS" "$WORK/logs"

CHR_BAM="${SPLIT}/${PREFIX}.bam"
CHR_BED="${SPLIT}/${PREFIX}.bed"
EVENTALIGN="${EVENT}/${PREFIX}_eventalign.txt"
SUMMARY="${EVENT}/${PREFIX}_summary.txt"

echo "========================================"
echo "SingleMod chromosome run"
echo "CHR=$CHR"
echo "THREADS=$THREADS"
echo "WORK=$WORK"
echo "Start: $(date)"
echo "========================================"

for f in "$FASTQ" "$BLOW5" "$REF" "$BAM" "${BAM}.bai"; do
    [[ -s "$f" ]] || { echo "ERROR: missing or empty input: $f" >&2; exit 1; }
done

for f in \
    "$SINGLEMOD_REPO/SingleMod/organize_from_eventalign.py" \
    "$SINGLEMOD_REPO/SingleMod/merge_motif_npy.py" \
    "$SINGLEMOD_REPO/SingleMod/SingleMod_m6A_prediction.py"
do
    [[ -s "$f" ]] || { echo "ERROR: SingleMod script missing: $f" >&2; exit 1; }
done

[[ -d "$SINGLEMOD_REPO/models/RNA004" ]] || {
    echo "ERROR: RNA004 model directory missing: $SINGLEMOD_REPO/models/RNA004" >&2
    exit 1
}

echo "[1/8] Create chromosome BAM..."
if [[ ! -s "$CHR_BAM" ]]; then
    "$SAMTOOLS" view -@ "$THREADS" -bh "$BAM" "$CHR" \
        | "$SAMTOOLS" sort -@ "$THREADS" -o "$CHR_BAM"
    "$SAMTOOLS" index -@ "$THREADS" "$CHR_BAM"
fi

N_BAM=$("$SAMTOOLS" view -c "$CHR_BAM")
echo "BAM records for $CHR: $N_BAM"

if [[ "$N_BAM" -eq 0 ]]; then
    echo "No BAM records for $CHR. Mark as empty."
    touch "$STATUS/${CHR}.empty.ok"
    exit 0
fi

echo "[2/8] BAM to BED..."
"$BEDTOOLS" bamtobed -i "$CHR_BAM" > "$CHR_BED"
[[ -s "$CHR_BED" ]] || { echo "ERROR: BED is empty for $CHR" >&2; exit 1; }

echo "[3/8] f5c index..."
"$F5C" index \
    -t "$THREADS" \
    --slow5 "$BLOW5" \
    "$FASTQ"

echo "[4/8] f5c eventalign..."
rm -f "${EVENTALIGN}.tmp" "${SUMMARY}.tmp"

"$F5C" eventalign \
    -r "$FASTQ" \
    -b "$CHR_BAM" \
    -g "$REF" \
    -t "$THREADS" \
    --slow5 "$BLOW5" \
    --pore rna004 \
    --rna \
    --scale-events \
    --samples \
    --signal-index \
    --summary "${SUMMARY}.tmp" \
    --print-read-names \
    > "${EVENTALIGN}.tmp"

[[ -s "${EVENTALIGN}.tmp" ]] || { echo "ERROR: eventalign output is empty for $CHR" >&2; exit 1; }

mv "${EVENTALIGN}.tmp" "$EVENTALIGN"
mv "${SUMMARY}.tmp" "$SUMMARY"

echo "[5/8] organize_from_eventalign..."
rm -rf "$TMP"
mkdir -p "$TMP"

"$PYTHON" "$SINGLEMOD_REPO/SingleMod/organize_from_eventalign.py" \
    -v 004 \
    -b "$CHR_BED" \
    -e "$EVENTALIGN" \
    -o "$TMP" \
    -p "$PREFIX" \
    -s "$S"

rm -f "$EVENTALIGN"

echo "[6/8] Filter empty motifs..."
cd "$TMP"

rm -f extra_info.txt
mkdir -p zero_empty_motifs

VALID=0
EMPTY=0

shopt -s nullglob
for info in ${PREFIX}-*-extra_info.txt; do
    n=$(wc -l < "$info")
    motif=${info#${PREFIX}-}
    motif=${motif%-extra_info.txt}

    if [[ "$n" -gt 0 ]]; then
        printf "%s\t%s\n" "$n" "$info" >> extra_info.txt
        VALID=$((VALID + 1))
    else
        EMPTY=$((EMPTY + 1))
        mv -f "$info" zero_empty_motifs/ || true
        for suffix in sequence signal extra; do
            f="${PREFIX}-${motif}-${suffix}.npy"
            [[ -e "$f" ]] && mv -f "$f" zero_empty_motifs/ || true
        done
    fi
done

echo "Valid motifs: $VALID"
echo "Empty motifs moved: $EMPTY"

if [[ "$VALID" -eq 0 ]]; then
    echo "No valid motifs for $CHR. Mark as no_valid_motifs."
    touch "$STATUS/${CHR}.no_valid_motifs"
    exit 0
fi

cd "$WORK"

echo "[7/8] merge_motif_npy..."
rm -rf "$FEAT"
mkdir -p "$FEAT"

"$PYTHON" "$SINGLEMOD_REPO/SingleMod/merge_motif_npy.py" \
    -v 004 \
    -d "$TMP" \
    -s "$S" \
    -o "$FEAT"

N_FEAT=$(find "$FEAT" -type f | wc -l)
echo "Feature files: $N_FEAT"

[[ "$N_FEAT" -gt 0 ]] || { echo "ERROR: no feature files generated for $CHR" >&2; exit 1; }

echo "[8/8] SingleMod RNA004 prediction..."
rm -rf "$PRED"
mkdir -p "$PRED"

MOTIFS=(
AAACA AAACC AAACG AAACT AAATA AAATT
AGACA AGACC AGACG AGACT AGATT ATACT
CAACT CGACT CTACT
GAACA GAACC GAACT GAATA GAATG GAATT
GGACA GGACC GGACG GGACT GGATA GGATC GGATG GGATT
GTACT TAACA TAACT TGACA TGACC TGACT TTACT
)

PRED_COUNT=0

for motif in "${MOTIFS[@]}"; do
    MODEL=$(find "$SINGLEMOD_REPO/models/RNA004" "$SINGLEMOD_REPO/models" \
        -type f \( -name "model_${motif}.pth.tar" -o -name "*${motif}*.pth.tar" -o -name "*${motif}*.pt" \) \
        2>/dev/null | head -n 1 || true)

    if [[ -z "$MODEL" || ! -s "$MODEL" ]]; then
        echo "WARNING: model not found for motif $motif, skip"
        continue
    fi

    if ! find "$FEAT" -type f -name "*${motif}*" | grep -q .; then
        echo "No feature file found for motif $motif, skip"
        continue
    fi

    echo "Predicting $CHR motif $motif"

    "$PYTHON" "$SINGLEMOD_REPO/SingleMod/SingleMod_m6A_prediction.py" \
        -v 004 \
        -d "$FEAT" \
        -k "$motif" \
        -m "$MODEL" \
        -b 3000 \
        -o "$PRED/${motif}_prediction.txt"

    [[ -s "$PRED/${motif}_prediction.txt" ]] && PRED_COUNT=$((PRED_COUNT + 1))
done

echo "Prediction files generated: $PRED_COUNT"
[[ "$PRED_COUNT" -gt 0 ]] || { echo "ERROR: no prediction files generated for $CHR" >&2; exit 1; }

echo "Generate chromosome-level prediction.txt and mr.bed..."
cat "$PRED"/*_prediction.txt > "$SITE/${CHR}.prediction.txt"

awk 'BEGIN{OFS=FS="\t"}{
    split($1,info,"|");
    s=info[1]"|"info[2]"|"info[3]"|"info[5];
    t[s]=t[s]+1;
    if($2 > 0.5){m[s]=m[s]+1}
}
END{
    for(i in t){
        split(i,info,"|");
        if(i in m){
            print info[1],info[2]-1,info[2],i,m[i]/t[i],info[3],t[i],m[i],info[4]
        }else{
            print info[1],info[2]-1,info[2],i,0,info[3],t[i],0,info[4]
        }
    }
}' "$SITE/${CHR}.prediction.txt" | sort -k1,1 -k2,2n > "$SITE/${CHR}.mr.bed"

[[ -s "$SITE/${CHR}.mr.bed" ]] || { echo "ERROR: mr.bed is empty for $CHR" >&2; exit 1; }

touch "$STATUS/${CHR}.sitelevel.ok"

rm -rf "$TMP" "$FEAT"

echo "Finished $CHR: $(date)"
echo "Output:"
ls -lh "$SITE/${CHR}.prediction.txt" "$SITE/${CHR}.mr.bed"
