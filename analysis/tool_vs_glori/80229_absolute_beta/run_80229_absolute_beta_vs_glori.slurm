#!/usr/bin/env bash
#SBATCH --job-name=80229_abs_glori
#SBATCH -p cpu8358
#SBATCH --qos=26cores
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --output=run.%j.out
#SBATCH --error=run.%j.err

export RNA004_STAT_ROOT="${RNA004_STAT_ROOT:-${SLURM_SUBMIT_DIR:-$PWD}/analysis}"


set -euo pipefail

BASE="${RNA004_STAT_ROOT}/tool_vs_glori/80229_absolute_beta"

STAT="${RNA004_STAT_ROOT}/80229"

GLORI="${RNA004_STAT_ROOT}/ground_truth/m6aconquer/processed/GLORI_HEK_minus_HeLa.delta_beta.tsv.gz"

DORADO_HEK="$STAT/dorado_Hek293T.sitelevel.genomic.standardized.collapsed.tsv"
DORADO_HELA="$STAT/dorado_HeLa.sitelevel.genomic.standardized.recollapsed.tsv"

M6ANET_HEK="$STAT/m6anet_Hek293T.sitelevel.genomic.standardized.collapsed.tsv"
M6ANET_HELA="$STAT/m6anet_HeLa.sitelevel.genomic.standardized.collapsed.tsv"

SINGLEMOD_HEK="$STAT/singlemod_Hek293T.sitelevel.genomic.standardized.collapsed.tsv"
SINGLEMOD_HELA="$STAT/singlemod_HeLa.sitelevel.genomic.standardized.collapsed.tsv"

OUT="$BASE/results"


echo "============================================================"
echo "PRJEB80229 absolute beta vs GLORI"
echo "GLORI: $GLORI"
echo
echo "Dorado HEK:    $DORADO_HEK"
echo "Dorado HeLa:   $DORADO_HELA"
echo
echo "m6Anet HEK:    $M6ANET_HEK"
echo "m6Anet HeLa:   $M6ANET_HELA"
echo
echo "SingleMod HEK:  $SINGLEMOD_HEK"
echo "SingleMod HeLa: $SINGLEMOD_HELA"
echo "============================================================"


for file in \
    "$GLORI" \
    "$DORADO_HEK" \
    "$DORADO_HELA" \
    "$M6ANET_HEK" \
    "$M6ANET_HELA" \
    "$SINGLEMOD_HEK" \
    "$SINGLEMOD_HELA"
do
    [[ -s "$file" ]] || {
        echo "ERROR: missing input:" >&2
        echo "$file" >&2
        exit 1
    }
done


mkdir -p "$OUT"

if find "$OUT" -mindepth 1 -print -quit | grep -q .; then
    echo "ERROR: output directory is not empty:" >&2
    echo "$OUT" >&2
    exit 1
fi


python3 \
    "$BASE/compare_80229_absolute_beta_vs_glori.py" \
    --glori "$GLORI" \
    --dorado-hek "$DORADO_HEK" \
    --dorado-hela "$DORADO_HELA" \
    --m6anet-hek "$M6ANET_HEK" \
    --m6anet-hela "$M6ANET_HELA" \
    --singlemod-hek "$SINGLEMOD_HEK" \
    --singlemod-hela "$SINGLEMOD_HELA" \
    --output-dir "$OUT" \
    --chunk-size 500000 \
    --sample-size 300000


echo
echo "Generated files:"

find "$OUT" \
    -maxdepth 1 \
    -type f \
    -printf '%s\t%f\n' \
    | sort -nr
