#!/usr/bin/env bash
#SBATCH --job-name=80229_vs_glori
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

BASE="${RNA004_STAT_ROOT}/tool_vs_glori/80229"

DELTA_DIR="${RNA004_STAT_ROOT}/80229/delta_beta"

GLORI="${RNA004_STAT_ROOT}/ground_truth/m6aconquer/processed/GLORI_HEK_minus_HeLa.delta_beta.tsv.gz"

OUT="$BASE/results"


find_delta_file() {

    local tool="$1"

    find "$DELTA_DIR" \
        -maxdepth 1 \
        -type f \
        \( \
            -iname "*${tool}*delta_beta*.tsv" \
            -o \
            -iname "*${tool}*delta_beta*.tsv.gz" \
        \) \
        ! -iname "*summary*" \
        ! -iname "*merged*" \
        ! -iname "*matched*" \
        | sort \
        | head -1
}


DORADO="$(find_delta_file dorado)"
M6ANET="$(find_delta_file m6anet)"
SINGLEMOD="$(find_delta_file singlemod)"


echo "============================================================"
echo "PRJEB80229 vs GLORI"
echo "GLORI:     $GLORI"
echo "Dorado:    $DORADO"
echo "m6Anet:    $M6ANET"
echo "SingleMod: $SINGLEMOD"
echo "Output:    $OUT"
echo "============================================================"


for file in \
    "$GLORI" \
    "$DORADO" \
    "$M6ANET" \
    "$SINGLEMOD"
do
    [[ -s "$file" ]] || {
        echo "ERROR: missing input: $file" >&2
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
    "$BASE/compare_80229_tools_vs_glori.py" \
    --glori "$GLORI" \
    --dorado "$DORADO" \
    --m6anet "$M6ANET" \
    --singlemod "$SINGLEMOD" \
    --output-dir "$OUT" \
    --sample-size 300000


echo
echo "Generated files:"

find "$OUT" \
    -maxdepth 1 \
    -type f \
    -printf '%s\t%f\n' \
    | sort -nr
