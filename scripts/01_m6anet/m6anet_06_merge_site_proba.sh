#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/../utilities/common.sh"

require_var INFERENCE_DIR
require_var OUTPUT

SITE_PROBA_NAME="${SITE_PROBA_NAME:-data.site_proba.csv}"
EXPECTED_FILES="${EXPECTED_FILES:-}"

mkdir -p "$(dirname "$OUTPUT")"

echo "========================================"
echo "m6Anet merge site_proba"
echo "INFERENCE_DIR=$INFERENCE_DIR"
echo "SITE_PROBA_NAME=$SITE_PROBA_NAME"
echo "OUTPUT=$OUTPUT"
echo "EXPECTED_FILES=${EXPECTED_FILES:-not_set}"
echo "Start: $(date)"
echo "========================================"

mapfile -t FILES < <(
    find "$INFERENCE_DIR" -maxdepth 2 -type f -name "$SITE_PROBA_NAME" | sort
)

N_FILES="${#FILES[@]}"
echo "Found $N_FILES site_proba files."

if [[ "$N_FILES" -eq 0 ]]; then
    echo "ERROR: no $SITE_PROBA_NAME files found under $INFERENCE_DIR" >&2
    exit 1
fi

if [[ -n "$EXPECTED_FILES" && "$N_FILES" -ne "$EXPECTED_FILES" ]]; then
    echo "ERROR: expected $EXPECTED_FILES files, found $N_FILES" >&2
    printf '%s\n' "${FILES[@]}" >&2
    exit 1
fi

FIRST="${FILES[0]}"
HEADER_FIRST="$(head -n 1 "$FIRST")"

TMP="${OUTPUT}.tmp"
rm -f "$TMP"

echo "$HEADER_FIRST" > "$TMP"

for f in "${FILES[@]}"; do
    echo "[MERGE] $f"

    HEADER="$(head -n 1 "$f")"
    if [[ "$HEADER" != "$HEADER_FIRST" ]]; then
        echo "ERROR: header mismatch in $f" >&2
        echo "Expected: $HEADER_FIRST" >&2
        echo "Observed: $HEADER" >&2
        exit 1
    fi

    tail -n +2 "$f" >> "$TMP"
done

mv "$TMP" "$OUTPUT"

echo "Finished: $(date)"
echo "Output:"
ls -lh "$OUTPUT"
wc -l "$OUTPUT"
