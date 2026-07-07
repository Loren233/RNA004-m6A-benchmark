#!/usr/bin/env bash
set -euo pipefail

source "$(dirname "$0")/../utilities/common.sh"

require_var CALL_DIR

if [[ -n "${DEEPRM_SITE_DIR:-}" ]]; then
    SITE_DIR="$DEEPRM_SITE_DIR"
elif [[ -d "${CALL_DIR}/prep/site-level" ]]; then
    SITE_DIR="${CALL_DIR}/prep/site-level"
elif [[ -d "${CALL_DIR}/site-level" ]]; then
    SITE_DIR="${CALL_DIR}/site-level"
else
    echo "ERROR: Cannot find DeepRM site-level directory under CALL_DIR=${CALL_DIR}" >&2
    echo "Expected one of:" >&2
    echo "  ${CALL_DIR}/prep/site-level" >&2
    echo "  ${CALL_DIR}/site-level" >&2
    exit 1
fi

PILEUP_BED="${DEEPRM_PILEUP_BED:-${SITE_DIR}/pileup.bed}"
OUT_TSV="${DEEPRM_SITELEVEL_TSV:-${CALL_DIR}/deeprm_sitelevel.tsv}"

require_file "$PILEUP_BED"

mkdir -p "$(dirname "$OUT_TSV")"

echo "[INFO] CALL_DIR=${CALL_DIR}" >&2
echo "[INFO] SITE_DIR=${SITE_DIR}" >&2
echo "[INFO] PILEUP_BED=${PILEUP_BED}" >&2
echo "[INFO] OUT_TSV=${OUT_TSV}" >&2

echo "[INFO] pileup.bed preview:" >&2
head -n 3 "$PILEUP_BED" >&2 || true

TMP_OUT="${OUT_TSV}.tmp"

LC_ALL=C awk '
BEGIN {
    FS = OFS = "\t";
    print "transcript_id", "pos_1based", "strand", "depth", "mod_count", "stoich", \
          "percent_modified", "modified_base", "score", \
          "bed_start_0based", "bed_end", \
          "canonical_count", "other_mod_count", "delete_count", \
          "fail_count", "diff_count", "nocall_count";
}

# skip comments / browser / track lines
/^#/ { next }
/^track/ { next }
/^browser/ { next }

NF < 11 {
    bad++;
    next;
}

{
    transcript_id = $1;
    start0 = $2;
    endpos = $3;
    modified_base = $4;
    score = $5;
    strand = $6;

    # bedMethyl / modkit convention:
    # col10 = coverage
    # col11 = percent modified
    depth = $10;
    percent_modified = $11;

    mod_count = (NF >= 12 ? $12 : "NA");
    canonical_count = (NF >= 13 ? $13 : "NA");
    other_mod_count = (NF >= 14 ? $14 : "NA");
    delete_count = (NF >= 15 ? $15 : "NA");
    fail_count = (NF >= 16 ? $16 : "NA");
    diff_count = (NF >= 17 ? $17 : "NA");
    nocall_count = (NF >= 18 ? $18 : "NA");

    pos_1based = start0 + 1;

    gsub("%", "", percent_modified);

    # Use mod_count / depth as the primary stoichiometry estimate.
    # Some pileup.bed rows can have non-numeric percent_modified, which would
    # otherwise create NA stoich and break downstream standardization.
    if (depth ~ /^[0-9.]+$/ && mod_count ~ /^[0-9.]+$/ && depth + 0 > 0) {
        stoich = sprintf("%.6f", (mod_count + 0) / (depth + 0));
    } else if (percent_modified ~ /^[0-9.]+$/) {
        stoich = sprintf("%.6f", percent_modified / 100.0);
    } else {
        stoich = "0.000000";
    }

    print transcript_id, pos_1based, strand, depth, mod_count, stoich, \
          percent_modified, modified_base, score, \
          start0, endpos, \
          canonical_count, other_mod_count, delete_count, \
          fail_count, diff_count, nocall_count;
}

END {
    if (bad > 0) {
        print "[WARN] skipped lines with fewer than 11 columns: " bad > "/dev/stderr";
    }
}
' "$PILEUP_BED" > "$TMP_OUT"

mv "$TMP_OUT" "$OUT_TSV"

echo "[DONE] DeepRM site-level table created: ${OUT_TSV}" >&2
echo "[INFO] Output preview:" >&2
head -n 5 "$OUT_TSV" >&2
echo "[INFO] Number of output lines:" >&2
wc -l "$OUT_TSV" >&2
