#!/usr/bin/env bash
set -euo pipefail

cat >&2 <<'MESSAGE'
TO VERIFY before enabling DeepRM transcript-to-genome conversion:
- transcript ID column
- transcript coordinate column
- whether the coordinate is 0-based or 1-based
- depth column
- stoichiometry column

Once the DeepRM site-level transcript table schema is verified, call:
scripts/09_standardization/convert_transcript_to_genome.py
with the verified column mappings.
MESSAGE

exit 2
