#!/usr/bin/env bash
set -euo pipefail

INPUT_MANIFEST=${INPUT_MANIFEST:?Set INPUT_MANIFEST to narrow_dbtt_long_growth_promotion_manifest.csv}
OUT=${OUT:-runs/mpz_v9_10_4_2d_validation_manifest_v1.csv}
CANDIDATE_COUNT=${CANDIDATE_COUNT:-3}

python prepare_mpz_v9_10_4_2d_validation.py \
  --input-manifest "$INPUT_MANIFEST" \
  --candidate-count "$CANDIDATE_COUNT" \
  --out "$OUT"
