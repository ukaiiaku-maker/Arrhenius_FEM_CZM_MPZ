#!/usr/bin/env bash
set -euo pipefail

INPUT_MANIFEST=${INPUT_MANIFEST:?Set INPUT_MANIFEST to narrow_dbtt_short_growth_promotion_manifest.csv}
OUT=${OUT:-runs/mpz_v9_10_4_narrow_dbtt_long_growth_500um_v1}
CLEAVAGE_SLOPE_MODE=${CLEAVAGE_SLOPE_MODE:-fixed_zero}
MAX_CANDIDATES=${MAX_CANDIDATES:-8}
LOCAL_MAXITER=${LOCAL_MAXITER:-60}

python refine_mpz_v9_10_4_growth.py \
  --input-manifest "$INPUT_MANIFEST" \
  --stage long \
  --target-extension-um 500 \
  --cleavage-slope-mode "$CLEAVAGE_SLOPE_MODE" \
  --max-candidates "$MAX_CANDIDATES" \
  --local-maxiter "$LOCAL_MAXITER" \
  --out "$OUT"
