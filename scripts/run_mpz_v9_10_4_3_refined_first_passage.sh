#!/usr/bin/env bash
set -euo pipefail

INPUT_MANIFEST=${INPUT_MANIFEST:?Set INPUT_MANIFEST to the dynamic refinement manifest}
OUT=${OUT:-runs/mpz_v9_10_4_3_refined_first_passage_v1}
COARSE_TEMPS=${COARSE_TEMPS:-"300 400 500 600 700 800 900 1000 1100"}
CLEAVAGE_SLOPE_MODE=${CLEAVAGE_SLOPE_MODE:-fixed_zero}
MAX_CANDIDATES=${MAX_CANDIDATES:-24}
LOCAL_MAXITER=${LOCAL_MAXITER:-150}
KDOT=${KDOT:-0.005}
KMAX=${KMAX:-80}
TARGET_EXT_UM=${TARGET_EXT_UM:-5}
MAX_DK_SUBSTEP=${MAX_DK_SUBSTEP:-0.05}
MAX_K_SHIELD=${MAX_K_SHIELD:-1.0}

mkdir -p "$OUT"

PYTHONUNBUFFERED=1 python refine_mpz_v9_10_4_3_first_passage.py \
  --input-manifest "$INPUT_MANIFEST" \
  --coarse-temperatures "$COARSE_TEMPS" \
  --cleavage-slope-mode "$CLEAVAGE_SLOPE_MODE" \
  --max-candidates "$MAX_CANDIDATES" \
  --local-maxiter "$LOCAL_MAXITER" \
  --Kdot "$KDOT" \
  --Kmax "$KMAX" \
  --target-extension-um "$TARGET_EXT_UM" \
  --max-dK-substep "$MAX_DK_SUBSTEP" \
  --max-K-shield "$MAX_K_SHIELD" \
  --out "$OUT"
