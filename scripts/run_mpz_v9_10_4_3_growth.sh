#!/usr/bin/env bash
set -euo pipefail

INPUT_MANIFEST=${INPUT_MANIFEST:?Set INPUT_MANIFEST to the promoted candidate manifest}
STAGE=${STAGE:?Set STAGE to short or long}
OUT=${OUT:?Set OUT}
COARSE_TEMPS=${COARSE_TEMPS:-"300 400 500 600 700 800 900 1000 1100"}
CLEAVAGE_SLOPE_MODE=${CLEAVAGE_SLOPE_MODE:-fixed_zero}
MAX_CANDIDATES=${MAX_CANDIDATES:-12}
LOCAL_MAXITER=${LOCAL_MAXITER:-100}
KDOT=${KDOT:-0.005}
KMAX=${KMAX:-80}
MAX_DK_SUBSTEP=${MAX_DK_SUBSTEP:-0.05}
MAX_K_SHIELD=${MAX_K_SHIELD:-1.0}
TARGET_EXT_UM=${TARGET_EXT_UM:-}

mkdir -p "$OUT"

ARGS=(
  --input-manifest "$INPUT_MANIFEST"
  --stage "$STAGE"
  --coarse-temperatures "$COARSE_TEMPS"
  --cleavage-slope-mode "$CLEAVAGE_SLOPE_MODE"
  --max-candidates "$MAX_CANDIDATES"
  --local-maxiter "$LOCAL_MAXITER"
  --Kdot "$KDOT"
  --Kmax "$KMAX"
  --max-dK-substep "$MAX_DK_SUBSTEP"
  --max-K-shield "$MAX_K_SHIELD"
  --out "$OUT"
)
if [[ -n "$TARGET_EXT_UM" ]]; then
  ARGS+=(--target-extension-um "$TARGET_EXT_UM")
fi

PYTHONUNBUFFERED=1 python refine_mpz_v9_10_4_3_growth.py "${ARGS[@]}"
