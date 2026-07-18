#!/usr/bin/env bash
set -euo pipefail

MANIFEST=${MANIFEST:?set MANIFEST to an analytical promotion manifest}
OUT=${OUT:-runs/mpz_v9_10_4_8_dynamic_1d_v1}
TARGET_EXT_UM=${TARGET_EXT_UM:-5}
PER_BRACKET_KEEP=${PER_BRACKET_KEEP:-2}

mkdir -p "$OUT"

echo "========================================================================"
echo "v9.10.4.8 corrected four-temperature 1-D evaluation"
echo "manifest=$MANIFEST"
echo "target_extension_um=$TARGET_EXT_UM"
echo "per_bracket_keep=$PER_BRACKET_KEEP"
echo "out=$OUT"
echo "========================================================================"

python evaluate_dynamic_1d_mpz_v9_10_4_8.py \
  --manifest "$MANIFEST" \
  --out "$OUT" \
  --target-extension-um "$TARGET_EXT_UM" \
  --per-bracket-keep "$PER_BRACKET_KEEP"
