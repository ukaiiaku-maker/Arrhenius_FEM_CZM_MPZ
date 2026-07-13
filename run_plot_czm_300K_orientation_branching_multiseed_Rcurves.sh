#!/usr/bin/env bash
set -euo pipefail
PYTHON_BIN=${PYTHON_BIN:-python}
ROOT=${ROOT:-runs/czm_Rcurve_300K_orientation_branching_ceramic_multiseed_v3}
OUT=${OUT:-${ROOT}/plots}
BIN_UM=${BIN_UM:-25}
MAX_EXT_UM=${MAX_EXT_UM:-1000}
TARGET_EXT_UM=${TARGET_EXT_UM:-1000}
"$PYTHON_BIN" plot_czm_300K_orientation_branching_multiseed_Rcurves.py \
  --root "$ROOT" \
  --out "$OUT" \
  --bin-um "$BIN_UM" \
  --max-ext-um "$MAX_EXT_UM" \
  --target-ext-um "$TARGET_EXT_UM"
