#!/usr/bin/env bash
set -euo pipefail

PYTHON_BIN=${PYTHON_BIN:-python}
ROOT=${ROOT:-runs/four_class_exp_floor_CZM_500K_5rep_1000um_theta45}
OUT=${OUT:-runs/four_class_exp_floor_CZM_500K_5rep_1000um_theta45/Rcurve_analysis}

"$PYTHON_BIN" analyze_saved_Rcurves.py \
  --root "$ROOT" \
  --out "$OUT" \
  --bin-um "${BIN_UM:-25}" \
  --max-ext-um "${MAX_EXT_UM:-1000}" \
  --target-ext-um "${TARGET_EXT_UM:-1000}" \
  --late-window "${LATE_WINDOW:-700 1000}"
