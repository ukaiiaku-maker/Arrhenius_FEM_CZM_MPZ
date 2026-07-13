#!/usr/bin/env bash
set -euo pipefail
PYTHON_BIN=${PYTHON_BIN:-python}
ROOT=${ROOT:-runs/czm_Rcurve_300K_orientation_branching_ceramic_v1}
OUT=${OUT:-${ROOT}/plots}
"$PYTHON_BIN" plot_czm_300K_orientation_branching_Rcurves.py --root "$ROOT" --out "$OUT" --bin-um "${BIN_UM:-25}" --max-ext-um "${MAX_EXT_UM:-1000}"
