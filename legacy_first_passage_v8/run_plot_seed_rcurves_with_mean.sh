#!/usr/bin/env bash
set -euo pipefail
PYTHON_BIN=${PYTHON_BIN:-python}
ROOT=${ROOT:-runs/four_class_exp_floor_CZM_500K_5rep_1000um_theta45}
OUT=${OUT:-runs/four_class_exp_floor_CZM_500K_5rep_1000um_theta45/seed_rcurve_plots}
"$PYTHON_BIN" plot_seed_rcurves_with_mean.py --root "$ROOT" --out "$OUT"
