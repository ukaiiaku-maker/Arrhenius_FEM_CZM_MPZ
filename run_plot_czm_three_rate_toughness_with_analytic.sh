#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV=${CONDA_ENV:-arrhenius-fem-czm}
if [[ -z "${PYTHON_BIN:-}" ]]; then
  if command -v conda >/dev/null 2>&1; then
    PYTHON_BIN="$(conda run -n "$CONDA_ENV" python -c 'import sys; print(sys.executable)' 2>&1 | tr -d '\r' | awk 'NF {last=$0} END {print last}')"
  else
    PYTHON_BIN=python
  fi
fi

ROOT=${ROOT:-runs/four_class_exp_floor_CZM_rates_no_branch_500um_theta45}
OUT=${OUT:-runs/CZM_three_rate_temperature_comparison}
RATES=${RATES:-"1 10 100"}
TEMPS=${TEMPS:-"300 400 500 600 700 800 900 1000 1100 1200"}
CLASSES=${CLASSES:-"ceramic peak weakT DBTT"}

printf 'python: %s\n' "$PYTHON_BIN"
"$PYTHON_BIN" plot_czm_three_rate_toughness_with_analytic.py \
  --root "$ROOT" \
  --out "$OUT" \
  --rates "$RATES" \
  --temps "$TEMPS" \
  --classes "$CLASSES" \
  --v1-script "${V1_SCRIPT:-run_v1_exp_floor_four_class_tuning.py}" \
  --params-csv "${PARAMS_CSV:-four_class_exp_floor_exact_model_inputs.csv}" \
  --analytic-csv "${ANALYTIC_CSV:-four_class_analytical_prediction_final.csv}" \
  --base-kdot "${BASE_KDOT:-0.005}" \
  --analytic-dK "${ANALYTIC_DK:-0.02}" \
  --analytic-Kmax "${ANALYTIC_KMAX:-100}"
