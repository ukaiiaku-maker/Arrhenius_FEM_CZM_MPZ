#!/usr/bin/env bash
set -euo pipefail

PF_ROOT=${PF_ROOT:-PF-four_class_exp_floor_PF_sharp_no_branch_500um_theta45}
CZM_ROOT=${CZM_ROOT:-four_class_exp_floor_CZM_rates_no_branch_500um_theta45/rate_1x}
ANALYTIC_CSV=${ANALYTIC_CSV:-four_class_analytical_prediction_final.csv}
OUT=${OUT:-PF_vs_CZM_first_passage_with_analytic}
PYTHON_BIN=${PYTHON_BIN:-python}

"$PYTHON_BIN" compare_pf_czm_first_passage_with_analytic.py \
  --pf-root "$PF_ROOT" \
  --czm-root "$CZM_ROOT" \
  --analytic-csv "$ANALYTIC_CSV" \
  --out "$OUT"
