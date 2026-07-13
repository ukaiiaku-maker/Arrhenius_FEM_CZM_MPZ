#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

CONDA_ENV="${CONDA_ENV:-arrhenius-fem-czm}"
OUTROOT="${OUTROOT:-runs/mixed_mode_fem_czm_v3_5_matrix_authoritative_500K}"
CALROOT="${CALROOT:-runs/mixed_mode_fem_czm_v3_5_matrix_authoritative_calibration}"
TARGET_PSI="${TARGET_PSI:--60 -45 -30 -15 0 15 30 45 60}"
CLASSES="${CLASSES:-ceramic DBTT}"

CALCSV="$CALROOT/mixed_mode_loading_calibration_v3_5.csv"
if [[ "${RECALIBRATE:-0}" == "1" || ! -f "$CALCSV" ]]; then
  conda run -n "$CONDA_ENV" python calibrate_mixed_mode_loading_v3_5.py \
    --out "$CALROOT" \
    --target-psi-deg="$TARGET_PSI" \
    --psi-tol-deg "${PSI_TOL_DEG:-0.75}" \
    --phase-spread-tol-deg "${PHASE_SPREAD_TOL_DEG:-20}"
else
  echo "Using existing calibration: $CALCSV"
fi

conda run -n "$CONDA_ENV" python run_mixed_mode_fem_czm_v3_5_campaign.py \
  --parameter-table "${PARAMETER_TABLE:-four_class_exp_floor_exact_model_inputs.csv}" \
  --calibration-csv "$CALCSV" \
  --classes "$CLASSES" \
  --target-psi-deg="$TARGET_PSI" \
  --T-K "${T_K:-500}" \
  --max-jobs "${MAX_JOBS:-1}" \
  --outroot "$OUTROOT"

conda run -n "$CONDA_ENV" python plot_mixed_mode_fem_czm_v3_5_results.py --root "$OUTROOT"
