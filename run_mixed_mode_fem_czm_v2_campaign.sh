#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"
CONDA_ENV="${CONDA_ENV:-arrhenius-fem-czm}"
OUTROOT="${OUTROOT:-runs/mixed_mode_fem_czm_v2_event_controlled_500K}"
CALROOT="${CALROOT:-runs/mixed_mode_fem_czm_v2_elastic_calibration}"
TARGET_PSI="${TARGET_PSI:--60 -45 -30 -15 0 15 30 45 60}"
CLASSES="${CLASSES:-ceramic DBTT}"
THETA_DEG="${THETA_DEG:-45}"

conda run -n "$CONDA_ENV" python calibrate_mixed_mode_loading_v2.py \
  --out "$CALROOT" --target-psi-deg="$TARGET_PSI" \
  --crystal-aniso --crystal-theta-deg "$THETA_DEG"

conda run -n "$CONDA_ENV" python run_mixed_mode_fem_czm_v2_campaign.py \
  --parameter-table "${PARAMETER_TABLE:-four_class_exp_floor_exact_model_inputs.csv}" \
  --calibration-csv "$CALROOT/mixed_mode_loading_calibration_v2.csv" \
  --classes "$CLASSES" --target-psi-deg="$TARGET_PSI" \
  --T-K "${T_K:-500}" --theta-deg "$THETA_DEG" \
  --psi-tol-deg "${PSI_TOL_DEG:-2}" --max-control-iters "${MAX_CONTROL_ITERS:-5}" \
  --max-jobs "${MAX_JOBS:-1}" --outroot "$OUTROOT"

conda run -n "$CONDA_ENV" python plot_mixed_mode_fem_czm_v2_results.py \
  --root "$OUTROOT"
