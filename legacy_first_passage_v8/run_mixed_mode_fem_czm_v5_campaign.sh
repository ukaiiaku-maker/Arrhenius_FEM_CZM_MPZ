#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"
ENV="${CONDA_ENV:-arrhenius-fem-czm}"
CALROOT="${CALROOT:-runs/mixed_mode_fem_czm_v5_anisotropic_calibration}"
OUTROOT="${OUTROOT:-runs/mixed_mode_fem_czm_v5_anisotropic_calibrated_tip_500K}"
PSI="${TARGET_PSI:--60 -45 -30 -15 0 15 30 45 60}"
THETA="${CRYSTAL_THETA_DEG:-45}"
if [[ "${RECALIBRATE:-0}" == "1" || ! -f "$CALROOT/mixed_mode_loading_calibration_v5.csv" ]]; then
  conda run -n "$ENV" python calibrate_mixed_mode_loading_v5.py \
    --out "$CALROOT" \
    --target-psi-deg="$PSI" \
    --crystal-theta-deg "$THETA" \
    --cleave-gamma-aniso "${CLEAVE_GAMMA_ANISO:-0.3}" \
    --traction-probe-radius-m "${TRACTION_PROBE_RADIUS_M:-1e-5}" \
    --shear-emission-weight "${SHEAR_EMISSION_WEIGHT:-1.0}" \
    --directional-factor-max "${DIRECTIONAL_FACTOR_MAX:-5.0}"
else
  echo "Using existing v5 anisotropic calibration: $CALROOT/mixed_mode_loading_calibration_v5.csv"
fi
conda run -n "$ENV" python run_mixed_mode_fem_czm_v5_campaign.py \
  --parameter-table "${PARAMETER_TABLE:-four_class_exp_floor_exact_model_inputs.csv}" \
  --calibration-csv "$CALROOT/mixed_mode_loading_calibration_v5.csv" \
  --classes "${CLASSES:-ceramic DBTT}" \
  --target-psi-deg="$PSI" \
  --T-K "${T_K:-500}" \
  --outroot "$OUTROOT" \
  --conda-env "$ENV" \
  --max-jobs "${MAX_JOBS:-1}" \
  --crystal-theta-deg "$THETA" \
  --cleave-gamma-aniso "${CLEAVE_GAMMA_ANISO:-0.3}" \
  --traction-probe-radius-m "${TRACTION_PROBE_RADIUS_M:-1e-5}" \
  --shear-emission-weight "${SHEAR_EMISSION_WEIGHT:-1.0}" \
  --directional-factor-max "${DIRECTIONAL_FACTOR_MAX:-5.0}" \
  --event-psi-tol-deg "${EVENT_PSI_TOL_DEG:-2.0}" \
  --max-control-iters "${MAX_CONTROL_ITERS:-4}" \
  --print-every "${PRINT_EVERY:-50}"
conda run -n "$ENV" python plot_mixed_mode_fem_czm_v5_results.py --root "$OUTROOT"
