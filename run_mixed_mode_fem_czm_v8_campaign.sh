#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"
ENV="${CONDA_ENV:-arrhenius-fem-czm}"
CALROOT="${CALROOT:-runs/mixed_mode_fem_czm_v8_full_circle_calibration}"
OUTROOT="${OUTROOT:-runs/mixed_mode_fem_czm_v8_full_circle_500K}"
PSI="${TARGET_PSI:--30 0 30}"
THETA="${CRYSTAL_THETA_DEG:-45}"

if [[ "${RECALIBRATE:-0}" == "1" || ! -f "$CALROOT/mixed_mode_loading_calibration_v8.csv" ]]; then
  conda run -n "$ENV" python calibrate_mixed_mode_loading_v8.py \
    --out "$CALROOT" \
    --target-psi-deg="$PSI" \
    --T-K "${T_K:-500}" \
    --nx "${NX:-24}" --ny "${NY:-48}" \
    --tip-h-fine "${TIP_H_FINE:-3e-6}" --tip-ratio "${TIP_RATIO:-1.25}" \
    --crystal-theta-deg "$THETA" \
    --cleave-gamma-aniso "${CLEAVE_GAMMA_ANISO:-0.3}" \
    --traction-probe-radius-m "${TRACTION_PROBE_RADIUS_M:-1e-5}" \
    --shear-emission-weight "${SHEAR_EMISSION_WEIGHT:-1.0}" \
    --directional-factor-max "${DIRECTIONAL_FACTOR_MAX:-5.0}" \
    --psi-tol-deg "${CAL_PSI_TOL_DEG:-0.75}" \
    --max-root-iters "${MAX_CAL_ROOT_ITERS:-20}" \
    --max-alpha-step-deg "${MAX_CAL_ALPHA_STEP_DEG:-20}"
else
  echo "Using existing v8 exact-backend full-circle calibration: $CALROOT/mixed_mode_loading_calibration_v8.csv"
fi

conda run -n "$ENV" python run_mixed_mode_fem_czm_v8_campaign.py \
  --parameter-table "${PARAMETER_TABLE:-four_class_exp_floor_exact_model_inputs.csv}" \
  --calibration-csv "$CALROOT/mixed_mode_loading_calibration_v8.csv" \
  --classes "${CLASSES:-ceramic DBTT}" \
  --target-psi-deg="$PSI" \
  --T-K "${T_K:-500}" \
  --outroot "$OUTROOT" \
  --conda-env "$ENV" \
  --max-jobs "${MAX_JOBS:-1}" \
  --nx "${NX:-24}" --ny "${NY:-48}" \
  --tip-h-fine "${TIP_H_FINE:-3e-6}" --tip-ratio "${TIP_RATIO:-1.25}" \
  --crystal-theta-deg "$THETA" \
  --cleave-gamma-aniso "${CLEAVE_GAMMA_ANISO:-0.3}" \
  --traction-probe-radius-m "${TRACTION_PROBE_RADIUS_M:-1e-5}" \
  --shear-emission-weight "${SHEAR_EMISSION_WEIGHT:-1.0}" \
  --directional-factor-max "${DIRECTIONAL_FACTOR_MAX:-5.0}" \
  --event-psi-tol-deg "${EVENT_PSI_TOL_DEG:-2.0}" \
  --max-control-iters "${MAX_CONTROL_ITERS:-8}" \
  --max-alpha-step-deg "${MAX_EVENT_ALPHA_STEP_DEG:-20}" \
  --print-every "${PRINT_EVERY:-50}"

conda run -n "$ENV" python plot_mixed_mode_fem_czm_v8_results.py --root "$OUTROOT"
