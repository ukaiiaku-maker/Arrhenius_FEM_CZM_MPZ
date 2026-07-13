#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)";cd "$ROOT"
ENV="${CONDA_ENV:-arrhenius-fem-czm}"
CALROOT="${CALROOT:-runs/mixed_mode_fem_czm_v4_anisotropic_calibration}"
OUTROOT="${OUTROOT:-runs/mixed_mode_fem_czm_v4_anisotropic_500K}"
PSI="${TARGET_PSI:--60 -45 -30 -15 0 15 30 45 60}"
THETA="${CRYSTAL_THETA_DEG:-0}"
if [[ "${RECALIBRATE:-0}" == "1" || ! -f "$CALROOT/mixed_mode_loading_calibration_v4.csv" ]]; then
 conda run -n "$ENV" python calibrate_mixed_mode_loading_v4.py --out "$CALROOT" --target-psi-deg="$PSI" --crystal-theta-deg "$THETA" --traction-probe-radius-m "${TRACTION_PROBE_RADIUS_M:-1e-5}"
else echo "Using existing anisotropic calibration: $CALROOT/mixed_mode_loading_calibration_v4.csv";fi
conda run -n "$ENV" python run_mixed_mode_fem_czm_v4_campaign.py --parameter-table "${PARAMETER_TABLE:-four_class_exp_floor_exact_model_inputs.csv}" --calibration-csv "$CALROOT/mixed_mode_loading_calibration_v4.csv" --classes "${CLASSES:-ceramic DBTT}" --target-psi-deg="$PSI" --T-K "${T_K:-500}" --outroot "$OUTROOT" --conda-env "$ENV" --max-jobs "${MAX_JOBS:-1}" --crystal-theta-deg "$THETA" --traction-probe-radius-m "${TRACTION_PROBE_RADIUS_M:-1e-5}"
conda run -n "$ENV" python plot_mixed_mode_fem_czm_v4_results.py --root "$OUTROOT"
