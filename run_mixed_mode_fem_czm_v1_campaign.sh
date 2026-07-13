#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"
CONDA_ENV="${CONDA_ENV:-arrhenius-fem-czm}"
PYTHON_EXE="${PYTHON_BIN:-}"
if [[ -z "$PYTHON_EXE" ]]; then
  PYTHON_EXE="$(conda run -n "$CONDA_ENV" python -c 'import sys;print(sys.executable)' | awk 'NF{a=$0}END{print a}')"
fi
TARGET_PSI="${TARGET_PSI:--60 -45 -30 -15 0 15 30 45 60}"
CAL_OUT="${CAL_OUT:-runs/mixed_mode_fem_czm_v1_elastic_calibration}"
OUTROOT="${OUTROOT:-runs/mixed_mode_fem_czm_v1_first_passage_500K}"
"$PYTHON_EXE" calibrate_mixed_mode_loading_v1.py --out "$CAL_OUT" --target-psi-deg "$TARGET_PSI" ${CALIBRATION_EXTRA_ARGS:-}
"$PYTHON_EXE" run_mixed_mode_fem_czm_v1_campaign.py \
  --parameter-table "${PARAMETER_TABLE:-four_class_exp_floor_exact_model_inputs.csv}" \
  --calibration-csv "$CAL_OUT/mixed_mode_loading_calibration.csv" \
  --classes "${CLASSES:-ceramic DBTT}" --target-psi-deg "$TARGET_PSI" \
  --seeds "${SEEDS:-1101 1102 1103}" --T-K "${T_K:-500}" --outroot "$OUTROOT" \
  --python-bin "$PYTHON_EXE" --max-jobs "${MAX_JOBS:-1}" ${CAMPAIGN_EXTRA_ARGS:-}
"$PYTHON_EXE" plot_mixed_mode_fem_czm_v1_results.py --root "$OUTROOT"
