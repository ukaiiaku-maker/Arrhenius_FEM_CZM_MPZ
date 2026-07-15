#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV=${CONDA_ENV:-arrhenius-fem-czm}
PARAMETER_ROOT=${PARAMETER_ROOT:-mpz_v9_11_parameters}
CALIBRATION_CSV=${CALIBRATION_CSV:?Set CALIBRATION_CSV to the verified v8 production-backend calibration CSV}
CLASSES=${CLASSES:-"ceramic weakT DBTT"}
TARGET_PSI_DEG=${TARGET_PSI_DEG:-"-60 -45 -30 -15 0 15 30 45 60"}
T_K=${T_K:-500}
OUTROOT=${OUTROOT:-runs/mixed_mode_fem_czm_v9_11_MPZ_${T_K}K}
MAX_JOBS=${MAX_JOBS:-1}

python run_mixed_mode_fem_czm_v9_11_campaign.py \
  --conda-env "$CONDA_ENV" \
  --parameter-root "$PARAMETER_ROOT" \
  --calibration-csv "$CALIBRATION_CSV" \
  --classes "$CLASSES" \
  --target-psi-deg "$TARGET_PSI_DEG" \
  --T-K "$T_K" \
  --outroot "$OUTROOT" \
  --max-jobs "$MAX_JOBS" \
  "$@"
