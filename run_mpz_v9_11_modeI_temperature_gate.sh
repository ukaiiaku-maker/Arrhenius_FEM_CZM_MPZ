#!/usr/bin/env bash
set -euo pipefail

CALIBRATION_CSV=${CALIBRATION_CSV:?Set CALIBRATION_CSV to the verified v8 calibration CSV containing psi=0}
CONDA_ENV=${CONDA_ENV:-arrhenius-fem-czm}
PARAMETER_ROOT=${PARAMETER_ROOT:-mpz_v9_11_parameters}
CLASSES=${CLASSES:-"ceramic weakT DBTT"}
TEMPS=${TEMPS:-"300 700 900 1200"}
ROOT=${ROOT:-runs/mpz_v9_11_modeI_temperature_gate_v1}
MAX_JOBS=${MAX_JOBS:-1}

for T in $TEMPS; do
  echo "=== v9.11 Mode-I gate T=${T} K ==="
  python run_mixed_mode_fem_czm_v9_11_campaign.py \
    --conda-env "$CONDA_ENV" \
    --parameter-root "$PARAMETER_ROOT" \
    --calibration-csv "$CALIBRATION_CSV" \
    --classes "$CLASSES" \
    --target-psi-deg "0" \
    --T-K "$T" \
    --outroot "$ROOT/T${T}" \
    --max-jobs "$MAX_JOBS" \
    "$@"
done
