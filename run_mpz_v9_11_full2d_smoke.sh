#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV=${CONDA_ENV:-arrhenius-fem-czm}
PARAMETER_ROOT=${PARAMETER_ROOT:-mpz_v9_11_parameters}
CALIBRATION_CSV=${CALIBRATION_CSV:?Set CALIBRATION_CSV to a verified v8 production-backend calibration CSV containing psi=0}
CLASS=${CLASS:-DBTT}
T_K=${T_K:-700}
OUTROOT=${OUTROOT:-runs/mpz_v9_11_full2d_${CLASS}_${T_K}K_smoke_v1}
RUN_SOLVER=${RUN_SOLVER:-1}

python verify_mpz_v9_11_install.py .
python verify_mpz_v9_11_physics.py --parameter-root "$PARAMETER_ROOT"
python -m pytest -q \
  tests/test_mpz_v9_11_2d_coupling.py \
  tests/test_mpz_v9_10_2_independent_shapes.py \
  tests/test_bulk_pt_plasticity.py

if [[ "$RUN_SOLVER" != "1" ]]; then
  echo "Preflight passed; RUN_SOLVER=$RUN_SOLVER so the FEM/CZM solve was skipped."
  exit 0
fi

python run_mixed_mode_fem_czm_v9_11_campaign.py \
  --conda-env "$CONDA_ENV" \
  --parameter-root "$PARAMETER_ROOT" \
  --calibration-csv "$CALIBRATION_CSV" \
  --classes "$CLASS" \
  --target-psi-deg "0" \
  --T-K "$T_K" \
  --outroot "$OUTROOT" \
  --max-jobs 1 \
  --nx 24 --ny 48 --steps 1200
