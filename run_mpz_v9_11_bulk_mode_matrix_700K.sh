#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV=${CONDA_ENV:-arrhenius-fem-czm}
PYTHON_BIN=${PYTHON_BIN:-}
PARAMETER_ROOT=${PARAMETER_ROOT:-mpz_v9_11_parameters}
OUTROOT=${OUTROOT:-runs/mpz_v9_11_bulk_mode_matrix_700K_v1}
T_K=${T_K:-700}
TARGET_EXT_UM=${TARGET_EXT_UM:-500}
STEPS=${STEPS:-6000}
NX=${NX:-36}
NY=${NY:-72}
TIP_H_FINE=${TIP_H_FINE:-1e-6}
TIP_RATIO=${TIP_RATIO:-1.20}
DU=${DU:-2e-7}
DT=${DT:-8.4}
N_STAGGER=${N_STAGGER:-2}
PRINT_EVERY=${PRINT_EVERY:-25}
ADAPTIVE_EVENT_TARGET=${ADAPTIVE_EVENT_TARGET:-0.15}
DA_PHYS_UM=${DA_PHYS_UM:-5}
MPZ_LENGTH_UM=${MPZ_LENGTH_UM:-100}
MPZ_N_BINS=${MPZ_N_BINS:-200}
CRYSTAL_THETA_DEG=${CRYSTAL_THETA_DEG:-45}
SAVE_SNAPSHOTS=${SAVE_SNAPSHOTS:-12}
SNAPSHOT_COLS=${SNAPSHOT_COLS:-4}
SNAPSHOT_BY_EXT_UM=${SNAPSHOT_BY_EXT_UM:-50}
SKIP_EXISTING=${SKIP_EXISTING:-1}

if [[ -z "$PYTHON_BIN" ]]; then
  if [[ "${CONDA_DEFAULT_ENV:-}" == "$CONDA_ENV" ]]; then
    PYTHON_BIN=$(command -v python)
  else
    PYTHON_BIN=$(conda run -n "$CONDA_ENV" python -c 'import sys; print(sys.executable)' | tail -n 1)
  fi
fi

mkdir -p "$OUTROOT"

"$PYTHON_BIN" verify_mpz_v9_11_install.py .
"$PYTHON_BIN" verify_mpz_v9_11_physics.py --parameter-root "$PARAMETER_ROOT"
"$PYTHON_BIN" -m pytest -q \
  tests/test_mode_i_first_passage_v9_11.py \
  tests/test_bulk_state_v911.py \
  tests/test_mpz_v9_11_2d_coupling.py \
  tests/test_mpz_v9_10_2_independent_shapes.py \
  tests/test_bulk_pt_plasticity.py

EXTRA=()
if [[ "$SKIP_EXISTING" == "1" ]]; then
  EXTRA+=(--skip-existing)
else
  EXTRA+=(--no-skip-existing)
fi

"$PYTHON_BIN" run_mpz_v9_11_bulk_mode_matrix_700K.py \
  --parameter-root "$PARAMETER_ROOT" \
  --outroot "$OUTROOT" \
  --T-K "$T_K" \
  --target-extension-um "$TARGET_EXT_UM" \
  --steps "$STEPS" \
  --nx "$NX" --ny "$NY" \
  --tip-h-fine "$TIP_H_FINE" --tip-ratio "$TIP_RATIO" \
  --dU "$DU" --dt "$DT" \
  --n-stagger "$N_STAGGER" \
  --print-every "$PRINT_EVERY" \
  --adaptive-event-target "$ADAPTIVE_EVENT_TARGET" \
  --da-phys-um "$DA_PHYS_UM" \
  --mpz-length-um "$MPZ_LENGTH_UM" \
  --mpz-n-bins "$MPZ_N_BINS" \
  --crystal-theta-deg "$CRYSTAL_THETA_DEG" \
  --save-snapshots "$SAVE_SNAPSHOTS" \
  --snapshot-cols "$SNAPSHOT_COLS" \
  --snapshot-by-extension-um "$SNAPSHOT_BY_EXT_UM" \
  "${EXTRA[@]}" \
  2>&1 | tee "$OUTROOT/driver.log"
