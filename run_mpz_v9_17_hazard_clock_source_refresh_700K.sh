#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV=${CONDA_ENV:-arrhenius-fem-czm}
PYTHON_BIN=${PYTHON_BIN:-}
PARAMETER_ROOT=${PARAMETER_ROOT:-mpz_v9_11_parameters}
OUTROOT=${OUTROOT:-runs/mpz_v9_17_hazard_clock_source_refresh_700K_v1}
AUDIT_OUT=${AUDIT_OUT:-runs/mpz_v9_17_matched_stress_audit_700K_v1}
SEEDS=${SEEDS:-"1"}
CLASSES=${CLASSES:-"ceramic weakT DBTT"}
T_K=${T_K:-700}
TARGET_EXT_UM=${TARGET_EXT_UM:-10}
STEPS=${STEPS:-10000}
NX=${NX:-36}
NY=${NY:-72}
TIP_H_FINE=${TIP_H_FINE:-1e-6}
TIP_RATIO=${TIP_RATIO:-1.20}
DU=${DU:-2e-7}
DT=${DT:-8.4}
MPZ_LENGTH_UM=${MPZ_LENGTH_UM:-100}
MPZ_N_BINS=${MPZ_N_BINS:-200}
SAVE_SNAPSHOTS=${SAVE_SNAPSHOTS:-3}
SNAPSHOT_COLS=${SNAPSHOT_COLS:-3}
SNAPSHOT_BY_EXT_UM=${SNAPSHOT_BY_EXT_UM:-5}
BULK_PLASTICITY_MODE=${BULK_PLASTICITY_MODE:-tip_only}
EVENT_TARGET_DQ=${EVENT_TARGET_DQ:-0.05}
EVENT_MIN_DT_S=${EVENT_MIN_DT_S:-1e-12}
EVENT_MAX_FIXED_HOLD_S=${EVENT_MAX_FIXED_HOLD_S:-inf}

if [[ -z "$PYTHON_BIN" ]]; then
  if [[ "${CONDA_DEFAULT_ENV:-}" == "$CONDA_ENV" ]]; then
    PYTHON_BIN=$(command -v python)
  else
    PYTHON_BIN=$(conda run -n "$CONDA_ENV" python -c 'import sys; print(sys.executable)' | tail -n 1)
  fi
fi

"$PYTHON_BIN" -m pytest -q \
  tests/test_hazard_clock_source_refresh_v917.py \
  tests/test_kinetic_trial_opening_v916.py \
  tests/test_coupled_event_relaxation_v915.py \
  tests/test_event_driven_remesh_v914.py \
  tests/test_event_hook_v914.py \
  tests/test_material_rcurve_audit_v913.py \
  tests/test_field_snapshots_v913.py

"$PYTHON_BIN" audit_matched_stress_classes_v917.py \
  --parameter-root "$PARAMETER_ROOT" \
  --T-K "$T_K" \
  --classes "$CLASSES" \
  --mpz-length-um "$MPZ_LENGTH_UM" \
  --mpz-n-bins "$MPZ_N_BINS" \
  --out "$AUDIT_OUT"

mkdir -p "$OUTROOT"
ARRHENIUS_BULK_PLASTICITY_MODE="$BULK_PLASTICITY_MODE" \
ARRHENIUS_EVENT_TARGET_DQ="$EVENT_TARGET_DQ" \
ARRHENIUS_EVENT_MIN_DT_S="$EVENT_MIN_DT_S" \
ARRHENIUS_EVENT_MAX_FIXED_HOLD_S="$EVENT_MAX_FIXED_HOLD_S" \
ARRHENIUS_EVENT_STATISTICS=deterministic \
ARRHENIUS_STOCHASTIC_EMISSION=0 \
ARRHENIUS_PROPAGATION_CONTROL=raw \
"$PYTHON_BIN" run_mpz_v9_17_hazard_clock_source_refresh.py \
  --parameter-root "$PARAMETER_ROOT" \
  --outroot "$OUTROOT" \
  --seeds "$SEEDS" \
  --classes "$CLASSES" \
  --T-K "$T_K" \
  --target-extension-um "$TARGET_EXT_UM" \
  --steps "$STEPS" \
  --nx "$NX" --ny "$NY" \
  --tip-h-fine "$TIP_H_FINE" --tip-ratio "$TIP_RATIO" \
  --dU "$DU" --dt "$DT" \
  --mpz-length-um "$MPZ_LENGTH_UM" --mpz-n-bins "$MPZ_N_BINS" \
  --save-snapshots "$SAVE_SNAPSHOTS" --snapshot-cols "$SNAPSHOT_COLS" \
  --snapshot-by-extension-um "$SNAPSHOT_BY_EXT_UM" \
  --event-statistics deterministic --no-stochastic-emission \
  --propagation-control raw --rng-coupling common \
  2>&1 | tee "$OUTROOT/driver.log"
