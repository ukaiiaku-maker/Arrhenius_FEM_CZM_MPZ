#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV=${CONDA_ENV:-arrhenius-fem-czm}
PYTHON_BIN=${PYTHON_BIN:-}
PARAMETER_ROOT=${PARAMETER_ROOT:-mpz_v9_11_parameters}
OUTROOT_BASE=${OUTROOT_BASE:-runs/mpz_v9_18_5_target_stop_quality_corridor_v1}
TEMPS=${TEMPS:-"700"}
SEEDS=${SEEDS:-"1"}
CLASSES=${CLASSES:-"ceramic weakT DBTT"}
TARGET_EXT_UM=${TARGET_EXT_UM:-60}
STEPS=${STEPS:-15000}
NX=${NX:-36}
NY=${NY:-72}
TIP_H_FINE=${TIP_H_FINE:-1e-6}
TIP_RATIO=${TIP_RATIO:-1.20}
DU=${DU:-2e-7}
DT=${DT:-8.4}
MPZ_LENGTH_UM=${MPZ_LENGTH_UM:-100}
MPZ_N_BINS=${MPZ_N_BINS:-200}
WAKE_LENGTH_UM=${WAKE_LENGTH_UM:-100}
WAKE_N_BINS=${WAKE_N_BINS:-0}
WAKE_SHIELDING=${WAKE_SHIELDING:-1}
WAKE_SHIELD_PROJECTION=${WAKE_SHIELD_PROJECTION:-1}
SAVE_SNAPSHOTS=${SAVE_SNAPSHOTS:-6}
SNAPSHOT_COLS=${SNAPSHOT_COLS:-3}
SNAPSHOT_BY_EXT_UM=${SNAPSHOT_BY_EXT_UM:-10}
BULK_PLASTICITY_MODE=${BULK_PLASTICITY_MODE:-tip_only}
EVENT_TARGET_DQ=${EVENT_TARGET_DQ:-0.05}
EVENT_MIN_DT_S=${EVENT_MIN_DT_S:-1e-12}
EVENT_MAX_FIXED_HOLD_S=${EVENT_MAX_FIXED_HOLD_S:-inf}
PREFINED_MODE_I_CORRIDOR=${PREFINED_MODE_I_CORRIDOR:-1}
CORRIDOR_CENTER_SPACING_UM=${CORRIDOR_CENTER_SPACING_UM:-25}
CORRIDOR_GUARD_UM=${CORRIDOR_GUARD_UM:-10}
MIN_ACCEPTED_TRIANGLE_QUALITY=${MIN_ACCEPTED_TRIANGLE_QUALITY:-0.035}
MIN_ACCEPTED_CHILD_AREA_RATIO=${MIN_ACCEPTED_CHILD_AREA_RATIO:-0.08}
MAX_TIP_H_OVER_DA=${MAX_TIP_H_OVER_DA:-0.75}
MAX_IDENTICAL_GEOMETRY_VETOES=${MAX_IDENTICAL_GEOMETRY_VETOES:-12}

if [[ -z "$PYTHON_BIN" ]]; then
  if [[ "${CONDA_DEFAULT_ENV:-}" == "$CONDA_ENV" ]]; then
    PYTHON_BIN=$(command -v python)
  else
    PYTHON_BIN=$(conda run -n "$CONDA_ENV" python -c 'import sys; print(sys.executable)' | tail -n 1)
  fi
fi

"$PYTHON_BIN" -m pytest -q \
  tests/test_target_stop_quality_corridor_v9185.py \
  tests/test_v9185_routing.py \
  tests/test_mechanically_valid_front_regularization_v9184.py \
  tests/test_v9184_bounded_routing.py \
  tests/test_edge_aware_geometry_recovery_v9183.py \
  tests/test_committed_completion_handshake_v9182.py \
  tests/test_active_event_renewal_rollback_v9181.py \
  tests/test_persistent_wake_v918.py \
  tests/test_persistent_wake_audit_v918.py \
  tests/test_one_fire_routing_v9171.py \
  tests/test_matched_stress_audit_routing_v9171.py \
  tests/test_hazard_clock_source_refresh_v917.py \
  tests/test_kinetic_trial_opening_v916.py \
  tests/test_coupled_event_relaxation_v915.py \
  tests/test_event_driven_remesh_v914.py \
  tests/test_event_hook_v914.py \
  tests/test_material_rcurve_audit_v913.py \
  tests/test_field_snapshots_v913.py

mkdir -p "$OUTROOT_BASE"
overall_rc=0
for T in $TEMPS; do
  case_root="$OUTROOT_BASE/T${T}K"
  mkdir -p "$case_root"
  echo "=== v9.18.5 target-stop quality corridor: T=${T} K, committed target=${TARGET_EXT_UM} um ==="
  set +e
  ARRHENIUS_BULK_PLASTICITY_MODE="$BULK_PLASTICITY_MODE" \
  ARRHENIUS_EVENT_TARGET_DQ="$EVENT_TARGET_DQ" \
  ARRHENIUS_EVENT_MIN_DT_S="$EVENT_MIN_DT_S" \
  ARRHENIUS_EVENT_MAX_FIXED_HOLD_S="$EVENT_MAX_FIXED_HOLD_S" \
  ARRHENIUS_NOMINAL_LOADING_DT_S="$DT" \
  ARRHENIUS_COMMITTED_TARGET_EXTENSION_UM="$TARGET_EXT_UM" \
  ARRHENIUS_WAKE_LENGTH_UM="$WAKE_LENGTH_UM" \
  ARRHENIUS_WAKE_N_BINS="$WAKE_N_BINS" \
  ARRHENIUS_WAKE_SHIELDING="$WAKE_SHIELDING" \
  ARRHENIUS_WAKE_SHIELD_PROJECTION="$WAKE_SHIELD_PROJECTION" \
  ARRHENIUS_PREFINED_MODE_I_CORRIDOR="$PREFINED_MODE_I_CORRIDOR" \
  ARRHENIUS_CORRIDOR_CENTER_SPACING_UM="$CORRIDOR_CENTER_SPACING_UM" \
  ARRHENIUS_CORRIDOR_GUARD_UM="$CORRIDOR_GUARD_UM" \
  ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY="$MIN_ACCEPTED_TRIANGLE_QUALITY" \
  ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO="$MIN_ACCEPTED_CHILD_AREA_RATIO" \
  ARRHENIUS_MAX_TIP_H_OVER_DA="$MAX_TIP_H_OVER_DA" \
  ARRHENIUS_MAX_IDENTICAL_GEOMETRY_VETOES="$MAX_IDENTICAL_GEOMETRY_VETOES" \
  ARRHENIUS_EVENT_STATISTICS=deterministic \
  ARRHENIUS_STOCHASTIC_EMISSION=0 \
  ARRHENIUS_PROPAGATION_CONTROL=raw \
  "$PYTHON_BIN" run_mpz_v9_18_5_persistent_plastic_wake.py \
    --parameter-root "$PARAMETER_ROOT" \
    --outroot "$case_root" \
    --seeds "$SEEDS" \
    --classes "$CLASSES" \
    --T-K "$T" \
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
    2>&1 | tee "$case_root/driver.log"
  rc=${PIPESTATUS[0]}
  set -e
  if [[ $rc -ne 0 ]]; then
    overall_rc=$rc
    echo "FAILED T=${T} K with return code ${rc}" >&2
  fi
done

exit "$overall_rc"
