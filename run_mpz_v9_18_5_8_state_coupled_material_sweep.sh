#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV=${CONDA_ENV:-arrhenius-fem-czm}
PYTHON_BIN=${PYTHON_BIN:-}
PARAMETER_ROOT=${PARAMETER_ROOT:-mpz_v9_11_parameters}
OUTROOT_BASE=${OUTROOT_BASE:-runs/mpz_v9_18_5_8_state_coupled_material_v1}
TEMPS=${TEMPS:-"700"}
SEEDS=${SEEDS:-"1"}
CLASSES=${CLASSES:-"weakT DBTT"}
TARGET_EXT_UM=${TARGET_EXT_UM:-40}
STEPS=${STEPS:-50000}
NX=${NX:-36}
NY=${NY:-72}
TIP_H_FINE=${TIP_H_FINE:-1e-6}
TIP_RATIO=${TIP_RATIO:-1.20}
DU=${DU:-2e-7}
DT=${DT:-8.4}
PHYSICAL_DA_UM=${PHYSICAL_DA_UM:-5}
MPZ_LENGTH_UM=${MPZ_LENGTH_UM:-100}
MPZ_N_BINS=${MPZ_N_BINS:-200}
WAKE_LENGTH_UM=${WAKE_LENGTH_UM:-100}
WAKE_N_BINS=${WAKE_N_BINS:-0}
WAKE_SHIELDING=${WAKE_SHIELDING:-1}
WAKE_SHIELD_PROJECTION=${WAKE_SHIELD_PROJECTION:-1}
SAVE_SNAPSHOTS=${SAVE_SNAPSHOTS:-8}
SNAPSHOT_COLS=${SNAPSHOT_COLS:-4}
SNAPSHOT_BY_EXT_UM=${SNAPSHOT_BY_EXT_UM:-5}
BULK_PLASTICITY_MODE=${BULK_PLASTICITY_MODE:-bulk_same_pt_km}
ALLOW_ELASTIC_COLLAPSE=${ALLOW_ELASTIC_COLLAPSE:-0}
MIN_NORMALIZED_K_SEPARATION=${MIN_NORMALIZED_K_SEPARATION:-0.02}
MIN_GEOMETRY_FACTOR_SEPARATION=${MIN_GEOMETRY_FACTOR_SEPARATION:-0.01}
EVENT_TARGET_DQ=${EVENT_TARGET_DQ:-0.05}
EVENT_MIN_DT_S=${EVENT_MIN_DT_S:-1e-12}
EVENT_MAX_FIXED_HOLD_S=${EVENT_MAX_FIXED_HOLD_S:-inf}
PREFINED_MODE_I_CORRIDOR=${PREFINED_MODE_I_CORRIDOR:-1}
CORRIDOR_MAX_CENTER_GAP_UM=${CORRIDOR_MAX_CENTER_GAP_UM:-35}
CORRIDOR_GUARD_UM=${CORRIDOR_GUARD_UM:-10}
MIN_INITIAL_TRIANGLE_QUALITY=${MIN_INITIAL_TRIANGLE_QUALITY:-0.035}
MIN_ACCEPTED_TRIANGLE_QUALITY=${MIN_ACCEPTED_TRIANGLE_QUALITY:-0.035}
MIN_ACCEPTED_CHILD_AREA_RATIO=${MIN_ACCEPTED_CHILD_AREA_RATIO:-0.08}
MAX_TIP_H_OVER_DA=${MAX_TIP_H_OVER_DA:-0.75}
MAX_IDENTICAL_GEOMETRY_VETOES=${MAX_IDENTICAL_GEOMETRY_VETOES:-12}

if [[ "$BULK_PLASTICITY_MODE" != "bulk_same_pt_km" && "$ALLOW_ELASTIC_COLLAPSE" != "1" ]]; then
  echo "ERROR: v9.18.5.8 material differentiation requires BULK_PLASTICITY_MODE=bulk_same_pt_km." >&2
  echo "tip_only leaves the FEM bulk elastic and reproduces the geometry-controlled normalized R-curve." >&2
  exit 2
fi

class_count=$(wc -w <<<"$CLASSES" | tr -d ' ')
if [[ "$class_count" -lt 2 ]]; then
  echo "ERROR: v9.18.5.8 differentiation gate requires at least two classes in CLASSES." >&2
  exit 2
fi

if [[ -z "$PYTHON_BIN" ]]; then
  if [[ "${CONDA_DEFAULT_ENV:-}" == "$CONDA_ENV" ]]; then
    PYTHON_BIN=$(command -v python)
  else
    PYTHON_BIN=$(conda run -n "$CONDA_ENV" python -c 'import sys; print(sys.executable)' | tail -n 1)
  fi
fi

"$PYTHON_BIN" -m pytest -q \
  tests/test_state_coupled_material_differentiation_v91858.py \
  tests/test_v91858_routing.py \
  tests/test_subsegment_aware_quality_certification_v91857.py \
  tests/test_v91857_routing.py \
  tests/test_explicit_quality_wrapper_chain_v91856.py \
  tests/test_v91856_routing.py

mkdir -p "$OUTROOT_BASE"
overall_rc=0
for T in $TEMPS; do
  case_root="$OUTROOT_BASE/T${T}K"
  mkdir -p "$case_root"
  echo "=== v9.18.5.8 state-coupled material differentiation: T=${T} K, target=${TARGET_EXT_UM} um ==="
  set +e
  ARRHENIUS_BULK_PLASTICITY_MODE="$BULK_PLASTICITY_MODE" \
  ARRHENIUS_EVENT_TARGET_DQ="$EVENT_TARGET_DQ" \
  ARRHENIUS_EVENT_MIN_DT_S="$EVENT_MIN_DT_S" \
  ARRHENIUS_EVENT_MAX_FIXED_HOLD_S="$EVENT_MAX_FIXED_HOLD_S" \
  ARRHENIUS_NOMINAL_LOADING_DT_S="$DT" \
  ARRHENIUS_COMMITTED_TARGET_EXTENSION_UM="$TARGET_EXT_UM" \
  ARRHENIUS_PHYSICAL_DA_UM="$PHYSICAL_DA_UM" \
  ARRHENIUS_WAKE_LENGTH_UM="$WAKE_LENGTH_UM" \
  ARRHENIUS_WAKE_N_BINS="$WAKE_N_BINS" \
  ARRHENIUS_WAKE_SHIELDING="$WAKE_SHIELDING" \
  ARRHENIUS_WAKE_SHIELD_PROJECTION="$WAKE_SHIELD_PROJECTION" \
  ARRHENIUS_PREFINED_MODE_I_CORRIDOR="$PREFINED_MODE_I_CORRIDOR" \
  ARRHENIUS_CORRIDOR_MAX_CENTER_GAP_UM="$CORRIDOR_MAX_CENTER_GAP_UM" \
  ARRHENIUS_CORRIDOR_GUARD_UM="$CORRIDOR_GUARD_UM" \
  ARRHENIUS_MIN_INITIAL_TRIANGLE_QUALITY="$MIN_INITIAL_TRIANGLE_QUALITY" \
  ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY="$MIN_ACCEPTED_TRIANGLE_QUALITY" \
  ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO="$MIN_ACCEPTED_CHILD_AREA_RATIO" \
  ARRHENIUS_MAX_TIP_H_OVER_DA="$MAX_TIP_H_OVER_DA" \
  ARRHENIUS_MAX_IDENTICAL_GEOMETRY_VETOES="$MAX_IDENTICAL_GEOMETRY_VETOES" \
  ARRHENIUS_EVENT_STATISTICS=deterministic \
  ARRHENIUS_STOCHASTIC_EMISSION=0 \
  ARRHENIUS_PROPAGATION_CONTROL=raw \
  "$PYTHON_BIN" run_mpz_v9_18_5_8_persistent_plastic_wake.py \
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
    --da-phys-um "$PHYSICAL_DA_UM" \
    --mpz-length-um "$MPZ_LENGTH_UM" --mpz-n-bins "$MPZ_N_BINS" \
    --save-snapshots "$SAVE_SNAPSHOTS" --snapshot-cols "$SNAPSHOT_COLS" \
    --snapshot-by-extension-um "$SNAPSHOT_BY_EXT_UM" \
    --event-statistics deterministic --no-stochastic-emission \
    --propagation-control raw --rng-coupling common \
    2>&1 | tee "$case_root/driver.log"
  outer_rc=${PIPESTATUS[0]}

  "$PYTHON_BIN" audit_v91857_subsegment_quality.py \
    "$case_root" "$TARGET_EXT_UM" "$PHYSICAL_DA_UM"
  quality_rc=$?

  "$PYTHON_BIN" audit_v91858_state_coupled_differentiation.py \
    "$case_root" \
    --min-normalized-k-separation "$MIN_NORMALIZED_K_SEPARATION" \
    --min-geometry-factor-separation "$MIN_GEOMETRY_FACTOR_SEPARATION"
  differentiation_rc=$?
  set -e

  rc=$outer_rc
  if [[ $quality_rc -ne 0 ]]; then rc=$quality_rc; fi
  if [[ $differentiation_rc -ne 0 ]]; then rc=$differentiation_rc; fi
  if [[ $rc -ne 0 ]]; then
    overall_rc=$rc
    echo "FAILED T=${T} K with effective return code ${rc}" >&2
    find "$case_root" -path '*/matrix_logs/*.log' -type f -print0 2>/dev/null | \
      while IFS= read -r -d '' log; do
        echo "### $log" >&2
        tail -n 120 "$log" >&2
      done
  fi
done

exit "$overall_rc"
