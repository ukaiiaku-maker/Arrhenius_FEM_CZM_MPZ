#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV=${CONDA_ENV:-arrhenius-fem-czm}
PYTHON_BIN=${PYTHON_BIN:-}
PARAMETER_ROOT=${PARAMETER_ROOT:-mpz_v9_11_parameters}
OUTROOT_BASE=${OUTROOT_BASE:-runs/mpz_v9_18_5_5_resolution_audit_only_v1}
TEMPS=${TEMPS:-"700"}
SEEDS=${SEEDS:-"1"}
CLASSES=${CLASSES:-"ceramic"}
TARGET_EXT_UM=${TARGET_EXT_UM:-60}
STEPS=${STEPS:-15000}
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
SAVE_SNAPSHOTS=${SAVE_SNAPSHOTS:-6}
SNAPSHOT_COLS=${SNAPSHOT_COLS:-3}
SNAPSHOT_BY_EXT_UM=${SNAPSHOT_BY_EXT_UM:-5}
BULK_PLASTICITY_MODE=${BULK_PLASTICITY_MODE:-tip_only}
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

if [[ -z "$PYTHON_BIN" ]]; then
  if [[ "${CONDA_DEFAULT_ENV:-}" == "$CONDA_ENV" ]]; then
    PYTHON_BIN=$(command -v python)
  else
    PYTHON_BIN=$(conda run -n "$CONDA_ENV" python -c 'import sys; print(sys.executable)' | tail -n 1)
  fi
fi

"$PYTHON_BIN" -m pytest -q \
  tests/test_resolution_audit_only_consecutive_veto_v91855.py \
  tests/test_v91855_routing.py \
  tests/test_active_tip_resolution_veto_guard_v91854.py \
  tests/test_v91854_routing.py \
  tests/test_quality_selected_corridor_v91853.py \
  tests/test_v91853_routing.py \
  tests/test_compact_corridor_mesh_v91852.py \
  tests/test_v91852_routing.py \
  tests/test_safe_target_stop_horizon_v91851.py \
  tests/test_v91851_routing.py \
  tests/test_target_stop_quality_corridor_v9185.py \
  tests/test_time_aware_sequence_audit_v9185.py \
  tests/test_v9185_routing.py

mkdir -p "$OUTROOT_BASE"
overall_rc=0
for T in $TEMPS; do
  case_root="$OUTROOT_BASE/T${T}K"
  mkdir -p "$case_root"
  echo "=== v9.18.5.5 resolution audit-only: T=${T} K, committed target=${TARGET_EXT_UM} um ==="
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
  "$PYTHON_BIN" run_mpz_v9_18_5_5_persistent_plastic_wake.py \
    --parameter-root "$PARAMETER_ROOT" \
    --outroot "$case_root" \
    --seeds "$SEEDS" --classes "$CLASSES" --T-K "$T" \
    --target-extension-um "$TARGET_EXT_UM" --steps "$STEPS" \
    --nx "$NX" --ny "$NY" \
    --tip-h-fine "$TIP_H_FINE" --tip-ratio "$TIP_RATIO" \
    --dU "$DU" --dt "$DT" --da-phys-um "$PHYSICAL_DA_UM" \
    --mpz-length-um "$MPZ_LENGTH_UM" --mpz-n-bins "$MPZ_N_BINS" \
    --save-snapshots "$SAVE_SNAPSHOTS" --snapshot-cols "$SNAPSHOT_COLS" \
    --snapshot-by-extension-um "$SNAPSHOT_BY_EXT_UM" \
    --event-statistics deterministic --no-stochastic-emission \
    --propagation-control raw --rng-coupling common \
    2>&1 | tee "$case_root/driver.log"
  outer_rc=${PIPESTATUS[0]}

  summary="$case_root/v9_13_campaign_summary.json"
  "$PYTHON_BIN" - "$summary" <<'PY'
import json, sys
from pathlib import Path
p = Path(sys.argv[1])
if not p.exists():
    print(f"campaign summary missing: {p}", file=sys.stderr)
    raise SystemExit(1)
rows = json.loads(p.read_text())
failed = [r for r in rows if int(r.get("subprocess_returncode", r.get("returncode", 1)) or 0) != 0]
if failed:
    print("inner solver failures:", [(r.get("class"), r.get("subprocess_returncode"), r.get("log")) for r in failed], file=sys.stderr)
    raise SystemExit(1)
PY
  summary_rc=$?
  set -e

  rc=$outer_rc
  if [[ $summary_rc -ne 0 ]]; then rc=$summary_rc; fi
  if [[ $rc -ne 0 ]]; then
    overall_rc=$rc
    echo "FAILED T=${T} K with effective return code ${rc}" >&2
    find "$case_root" -path '*/matrix_logs/*.log' -type f -print0 2>/dev/null | \
      while IFS= read -r -d '' log; do
        echo "### $log" >&2
        tail -n 180 "$log" >&2
      done
  fi
done

exit "$overall_rc"
