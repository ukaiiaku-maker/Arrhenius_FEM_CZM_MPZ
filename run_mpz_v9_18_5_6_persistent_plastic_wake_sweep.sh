#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV=${CONDA_ENV:-arrhenius-fem-czm}
PYTHON_BIN=${PYTHON_BIN:-}
PARAMETER_ROOT=${PARAMETER_ROOT:-mpz_v9_11_parameters}
OUTROOT_BASE=${OUTROOT_BASE:-runs/mpz_v9_18_5_6_explicit_quality_gate_v1}
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
  tests/test_explicit_quality_wrapper_chain_v91856.py \
  tests/test_v91856_routing.py \
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
  echo "=== v9.18.5.6 explicit quality chain: T=${T} K, committed target=${TARGET_EXT_UM} um ==="
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
  "$PYTHON_BIN" run_mpz_v9_18_5_6_persistent_plastic_wake.py \
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

  "$PYTHON_BIN" - "$case_root" "$TARGET_EXT_UM" "$PHYSICAL_DA_UM" <<'PY'
import json, math, sys
from pathlib import Path
root = Path(sys.argv[1])
target = float(sys.argv[2])
da = float(sys.argv[3])
summary = root / "v9_13_campaign_summary.json"
if not summary.exists():
    raise SystemExit(f"campaign summary missing: {summary}")
rows = json.loads(summary.read_text())
failed = [r for r in rows if int(r.get("subprocess_returncode", r.get("returncode", 1)) or 0) != 0]
if failed:
    raise SystemExit(f"inner solver failures: {[(r.get('class'), r.get('subprocess_returncode'), r.get('log')) for r in failed]}")
expected = int(round(target / da))
for row in rows:
    case_dir = Path(row["case_dir"])
    audit_path = case_dir / "explicit_quality_wrapper_chain_v91856.json"
    if not audit_path.exists():
        raise SystemExit(f"v9.18.5.6 quality audit missing: {audit_path}")
    audit = json.loads(audit_path.read_text())
    accepted = audit.get("accepted_events", [])
    vetoes = audit.get("quality_vetoes", [])
    if not audit.get("run_completed_without_exception", False):
        raise SystemExit(f"quality audit reports runtime failure: {audit_path}")
    if len(accepted) != expected:
        raise SystemExit(f"quality gate count mismatch: accepted={len(accepted)} expected={expected} path={audit_path}")
    if vetoes:
        raise SystemExit(f"unexpected production quality vetoes: {len(vetoes)} path={audit_path}")
    if not all(bool(x.get("accepted", False)) for x in accepted):
        raise SystemExit(f"nonaccepted row in quality audit: {audit_path}")
    print(
        f"QUALITY AUDIT {case_dir.name}: accepted={len(accepted)} "
        f"resolution_warnings={len(audit.get('resolution_warnings', []))} "
        f"qmin={min(x['min_triangle_quality'] for x in accepted):.6g} "
        f"area_ratio_min={min(x['min_child_area_ratio'] for x in accepted):.6g}"
    )
PY
  audit_rc=$?
  set -e

  rc=$outer_rc
  if [[ $audit_rc -ne 0 ]]; then rc=$audit_rc; fi
  if [[ $rc -ne 0 ]]; then
    overall_rc=$rc
    echo "FAILED T=${T} K with effective return code ${rc}" >&2
    echo "--- per-case traceback tail ---" >&2
    find "$case_root" -path '*/matrix_logs/*.log' -type f -print0 2>/dev/null | \
      while IFS= read -r -d '' log; do
        echo "### $log" >&2
        tail -n 160 "$log" >&2
      done
  fi
done

exit "$overall_rc"
