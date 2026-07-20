#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
cd "$ROOT_DIR"

CONDA_ENV=${CONDA_ENV:-arrhenius-fem-czm}
PYTHON_BIN=${PYTHON_BIN:-}
MODE=${MODE:-full}
MAX_JOBS=${MAX_JOBS:-2}
SKIP_FINISHED=${SKIP_FINISHED:-1}
RUN_TESTS=${RUN_TESTS:-1}
OPTIONS=${OPTIONS:-}
TEMPERATURES=${TEMPERATURES:-}
TARGET_EXT_UM=${TARGET_EXT_UM:-}
SAVE_SNAPSHOTS=${SAVE_SNAPSHOTS:-3}
SNAPSHOT_COLS=${SNAPSHOT_COLS:-3}
SNAPSHOT_INTERVAL_UM=${SNAPSHOT_INTERVAL_UM:-50}
PRINT_EVERY=${PRINT_EVERY:-25}

# Rate-preserving adaptive macro-step. The historical physical loading rate is
# BASE_DU/BASE_DT. DU and DT may be coarsened together; the existing adaptive
# event controller reduces both by the same accepted trial fraction near events.
BASE_DU=${BASE_DU:-2e-7}
BASE_DT=${BASE_DT:-8.4}
DU=${DU:-2e-5}
DT=${DT:-840}
ALLOW_RATE_CHANGE=${ALLOW_RATE_CHANGE:-0}

case "$MODE" in
  smoke)
    DEFAULT_OUTROOT="$ROOT_DIR/runs/v10_0_5_13_3_tip_only_DBTT_700K_20um_macro100_v1"
    ;;
  full)
    DEFAULT_OUTROOT="$ROOT_DIR/runs/v10_0_5_13_3_tip_only_4class_300_1200K_100um_macro100_v1"
    ;;
  *)
    echo "ERROR: MODE must be smoke or full" >&2
    exit 2
    ;;
esac
OUTROOT=${OUTROOT:-$DEFAULT_OUTROOT}
mkdir -p "$OUTROOT"

if [[ -z "$PYTHON_BIN" ]]; then
  if [[ "${CONDA_DEFAULT_ENV:-}" == "$CONDA_ENV" ]]; then
    PYTHON_BIN=$(command -v python)
  else
    PYTHON_BIN=$(conda run -n "$CONDA_ENV" python -c 'import sys; print(sys.executable)' | tail -n 1)
  fi
fi

ACTUAL_PACKAGE=$(
  "$PYTHON_BIN" -c 'from pathlib import Path; import arrhenius_fracture; print(Path(arrhenius_fracture.__file__).resolve().parent)'
)
EXPECTED_PACKAGE=$(cd "$ROOT_DIR/arrhenius_fracture" && pwd)
if [[ "$ACTUAL_PACKAGE" != "$EXPECTED_PACKAGE" ]]; then
  echo "ERROR: arrhenius_fracture resolves outside this FEM/CZM installation." >&2
  echo "  expected: $EXPECTED_PACKAGE" >&2
  echo "  actual:   $ACTUAL_PACKAGE" >&2
  echo "Reinstall with: $PYTHON_BIN -m pip install -e $ROOT_DIR --no-deps" >&2
  exit 2
fi

RATE_AUDIT=$(
  "$PYTHON_BIN" - "$BASE_DU" "$BASE_DT" "$DU" "$DT" "$ALLOW_RATE_CHANGE" <<'PY'
import math
import sys
base_du, base_dt, du, dt = map(float, sys.argv[1:5])
allow = sys.argv[5] == "1"
if base_du <= 0 or base_dt <= 0 or du <= 0 or dt <= 0:
    raise SystemExit("dU and dt must be positive")
base_rate = base_du / base_dt
rate = du / dt
rel = abs(rate - base_rate) / max(abs(base_rate), 1e-300)
if rel > 1e-12 and not allow:
    raise SystemExit(
        f"loading-rate change rejected: base={base_rate:.16g} requested={rate:.16g}; "
        "set ALLOW_RATE_CHANGE=1 only for an intentional rate study"
    )
print(
    f"base_rate_m_per_s={base_rate:.16g} requested_rate_m_per_s={rate:.16g} "
    f"macro_step_factor={du/base_du:.16g} time_step_factor={dt/base_dt:.16g}"
)
PY
)

args=(
  "$PYTHON_BIN" "$ROOT_DIR/run_v10_0_5_13_3_barrier_only_monotonic.py"
  --mode "$MODE"
  --python-bin "$PYTHON_BIN"
  --conda-env "$CONDA_ENV"
  --outroot "$OUTROOT"
  --max-jobs "$MAX_JOBS"
  --save-snapshots "$SAVE_SNAPSHOTS"
  --snapshot-cols "$SNAPSHOT_COLS"
  --snapshot-interval-um "$SNAPSHOT_INTERVAL_UM"
  --print-every "$PRINT_EVERY"
  --dU "$DU"
  --dt "$DT"
)

[[ -n "$OPTIONS" ]] && args+=(--options "$OPTIONS")
[[ -n "$TEMPERATURES" ]] && args+=(--temperatures "$TEMPERATURES")
[[ -n "$TARGET_EXT_UM" ]] && args+=(--target-extension-um "$TARGET_EXT_UM")

if [[ "$SKIP_FINISHED" == "1" ]]; then
  args+=(--skip-finished)
else
  args+=(--no-skip-finished)
fi
if [[ "$RUN_TESTS" == "1" ]]; then
  args+=(--run-tests)
else
  args+=(--no-run-tests)
fi

printf 'FEM/CZM installation: %s\n' "$ROOT_DIR"
printf 'Python: %s\n' "$PYTHON_BIN"
printf 'Output: %s\n' "$OUTROOT"
printf 'Point release: 10.0.5.13.3\n'
printf 'Plasticity scope: tip_only (elastic bulk + moving crack-tip MPZ)\n'
printf 'Ramp audit: %s\n' "$RATE_AUDIT"
printf 'Launching:'
printf ' %q' "${args[@]}"
printf '\n'

exec "${args[@]}"
