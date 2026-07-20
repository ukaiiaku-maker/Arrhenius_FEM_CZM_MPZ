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

case "$MODE" in
  smoke)
    DEFAULT_OUTROOT="$ROOT_DIR/runs/v10_0_5_13_barrier_only_DBTT_700K_20um_smoke_v1"
    ;;
  full)
    DEFAULT_OUTROOT="$ROOT_DIR/runs/v10_0_5_13_barrier_only_4class_300_1200K_100um_v1"
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

# Fail closed against the historical shared arrhenius_fracture namespace.  The
# active editable package must resolve to this installation directory.
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

args=(
  "$PYTHON_BIN" "$ROOT_DIR/run_v10_0_5_13_barrier_only_monotonic.py"
  --mode "$MODE"
  --python-bin "$PYTHON_BIN"
  --conda-env "$CONDA_ENV"
  --outroot "$OUTROOT"
  --max-jobs "$MAX_JOBS"
  --save-snapshots "$SAVE_SNAPSHOTS"
  --snapshot-cols "$SNAPSHOT_COLS"
  --snapshot-interval-um "$SNAPSHOT_INTERVAL_UM"
  --print-every "$PRINT_EVERY"
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
printf 'Launching:'
printf ' %q' "${args[@]}"
printf '\n'

exec "${args[@]}"
