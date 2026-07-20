#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV=${CONDA_ENV:-arrhenius-fem-czm}
PYTHON_BIN=${PYTHON_BIN:-}
MODE=${MODE:-smoke}
MAX_JOBS=${MAX_JOBS:-2}
SKIP_EXISTING=${SKIP_EXISTING:-1}
RUN_TESTS=${RUN_TESTS:-1}
FAIL_FAST=${FAIL_FAST:-0}
OPTIONS=${OPTIONS:-}
TEMPERATURES=${TEMPERATURES:-}
TARGET_EXT_UM=${TARGET_EXT_UM:-}

case "$MODE" in
  smoke)
    DEFAULT_OUTROOT=runs/v10_0_5_12_2_phase_c_smoke_DBTT_700K_50um_v1
    ;;
  anchors)
    DEFAULT_OUTROOT=runs/v10_0_5_12_2_phase_c_500um_theta45_v1
    ;;
  full)
    DEFAULT_OUTROOT=runs/v10_0_5_12_2_phase_c_500um_theta45_v1
    ;;
  *)
    echo "ERROR: MODE must be smoke, anchors, or full" >&2
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

args=(
  "$PYTHON_BIN" run_v10_0_5_12_2_phase_c_monotonic.py
  --mode "$MODE"
  --python-bin "$PYTHON_BIN"
  --conda-env "$CONDA_ENV"
  --outroot "$OUTROOT"
  --max-jobs "$MAX_JOBS"
)

[[ -n "$OPTIONS" ]] && args+=(--options "$OPTIONS")
[[ -n "$TEMPERATURES" ]] && args+=(--temperatures "$TEMPERATURES")
[[ -n "$TARGET_EXT_UM" ]] && args+=(--target-extension-um "$TARGET_EXT_UM")

if [[ "$SKIP_EXISTING" == "1" ]]; then
  args+=(--skip-existing)
else
  args+=(--no-skip-existing)
fi
if [[ "$RUN_TESTS" == "1" ]]; then
  args+=(--run-tests)
else
  args+=(--no-run-tests)
fi
if [[ "$FAIL_FAST" == "1" ]]; then
  args+=(--fail-fast)
fi

printf 'Launching:'
printf ' %q' "${args[@]}"
printf '\n'
exec "${args[@]}"
