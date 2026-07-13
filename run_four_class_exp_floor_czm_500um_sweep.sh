#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV=${CONDA_ENV:-arrhenius-fem-czm}
OUTROOT=${OUTROOT:-runs/four_class_exp_floor_CZM_no_branch_500um_theta45}
PARAMETERS=${PARAMETERS:-four_class_exp_floor_exact_model_inputs.csv}
CLASSES=${CLASSES:-"ceramic peak weakT DBTT"}
TEMPS=${TEMPS:-"300 400 500 600 700 800 900 1000 1100 1200"}
THETA=${THETA:-45}
TARGET_EXT_UM=${TARGET_EXT_UM:-500}
LONG_STEPS=${LONG_STEPS:-20000}
MAX_JOBS=${MAX_JOBS:-1}
SAVE_SNAPSHOTS=${SAVE_SNAPSHOTS:-0}
SNAPSHOT_COLS=${SNAPSHOT_COLS:-5}
SNAPSHOT_BY_EXT_UM=${SNAPSHOT_BY_EXT_UM:-50}
FORCE=${FORCE:-0}

ARGS=(
  --parameters "$PARAMETERS"
  --outroot "$OUTROOT"
  --classes "$CLASSES"
  --temps "$TEMPS"
  --theta "$THETA"
  --target-ext-um "$TARGET_EXT_UM"
  --long-steps "$LONG_STEPS"
  --max-jobs "$MAX_JOBS"
  --conda-env "$CONDA_ENV"
  --save-snapshots "$SAVE_SNAPSHOTS"
  --snapshot-cols "$SNAPSHOT_COLS"
  --snapshot-by-ext-um "$SNAPSHOT_BY_EXT_UM"
)
[[ -n "${PYTHON_BIN:-}" ]] && ARGS+=(--python-bin "$PYTHON_BIN")
[[ "$FORCE" == "1" ]] && ARGS+=(--force)

# Resolve orchestration Python from the isolated environment unless explicitly supplied.
if [[ -n "${PYTHON_BIN:-}" ]]; then
  PY="$PYTHON_BIN"
else
  if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: conda not found; set PYTHON_BIN explicitly" >&2
    exit 2
  fi
  PY="$(conda run -n "$CONDA_ENV" python -c 'import sys; print(sys.executable)' 2>&1 | tr -d '\r' | awk 'NF {last=$0} END {print last}')"
fi
if [[ -z "$PY" || ! -x "$PY" ]]; then
  echo "ERROR: could not resolve Python interpreter" >&2
  exit 2
fi

exec "$PY" run_four_class_exp_floor_czm_500um_sweep.py "${ARGS[@]}"
