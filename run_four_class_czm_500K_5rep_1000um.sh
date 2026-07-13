#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV=${CONDA_ENV:-arrhenius-fem-czm}
PARAMETERS=${PARAMETERS:-four_class_exp_floor_exact_model_inputs.csv}
OUTROOT=${OUTROOT:-runs/four_class_exp_floor_CZM_500K_5rep_1000um_theta45}
CLASSES=${CLASSES:-"ceramic peak weakT DBTT"}
SOLVER_SEEDS=${SOLVER_SEEDS:-"1101 1102 1103 1104 1105"}
THETA=${THETA:-45}
TARGET_EXT_UM=${TARGET_EXT_UM:-1000}
LONG_STEPS=${LONG_STEPS:-50000}
MAX_JOBS=${MAX_JOBS:-1}
FORCE=${FORCE:-0}
SAVE_SNAPSHOTS=${SAVE_SNAPSHOTS:-0}
SNAPSHOT_COLS=${SNAPSHOT_COLS:-5}
SNAPSHOT_BY_EXT_UM=${SNAPSHOT_BY_EXT_UM:-100}
GRID_STEP_UM=${GRID_STEP_UM:-5}
BASE_DU=${BASE_DU:-2e-7}
BASE_DT=${BASE_DT:-8.4}
BURGERS_VECTOR_M=${BURGERS_VECTOR_M:-2.74e-10}
PLOT_ONLY=${PLOT_ONLY:-0}

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
  echo "ERROR: could not resolve Python for environment $CONDA_ENV" >&2
  exit 2
fi

for f in \
  run_four_class_exp_floor_czm_500um_sweep.py \
  run_four_class_czm_500K_seeded_replicates.py \
  run_seeded_sharp_front.py \
  plot_four_class_500K_seeded_rcurves.py \
  "$PARAMETERS"; do
  if [[ ! -f "$f" ]]; then
    echo "ERROR: required file not found: $f" >&2
    exit 2
  fi
done

if [[ "$PLOT_ONLY" == "1" ]]; then
  exec "$PY" plot_four_class_500K_seeded_rcurves.py \
    --root "$OUTROOT" \
    --classes "$CLASSES" \
    --temperature 500 \
    --target-ext-um "$TARGET_EXT_UM" \
    --grid-step-um "$GRID_STEP_UM" \
    --burgers-vector-m "$BURGERS_VECTOR_M"
fi

ARGS=(
  --parameters "$PARAMETERS"
  --outroot "$OUTROOT"
  --classes "$CLASSES"
  --temperature 500
  --solver-seeds "$SOLVER_SEEDS"
  --theta "$THETA"
  --target-ext-um "$TARGET_EXT_UM"
  --long-steps "$LONG_STEPS"
  --max-jobs "$MAX_JOBS"
  --python-bin "$PY"
  --dU "$BASE_DU"
  --dt "$BASE_DT"
  --save-snapshots "$SAVE_SNAPSHOTS"
  --snapshot-cols "$SNAPSHOT_COLS"
  --snapshot-by-ext-um "$SNAPSHOT_BY_EXT_UM"
  --grid-step-um "$GRID_STEP_UM"
  --burgers-vector-m "$BURGERS_VECTOR_M"
)
[[ "$FORCE" == "1" ]] && ARGS+=(--force)

exec "$PY" run_four_class_czm_500K_seeded_replicates.py "${ARGS[@]}"
