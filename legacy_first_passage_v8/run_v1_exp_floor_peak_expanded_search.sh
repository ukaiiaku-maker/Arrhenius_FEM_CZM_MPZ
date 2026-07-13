#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV=${CONDA_ENV:-arrhenius-fem-czm}
SOURCE_RUN=${SOURCE_RUN:-runs/v1_exp_floor_four_class_tuning}
FOCUSED_RUN=${FOCUSED_RUN:-runs/v1_exp_floor_peak_focused_refinement}
OUT=${OUT:-runs/v1_exp_floor_peak_expanded_search}
SMOKE=${SMOKE:-0}
RESUME=${RESUME:-1}

if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda not found" >&2
  exit 2
fi
PYTHON_EXE="$(conda run -n "$CONDA_ENV" python -c 'import sys; print(sys.executable)' 2>&1 | tr -d '\r' | awk 'NF {last=$0} END {print last}')"
if [[ -z "$PYTHON_EXE" || ! -x "$PYTHON_EXE" ]]; then
  echo "ERROR: could not resolve Python in Conda env '$CONDA_ENV'" >&2
  exit 2
fi

ARGS=(
  --source-run "$SOURCE_RUN"
  --focused-run "$FOCUSED_RUN"
  --out "$OUT"
  --Kdot "${KDOT:-0.005}"
  --Kmax "${KMAX:-100}"
  --global-n "${GLOBAL_N:-16384}"
  --global-dK "${GLOBAL_DK:-0.5}"
  --n-seeds "${N_SEEDS:-64}"
  --gen1-perturb "${GEN1_PERTURB:-72}"
  --gen2-seeds "${GEN2_SEEDS:-40}"
  --gen2-perturb "${GEN2_PERTURB:-80}"
  --gen3-seeds "${GEN3_SEEDS:-24}"
  --gen3-perturb "${GEN3_PERTURB:-64}"
  --finalists "${FINALISTS:-12}"
  --gen1-dK "${GEN1_DK:-0.075}"
  --gen2-dK "${GEN2_DK:-0.04}"
  --gen3-dK "${GEN3_DK:-0.025}"
  --final-dK "${FINAL_DK:-0.02}"
  --min-prominence "${MIN_PROMINENCE:-2.0}"
  --seed "${SEED:-20260712}"
)
[[ "$RESUME" == "1" ]] && ARGS+=(--resume)
[[ "$SMOKE" == "1" ]] && ARGS+=(--smoke)

printf '=== Expanded EXP-floor peak search ===\n'
printf 'python:       %s\n' "$PYTHON_EXE"
printf 'source_run:   %s\n' "$SOURCE_RUN"
printf 'focused_run:  %s\n' "$FOCUSED_RUN"
printf 'out:          %s\n' "$OUT"
printf 'Kdot:         %s\n' "${KDOT:-0.005}"
printf 'smoke:        %s\n' "$SMOKE"

"$PYTHON_EXE" run_v1_exp_floor_peak_expanded_search.py "${ARGS[@]}"
