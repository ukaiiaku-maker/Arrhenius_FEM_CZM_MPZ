#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV=${CONDA_ENV:-arrhenius-fem-czm}
SOURCE_RUN=${SOURCE_RUN:-runs/v1_exp_floor_four_class_tuning}
PEAK_RUN=${PEAK_RUN:-runs/v1_exp_floor_peak_expanded_search}
OUT=${OUT:-runs/v1_exp_floor_weakT_convex_refinement}
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
  --peak-run "$PEAK_RUN"
  --out "$OUT"
  --Kdot "${KDOT:-0.005}"
  --Kmax "${KMAX:-80}"
  --tau-K "${TAU_K:-450}"
  --convex-tol "${CONVEX_TOL:-0.002}"
  --n-seeds "${N_SEEDS:-64}"
  --gen1-perturb "${GEN1_PERTURB:-64}"
  --gen2-seeds "${GEN2_SEEDS:-36}"
  --gen2-perturb "${GEN2_PERTURB:-72}"
  --gen3-seeds "${GEN3_SEEDS:-20}"
  --gen3-perturb "${GEN3_PERTURB:-64}"
  --finalists "${FINALISTS:-12}"
  --gen1-dK "${GEN1_DK:-0.075}"
  --gen2-dK "${GEN2_DK:-0.04}"
  --gen3-dK "${GEN3_DK:-0.025}"
  --final-dK "${FINAL_DK:-0.02}"
  --seed "${SEED:-20260713}"
)
[[ "$RESUME" == "1" ]] && ARGS+=(--resume)
[[ "$SMOKE" == "1" ]] && ARGS+=(--smoke)

printf '=== V1 EXP-floor weak-T convex refinement ===\n'
printf 'python:       %s\n' "$PYTHON_EXE"
printf 'source_run:   %s\n' "$SOURCE_RUN"
printf 'peak_run:     %s\n' "$PEAK_RUN"
printf 'out:          %s\n' "$OUT"
printf 'Kdot:         %s\n' "${KDOT:-0.005}"
printf 'tau_K:        %s\n' "${TAU_K:-450}"
printf 'smoke:        %s\n' "$SMOKE"

"$PYTHON_EXE" run_v1_exp_floor_weakT_convex_refinement.py "${ARGS[@]}"
