#!/usr/bin/env bash
set -euo pipefail
CONDA_ENV=${CONDA_ENV:-arrhenius-fem-czm}
SOURCE_RUN=${SOURCE_RUN:-runs/v1_exp_floor_four_class_tuning}
OUT=${OUT:-runs/v1_exp_floor_peak_focused_refinement}
TARGETS=${TARGETS:-exp_floor_four_class_target_curves.csv}

if [[ -n "${PYTHON_BIN:-}" ]]; then
  PYTHON_EXE="$PYTHON_BIN"
elif command -v conda >/dev/null 2>&1; then
  PYTHON_EXE="$(conda run -n "$CONDA_ENV" python -c 'import sys; print(sys.executable)' 2>&1 | tr -d '\r' | awk 'NF {last=$0} END {print last}')"
else
  echo "ERROR: conda not found and PYTHON_BIN not set" >&2; exit 2
fi
[[ -n "$PYTHON_EXE" && -x "$PYTHON_EXE" ]] || { echo "ERROR: no usable Python for '$CONDA_ENV'" >&2; exit 2; }

EXTRA=()
[[ "${RESUME:-1}" == "1" ]] && EXTRA+=(--resume)
[[ "${SMOKE:-0}" == "1" ]] && EXTRA+=(--smoke)

echo "=== focused EXP-floor peak refinement ==="
echo "source run: $SOURCE_RUN"
echo "out:        $OUT"
echo "Kdot:       ${KDOT:-0.005}"
echo "seeds:      ${N_SEEDS:-24}"

"$PYTHON_EXE" run_v1_exp_floor_peak_focused_refinement.py \
  --source-run "$SOURCE_RUN" --targets "$TARGETS" --out "$OUT" \
  --Kdot "${KDOT:-0.005}" --Kmax "${KMAX:-80}" \
  --n-seeds "${N_SEEDS:-24}" \
  --gen1-perturb "${GEN1_PERTURB:-48}" --gen2-seeds "${GEN2_SEEDS:-24}" \
  --gen2-perturb "${GEN2_PERTURB:-64}" --finalists "${FINALISTS:-8}" \
  --gen1-dK "${GEN1_DK:-0.05}" --gen2-dK "${GEN2_DK:-0.025}" --final-dK "${FINAL_DK:-0.02}" \
  --min-prominence "${MIN_PROMINENCE:-1.0}" \
  --rates "${RATES:-0.00125 0.0025 0.005 0.01 0.02}" \
  --seed "${SEED:-20260711}" \
  "${EXTRA[@]}"
