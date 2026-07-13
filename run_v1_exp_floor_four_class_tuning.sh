#!/usr/bin/env bash
set -euo pipefail
CONDA_ENV=${CONDA_ENV:-arrhenius-fem-czm}
OUT=${OUT:-runs/v1_exp_floor_four_class_tuning}
DESIGN=${DESIGN:-exp_floor_emission_design_seed.csv}
TARGETS=${TARGETS:-exp_floor_four_class_target_curves.csv}
VALIDATION_2D=${VALIDATION_2D:-exp_floor_four_class_2D_validation.csv}

if [[ -n "${PYTHON_BIN:-}" ]]; then
  PYTHON_EXE="$PYTHON_BIN"
elif command -v conda >/dev/null 2>&1; then
  PYTHON_EXE="$(conda run -n "$CONDA_ENV" python -c 'import sys; print(sys.executable)' 2>&1 | tr -d '\r' | awk 'NF {last=$0} END {print last}')"
else
  echo "ERROR: conda not found and PYTHON_BIN not set" >&2; exit 2
fi
[[ -n "$PYTHON_EXE" && -x "$PYTHON_EXE" ]] || { echo "ERROR: no usable Python for '$CONDA_ENV'" >&2; exit 2; }

"$PYTHON_EXE" - <<'PY'
import numpy,pandas,scipy,arrhenius_fracture
print(f"preflight OK: numpy {numpy.__version__}, pandas {pandas.__version__}, scipy {scipy.__version__}")
PY

echo "=== fully EXP-floor four-class V1 tuning ==="
echo "python:      $PYTHON_EXE"
echo "design:      $DESIGN"
echo "targets:     $TARGETS"
echo "validation:  $VALIDATION_2D"
echo "out:         $OUT"
echo "contexts:    ${N_CONTEXTS:-256}"

EXTRA=()
[[ "${RESUME:-1}" == "1" ]] && EXTRA+=(--resume)
[[ "${SMOKE:-0}" == "1" ]] && EXTRA+=(--smoke)
"$PYTHON_EXE" run_v1_exp_floor_four_class_tuning.py \
  --design "$DESIGN" --targets "$TARGETS" --validation-2d "$VALIDATION_2D" --out "$OUT" \
  --n-contexts "${N_CONTEXTS:-256}" --seed "${SEED:-20260710}" \
  --Kmax "${KMAX:-80}" --Kdot "${KDOT:-0.005}" \
  --broad-dK "${BROAD_DK:-0.25}" --local-dK "${LOCAL_DK:-0.05}" --final-dK "${FINAL_DK:-0.02}" \
  --broad-top-per-class "${BROAD_TOP_PER_CLASS:-24}" \
  --local-perturb-per-seed "${LOCAL_PERTURB_PER_SEED:-16}" \
  --final-top-per-class "${FINAL_TOP_PER_CLASS:-6}" \
  "${EXTRA[@]}"
