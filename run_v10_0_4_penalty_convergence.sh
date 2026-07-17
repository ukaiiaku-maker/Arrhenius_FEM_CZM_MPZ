#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV=${CONDA_ENV:-arrhenius-fem-czm}
PYTHON_BIN=${PYTHON_BIN:-}
OUTROOT=${OUTROOT:-runs/v10_0_4_penalty_convergence_weakT_700K_v1}
MATERIAL=${MATERIAL:-weakT}
T_K=${T_K:-700}
THETA_DEG=${THETA_DEG:-45}
NORMAL_PENALTIES=${NORMAL_PENALTIES:-"5e17 1e18 2e18"}
TANGENT_PENALTIES=${TANGENT_PENALTIES:-"5e17 1e18 2e18"}
REFERENCE_PENALTY=${REFERENCE_PENALTY:-1e18}
NO_PLOTS=${NO_PLOTS:-1}

if [[ -z "$PYTHON_BIN" ]]; then
  if [[ "${CONDA_DEFAULT_ENV:-}" == "$CONDA_ENV" ]]; then
    PYTHON_BIN=$(command -v python)
  else
    PYTHON_BIN=$(conda run -n "$CONDA_ENV" python -c 'import sys; print(sys.executable)' | tail -n 1)
  fi
fi

if [[ -e "$OUTROOT" ]]; then
  echo "ERROR: output path already exists: $OUTROOT"
  exit 1
fi

CONDA_ENV="$CONDA_ENV" PYTHON_BIN="$PYTHON_BIN" \
  bash run_v10_0_3_integration_tests.sh

"$PYTHON_BIN" prepare_v10_0_4_penalty_plan.py \
  --normal-penalties "$NORMAL_PENALTIES" \
  --tangent-penalties "$TANGENT_PENALTIES" \
  --reference-penalty "$REFERENCE_PENALTY" \
  --out "$OUTROOT"

PLAN="$OUTROOT/v10_0_4_penalty_plan.csv"
tail -n +2 "$PLAN" | while IFS=, read -r case_id axis normal_penalty tangent_penalty; do
  echo "========================================================================"
  echo "v10.0.4 penalty case: $case_id"
  echo "axis=$axis normal=$normal_penalty tangent=$tangent_penalty"
  echo "========================================================================"
  CONDA_ENV="$CONDA_ENV" PYTHON_BIN="$PYTHON_BIN" \
  MATERIAL="$MATERIAL" T_K="$T_K" THETA_DEG="$THETA_DEG" \
  CZM_PENALTY_NORMAL="$normal_penalty" \
  CZM_PENALTY_TANGENT="$tangent_penalty" \
  RUN_TESTS=0 NO_PLOTS="$NO_PLOTS" \
  OUTROOT="$OUTROOT/$case_id" \
    bash run_v10_0_4_single_case.sh
done

"$PYTHON_BIN" analyze_v10_0_4_penalty_convergence.py \
  "$OUTROOT" --reference-penalty "$REFERENCE_PENALTY"

cat <<EOF
V10.0.4 PENALTY CONVERGENCE GATE PASSED
out=$OUTROOT
material=$MATERIAL
T_K=$T_K
parameterization_matrix_authorized=true
long_growth_authorized=false
EOF
