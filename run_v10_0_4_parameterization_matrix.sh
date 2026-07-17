#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV=${CONDA_ENV:-arrhenius-fem-czm}
PYTHON_BIN=${PYTHON_BIN:-}
PENALTY_CERT=${PENALTY_CERT:?Set PENALTY_CERT to v10_0_4_penalty_convergence_certification.json}
OUTROOT=${OUTROOT:-runs/v10_0_4_parameterization_matrix_3class_3temp_v1}
MATERIALS=${MATERIALS:-"ceramic weakT DBTT"}
TEMPERATURES=${TEMPERATURES:-"300 700 1100"}
THETA_DEG=${THETA_DEG:-45}
NO_PLOTS=${NO_PLOTS:-1}

if [[ -z "$PYTHON_BIN" ]]; then
  if [[ "${CONDA_DEFAULT_ENV:-}" == "$CONDA_ENV" ]]; then
    PYTHON_BIN=$(command -v python)
  else
    PYTHON_BIN=$(conda run -n "$CONDA_ENV" python -c 'import sys; print(sys.executable)' | tail -n 1)
  fi
fi

if [[ ! -f "$PENALTY_CERT" ]]; then
  echo "ERROR: penalty certification not found: $PENALTY_CERT"
  exit 1
fi
if [[ -e "$OUTROOT" ]]; then
  echo "ERROR: output path already exists: $OUTROOT"
  exit 1
fi

CONDA_ENV="$CONDA_ENV" PYTHON_BIN="$PYTHON_BIN" \
  bash run_v10_0_3_integration_tests.sh

"$PYTHON_BIN" prepare_v10_0_4_parameter_matrix.py \
  --penalty-certification "$PENALTY_CERT" \
  --materials "$MATERIALS" \
  --temperatures "$TEMPERATURES" \
  --out "$OUTROOT"

PLAN="$OUTROOT/v10_0_4_parameter_matrix_plan.csv"
tail -n +2 "$PLAN" | while IFS=, read -r case_id material temperature normal_penalty tangent_penalty; do
  echo "========================================================================"
  echo "v10.0.4 parameterization case: $case_id"
  echo "material=$material T_K=$temperature"
  echo "normal_penalty=$normal_penalty tangent_penalty=$tangent_penalty"
  echo "========================================================================"
  CONDA_ENV="$CONDA_ENV" PYTHON_BIN="$PYTHON_BIN" \
  MATERIAL="$material" T_K="$temperature" THETA_DEG="$THETA_DEG" \
  CZM_PENALTY_NORMAL="$normal_penalty" \
  CZM_PENALTY_TANGENT="$tangent_penalty" \
  RUN_TESTS=0 NO_PLOTS="$NO_PLOTS" \
  OUTROOT="$OUTROOT/$case_id" \
    bash run_v10_0_4_single_case.sh
done

"$PYTHON_BIN" analyze_v10_0_4_parameter_matrix.py "$OUTROOT"

cat <<EOF
V10.0.4 THREE-PARAMETERIZATION VALIDATION GATE PASSED
out=$OUTROOT
materials=$MATERIALS
temperatures=$TEMPERATURES
short_growth_authorized=true
long_growth_authorized=false
temperature_sweep_authorized=false
EOF
