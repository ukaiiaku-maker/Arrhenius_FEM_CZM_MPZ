#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV=${CONDA_ENV:-arrhenius-fem-czm}
PYTHON_BIN=${PYTHON_BIN:-}
OUTROOT=${OUTROOT:-}

if [[ -z "$OUTROOT" ]]; then
  echo "ERROR: set OUTROOT to a certified v10.0.3 run directory"
  exit 1
fi
if [[ ! -d "$OUTROOT" ]]; then
  echo "ERROR: output directory does not exist: $OUTROOT"
  exit 1
fi

if [[ -z "$PYTHON_BIN" ]]; then
  if [[ "${CONDA_DEFAULT_ENV:-}" == "$CONDA_ENV" ]]; then
    PYTHON_BIN=$(command -v python)
  else
    PYTHON_BIN=$(conda run -n "$CONDA_ENV" python -c 'import sys; print(sys.executable)' | tail -n 1)
  fi
fi

CONDA_ENV="$CONDA_ENV" PYTHON_BIN="$PYTHON_BIN" \
  bash run_v10_0_3_integration_tests.sh

"$PYTHON_BIN" audit_v10_0_3_progressive_integration.py \
  "$OUTROOT" --target-um 5

"$PYTHON_BIN" normalize_v10_0_3_1_reporting.py "$OUTROOT"

cat <<EOF
V10.0.3.1 EXISTING OUTPUT REPORTING NORMALIZED
out=$OUTROOT
No FEM solve was launched.
No mechanics, kinetics, cohesive state, or material parameters were changed.
EOF
