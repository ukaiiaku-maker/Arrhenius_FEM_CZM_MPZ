#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV=${CONDA_ENV:-arrhenius-fem-czm}
PYTHON_BIN=${PYTHON_BIN:-}
OUTROOT=${OUTROOT:-}

if [[ -z "$OUTROOT" ]]; then
  echo "ERROR: set OUTROOT to a completed v10.0.5 output directory"
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

"$PYTHON_BIN" -m pip install -e . --no-deps
"$PYTHON_BIN" normalize_v10_0_5_1_slip_trace_reporting.py "$OUTROOT"

cat <<EOF
V10.0.5.1 EXISTING-OUTPUT NORMALIZATION COMPLETE
out=$OUTROOT
No FEM solve was launched.
No mechanics, kinetics, barriers, sources, or material parameters were changed.
EOF
