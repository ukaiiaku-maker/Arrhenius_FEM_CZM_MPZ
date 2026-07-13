#!/usr/bin/env bash
set -euo pipefail

# Create a dedicated, isolated environment for Figure 1 V1 workflows.
# Default strategy: clone the already-tested fatigue-pf environment so no
# package resolution or changes to base are required.

SOURCE_ENV="${SOURCE_ENV:-fatigue-pf}"
TARGET_ENV="${TARGET_ENV:-fatigue-pf-fig1}"

if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda is not on PATH." >&2
  exit 2
fi

if conda env list | awk '{print $1}' | grep -Fxq "$TARGET_ENV"; then
  echo "Environment '$TARGET_ENV' already exists; leaving it unchanged."
else
  if ! conda env list | awk '{print $1}' | grep -Fxq "$SOURCE_ENV"; then
    echo "ERROR: source environment '$SOURCE_ENV' does not exist." >&2
    echo "Expected the previously tested clean environment named '$SOURCE_ENV'." >&2
    exit 2
  fi

  echo "Cloning '$SOURCE_ENV' -> '$TARGET_ENV' ..."
  conda create -y -n "$TARGET_ENV" --clone "$SOURCE_ENV"
fi

echo
echo "Preflight for '$TARGET_ENV':"
conda run -n "$TARGET_ENV" python - <<'PY'
import sys
import numpy
import scipy
import scipy.sparse.linalg
import pandas
import matplotlib
from scipy.sparse.linalg import spsolve

print('python    :', sys.executable)
print('numpy     :', numpy.__version__)
print('scipy     :', scipy.__version__)
print('pandas    :', pandas.__version__)
print('matplotlib:', matplotlib.__version__)
print('scipy sparse solver import: OK')
PY

echo
echo "Created/verified isolated Figure 1 environment: $TARGET_ENV"
echo "The base environment was not modified."
