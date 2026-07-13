#!/usr/bin/env bash
set -euo pipefail

ENV_NAME="${CONDA_ENV:-arrhenius-fem-czm}"
PYTHON_VERSION="${PYTHON_VERSION:-3.12}"
CHANNEL="${CONDA_CHANNEL:-conda-forge}"

if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda is not on PATH." >&2
  exit 2
fi

if conda env list | awk '{print $1}' | grep -Fxq "$ENV_NAME"; then
  echo "Environment '$ENV_NAME' already exists; leaving packages unchanged."
else
  echo "Creating isolated environment '$ENV_NAME' with Python $PYTHON_VERSION ..."
  conda create -y -n "$ENV_NAME" -c "$CHANNEL" --strict-channel-priority \
    "python=$PYTHON_VERSION" numpy scipy pandas matplotlib pip
fi

ENV_PYTHON="$(conda run -n "$ENV_NAME" python -c 'import sys; print(sys.executable)' 2>&1 \
  | tr -d '\r' \
  | awk 'NF {last=$0} END {print last}')"

if [[ -z "$ENV_PYTHON" || ! -x "$ENV_PYTHON" ]]; then
  echo "ERROR: could not resolve Python for Conda environment '$ENV_NAME'." >&2
  exit 2
fi

echo
echo "Installing the local package in editable mode without changing dependency resolution ..."
"$ENV_PYTHON" -m pip install -e . --no-deps

echo
echo "Preflight for '$ENV_NAME':"
"$ENV_PYTHON" - <<'PY'
import sys
import numpy
import scipy
import scipy.sparse.linalg
import pandas
import matplotlib
import arrhenius_fracture
from scipy.sparse.linalg import spsolve

print('python    :', sys.executable)
print('numpy     :', numpy.__version__)
print('scipy     :', scipy.__version__)
print('pandas    :', pandas.__version__)
print('matplotlib:', matplotlib.__version__)
print('arrhenius_fracture import: OK')
print('scipy sparse solver import: OK')
PY

echo
echo "Environment ready. The sweep runner defaults to CONDA_ENV=$ENV_NAME."
