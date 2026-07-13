#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV="${CONDA_ENV:-arrhenius-stateful-pd}"
PYTHON_VERSION="${PYTHON_VERSION:-3.12}"

if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda is not available on PATH" >&2
  exit 2
fi

if conda env list | awk '{print $1}' | grep -Fxq "$CONDA_ENV"; then
  echo "Conda environment '$CONDA_ENV' already exists."
else
  conda create -y -n "$CONDA_ENV" "python=$PYTHON_VERSION" numpy scipy matplotlib pandas
fi

conda run -n "$CONDA_ENV" python -m compileall -q arrhenius_fracture
conda run -n "$CONDA_ENV" python -m unittest tests.test_stateful_pd_core

echo "Environment ready: $CONDA_ENV"
