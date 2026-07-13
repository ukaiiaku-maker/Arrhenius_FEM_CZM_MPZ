#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"
ENV_PREFIX="${PANELS_CD_ENV_PREFIX:-$ROOT_DIR/.conda-envs/panels-cd-entropy}"

if [[ -x "$ENV_PREFIX/bin/python" ]]; then
  echo "Local environment already exists: $ENV_PREFIX"
else
  mkdir -p "$(dirname "$ENV_PREFIX")"
  echo "Creating local Conda environment: $ENV_PREFIX"
  conda create -y -p "$ENV_PREFIX" \
    -c conda-forge \
    --strict-channel-priority \
    python=3.13 \
    'numpy<3' \
    'scipy<2' \
    'pandas<4' \
    'matplotlib<4'
fi

PY="$ENV_PREFIX/bin/python"
export PYTHONNOUSERSITE=1

"$PY" - <<'PY'
import sys
import numpy
import scipy
import scipy.optimize
import scipy.special
import pandas
import matplotlib
print("python:", sys.executable)
print("numpy:", numpy.__version__)
print("scipy:", scipy.__version__)
print("pandas:", pandas.__version__)
print("matplotlib:", matplotlib.__version__)
print("Panels C/D numerical stack OK")
PY

conda list -p "$ENV_PREFIX" --explicit > "$ROOT_DIR/panels_CD_local_environment_explicit.txt"
echo "Wrote exact environment export: $ROOT_DIR/panels_CD_local_environment_explicit.txt"
echo "Local environment ready: $ENV_PREFIX"
