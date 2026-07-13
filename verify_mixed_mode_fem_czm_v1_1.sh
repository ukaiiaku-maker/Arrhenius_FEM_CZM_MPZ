#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${PROJECT_ROOT:-$PWD}"
CONDA_ENV="${CONDA_ENV:-arrhenius-fem-czm}"

cd "$PROJECT_ROOT"

if ! command -v conda >/dev/null 2>&1; then
  CONDA_SH="/opt/homebrew/Caskroom/miniconda/base/etc/profile.d/conda.sh"
  if [[ -f "$CONDA_SH" ]]; then
    # shellcheck disable=SC1090
    source "$CONDA_SH"
  else
    echo "ERROR: conda is not available on PATH and $CONDA_SH was not found." >&2
    exit 2
  fi
fi

if ! conda env list | awk '{print $1}' | grep -Fxq "$CONDA_ENV"; then
  echo "ERROR: Conda environment '$CONDA_ENV' does not exist." >&2
  echo "Available environments:" >&2
  conda env list >&2
  exit 2
fi

if [[ ! -f verify_mixed_mode_fem_czm_v1.py ]]; then
  echo "ERROR: verify_mixed_mode_fem_czm_v1.py is not in project root: $PROJECT_ROOT" >&2
  exit 2
fi

PYTHON_EXE="$(conda run -n "$CONDA_ENV" python -c 'import sys; print(sys.executable)' | awk 'NF{a=$0}END{print a}')"

echo "MIXED_MODE_V1_1 project root: $PROJECT_ROOT"
echo "MIXED_MODE_V1_1 conda env:    $CONDA_ENV"
echo "MIXED_MODE_V1_1 python:       $PYTHON_EXE"

conda run -n "$CONDA_ENV" python -c \
  'import sys, numpy, scipy; print("Python", sys.version.split()[0]); print("NumPy", numpy.__version__); print("SciPy", scipy.__version__)'

conda run -n "$CONDA_ENV" python verify_mixed_mode_fem_czm_v1.py
