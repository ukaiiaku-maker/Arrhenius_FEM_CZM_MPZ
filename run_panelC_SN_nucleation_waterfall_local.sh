#!/usr/bin/env bash
set -euo pipefail

# Figure 1C: project-local Conda environment + restartable V1 S-N waterfall.
# Run from the Fatigue-PF project root.
# No `conda activate` is required.

PROJECT_ROOT="${PROJECT_ROOT:-$PWD}"
ENV_PREFIX="${PANELC_ENV_PREFIX:-$PROJECT_ROOT/.conda-envs/panelC-sn}"
PYTHON_BIN="$ENV_PREFIX/bin/python"
OUT="${OUT:-runs/panelC_SN_nucleation_waterfall}"

ENTROPY_MAG_KB="${ENTROPY_MAG_KB:-0 4 8 12 16 20 24}"
LAMBDA_SH="${LAMBDA_SH:-0 0.15 0.30 0.45 0.60}"
STRESSES="${STRESSES:-150 200 250 300 350 400 450 500 550 600 700 800 900 1050 1200 1400}"
T_K="${T_K:-300}"
R="${R:-0.1}"
FREQUENCY_HZ="${FREQUENCY_HZ:-1000}"
CHI_BACK="${CHI_BACK:-0.60}"
CYCLES_MAX="${CYCLES_MAX:-1e12}"
MAX_BLOCKS="${MAX_BLOCKS:-20000}"
N_PHASE="${N_PHASE:-64}"

# Project-local env creation. The environment lives inside this repository tree
# and does not modify base or any named Conda environment.
if [[ ! -x "$PYTHON_BIN" ]]; then
  if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: conda was not found on PATH." >&2
    exit 2
  fi
  echo "=== Creating project-local Panel C environment ==="
  echo "prefix: $ENV_PREFIX"
  mkdir -p "$(dirname "$ENV_PREFIX")"
  conda create -y \
    --prefix "$ENV_PREFIX" \
    -c conda-forge \
    --strict-channel-priority \
    python=3.13 \
    'numpy<3' \
    'scipy<2' \
    'pandas<4' \
    'matplotlib<4'

  conda list --prefix "$ENV_PREFIX" --explicit \
    > "$PROJECT_ROOT/panelC_local_environment_explicit.txt"
fi

# Prevent packages from ~/.local or other user-site locations from leaking in.
export PYTHONNOUSERSITE=1
export PYTHONPATH="$PROJECT_ROOT${PYTHONPATH:+:$PYTHONPATH}"

cd "$PROJECT_ROOT"

# Safe optional scalar flags; compatible with macOS Bash 3.2 + set -u.
RESUME_FLAG=""
if [[ "${RESUME:-1}" == "1" ]]; then
  RESUME_FLAG="--resume"
fi
PLOT_ONLY_FLAG=""
if [[ "${PLOT_ONLY:-0}" == "1" ]]; then
  PLOT_ONLY_FLAG="--plot-only"
fi
PRINT_FLAG=""
if [[ "${PRINT_EVERY_POINT:-0}" == "1" ]]; then
  PRINT_FLAG="--print-every-point"
fi

printf 'Project root: %s\n' "$PROJECT_ROOT"
printf 'Local env:    %s\n' "$ENV_PREFIX"
printf 'Python:       %s\n' "$PYTHON_BIN"

"$PYTHON_BIN" - <<'PY'
import sys
import numpy
import scipy
import scipy.sparse.linalg
import pandas
import matplotlib
print("python", sys.version.split()[0])
print("numpy", numpy.__version__)
print("scipy", scipy.__version__)
print("pandas", pandas.__version__)
print("matplotlib", matplotlib.__version__)
print("scipy sparse import OK")
PY

"$PYTHON_BIN" -m compileall -q arrhenius_fracture

"$PYTHON_BIN" build_panelC_SN_nucleation_waterfall_3d.py \
  --project-root "$PROJECT_ROOT" \
  --out "$OUT" \
  --entropy-mag-kB "$ENTROPY_MAG_KB" \
  --lambda-sh "$LAMBDA_SH" \
  --stresses-MPa "$STRESSES" \
  --T-K "$T_K" \
  --R "$R" \
  --frequency-Hz "$FREQUENCY_HZ" \
  --chi-back "$CHI_BACK" \
  --cycles-max "$CYCLES_MAX" \
  --max-blocks "$MAX_BLOCKS" \
  --n-phase "$N_PHASE" \
  ${RESUME_FLAG:+$RESUME_FLAG} \
  ${PLOT_ONLY_FLAG:+$PLOT_ONLY_FLAG} \
  ${PRINT_FLAG:+$PRINT_FLAG}
