#!/usr/bin/env bash
set -euo pipefail

# Figure 1C: fixed-entropy Lambda_sh sweep using a project-local Conda prefix.
# Run from the Fatigue-PF project root. No `conda activate` is required.

PROJECT_ROOT="${PROJECT_ROOT:-$PWD}"
ENV_PREFIX="${PANELC_ENV_PREFIX:-$PROJECT_ROOT/.conda-envs/panelC-sn}"
PYTHON_BIN="$ENV_PREFIX/bin/python"
OUT="${OUT:-runs/panelC_SN_lambda_axis}"

ENTROPY_MAG_KB="${ENTROPY_MAG_KB:-0}"
LAMBDA_SH="${LAMBDA_SH:-0 0.10 0.20 0.30 0.40 0.50 0.60 0.80 1.00}"
STRESSES="${STRESSES:-250 350 450 550 650 750 850 1000 1150 1300}"
T_K="${T_K:-300}"
R="${R:-0.1}"
FREQUENCY_HZ="${FREQUENCY_HZ:-1000}"
CHI_BACK="${CHI_BACK:-0.60}"
CYCLES_MAX="${CYCLES_MAX:-1e12}"
MAX_BLOCKS="${MAX_BLOCKS:-10000}"
N_PHASE="${N_PHASE:-64}"

cd "$PROJECT_ROOT"

# Create the repository-local environment if it does not exist.
if [[ ! -x "$PYTHON_BIN" ]]; then
  if [[ -x "$PROJECT_ROOT/setup_panelC_local_env.sh" ]]; then
    PROJECT_ROOT="$PROJECT_ROOT" PANELC_ENV_PREFIX="$ENV_PREFIX" \
      bash "$PROJECT_ROOT/setup_panelC_local_env.sh"
  else
    echo "ERROR: local env is missing and setup_panelC_local_env.sh was not found." >&2
    exit 2
  fi
fi

# Use the corrected builder if the canonical file has not yet been replaced.
BUILDER="$PROJECT_ROOT/build_panelC_SN_nucleation_waterfall_3d.py"
if [[ -f "$PROJECT_ROOT/build_panelC_SN_nucleation_waterfall_3d_fixed.py" ]]; then
  BUILDER="$PROJECT_ROOT/build_panelC_SN_nucleation_waterfall_3d_fixed.py"
fi
PLOTTER="$PROJECT_ROOT/plot_panelC_SN_lambda_axis.py"

if [[ ! -f "$BUILDER" ]]; then
  echo "ERROR: Panel C builder not found: $BUILDER" >&2
  exit 2
fi
if [[ ! -f "$PLOTTER" ]]; then
  echo "ERROR: Lambda-axis plotter not found: $PLOTTER" >&2
  exit 2
fi

export PYTHONNOUSERSITE=1
export PYTHONPATH="$PROJECT_ROOT${PYTHONPATH:+:$PYTHONPATH}"

RESUME_FLAG=""
if [[ "${RESUME:-1}" == "1" ]]; then
  RESUME_FLAG="--resume"
fi
PRINT_FLAG=""
if [[ "${PRINT_EVERY_POINT:-0}" == "1" ]]; then
  PRINT_FLAG="--print-every-point"
fi

printf 'Project root: %s\n' "$PROJECT_ROOT"
printf 'Local env:    %s\n' "$ENV_PREFIX"
printf 'Python:       %s\n' "$PYTHON_BIN"
printf 'Output:       %s\n' "$OUT"

"$PYTHON_BIN" - <<'PY'
import sys, numpy, scipy, pandas, matplotlib
import scipy.sparse.linalg
print("python", sys.version.split()[0])
print("numpy", numpy.__version__)
print("scipy", scipy.__version__)
print("pandas", pandas.__version__)
print("matplotlib", matplotlib.__version__)
print("scipy sparse import OK")
PY

"$PYTHON_BIN" -m compileall -q arrhenius_fracture

"$PYTHON_BIN" "$BUILDER" \
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
  ${PRINT_FLAG:+$PRINT_FLAG}

RAW_CSV="$OUT/panelC_SN_nucleation_waterfall_raw.csv"
"$PYTHON_BIN" "$PLOTTER" \
  --raw-csv "$RAW_CSV" \
  --outdir "$OUT" \
  --entropy-mag-kB "$ENTROPY_MAG_KB" \
  --cycles-max "$CYCLES_MAX"

printf '\nPanel C Lambda-axis workflow complete.\n'
printf 'Main figure: %s/panelC_SN_lambda_waterfall_3d.png\n' "$OUT"
