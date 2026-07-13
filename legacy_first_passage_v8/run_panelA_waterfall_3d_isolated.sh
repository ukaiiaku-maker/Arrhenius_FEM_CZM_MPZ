#!/usr/bin/env bash
set -euo pipefail

# Run the Panel A waterfall inside a dedicated Conda environment without
# activating it and without modifying the caller's shell or base environment.

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONDA_ENV="${CONDA_ENV:-fatigue-pf-fig1}"
PROJECT_ROOT="${PROJECT_ROOT:-$PWD}"
OUTROOT="${OUTROOT:-runs/panelA_waterfall_3d_dense}"
SCRIPT="${SCRIPT:-$PROJECT_ROOT/build_panelA_waterfall_3d.py}"

N_LINES="${N_LINES:-20}"
H0_MIN_EV="${H0_MIN_EV:-2.6}"
H0_MAX_EV="${H0_MAX_EV:-6.0}"
T_MIN="${T_MIN:-300}"
T_MAX="${T_MAX:-1200}"
T_STEP="${T_STEP:-10}"
KDOT="${KDOT:-0.02}"
KMAX="${KMAX:-100}"
DT="${DT:-1.0}"
VIEW_ELEV="${VIEW_ELEV:-23}"
VIEW_AZIM="${VIEW_AZIM:--63}"

if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda is not on PATH." >&2
  exit 2
fi

if ! conda env list | awk '{print $1}' | grep -Fxq "$CONDA_ENV"; then
  echo "ERROR: Conda environment '$CONDA_ENV' does not exist." >&2
  echo "Create it with: bash $HERE/setup_fig1_isolated_env.sh" >&2
  exit 2
fi

if [[ ! -f "$SCRIPT" ]]; then
  echo "ERROR: builder not found: $SCRIPT" >&2
  exit 2
fi

if [[ ! -d "$PROJECT_ROOT/arrhenius_fracture" ]]; then
  echo "ERROR: arrhenius_fracture package not found under PROJECT_ROOT=$PROJECT_ROOT" >&2
  exit 2
fi

mkdir -p "$OUTROOT"
LOG="$OUTROOT/panelA_waterfall_run.log"

# Verify the exact environment before expensive work.
echo "=== isolated environment preflight: $CONDA_ENV ===" | tee "$LOG"
conda run -n "$CONDA_ENV" python - <<'PY' | tee -a "$LOG"
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

# Verify project import from the current checkout, not an installed stale copy.
(
  cd "$PROJECT_ROOT"
  conda run -n "$CONDA_ENV" python - <<'PY'
import pathlib
import arrhenius_fracture
print('arrhenius_fracture import: OK')
print('package path:', pathlib.Path(arrhenius_fracture.__file__).resolve())
PY
) | tee -a "$LOG"

echo "=== running Panel A waterfall ===" | tee -a "$LOG"
(
  cd "$PROJECT_ROOT"
  conda run -n "$CONDA_ENV" python "$SCRIPT" \
    --project-root "$PROJECT_ROOT" \
    --out "$OUTROOT" \
    --n-lines "$N_LINES" \
    --H0-min-eV "$H0_MIN_EV" \
    --H0-max-eV "$H0_MAX_EV" \
    --T-min "$T_MIN" \
    --T-max "$T_MAX" \
    --T-step "$T_STEP" \
    --Kdot-MPa-sqrtm-per-s "$KDOT" \
    --Kmax-MPa-sqrtm "$KMAX" \
    --dt "$DT" \
    --view-elev "$VIEW_ELEV" \
    --view-azim "$VIEW_AZIM"
) 2>&1 | tee -a "$LOG"

echo
echo "Panel A waterfall complete."
echo "Environment: $CONDA_ENV"
echo "Outputs:     $OUTROOT"
echo "Log:         $LOG"
