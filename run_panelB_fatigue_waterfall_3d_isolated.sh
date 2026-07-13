#!/usr/bin/env bash
set -euo pipefail

# Isolated Figure 1 Panel B fatigue-waterfall workflow.
# Runs through a named conda environment without activating or modifying base.

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONDA_ENV="${CONDA_ENV:-fatigue-pf-fig1}"
PROJECT_ROOT="${PROJECT_ROOT:-$PWD}"
OUTROOT="${OUTROOT:-runs/panelB_fatigue_waterfall_3d}"
SCRIPT="${SCRIPT:-$HERE/build_panelB_fatigue_waterfall_3d.py}"

N_LINES="${N_LINES:-10}"
KMAX_GRID="${KMAX_GRID:-6 8 10 12 14 16 18 20 22 24}"
T_K="${T_K:-300}"
R_RATIO="${R_RATIO:-0.1}"
FREQUENCY_HZ="${FREQUENCY_HZ:-1000}"
CYCLES_MAX="${CYCLES_MAX:-1e12}"
MAX_BLOCKS="${MAX_BLOCKS:-20000}"
N_PHASE="${N_PHASE:-96}"
TARGET_DB="${TARGET_DB:-0.02}"
TARGET_DN_STORE="${TARGET_DN_STORE:-0.01}"
MIN_BLOCK_CYCLES="${MIN_BLOCK_CYCLES:-1e-6}"
VIEW_ELEV="${VIEW_ELEV:-24}"
VIEW_AZIM="${VIEW_AZIM:--62}"
PRINT_EVERY_CASE="${PRINT_EVERY_CASE:-0}"

if [[ ! -f "$SCRIPT" ]]; then
  echo "ERROR: workflow script not found: $SCRIPT" >&2
  exit 2
fi

if ! command -v conda >/dev/null 2>&1; then
  echo "ERROR: conda not found on PATH" >&2
  exit 2
fi

echo "=== isolated environment preflight: $CONDA_ENV ==="
conda run -n "$CONDA_ENV" python - <<'PY'
import sys
import numpy, scipy, pandas, matplotlib
from scipy.sparse.linalg import spsolve
print('python    :', sys.executable)
print('numpy     :', numpy.__version__)
print('scipy     :', scipy.__version__)
print('pandas    :', pandas.__version__)
print('matplotlib:', matplotlib.__version__)
print('scipy sparse solver import: OK')
PY

# Verify the current checkout, not a stale installed package.
PROJECT_ROOT="$PROJECT_ROOT" conda run -n "$CONDA_ENV" python - <<'PY'
import os, sys
root = os.path.abspath(os.environ['PROJECT_ROOT'])
sys.path.insert(0, root)
import arrhenius_fracture
print('arrhenius_fracture import: OK')
print('package path:', arrhenius_fracture.__file__)
if not os.path.abspath(arrhenius_fracture.__file__).startswith(root):
    raise SystemExit('ERROR: imported arrhenius_fracture is outside PROJECT_ROOT')
PY

EXTRA_ARGS=()
if [[ "$PRINT_EVERY_CASE" == "1" ]]; then
  EXTRA_ARGS+=(--print-every-case)
fi

conda run -n "$CONDA_ENV" python "$SCRIPT" \
  --project-root "$PROJECT_ROOT" \
  --out "$OUTROOT" \
  --n-lines "$N_LINES" \
  --Kmax-grid "$KMAX_GRID" \
  --T-K "$T_K" \
  --R "$R_RATIO" \
  --frequency-Hz "$FREQUENCY_HZ" \
  --cycles-max "$CYCLES_MAX" \
  --max-blocks "$MAX_BLOCKS" \
  --n-phase "$N_PHASE" \
  --target-dB "$TARGET_DB" \
  --target-dN-store "$TARGET_DN_STORE" \
  --min-block-cycles "$MIN_BLOCK_CYCLES" \
  --view-elev "$VIEW_ELEV" \
  --view-azim "$VIEW_AZIM" \
  "${EXTRA_ARGS[@]}"

echo
echo "Panel B fatigue waterfall complete."
echo "Outputs: $OUTROOT"
