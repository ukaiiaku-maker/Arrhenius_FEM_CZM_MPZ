#!/usr/bin/env bash
set -euo pipefail

# Figure 1C: systematic V1 S-N crack-initiation waterfall.
# Run from the Fatigue-PF project root after:
#   conda activate fatigue-pf

PYTHON_BIN="${PYTHON_BIN:-python}"
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

# Safe optional flags for macOS's older Bash as well as modern Bash.
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

printf 'Python: %s\n' "$(command -v "$PYTHON_BIN")"
"$PYTHON_BIN" - <<'PY'
import numpy
import pandas
import matplotlib
print("numpy", numpy.__version__)
print("pandas", pandas.__version__)
print("matplotlib", matplotlib.__version__)
PY

"$PYTHON_BIN" build_panelC_SN_nucleation_waterfall_3d.py \
  --project-root . \
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
