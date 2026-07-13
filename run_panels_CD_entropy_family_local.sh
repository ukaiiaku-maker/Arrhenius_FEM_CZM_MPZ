#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"
ENV_PREFIX="${PANELS_CD_ENV_PREFIX:-$ROOT_DIR/.conda-envs/panels-cd-entropy}"
SETUP="$ROOT_DIR/setup_panels_CD_local_env.sh"
DRIVER="$ROOT_DIR/build_panels_CD_entropy_family.py"
SHARED="$ROOT_DIR/shared_entropy_family.py"

[[ -f "$DRIVER" ]] || { echo "ERROR: missing $DRIVER" >&2; exit 2; }
[[ -f "$SHARED" ]] || { echo "ERROR: missing $SHARED" >&2; exit 2; }
[[ -f "$SETUP" ]] || { echo "ERROR: missing $SETUP" >&2; exit 2; }

if [[ ! -x "$ENV_PREFIX/bin/python" ]]; then
  bash "$SETUP"
fi

PY="$ENV_PREFIX/bin/python"
export PYTHONNOUSERSITE=1
export PYTHONPATH="$ROOT_DIR${PYTHONPATH:+:$PYTHONPATH}"

echo "Python: $PY"
"$PY" - <<'PY'
import numpy, scipy, pandas, matplotlib
from scipy.optimize import brentq
from scipy.special import gammainc
print("numpy", numpy.__version__)
print("scipy", scipy.__version__)
print("pandas", pandas.__version__)
print("matplotlib", matplotlib.__version__)
print("preflight OK")
PY

RESUME_ARG=""
if [[ "${RESUME:-1}" == "1" || "${RESUME:-1}" == "true" ]]; then
  RESUME_ARG="--resume"
fi

"$PY" "$DRIVER" \
  --out "${OUT:-runs/panels_CD_entropy_family}" \
  --A-T-kB "${A_T_KB:-0 2 4}" \
  --A-sigma-kB "${A_SIGMA_KB:-0 2 4 6 8}" \
  --T0-K "${T0_K:-300}" \
  --T-S-K "${T_S_K:-400}" \
  --T-gate-power "${T_GATE_POWER:-4}" \
  --sigma-S-GPa "${SIGMA_S_GPA:-3}" \
  --sigma-gate-power "${SIGMA_GATE_POWER:-1}" \
  --S-min-kB "${S_MIN_KB:--40}" \
  --H0-eV "${H0_EV:-0.8}" \
  --sigma0-H-GPa "${SIGMA0_H_GPA:-2.5}" \
  --v0-b3 "${V0_B3:-0.6}" \
  --sigma0-v-GPa "${SIGMA0_V_GPA:-2.5}" \
  --b-m "${B_M:-2.74e-10}" \
  --nu0-s "${NU0_S:-1e11}" \
  --C-T-K "${C_T_K:-300}" \
  --C-stresses-MPa "${C_STRESSES_MPA:-100 150 200 250 300 350 400 450 500 600 700 800}" \
  --R "${R_RATIO:-0.1}" \
  --frequency-Hz "${FREQUENCY_HZ:-1000}" \
  --Kt "${KT:-3}" \
  --cycles-max "${CYCLES_MAX:-1e12}" \
  --n-phase "${N_PHASE:-128}" \
  --multihit-m "${MULTIHIT_M:-3}" \
  --multihit-tau-s "${MULTIHIT_TAU_S:-1e-6}" \
  --D-temperatures-K "${D_TEMPERATURES_K:-250 300 350 400 450 500 550 600 650 700 750 800 850 900 950 1000 1050 1100 1150 1200}" \
  --yield-strain-rate-s "${YIELD_STRAIN_RATE_S:-1e-4}" \
  --yield-event-strain "${YIELD_EVENT_STRAIN:-1e-5}" \
  --yield-sigma-max-GPa "${YIELD_SIGMA_MAX_GPA:-100}" \
  ${RESUME_ARG:+$RESUME_ARG}
