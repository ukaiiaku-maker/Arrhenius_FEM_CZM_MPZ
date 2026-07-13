#!/usr/bin/env bash
set -euo pipefail

# Figure 1, Panel A systematic 3-D waterfall workflow.
#
# Axes:
#   x = T
#   y = chi_shield
#   z = Kc
# Line family parameter:
#   H0,c
#
# The workflow follows a continuation path through the canonical regimes and
# writes both the raw Kc(T) data and the final plot.

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON="${PYTHON:-python3}"
PROJECT_ROOT="${PROJECT_ROOT:-$PWD}"
OUTROOT="${OUTROOT:-runs/panelA_waterfall_3d}"
SCRIPT="${SCRIPT:-$HERE/build_panelA_waterfall_3d.py}"

N_LINES="${N_LINES:-10}"
H0_MIN_EV="${H0_MIN_EV:-2.6}"
H0_MAX_EV="${H0_MAX_EV:-6.0}"
T_MIN="${T_MIN:-300}"
T_MAX="${T_MAX:-1200}"
T_STEP="${T_STEP:-50}"
KDOT="${KDOT:-0.02}"
KMAX="${KMAX:-80}"
DT="${DT:-1.0}"
VIEW_ELEV="${VIEW_ELEV:-23}"
VIEW_AZIM="${VIEW_AZIM:--63}"

"$PYTHON" - <<'PY'
import sys, numpy, pandas, matplotlib
print('python    :', sys.executable)
print('numpy     :', numpy.__version__)
print('pandas    :', pandas.__version__)
print('matplotlib:', matplotlib.__version__)
PY

if [[ ! -f "$SCRIPT" ]]; then
  echo "ERROR: workflow script not found: $SCRIPT" >&2
  exit 2
fi

"$PYTHON" "$SCRIPT" \
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

echo
echo "Panel A waterfall complete."
echo "Outputs: $OUTROOT"
