#!/usr/bin/env bash
set -euo pipefail

# Targeted extension of the existing 2-D full-plastic S-N sweep.
# Existing stresses 500, 550, 600, 700 MPa are not repeated.
# New grid:
#   450       : low-stress censored reference
#   575-675   : transition refinement for unshielded response
#   750-1000  : shielded crack-onset and clear-cracking branch

NEW_STRESSES_DEFAULT="450 575 625 650 675 750 800 850 900 1000"
NEW_STRESSES="${NEW_STRESSES:-$NEW_STRESSES_DEFAULT}"
CYCLES_MAX="${CYCLES_MAX:-1e10}"
NX="${NX:-36}"
NY="${NY:-72}"
RESULT_ROOT="${RESULT_ROOT:-runs/sn_pf2d_arrhenius_fullplastic_production}"
ANALYSIS_OUT="${ANALYSIS_OUT:-${RESULT_ROOT}/staged_initiation_analysis_v2}"
PRODUCTION_SCRIPT="${PRODUCTION_SCRIPT:-run_sn_fullplastic_arrhenius_production.sh}"
ANALYSIS_SCRIPT="${ANALYSIS_SCRIPT:-analyze_sn_pf2d_staged_events_v2.py}"
LOG_DIR="${LOG_DIR:-${RESULT_ROOT}/extension_logs}"

mkdir -p "$LOG_DIR"

if [[ ! -f "$PRODUCTION_SCRIPT" ]]; then
  echo "ERROR: production script not found: $PRODUCTION_SCRIPT" >&2
  exit 2
fi
if [[ ! -f "$ANALYSIS_SCRIPT" ]]; then
  echo "ERROR: analysis script not found: $ANALYSIS_SCRIPT" >&2
  exit 2
fi

# Fail before launching expensive work if the active Python stack is not healthy.
python - <<'PY'
import sys
import numpy
import scipy
import scipy.sparse.linalg
print("python:", sys.executable)
print("numpy:", numpy.__version__)
print("scipy:", scipy.__version__)
print("SciPy sparse import OK")
PY

python -m compileall -q arrhenius_fracture

case_complete() {
  local case_name="$1"
  local stress="$2"
  local d="${RESULT_ROOT}/${case_name}/sigmaA_${stress}MPa"
  [[ -s "${d}/sn_pf2d_fullplastic_history.csv" && -s "${d}/summary.json" ]]
}

run_one_stress() {
  local stress="$1"
  local log="${LOG_DIR}/sigmaA_${stress}MPa.log"

  if case_complete no_shield "$stress" && case_complete shielded "$stress"; then
    echo "=== SKIP sigma_a=${stress} MPa: both cases already complete ==="
    return 0
  fi

  echo "========================================================================"
  echo "  EXTENSION RUN sigma_a=${stress} MPa"
  echo "  CYCLES_MAX=${CYCLES_MAX} NX=${NX} NY=${NY}"
  echo "========================================================================"

  # The existing production script remains the source of truth for all physics.
  # We only give it one stress at a time so this wrapper is restartable.
  STRESSES="$stress" \
  CYCLES_MAX="$CYCLES_MAX" \
  NX="$NX" NY="$NY" \
  bash "$PRODUCTION_SCRIPT" 2>&1 | tee "$log"
}

for stress in $NEW_STRESSES; do
  run_one_stress "$stress"
done

echo "========================================================================"
echo "  RUNNING COMBINED STAGED INITIATION ANALYSIS"
echo "  root: $RESULT_ROOT"
echo "========================================================================"
python "$ANALYSIS_SCRIPT" "$RESULT_ROOT" --out "$ANALYSIS_OUT"

echo
echo "Completed targeted extension and analysis."
echo "Results:  $RESULT_ROOT"
echo "Analysis: $ANALYSIS_OUT"
