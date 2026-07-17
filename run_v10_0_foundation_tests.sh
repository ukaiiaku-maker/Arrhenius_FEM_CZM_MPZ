#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV=${CONDA_ENV:-arrhenius-fem-czm}
PYTHON_BIN=${PYTHON_BIN:-}
OUTROOT=${OUTROOT:-runs/v10_0_foundation}

if [[ -z "$PYTHON_BIN" ]]; then
  if [[ "${CONDA_DEFAULT_ENV:-}" == "$CONDA_ENV" ]]; then
    PYTHON_BIN=$(command -v python)
  else
    PYTHON_BIN=$(conda run -n "$CONDA_ENV" python -c 'import sys; print(sys.executable)' | tail -n 1)
  fi
fi

"$PYTHON_BIN" -m pip install -e . --no-deps
"$PYTHON_BIN" -m compileall -q arrhenius_fracture
"$PYTHON_BIN" -m pytest -q \
  tests/test_pf_equivalent_manifest_v10.py \
  tests/test_kinetic_campaign_czm_v10.py \
  tests/test_cohesive_trial_state_v10.py \
  tests/test_kinetic_cohesive_stepper_v10.py \
  tests/test_progressive_run_2d_transform_v10.py

mkdir -p "$OUTROOT"
for material in ceramic weakT DBTT; do
  "$PYTHON_BIN" -m arrhenius_fracture.kinetic_campaign_parity \
    --material "$material" \
    --out "$OUTROOT/isolated_$material"
done

cat <<EOF
FOUNDATION COMPLETE
- unit tests passed
- actual sharp_front.run_2d progressive source transform compiled
- isolated CZM traces written under $OUTROOT
- PF parity is NOT certified until --reference-json traces from PF v10.1.7.1 are compared
- progressive 2-D production remains blocked until a one-segment smoke, rollback audit, and penalty convergence pass
EOF
