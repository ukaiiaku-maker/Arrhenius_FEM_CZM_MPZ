#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV=${CONDA_ENV:-arrhenius-fem-czm}
PYTHON_BIN=${PYTHON_BIN:-}
OUTROOT=${OUTROOT:-runs/v10_0_1_foundation}

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
  tests/test_kinetic_campaign_reset_v1001.py \
  tests/test_cohesive_trial_state_v10.py \
  tests/test_kinetic_cohesive_stepper_v10.py \
  tests/test_progressive_run_2d_transform_v10.py

"$PYTHON_BIN" - <<'PY'
from arrhenius_fracture import sharp_front
from arrhenius_fracture.kinetic_progressive_2d_v10 import build_progressive_run_2d
from arrhenius_fracture.mode_i_first_passage_v10_0_1 import main as abrupt_main
from arrhenius_fracture.mode_i_first_passage_v10_0_1_progressive import main as progressive_main

transformed = build_progressive_run_2d(sharp_front.run_2d)
assert transformed._v10_progressive_source_transform is True
assert callable(abrupt_main)
assert callable(progressive_main)
print("v10.0.1 entry points and guarded production-loop transform import successfully")
PY

mkdir -p "$OUTROOT"
for material in ceramic weakT DBTT; do
  "$PYTHON_BIN" -m arrhenius_fracture.kinetic_campaign_parity \
    --material "$material" \
    --out "$OUTROOT/isolated_$material"
done

cat <<EOF
V10.0.1 FOUNDATION COMPLETE
- focused unit tests passed
- reset-safe kinetic campaign state passed
- guarded sharp_front.run_2d transform compiled
- isolated CZM traces written under $OUTROOT
- PF parity is NOT certified until matching PF v10.1.7.1 reference traces are compared
- longer progressive runs remain blocked until rejected-step retry and unused-time carry are integrated
EOF
