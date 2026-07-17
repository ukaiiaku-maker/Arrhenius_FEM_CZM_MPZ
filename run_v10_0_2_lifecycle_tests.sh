#!/usr/bin/env bash
set -euo pipefail

CONDA_ENV=${CONDA_ENV:-arrhenius-fem-czm}
PYTHON_BIN=${PYTHON_BIN:-}

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
  tests/test_kinetic_campaign_reset_v1001.py \
  tests/test_kinetic_event_lifecycle_v1002.py \
  tests/test_progressive_run_2d_transform_v10.py \
  tests/test_progressive_event_lifecycle_transform_v1002.py

"$PYTHON_BIN" - <<'PY'
import arrhenius_fracture
from arrhenius_fracture import sharp_front
from arrhenius_fracture.kinetic_progressive_2d_v1002 import (
    build_progressive_run_2d_v1002,
)
from arrhenius_fracture.mode_i_first_passage_v10_0_2_progressive import main

assert arrhenius_fracture.__version__ == "10.0.2"
transformed = build_progressive_run_2d_v1002(sharp_front.run_2d)
assert transformed._v1002_event_lifecycle is True
assert callable(main)
print("v10.0.2 package version, entry point, and production-loop transform verified")
PY

cat <<'EOF'
V10.0.2 EVENT-LIFECYCLE FOUNDATION COMPLETE
- reduced-dt transactional retry tests passed
- unused-time carry and target-stop accounting tests passed
- dot_ep is included in full trial rollback history
- guarded sharp_front.run_2d lifecycle transform compiled
- no 2-D progressive result is certified by this foundation test
EOF
