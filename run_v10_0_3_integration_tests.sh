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
"$PYTHON_BIN" -m compileall -q arrhenius_fracture tests
"$PYTHON_BIN" -m py_compile \
  audit_v10_0_3_progressive_integration.py \
  normalize_v10_0_3_1_reporting.py

"$PYTHON_BIN" - <<'PY'
import importlib.metadata

from arrhenius_fracture import sharp_front
from arrhenius_fracture.kinetic_campaign_czm_v1003 import (
    CampaignAwareV1003TipEngineMixin,
    STATE_MODEL,
)
from arrhenius_fracture.kinetic_progressive_2d_v1003_source import (
    build_progressive_run_2d_v1003_source,
)

version = importlib.metadata.version("arrhenius-fem-czm")
assert version == "10.0.3.1", version
assert STATE_MODEL == "kinetic_campaign_czm"
assert CampaignAwareV1003TipEngineMixin.supports_progressive_kinetic_czm is True
transformed = build_progressive_run_2d_v1003_source(sharp_front.run_2d)
assert transformed._v1002_event_lifecycle is True
assert transformed._v1003_source_adapter is True
assert transformed._v1003_campaign_state_compatibility is True
assert transformed._v1003_nondeflect_summary_accounting is True
print("package version:", version)
print("v10.0.3 integration source transform preflight: PASS")
PY

"$PYTHON_BIN" -m pytest -q \
  tests/test_pf_equivalent_manifest_v10.py \
  tests/test_kinetic_campaign_czm_v10.py \
  tests/test_cohesive_trial_state_v10.py \
  tests/test_kinetic_cohesive_stepper_v10.py \
  tests/test_kinetic_campaign_reset_v1001.py \
  tests/test_kinetic_event_lifecycle_v1002.py \
  tests/test_progressive_run_2d_transform_v10.py \
  tests/test_progressive_event_lifecycle_transform_v1002.py \
  tests/test_v1003_campaign_dispatch.py \
  tests/test_v1003_live_binding_capture.py \
  tests/test_v1003_source_population_bound.py \
  tests/test_v10031_reporting_normalization.py

cat <<'EOF'
V10.0.3.1 TESTS-ONLY INTEGRATION/REPORTING GATE PASSED
No FEM solve was launched.
The certified integration kernel remains v10.0.3; v10.0.3.1 normalizes reporting only.
EOF
