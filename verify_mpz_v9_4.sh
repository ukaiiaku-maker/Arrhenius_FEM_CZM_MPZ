#!/usr/bin/env bash
set -euo pipefail

PYTHON_BIN="${PYTHON_BIN:-python}"

"$PYTHON_BIN" -m compileall -q \
  arrhenius_fracture \
  prepare_mpz_v9_3_pt_input.py \
  search_mpz_peierls_taylor_parameters.py \
  search_mpz_peierls_taylor_parameters_v94.py

PYTHONPATH=. "$PYTHON_BIN" -m pytest -q \
  tests/test_emission_derived_peierls_taylor.py \
  tests/test_bulk_pt_plasticity.py \
  tests/test_moving_process_zone.py \
  tests/test_prepare_mpz_v9_3_pt_input.py

"$PYTHON_BIN" - <<'PY'
import numpy as np
import arrhenius_fracture as af
from arrhenius_fracture.emission_derived_plasticity import (
    CorrelatedTaylorConfig,
    EmissionDerivedPeierlsTaylorConfig,
    EmissionDerivedPeierlsTaylorModel,
)

assert af.__version__ == "0.9.4"
assert EmissionDerivedPeierlsTaylorModel.__module__.endswith(
    "emission_derived_plasticity_v94"
)
model = EmissionDerivedPeierlsTaylorModel(
    EmissionDerivedPeierlsTaylorConfig(
        correlated_taylor=CorrelatedTaylorConfig(
            rho_c_m2=1.0e11,
            renewal_time_s=1.0e-10,
            m_cap=22.0,
        )
    )
)
rho = np.logspace(12.0, 17.0, 20)
zero = model.rates(0.0, rho, 700.0, 2.74e-10)
assert np.all(zero["peierls_rate_s"] == 0.0)
assert np.all(zero["taylor_completion_rate_s"] == 0.0)
assert np.all(zero["series_rate_s"] == 0.0)
assert np.all(zero["equivalent_plastic_rate_s"] == 0.0)
print("package version:", af.__version__)
print("active PT model:", EmissionDerivedPeierlsTaylorModel.__module__)
print("zero-stress maximum rate:", float(np.max(zero["equivalent_plastic_rate_s"])))
PY

echo "MPZ v9.4 focused verification passed."
