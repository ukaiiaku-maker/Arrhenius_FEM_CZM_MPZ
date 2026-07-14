#!/usr/bin/env bash
set -euo pipefail

PYTHON_BIN="${PYTHON_BIN:-python}"

"$PYTHON_BIN" -m compileall -q \
  arrhenius_fracture \
  prepare_mpz_v9_3_pt_input.py \
  prepare_mpz_v9_6_canonical_proxies.py \
  audit_mpz_v9_6_uncapped_pt.py \
  search_mpz_v9_6_broad_dbtt_map.py \
  search_mpz_v9_4_developed_state.py \
  audit_mpz_v9_5_state_continuation.py

PYTHONPATH=. "$PYTHON_BIN" -m pytest -q \
  tests/test_emission_derived_peierls_taylor_v96.py \
  tests/test_emission_derived_peierls_taylor.py \
  tests/test_bulk_pt_plasticity.py \
  tests/test_moving_process_zone.py \
  tests/test_moving_process_zone_v95.py \
  tests/test_prepare_mpz_v9_3_pt_input.py \
  tests/test_prepare_mpz_v9_6_canonical_proxies.py \
  tests/test_pt_search_v94_wrapper.py \
  tests/test_mpz_v9_4_developed_state_search.py

"$PYTHON_BIN" - <<'PY'
import numpy as np
import arrhenius_fracture as af
from arrhenius_fracture.emission_derived_plasticity import (
    CorrelatedTaylorConfig,
    EmissionDerivedPeierlsTaylorConfig,
    EmissionDerivedPeierlsTaylorModel,
)

assert af.__version__ == "0.9.6"
assert af.MovingProcessZoneState.__module__.endswith(
    "moving_process_zone_v95"
)
assert EmissionDerivedPeierlsTaylorModel.__module__.endswith(
    "emission_derived_plasticity_v96"
)
model = EmissionDerivedPeierlsTaylorModel(
    EmissionDerivedPeierlsTaylorConfig(
        correlated_taylor=CorrelatedTaylorConfig(
            rho_c_m2=1.0e14,
            renewal_time_s=1.0e-10,
            m_cap=2.0,
        ),
        taylor_phi_max=2.0,
        mobile_saturation_density_m2=1.0e12,
        jump_length_min_m=1.0e-6,
        rate_cap_s=1.0e-30,
    )
)
rho = np.logspace(12.0, 17.0, 20)
zero = model.rates(0.0, rho, 700.0, 2.74e-10)
driven = model.rates(2.0e9, rho, 700.0, 2.74e-10)
assert np.all(zero["peierls_rate_s"] == 0.0)
assert np.all(zero["taylor_completion_rate_s"] == 0.0)
assert np.all(zero["series_rate_s"] == 0.0)
assert np.all(zero["equivalent_plastic_rate_s"] == 0.0)
assert np.all(np.diff(driven["taylor_m_eff"]) > 0.0)
assert np.all(np.diff(driven["taylor_amplification"]) < 0.0)
assert not bool(np.asarray(driven["constitutive_caps_active"]))
print("package version:", af.__version__)
print("active MPZ state:", af.MovingProcessZoneState.__module__)
print("active PT model:", EmissionDerivedPeierlsTaylorModel.__module__)
print("zero-stress maximum rate:", float(np.max(zero["equivalent_plastic_rate_s"])))
print("constitutive caps active:", bool(np.asarray(driven["constitutive_caps_active"])))
PY

echo "MPZ v9.6 focused verification passed."
