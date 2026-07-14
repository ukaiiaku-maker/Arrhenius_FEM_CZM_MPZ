#!/usr/bin/env bash
set -euo pipefail

PYTHON_BIN="${PYTHON_BIN:-python}"

"$PYTHON_BIN" -m compileall -q \
  arrhenius_fracture \
  prepare_mpz_v9_3_pt_input.py \
  prepare_mpz_v9_6_canonical_proxies.py \
  audit_mpz_v9_6_uncapped_pt.py \
  search_mpz_v9_6_broad_dbtt_map.py \
  calibrate_mpz_v9_7_pt_entropy.py \
  optimize_mpz_v9_8_joint_response.py \
  optimize_mpz_v9_8_1_joint_response.py \
  continue_mpz_v9_9_barrier_scale.py \
  promote_mpz_v9_9_spatial.py \
  search_mpz_v9_4_developed_state.py \
  audit_mpz_v9_5_state_continuation.py

PYTHONPATH=. "$PYTHON_BIN" -m pytest -q \
  tests/test_emission_derived_peierls_taylor_v96.py \
  tests/test_emission_derived_peierls_taylor_v97.py \
  tests/test_emission_derived_peierls_taylor.py \
  tests/test_bulk_pt_plasticity.py \
  tests/test_moving_process_zone.py \
  tests/test_moving_process_zone_v95.py \
  tests/test_prepare_mpz_v9_3_pt_input.py \
  tests/test_prepare_mpz_v9_6_canonical_proxies.py \
  tests/test_pt_search_v94_wrapper.py \
  tests/test_mpz_v9_4_developed_state_search.py \
  tests/test_mpz_v9_8_joint_optimizer.py \
  tests/test_mpz_v9_9_barrier_continuation.py

"$PYTHON_BIN" - <<'PY'
import numpy as np
import arrhenius_fracture as af
from arrhenius_fracture.emission_derived_plasticity import (
    CorrelatedTaylorConfig,
    EmissionDerivedPeierlsTaylorConfig,
    EmissionDerivedPeierlsTaylorModel,
)
from arrhenius_fracture.emission_derived_plasticity_v97 import (
    EmissionDerivedPeierlsTaylorModel as EntropyCalibrationModel,
    IndependentEntropyMechanismScale,
)
from arrhenius_fracture.moving_process_zone_v99 import (
    MovingProcessZoneState as PromotionMPZState,
)
from optimize_mpz_v9_8_joint_response import PARAMETER_NAMES, bounds_array
from continue_mpz_v9_9_barrier_scale import LOCAL_NAMES, LOCAL_BOUNDS

assert af.__version__ == "0.9.6"
assert af.MovingProcessZoneState.__module__.endswith(
    "moving_process_zone_v95"
)
assert PromotionMPZState.__module__.endswith("moving_process_zone_v99")
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

calibration = EntropyCalibrationModel(
    EmissionDerivedPeierlsTaylorConfig(
        peierls=IndependentEntropyMechanismScale(0.5, -20.0),
        taylor=IndependentEntropyMechanismScale(0.5, -10.0, rate_prefactor_s=1.0e11),
        correlated_taylor=CorrelatedTaylorConfig(
            rho_c_m2=1.0e14,
            renewal_time_s=1.0,
            m_cap=float("inf"),
        ),
        mobile_saturation_density_m2=float("inf"),
        mobile_density_floor_m2=0.0,
        jump_length_min_m=0.0,
        taylor_phi_max=float("inf"),
        rate_cap_s=float("inf"),
    )
)
cal_zero = calibration.rates(0.0, rho, 700.0, 2.74e-10)
assert np.all(cal_zero["equivalent_plastic_rate_s"] == 0.0)
assert bool(np.asarray(cal_zero["entropy_decoupled_from_emission"]))
assert len(PARAMETER_NAMES) == len(bounds_array()) == 17
assert len(LOCAL_NAMES) == len(LOCAL_BOUNDS) == 11

print("package version:", af.__version__)
print("active MPZ state:", af.MovingProcessZoneState.__module__)
print("promotion MPZ state:", PromotionMPZState.__module__)
print("active PT model:", EmissionDerivedPeierlsTaylorModel.__module__)
print("entropy calibration model:", EntropyCalibrationModel.__module__)
print("joint optimizer parameters:", len(PARAMETER_NAMES))
print("continuation local parameters:", len(LOCAL_NAMES))
print("zero-stress maximum rate:", float(np.max(zero["equivalent_plastic_rate_s"])))
print("constitutive caps active:", bool(np.asarray(driven["constitutive_caps_active"])))
PY

echo "MPZ v9.6 production through v9.9 continuation/promotion verification passed."
