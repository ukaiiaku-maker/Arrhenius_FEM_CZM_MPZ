import numpy as np

from arrhenius_fracture.emission_derived_plasticity import (
    CorrelatedTaylorConfig,
    EmissionDerivedPeierlsTaylorConfig,
)
from arrhenius_fracture.emission_derived_plasticity_v96 import (
    EmissionDerivedPeierlsTaylorModel,
)


def _model():
    return EmissionDerivedPeierlsTaylorModel(
        EmissionDerivedPeierlsTaylorConfig(
            correlated_taylor=CorrelatedTaylorConfig(
                rho_c_m2=1.0e14,
                renewal_time_s=1.0e-9,
                m_exponent=1.0,
                m_scale=1.0,
                m_cap=2.0,
            ),
            taylor_phi_max=2.0,
            mobile_fraction_low_density=0.01,
            mobile_saturation_density_m2=1.0e12,
            mobile_density_floor_m2=1.0e9,
            jump_length_min_m=1.0e-6,
            rate_cap_s=1.0e-30,
        )
    )


def test_v96_ignores_exploratory_caps_and_saturations():
    model = _model()
    rho = np.logspace(10, 20, 101)
    b = 2.74e-10

    phi = model.taylor_amplification(rho, b)
    order = model.natural_hit_order(rho)
    mobile = model.mobile_density(rho)
    jump = model.jump_length(rho)

    assert phi[0] > 2.0
    assert phi[-1] < 2.0
    assert np.all(np.diff(phi) < 0.0)
    assert order[-1] > 2.0
    assert np.all(np.diff(order) > 0.0)
    assert np.all(np.diff(mobile) > 0.0)
    assert np.all(np.diff(jump) < 0.0)
    assert mobile[-1] > 1.0e12
    assert jump[-1] < 1.0e-6


def test_v96_hit_order_matches_natural_correlation_length_formula():
    model = _model()
    rho = np.array([1.0e12, 1.0e14, 1.0e16])
    L = model.correlation_length_m()
    expected = 1.0 + 2.0 * L * np.sqrt(rho)
    assert np.allclose(model.natural_hit_order(rho), expected)


def test_v96_mean_gamma_completion_and_detailed_balance():
    model = _model()
    rho = np.logspace(12, 17, 20)
    zero = model.rates(0.0, rho, 700.0, 2.74e-10)
    assert np.all(zero["series_rate_s"] == 0.0)
    assert np.all(zero["equivalent_plastic_rate_s"] == 0.0)

    driven = model.rates(2.0e9, rho, 700.0, 2.74e-10)
    order = driven["taylor_m_eff"]
    assert np.allclose(
        driven["taylor_completion_forward_rate_s"],
        driven["taylor_single_hit_forward_rate_s"] / order,
    )
    assert np.all(driven["series_rate_s"] >= 0.0)
    assert not bool(np.asarray(driven["constitutive_caps_active"]))


def test_v96_accepts_explicit_mobile_density_state():
    model = _model()
    rho_f = np.array([1.0e13, 1.0e14])
    rho_m = np.array([2.0e10, 8.0e10])
    out = model.rates(2.0e9, rho_f, 700.0, 2.74e-10, rho_mobile_m2=rho_m)
    assert np.allclose(out["rho_mobile_m2"], rho_m)
