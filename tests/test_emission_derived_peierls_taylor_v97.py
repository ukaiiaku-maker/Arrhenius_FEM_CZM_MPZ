import numpy as np

from arrhenius_fracture.emission_derived_plasticity import (
    CorrelatedTaylorConfig,
    EmissionDerivedPeierlsTaylorConfig,
    ExpFloorSurface,
)
from arrhenius_fracture.emission_derived_plasticity_v97 import (
    EmissionDerivedPeierlsTaylorModel,
    IndependentEntropyMechanismScale,
    KB_EV_PER_K,
)


def make_model(parent_gT, p_entropy, t_entropy, p_energy=0.5, t_energy=0.5):
    parent = ExpFloorSurface(
        G00_eV=2.0,
        gT_eV_per_K=parent_gT,
        sigc0_Pa=3.0e9,
        sT_Pa_per_K=0.0,
        Tref_K=500.0,
        a=0.2,
        n=1.0,
        floor_fraction=0.02,
    )
    cfg = EmissionDerivedPeierlsTaylorConfig(
        parent=parent,
        peierls=IndependentEntropyMechanismScale(
            p_energy, p_entropy, 1.0, 1.0e12
        ),
        taylor=IndependentEntropyMechanismScale(
            t_energy, t_entropy, 1.0, 1.0e11
        ),
        correlated_taylor=CorrelatedTaylorConfig(
            rho_c_m2=1.0e14,
            renewal_time_s=1.0,
            m_exponent=1.0,
            m_scale=1.0,
            m_cap=float("inf"),
        ),
        mobile_fraction_low_density=0.01,
        mobile_saturation_density_m2=float("inf"),
        mobile_density_floor_m2=0.0,
        jump_length_min_m=0.0,
        taylor_phi_max=float("inf"),
        rate_cap_s=float("inf"),
    )
    return EmissionDerivedPeierlsTaylorModel(cfg)


def test_independent_entropy_does_not_inherit_emission_temperature_slope():
    a = make_model(0.001, -20.0, 10.0)
    b = make_model(0.020, -20.0, 10.0)
    for T in (300.0, 500.0, 1200.0):
        assert np.isclose(
            a.raw_zero_stress_barrier_eV("peierls", T),
            b.raw_zero_stress_barrier_eV("peierls", T),
        )
        assert np.isclose(
            a.raw_zero_stress_barrier_eV("taylor", T),
            b.raw_zero_stress_barrier_eV("taylor", T),
        )


def test_activation_entropy_has_expected_free_energy_sign():
    positive = make_model(0.004, 20.0, 20.0)
    negative = make_model(0.004, -20.0, -20.0)
    assert np.isclose(
        positive.raw_zero_stress_barrier_eV("peierls", 500.0), 1.0
    )
    assert positive.raw_zero_stress_barrier_eV(
        "peierls", 1200.0
    ) < positive.raw_zero_stress_barrier_eV("peierls", 300.0)
    assert negative.raw_zero_stress_barrier_eV(
        "peierls", 1200.0
    ) > negative.raw_zero_stress_barrier_eV("peierls", 300.0)
    expected_change = -20.0 * KB_EV_PER_K * (1200.0 - 500.0)
    assert np.isclose(
        positive.raw_zero_stress_barrier_eV("peierls", 1200.0) - 1.0,
        expected_change,
    )


def test_v97_preserves_uncapped_detailed_balance():
    model = make_model(0.004, -20.0, -10.0)
    rho = np.logspace(12, 18, 20)
    out = model.rates(0.0, rho, 700.0, 2.74e-10)
    assert np.all(out["series_rate_s"] == 0.0)
    assert np.all(out["equivalent_plastic_rate_s"] == 0.0)
    assert bool(np.asarray(out["entropy_decoupled_from_emission"]))
    assert not bool(np.asarray(out["constitutive_caps_active"]))
