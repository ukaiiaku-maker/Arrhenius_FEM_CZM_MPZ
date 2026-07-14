import numpy as np

from arrhenius_fracture.emission_derived_plasticity import (
    CorrelatedTaylorConfig,
    EmissionDerivedPeierlsTaylorConfig,
    EmissionDerivedPeierlsTaylorModel,
)


def test_scaled_barriers_preserve_selected_ratios_at_tref_zero_stress():
    model = EmissionDerivedPeierlsTaylorModel(
        EmissionDerivedPeierlsTaylorConfig()
    )
    T = model.cfg.parent.Tref_K
    ge = float(model.barrier_eV("emission", 0, T))
    gp = float(model.barrier_eV("peierls", 0, T))
    gt = float(model.barrier_eV("taylor", 0, T))
    assert np.isclose(gp / ge, 0.005, rtol=1e-8)
    assert np.isclose(gt / ge, 0.02, rtol=1e-8)


def test_multihit_order_is_uncapped_and_density_monotone():
    cfg = CorrelatedTaylorConfig(
        rho_c_m2=1e12,
        renewal_time_s=1e-9,
        m_exponent=1,
        m_cap=float("inf"),
    )
    rho = np.logspace(10, 18, 81)
    order = cfg.hit_order(rho)
    assert np.all(np.diff(order) >= 0)
    assert order[-1] > 100


def test_correlated_completion_never_exceeds_single_hit_rate():
    model = EmissionDerivedPeierlsTaylorModel(
        EmissionDerivedPeierlsTaylorConfig(
            correlated_taylor=CorrelatedTaylorConfig(rho_c_m2=1e12)
        )
    )
    rho = np.logspace(12, 17, 20)
    rates = model.rates(2e9, rho, 700, 2.74e-10)
    assert np.all(
        rates["taylor_completion_rate_s"]
        <= rates["taylor_single_hit_rate_s"] * (1 + 1e-12)
    )
    assert np.all(
        rates["series_rate_s"]
        <= rates["peierls_rate_s"] * (1 + 1e-12)
    )
    assert np.all(
        rates["series_rate_s"]
        <= rates["taylor_completion_rate_s"] * (1 + 1e-12)
    )


def test_finite_correlation_domain_is_not_a_density_cap_and_can_remain_monotone():
    from arrhenius_fracture.emission_derived_plasticity import (
        ExpFloorSurface,
        MechanismScale,
    )

    parent = ExpFloorSurface(
        G00_eV=1.1173,
        gT_eV_per_K=0.0063971,
        sigc0_Pa=0.9506e9,
        sT_Pa_per_K=0.0009817e9,
        a=0.5055,
        n=0.8432,
        floor_fraction=0.0309,
    )
    cfg = EmissionDerivedPeierlsTaylorConfig(
        parent=parent,
        peierls=MechanismScale(0.005, 0.05, 1.0, 1.0e12),
        taylor=MechanismScale(0.02, 0.20, 1.0, 1.0e11),
        correlated_taylor=CorrelatedTaylorConfig(
            rho_c_m2=1.0e11,
            renewal_time_s=1.0e-10,
            m_exponent=1.0,
            m_scale=1.0,
            m_cap=15.0,
        ),
        mobile_fraction_low_density=0.01,
        mobile_saturation_density_m2=1.0e14,
    )
    model = EmissionDerivedPeierlsTaylorModel(cfg)
    rho = np.logspace(np.log10(5.0e12), 18.0, 65)
    sig = model.flow_stress(rho, 700.0, 1.0e-5, 2.74e-10)
    assert np.all(np.isfinite(sig))
    assert np.min(np.diff(sig)) >= -2.0e7
    assert rho[-1] == 1.0e18
