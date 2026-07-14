import numpy as np

from arrhenius_fracture.config import (
    DislocationConfig,
    ElasticProperties,
    PlasticityBarrier,
)
from arrhenius_fracture.materials import PlasticityModel
from arrhenius_fracture.plasticity import update_plasticity


def test_bulk_fem_uses_emission_derived_pt_and_preserves_high_density_state():
    mat = ElasticProperties()
    disl = DislocationConfig()
    disl.bulk_kinetics_model = (
        "emission_derived_peierls_taylor_multihit"
    )
    disl.thermo_consistency_mode = "time_cone"
    disl.freeze_rho = True
    disl.rho_cap = 1.0e20
    disl.pt_taylor_corr_rho_c = 1.0e11
    disl.pt_taylor_renewal_time_s = 1.0e-10
    disl.pt_taylor_m_scale = 1.0
    disl.pt_taylor_m_cap = 22.0

    model = PlasticityModel(PlasticityBarrier(), mat)
    ep = np.zeros((3, 2))
    rho = np.array([5.0e12, 1.0e18])
    stress = np.array([
        [2.0e9, 2.0e9],
        [0.0, 0.0],
        [0.0, 0.0],
    ])

    _, rho_out, rate, info = update_plasticity(
        ep,
        rho,
        stress,
        mat,
        700.0,
        1.0e-6,
        model,
        disl,
        return_info=True,
    )

    assert info["bulk_pt_active"] is True
    assert info["bulk_kinetics_model"] == (
        "emission_derived_peierls_taylor_multihit"
    )
    assert np.all(np.isfinite(info["pt_peierls_rate_gp"]))
    assert np.all(np.isfinite(info["pt_taylor_completion_rate_gp"]))
    assert np.all(
        info["pt_series_rate_gp"]
        <= info["pt_peierls_rate_gp"] * (1.0 + 1.0e-12)
    )
    assert np.all(
        info["pt_series_rate_gp"]
        <= info["pt_taylor_completion_rate_gp"] * (1.0 + 1.0e-12)
    )
    assert np.all(np.isfinite(rate))
    assert rho_out[1] == 1.0e18
