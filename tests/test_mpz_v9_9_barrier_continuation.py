import numpy as np
import pandas as pd

from arrhenius_fracture.moving_process_zone import MovingProcessZoneConfig
from arrhenius_fracture.moving_process_zone_v99 import MovingProcessZoneState
from continue_mpz_v9_9_barrier_scale import (
    build_model,
    fixed_parameters,
    initial_local_vector,
    local_to_parameters,
    shape_from_row,
)
from promote_mpz_v9_9_spatial import material_row


def candidate_row():
    return pd.Series(
        {
            "continuation_candidate_id": "weakT_rank00_scale0.6",
            "target_class": "weakT",
            "peierls_H0_eV": 4.0,
            "taylor_H0_eV": 10.0,
            "peierls_activation_entropy_kB": -20.0,
            "taylor_activation_entropy_kB": -10.0,
            "peierls_nu0_s": 1.0e12,
            "taylor_nu0_s": 1.0e11,
            "cleave_G00_eV": 2.0,
            "cleave_gT_eV_per_K": 0.0,
            "cleave_sigc0_GPa": 4.0,
            "emit_G00_eV": 1.5,
            "emit_gT_eV_per_K": 0.0,
            "emit_sigc0_GPa": 2.5,
            "shape_cleave_sT_GPa_per_K": 0.0,
            "shape_cleave_exp_a": 0.2,
            "shape_cleave_exp_n": 1.0,
            "shape_cleave_floor_frac": 0.02,
            "shape_emit_sT_GPa_per_K": 0.0,
            "shape_emit_exp_a": 0.2,
            "shape_emit_exp_n": 1.0,
            "shape_emit_floor_frac": 0.02,
            "log10_taylor_corr_rho_c_m2": 14.0,
            "log10_taylor_corr_scale": 0.0,
            "log10_mobile_fraction": -2.0,
            "log10_source_sites_per_system": 2.0,
            "log10_recovery_rate_s": -5.0,
            "log10_source_refresh_length_um": 1.0,
            "c_blunt": 1.0,
            "source_sites_per_system": 100.0,
            "source_refresh_length_um": 10.0,
            "recovery_rate_s": 1.0e-5,
            "taylor_corr_rho_c_m2": 1.0e14,
            "taylor_corr_scale": 1.0,
            "mobile_fraction": 0.01,
            "barrier_scale": 0.6,
        }
    )


def test_common_scaling_preserves_pt_barrier_ratio_and_order():
    row = candidate_row()
    fixed = fixed_parameters(row, 0.4)
    assert np.isclose(fixed["peierls_H0_eV"], 1.6)
    assert np.isclose(fixed["taylor_H0_eV"], 4.0)
    assert np.isclose(
        fixed["taylor_H0_eV"] / fixed["peierls_H0_eV"], 10.0 / 4.0
    )
    assert fixed["taylor_H0_eV"] >= fixed["peierls_H0_eV"]


def test_continuation_warm_start_and_model_are_finite():
    row = candidate_row()
    x = initial_local_vector(row, 0.6, 700.0)
    assert np.all(np.isfinite(x))
    local = local_to_parameters(x)
    fixed = fixed_parameters(row, 0.6)
    model = build_model(fixed, local, shape_from_row(row), 481.33)
    rates = model.rates(2.0e9, 1.0e14, 700.0, 2.74e-10)
    assert np.isfinite(float(np.asarray(rates["peierls_rate_s"])))
    assert np.isfinite(float(np.asarray(rates["taylor_completion_rate_s"])))
    assert not bool(np.asarray(rates["constitutive_caps_active"]))


def test_v99_spatial_adapter_uses_independent_entropy_and_no_caps():
    cfg = MovingProcessZoneConfig(
        n_bins=20,
        n_systems=2,
        pt_emit_G00_eV=2.0,
        pt_emit_gT_eV_per_K=0.004,
        pt_peierls_energy_ratio=0.5,
        pt_peierls_entropy_ratio=-20.0,
        pt_peierls_nu0_s=1.0e12,
        pt_taylor_energy_ratio=1.0,
        pt_taylor_entropy_ratio=-10.0,
        pt_taylor_nu0_s=1.0e11,
        pt_taylor_m_cap=float("inf"),
        pt_mobile_saturation_density_m2=float("inf"),
        pt_mobile_density_floor_m2=0.0,
        pt_jump_length_min_m=0.0,
        pt_taylor_phi_max=float("inf"),
    )
    state = MovingProcessZoneState(cfg)
    out = state.evolve(0.0, 700.0, 2.0e9, 2.74e-10, 0.0)
    assert out["pt_independent_entropy_active"] == 1.0
    assert out["pt_peierls_activation_entropy_kB"] == -20.0
    assert out["pt_taylor_activation_entropy_kB"] == -10.0
    assert np.isfinite(out["peierls_rate_s"])
    assert np.isfinite(out["taylor_completion_rate_s"])


def test_spatial_material_row_is_complete_and_ordered():
    row = material_row(candidate_row(), 100.0e-6, 200)
    assert row.mpz_n_bins == 200
    assert row.mpz_length_m == 100.0e-6
    assert row.pt_taylor_energy_ratio >= row.pt_peierls_energy_ratio
    assert row.pt_peierls_entropy_ratio == -20.0
    assert row.pt_taylor_entropy_ratio == -10.0
    assert row.mpz_source_sites_per_system == 100.0
