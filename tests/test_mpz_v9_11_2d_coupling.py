from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np

from arrhenius_fracture.mpz_parameterization_v911 import load_selected_row
from arrhenius_fracture.process_zone_2d_v911 import sample_process_zone_profile


def test_exact_dbtt_primary_row_is_selected():
    root = Path(__file__).resolve().parents[1] / "mpz_v9_11_parameters"
    row = load_selected_row(root / "DBTT" / "spatial_promotion_manifest.csv", "DBTT")
    assert row["candidate_id"] == "DBTT_restart01_candidate05"
    assert np.isclose(row["peierls_H0_eV"], 3.35812, rtol=2e-6)
    assert np.isclose(row["taylor_H0_eV"], 10.4361, rtol=2e-6)
    assert row["peierls_exp_a"] != row["emit_exp_a"]
    assert row["taylor_exp_n"] != row["emit_exp_n"]


def test_2d_profile_uses_density_but_only_dimensionless_stress_shape():
    nodes = []
    elems = []
    for i in range(4):
        x0, x1 = i * 1e-6, (i + 1) * 1e-6
        base = len(nodes)
        nodes += [(x0, -0.2e-6), (x1, -0.2e-6), (x1, 0.2e-6), (x0, 0.2e-6)]
        elems += [(base, base + 1, base + 2), (base, base + 2, base + 3)]
    mesh = SimpleNamespace(
        nodes=np.asarray(nodes, float),
        elems=np.asarray(elems, int),
        area_e=np.full(len(elems), 0.2e-12),
    )
    ne = len(elems)
    sigma = np.zeros((3, ne))
    sigma[0] = np.linspace(4e9, 1e9, ne)
    sigma[1] = 0.2 * sigma[0]
    rho = np.linspace(5e12, 2e14, ne)
    damage = np.zeros(len(nodes))
    profile = sample_process_zone_profile(
        mesh, sigma, rho, damage, [0.0, 0.0], [1.0, 0.0],
        length_m=4e-6, n_bins=4, min_elements=4,
    )
    assert profile.reliable
    assert profile.forest_density_m2[-1] > profile.forest_density_m2[0]
    assert profile.stress_shape[0] > profile.stress_shape[-1]
    assert np.max(profile.stress_shape) <= 1.0
    assert profile.diagnostics()["bulk_scalar_rho_used_for_signed_shielding"] is False


def test_bulk_pt_config_preserves_independent_shapes_and_uncapped_limits():
    from arrhenius_fracture.bulk_plasticity_v9102 import independent_config_from_dislocation_config
    from arrhenius_fracture.mpz_parameterization_v911 import apply_pt_dislocation_config

    root = Path(__file__).resolve().parents[1] / "mpz_v9_11_parameters"
    row = load_selected_row(root / "weakT" / "spatial_promotion_manifest.csv", "weakT")
    cfg_in = SimpleNamespace()
    apply_pt_dislocation_config(cfg_in, row)
    cfg = independent_config_from_dislocation_config(cfg_in)
    assert np.isclose(cfg.peierls.exp_a, row["peierls_exp_a"])
    assert np.isclose(cfg.peierls.exp_n, row["peierls_exp_n"])
    assert np.isclose(cfg.taylor.exp_a, row["taylor_exp_a"])
    assert np.isclose(cfg.taylor.exp_n, row["taylor_exp_n"])
    assert np.isinf(cfg.correlated_taylor.m_cap)
    assert np.isinf(cfg.mobile_saturation_density_m2)
    assert cfg.mobile_density_floor_m2 == 0.0
    assert cfg.jump_length_min_m == 0.0
    assert np.isinf(cfg.rate_cap_s)


def test_v911_state_split_and_advance_conserve_front_local_inventory():
    from arrhenius_fracture.mpz_parameterization_v911 import build_mpz_config
    from arrhenius_fracture.moving_process_zone_v911 import MovingProcessZoneState

    root = Path(__file__).resolve().parents[1] / "mpz_v9_11_parameters"
    row = load_selected_row(root / "DBTT" / "spatial_promotion_manifest.csv", "DBTT")
    args = SimpleNamespace(mpz_length_um=100.0, mpz_n_bins=40, r_pz=1.0e-6)
    state = MovingProcessZoneState(build_mpz_config(args, row))
    state.mobile[:, :3] = 2.0
    state.retained[:, :3] = 1.0
    state.accumulated_slip[:, :3] = 0.5
    total_before = state.mobile_count + state.retained_count
    child = state.split(0.25)
    assert np.isclose(
        state.mobile_count + state.retained_count + child.mobile_count + child.retained_count,
        total_before,
    )
    before_parent = state.mobile_count + state.retained_count
    wake = state.advance(2.0 * state.dx)
    after_parent = state.mobile_count + state.retained_count
    assert after_parent <= before_parent + 1.0e-12
    assert wake["wake_mobile"] + wake["wake_retained"] >= 0.0


def test_v911_bulk_density_changes_forest_resistance_not_direct_K_shield():
    from arrhenius_fracture.mpz_parameterization_v911 import build_mpz_config
    from arrhenius_fracture.moving_process_zone_v911 import MovingProcessZoneState
    from arrhenius_fracture.process_zone_2d_v911 import ProcessZone2DProfile

    root = Path(__file__).resolve().parents[1] / "mpz_v9_11_parameters"
    row = load_selected_row(root / "weakT" / "spatial_promotion_manifest.csv", "weakT")
    args = SimpleNamespace(mpz_length_um=100.0, mpz_n_bins=40, r_pz=1.0e-6)
    state = MovingProcessZoneState(build_mpz_config(args, row))
    k0 = state.shielding_K(160.0e9, 0.28, 2.74e-10)
    profile = ProcessZone2DProfile(
        x_m=state.x.copy(),
        forest_density_m2=np.full(state.n_bins, 2.0e14),
        stress_shape=np.linspace(1.0, 0.2, state.n_bins),
        reliable=True,
        coverage_fraction=1.0,
        selected_elements=state.n_bins,
        reason="test",
    )
    state.set_2d_profile(profile)
    assert np.max(state.local_forest_density_m2()) >= 2.0e14
    assert state.shielding_K(160.0e9, 0.28, 2.74e-10) == k0
