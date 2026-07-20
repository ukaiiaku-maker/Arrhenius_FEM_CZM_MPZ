from __future__ import annotations

from types import SimpleNamespace

import pytest

from arrhenius_fracture.barrier_only_response_registry_v100513 import (
    TWO_D_STATE_POLICY,
    load_barrier_option,
)
from arrhenius_fracture import mode_i_first_passage_v10_0_5_13_1_barrier_only as entry
import run_v10_0_5_13_1_barrier_only_monotonic as campaign


def test_pt_barrier_hook_does_not_reassign_nonbarrier_state_fields():
    row = load_barrier_option("peak_primary").legacy_row()
    cfg = SimpleNamespace(
        source_sites_per_system=731.0,
        source_recovery_rate_s=2.3,
        source_refresh_length_m=8.2e-6,
        pt_encounter_efficiency=4.7,
        pt_forest_density_floor_m2=9.1e13,
        pt_peierls_stress_fraction=0.41,
        pt_taylor_stress_fraction=0.37,
        pt_taylor_corr_rho_c=6.2e15,
        pt_taylor_renewal_time_s=4.5e-8,
        pt_taylor_m_exponent=2.4,
        pt_taylor_m_scale=7.1,
        pt_taylor_m_cap=13.0,
        mobile_shield_fraction=0.62,
        shielding_core_m=8.8e-10,
        retained_recovery_nu0_s=4.4e7,
        retained_recovery_barrier_eV=2.2,
        blunting_length_m=6.4e-6,
        blunting_slip_fraction=0.33,
    )
    before = dict(vars(cfg))
    entry._apply_barrier_pt_config(cfg, row)

    for name, value in before.items():
        assert getattr(cfg, name) == value, name

    assert cfg.pt_emit_G00_eV == pytest.approx(row["emit_G00_eV"])
    assert cfg.pt_peierls_energy_ratio == pytest.approx(
        row["peierls_H0_eV"] / row["emit_G00_eV"]
    )
    assert cfg.pt_taylor_energy_ratio == pytest.approx(
        row["taylor_H0_eV"] / row["emit_G00_eV"]
    )
    assert cfg.pt_peierls_exp_a == pytest.approx(row["peierls_exp_a"])
    assert cfg.pt_taylor_exp_n == pytest.approx(row["taylor_exp_n"])


def test_mpz_builder_preserves_existing_namespace_state_configuration():
    row = load_barrier_option("weakT_primary").legacy_row()
    args = SimpleNamespace(
        mpz_length_m=100.0e-6,
        mpz_n_bins=80,
        mpz_n_systems=2,
        mpz_source_sites_per_system=321.0,
        mpz_source_recovery_rate_s=0.017,
        mpz_source_refresh_length_m=9.0e-6,
        mpz_source_bin_count=7,
        mpz_shielding_factors="0.25 0.75",
        mpz_mobile_shield_fraction=0.42,
        mpz_shielding_core_m=7.5e-10,
        mpz_retained_recovery_nu0_s=6.0e8,
        mpz_retained_recovery_barrier_eV=1.75,
        mpz_retained_recovery_activation_volume_b3=3.0,
        mpz_mobile_recovery_rate_s=0.14,
        mpz_pair_annihilation_rate_per_count_s=0.002,
        mpz_blunting_length_m=3.2e-6,
        mpz_blunting_slip_fraction=0.58,
        pt_encounter_efficiency=2.6,
        pt_forest_density_floor_m2=8.0e13,
        pt_peierls_stress_fraction=0.45,
        pt_taylor_stress_fraction=0.35,
        pt_taylor_corr_rho_c=4.0e15,
        pt_taylor_renewal_time_s=2.0e-8,
        pt_taylor_m_exponent=1.8,
        pt_taylor_m_scale=3.3,
        pt_taylor_m_cap=12.0,
    )
    cfg = entry._build_barrier_only_mpz_config(args, row)

    assert cfg.length_m == pytest.approx(TWO_D_STATE_POLICY["mpz_length_um"] * 1.0e-6)
    assert cfg.n_bins == TWO_D_STATE_POLICY["mpz_n_bins"]
    assert cfg.source_sites_per_system == pytest.approx(321.0)
    assert cfg.source_recovery_rate_s == pytest.approx(0.017)
    assert cfg.source_refresh_length_m == pytest.approx(9.0e-6)
    assert cfg.source_bin_count == 7
    assert cfg.shielding_orientation_factors == pytest.approx((0.25, 0.75))
    assert cfg.mobile_shield_fraction == pytest.approx(0.42)
    assert cfg.retained_recovery_nu0_s == pytest.approx(6.0e8)
    assert cfg.retained_recovery_barrier_eV == pytest.approx(1.75)
    assert cfg.blunting_length_m == pytest.approx(3.2e-6)
    assert cfg.blunting_slip_fraction == pytest.approx(0.58)
    assert cfg.pt_encounter_efficiency == pytest.approx(2.6)
    assert cfg.pt_forest_density_floor_m2 == pytest.approx(8.0e13)
    assert cfg.pt_peierls_stress_fraction == pytest.approx(0.45)
    assert cfg.pt_taylor_stress_fraction == pytest.approx(0.35)
    assert cfg.pt_taylor_corr_rho_c == pytest.approx(4.0e15)
    assert cfg.pt_taylor_m_scale == pytest.approx(3.3)


def test_campaign_wrapper_routes_point_release_entry(tmp_path):
    from argparse import Namespace

    args = Namespace(
        registry=None,
        tip_refinement_radius_um=330.0,
        cluster_J_outer_um=240.0,
        local_J_outer_um=100.0,
        steps=50000,
        nx=36,
        ny=72,
        tip_h_fine=2.5e-6,
        tip_ratio=1.15,
        dU=2.0e-7,
        dt=8.4,
        n_stagger=2,
        print_every=25,
        adaptive_event_target=0.15,
        da_um=5.0,
        theta_deg=45.0,
        save_snapshots=3,
        snapshot_cols=3,
        snapshot_interval_um=50.0,
    )
    cmd = campaign._build_command(
        "/example/python", args, "ceramic_primary", 300, 100.0, tmp_path / "case"
    )
    assert campaign.ENTRY_MODULE in cmd
    assert "arrhenius_fracture.mode_i_first_passage_v10_0_5_13_barrier_only" not in cmd
