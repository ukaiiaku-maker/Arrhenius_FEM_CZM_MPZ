import numpy as np

import arrhenius_fracture as af
from arrhenius_fracture.moving_process_zone import MovingProcessZoneConfig


def test_v95_state_is_active_and_local_density_scales_with_bin_area():
    assert af.__version__ == "0.9.5"
    state = af.MovingProcessZoneState(
        MovingProcessZoneConfig(
            length_m=5.0e-5,
            n_bins=100,
            n_systems=2,
            blunting_length_m=1.0e-6,
            pt_forest_density_floor_m2=5.0e12,
        )
    )
    assert state.__class__.__module__.endswith("moving_process_zone_v95")
    base = state.local_forest_density_m2()
    state.retained[0, 0] = 1.0
    rho = state.local_forest_density_m2()
    expected_increment = 1.0 / (state.dx * 1.0e-6)
    assert np.isclose(rho[0] - base[0], expected_increment)
    assert rho[0] > rho[-1]


def test_v95_local_stress_decays_away_from_tip():
    state = af.MovingProcessZoneState(
        MovingProcessZoneConfig(length_m=5.0e-5, n_bins=100)
    )
    profile = state.local_stress_profile_Pa(4.0e9)
    assert profile.shape == (state.n_bins,)
    assert np.all(np.diff(profile) < 0.0)
    assert profile[0] < 4.0e9
    assert profile[-1] < profile[0]


def test_v95_seed_profile_has_requested_tip_density_scale():
    cfg = MovingProcessZoneConfig(
        length_m=5.0e-5,
        n_bins=100,
        n_systems=2,
        blunting_length_m=1.0e-6,
        pt_forest_density_floor_m2=5.0e12,
    )
    state = af.MovingProcessZoneState(cfg)
    state.initialize_forest_profile(
        1.0e14,
        decay_length_m=5.0e-6,
        available_site_fraction=0.5,
    )
    rho = state.local_forest_density_m2()
    assert 8.0e13 < rho[0] < 1.0e14
    assert rho[-1] < rho[0]
    assert np.isclose(state.available_site_fraction, 0.5)
    assert state.retained_count > 0.0


def test_v95_evolve_reports_local_density_diagnostics():
    cfg = MovingProcessZoneConfig(
        length_m=2.0e-5,
        n_bins=40,
        n_systems=2,
        source_sites_per_system=10.0,
        blunting_length_m=1.0e-6,
        pt_forest_density_floor_m2=5.0e12,
    )
    state = af.MovingProcessZoneState(cfg)
    state.initialize_forest_profile(1.0e14, decay_length_m=3.0e-6)
    out = state.evolve(
        1.0e-12,
        900.0,
        2.0e9,
        2.74e-10,
        emission_hazard_integral=0.0,
    )
    assert out["rho_forest_max_m2"] > out["rho_forest_min_m2"]
    assert out["local_stress_max_Pa"] > out["local_stress_min_Pa"]
    assert np.isfinite(out["peierls_rate_s"])
    assert np.isfinite(out["taylor_completion_rate_s"])
