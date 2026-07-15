from __future__ import annotations

import numpy as np

from arrhenius_fracture.moving_process_zone_v9102 import MovingProcessZoneState
from arrhenius_fracture.moving_process_zone import MovingProcessZoneConfig


def _state():
    cfg = MovingProcessZoneConfig(
        length_m=10e-6,
        n_bins=20,
        n_systems=1,
        source_sites_per_system=1.0,
        source_bin_count=2,
        source_recovery_rate_s=0.0,
        pt_forest_density_floor_m2=5e12,
        pt_mobile_fraction=0.0,
        retained_recovery_nu0_s=0.0,
        mobile_recovery_rate_s=0.0,
        pair_annihilation_rate_per_count_s=0.0,
    )
    return MovingProcessZoneState(cfg)


def test_one_emitted_line_adds_one_blunting_ledger_count():
    state = _state()
    emitted = state._source_commit_from_hazard(100.0)
    assert np.isclose(np.sum(emitted), 1.0)
    assert np.isclose(state.emitted_total, 1.0)
    assert np.isclose(np.sum(state.accumulated_slip), 1.0)


def test_peierls_taylor_transport_does_not_multiply_blunting_ledger():
    state = _state()
    state._source_commit_from_hazard(100.0)
    before = float(np.sum(state.accumulated_slip))
    state.evolve(
        dt_s=1.0e-6,
        T_K=700.0,
        stress_Pa=2.0e9,
        b=2.74e-10,
        emission_hazard_integral=0.0,
    )
    after = float(np.sum(state.accumulated_slip))
    assert np.isclose(before, 1.0)
    assert np.isclose(after, before)
