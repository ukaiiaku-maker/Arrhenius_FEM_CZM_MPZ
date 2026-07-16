from __future__ import annotations

import numpy as np

from arrhenius_fracture.moving_process_zone import MovingProcessZoneConfig
from arrhenius_fracture.moving_process_zone_v918 import MovingProcessZoneState
from arrhenius_fracture.mode_i_first_passage_v9_18 import (
    PersistentWakeHazardController,
)


def _cfg():
    return MovingProcessZoneConfig(
        length_m=10.0e-6,
        n_bins=10,
        n_systems=1,
        source_sites_per_system=10.0,
        source_refresh_length_m=5.0e-6,
        source_bin_count=1,
        shielding_orientation_factors=(1.0,),
        mobile_shield_fraction=0.0,
        blunting_length_m=1.0e-6,
    )


def test_retained_and_slip_cross_into_wake_conservatively(monkeypatch):
    monkeypatch.setenv("ARRHENIUS_WAKE_LENGTH_UM", "10")
    monkeypatch.setenv("ARRHENIUS_WAKE_N_BINS", "10")
    state = MovingProcessZoneState(_cfg())
    state.retained[0, 0] = 1.0
    state.accumulated_slip[0, 0] = 2.0

    before_retained = state.retained_count + state.wake_retained_count
    before_slip = float(np.sum(state.accumulated_slip)) + state.wake_slip_count
    active_K_before = state.active_K_shielding(160.0e9, 0.28, 2.74e-10)

    out = state.advance(2.0e-6)

    after_retained = state.retained_count + state.wake_retained_count
    after_slip = float(np.sum(state.accumulated_slip)) + state.wake_slip_count
    assert abs(after_retained - before_retained) < 1.0e-12
    assert abs(after_slip - before_slip) < 1.0e-12
    assert state.retained_count == 0.0
    assert state.wake_retained_count == 1.0
    assert out["wake_retained"] == 1.0
    assert out["wake_retained_discarded"] == 0.0
    assert active_K_before > 0.0
    assert state.active_K_shielding(160.0e9, 0.28, 2.74e-10) == 0.0
    assert state.wake_K_shielding(160.0e9, 0.28, 2.74e-10) > 0.0
    # Wake slip is stored but is not included in the current-tip blunting measure.
    assert state.local_slip_count() == 0.0
    assert state.wake_slip_count == 2.0


def test_wake_state_split_is_conservative(monkeypatch):
    monkeypatch.setenv("ARRHENIUS_WAKE_LENGTH_UM", "10")
    state = MovingProcessZoneState(_cfg())
    state.wake_retained[0, 0] = 4.0
    state.wake_mobile[0, 1] = 2.0
    child = state.split(0.25)
    assert abs(state.wake_retained_count - 3.0) < 1.0e-12
    assert abs(child.wake_retained_count - 1.0) < 1.0e-12
    assert abs(state.wake_mobile_count - 1.5) < 1.0e-12
    assert abs(child.wake_mobile_count - 0.5) < 1.0e-12


def test_nominal_loading_dt_not_adaptive_fraction_controls_hold(monkeypatch):
    monkeypatch.setenv("ARRHENIUS_NOMINAL_LOADING_DT_S", "8.4")
    monkeypatch.setenv("ARRHENIUS_EVENT_MAX_FIXED_HOLD_S", "inf")
    ctl = PersistentWakeHazardController()
    ctl.external_dt_s = 8.4e-8
    assert abs(ctl._hold_cap_s() - 8.4) < 1.0e-12
