from __future__ import annotations

import copy
from types import SimpleNamespace

import numpy as np

from arrhenius_fracture.emergent_gnd_campaign_v913 import candidate_from_registry_row
from arrhenius_fracture.emergent_gnd_state_v913 import EmergentGNDState
from reduced_fracture_v2.spatial_transfer import (
    ObservedEmergentGNDState,
    physical_onset_event_indices,
    state_snapshot,
)
from reduced_fracture_v2.taylor_peierls_contract import validate_candidate
from scripts.run_oneD_v2_predictive_campaign import inputs
from scripts.run_oneD_v2_taylor_peierls_spatial_transfer import (
    CONTROL_ID,
    _screen_pool,
)


def _states():
    physics, _, _ = inputs()
    row = _screen_pool().query("candidate_id == @CONTROL_ID").iloc[0]
    candidate = candidate_from_registry_row(row)
    return EmergentGNDState(candidate, physics), ObservedEmergentGNDState(candidate, physics)


def test_observational_wake_ledger_does_not_change_active_state() -> None:
    base, observed = _states()
    rng = np.random.default_rng(91)
    for name in ("mobile_m2", "retained_m2", "accumulated_slip_m2"):
        values = rng.uniform(0.0, 1.0e12, size=getattr(base, name).shape)
        setattr(base, name, values.copy())
        setattr(observed, name, values.copy())
    base.translate_tip(0.37 * base.dx)
    observed.translate_tip(0.37 * observed.dx)
    for name in ("mobile_m2", "retained_m2", "accumulated_slip_m2"):
        assert np.array_equal(getattr(base, name), getattr(observed, name))
    assert observed.observer_mobile_departed_line_content >= 0.0
    assert observed.observer_retained_departed_line_content >= 0.0


def test_state_snapshot_is_neutral_and_contains_required_spatial_features() -> None:
    _, state = _states()
    before = copy.deepcopy(state)
    record = state_snapshot(
        "INITIAL",
        {"event_index": -1, "projected_extension_m": 0.0, "K_MPa_sqrt_m": 0.0},
        state,
        temperature_K=1100.0,
    )
    assert record["observer_feedback"] is False
    for field in (
        "mobile_line_content", "retained_line_content", "retained_fraction",
        "mobile_centroid_tip_relative_m", "retained_width_m",
        "near_tip_retained_line_content", "wake_retained_departed_line_content",
        "tip_radius_m", "front_width_m", "backstress_Pa",
        "K_shield_MPa_sqrt_m", "source_multiplicity",
        "peierls_transport_velocity_abs_max_m_s",
        "encounter_retention_rate_max_s", "taylor_completion_rate_max_s",
    ):
        assert field in record and np.isfinite(record[field])
    assert np.array_equal(state.mobile_m2, before.mobile_m2)
    assert np.array_equal(state.retained_m2, before.retained_m2)


def test_physical_onset_latch_does_not_reopen_after_unresolved_pair() -> None:
    events = [
        SimpleNamespace(applied_displacement_m=value * 1.0e-6)
        for value in (1.0, 15.0, 16.0, 30.0)
    ]
    onsets, avalanches = physical_onset_event_indices(
        events, reload_gap_threshold_m=10.0e-6
    )
    assert onsets == [0, 1]
    assert avalanches == [0, 1, 1, 1]


def test_spatial_pool_is_fully_evaluated_dbtt_and_keeps_barriers_fixed() -> None:
    pool = _screen_pool()
    assert len(pool) == 513 and pool.candidate_id.is_unique
    control = pool[pool.candidate_id.eq(CONTROL_ID)].iloc[0]
    for _, row in pool.sample(32, random_state=83).iterrows():
        validate_candidate(row, control)
        assert row.cleavage_barrier_sha256 == control.cleavage_barrier_sha256
        assert row.emission_barrier_sha256 == control.emission_barrier_sha256
