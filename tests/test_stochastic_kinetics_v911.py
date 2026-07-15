from __future__ import annotations

import os

import numpy as np

from arrhenius_fracture.config import FractureBarrier
from arrhenius_fracture.moving_process_zone import MovingProcessZoneConfig
from arrhenius_fracture.moving_process_zone_v911 import MovingProcessZoneState
from arrhenius_fracture.mpz_front_engine_v911 import MovingProcessZone2DFrontEngine
from arrhenius_fracture.sharp_front import FrontConfig
from arrhenius_fracture.stochastic_kinetics_v911 import (
    HazardThresholdStream,
    sample_effective_binomial,
)


def test_deterministic_threshold_stream_matches_legacy_unit_action():
    s = HazardThresholdStream("deterministic", seed=1)
    assert s.target == 1.0
    residual, n, crossed = s.consume(2.4, max_events=1)
    assert n == 1
    assert crossed == [1.0]
    assert np.isclose(residual, 1.4)
    assert s.target == 1.0


def test_stochastic_threshold_stream_is_reproducible_and_restorable():
    a = HazardThresholdStream("stochastic", seed=17, stream=3)
    b = HazardThresholdStream("stochastic", seed=17, stream=3)
    c = HazardThresholdStream("stochastic", seed=18, stream=3)
    assert np.isclose(a.target, b.target)
    assert not np.isclose(a.target, c.target)

    snap = a.snapshot()
    residual1, n1, crossed1 = a.consume(a.target + 0.25, max_events=1)
    next1 = a.target
    a.restore(snap)
    residual2, n2, crossed2 = a.consume(a.target + 0.25, max_events=1)
    assert n1 == n2 == 1
    assert crossed1 == crossed2
    assert np.isclose(residual1, residual2)
    assert np.isclose(next1, a.target)


def test_effective_binomial_preserves_expected_site_count():
    rng = np.random.default_rng(123)
    vals = [sample_effective_binomial(rng, 12.5, 0.2) for _ in range(20000)]
    assert abs(np.mean(vals) - 2.5) < 0.05
    assert min(vals) >= 0.0
    assert max(vals) <= 12.5


def test_stochastic_finite_site_emission_reproducible_by_seed(monkeypatch):
    monkeypatch.setenv("ARRHENIUS_EVENT_STATISTICS", "stochastic")
    monkeypatch.setenv("ARRHENIUS_STOCHASTIC_EMISSION", "1")
    monkeypatch.setenv("ARRHENIUS_STOCHASTIC_SEED", "42")
    cfg = MovingProcessZoneConfig(
        n_bins=8,
        n_systems=2,
        source_sites_per_system=20.0,
        source_bin_count=2,
    )
    a = MovingProcessZoneState(cfg)
    b = MovingProcessZoneState(cfg)
    ea = a._source_commit_from_hazard(0.25)
    eb = b._source_commit_from_hazard(0.25)
    assert np.array_equal(ea, eb)
    assert np.sum(ea) == a.emitted_total
    assert np.all(a.available_sites <= a.site_capacity)


def test_event_reload_gate_holds_only_cleavage_until_load_increases(monkeypatch):
    monkeypatch.setenv("ARRHENIUS_EVENT_STATISTICS", "stochastic")
    monkeypatch.setenv("ARRHENIUS_STOCHASTIC_SEED", "5")
    monkeypatch.setenv("ARRHENIUS_PROPAGATION_CONTROL", "event_reload")
    f = FrontConfig()
    f.sigma_cap = 0.0
    cfg = MovingProcessZoneConfig(n_bins=8, source_sites_per_system=0.0)
    eng = MovingProcessZone2DFrontEngine(
        f, FractureBarrier(), FractureBarrier(), 160e9, 0.28, 2.74e-10, cfg
    )
    eng._reload_until_K = 10.0
    assert eng._reload_gate_active(5.0)
    assert eng.predict_clock_increment_drives(5.0, 5.0, 700.0, 1.0) == 0.0
    assert not eng._reload_gate_active(10.1)
