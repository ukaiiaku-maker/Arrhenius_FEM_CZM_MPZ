from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from arrhenius_fracture import pf_theta0_stochastic_fingerprint_v10051840 as fp


class _MPZState:
    def __init__(self, mobile, retained, emitted=0.0, escaped=0.0, recovered=0.0,
                 wake_mobile=0.0, wake_retained=0.0):
        self.mobile = np.asarray(mobile, dtype=float)
        self.retained = np.asarray(retained, dtype=float)
        self.emitted_total = emitted
        self.escaped_total = escaped
        self.recovered_total = recovered
        self.wake_mobile_total = wake_mobile
        self.wake_retained_total = wake_retained


def _minimal_engine():
    # Mirrors the base FrontEngine: only B and N_em exist.
    return SimpleNamespace(B=0.0, N_em=0.0)


def _full_engine(seed=8666, threshold=0.42, event_factor=1.0):
    rng = np.random.default_rng(seed)
    return SimpleNamespace(
        B=0.1,
        N_em=5.0,
        _hazard_rng=rng,
        _emission_rngs=[np.random.default_rng(seed + 1), np.random.default_rng(seed + 2)],
        hazard_threshold_action=threshold,
        stochastic_event_length_factor=event_factor,
        mpz_state=_MPZState(mobile=[[1.0, 2.0], [3.0, 4.0]], retained=[[0.5, 0.5], [0.5, 0.5]]),
    )


def test_minimal_engine_reports_not_materialized_for_absent_fields():
    fingerprint = fp.capture_stochastic_fingerprint(_minimal_engine())
    assert fingerprint.hazard_rng_digest == fp.NOT_MATERIALIZED
    assert fingerprint.emission_rng_digest == fp.NOT_MATERIALIZED
    assert fingerprint.cleave_threshold_action == fp.NOT_MATERIALIZED
    assert fingerprint.stochastic_event_length_factor == fp.NOT_MATERIALIZED
    assert fingerprint.mpz.digest == fp.NOT_MATERIALIZED
    # B and N_em DO exist on the minimal engine and must be captured, not
    # reported as absent.
    assert fingerprint.B == 0.0
    assert fingerprint.N_em == 0.0


def test_full_engine_captures_real_values_without_dumping_full_state():
    eng = _full_engine()
    fingerprint = fp.capture_stochastic_fingerprint(eng)
    assert fingerprint.B == 0.1
    assert fingerprint.N_em == 5.0
    assert fingerprint.cleave_threshold_action == 0.42
    assert fingerprint.stochastic_event_length_factor == 1.0
    assert fingerprint.hazard_rng_digest != fp.NOT_MATERIALIZED
    assert len(fingerprint.hazard_rng_digest) == 16  # short digest, not full state
    assert fingerprint.mpz.mobile_total == pytest.approx(10.0)
    assert fingerprint.mpz.retained_total == pytest.approx(2.0)
    serialized = fingerprint.as_dict()
    # The serialized form must never contain a raw bit_generator state blob
    # or a raw MPZ array -- only scalars and short digests.
    assert "bit_generator" not in str(serialized)
    assert "mobile" not in serialized["mpz"] or not isinstance(
        serialized["mpz"].get("mobile"), (list, np.ndarray)
    )


def test_identical_engine_state_produces_matching_fingerprints():
    eng_a = _full_engine(seed=8666)
    eng_b = _full_engine(seed=8666)
    fp_a = fp.capture_stochastic_fingerprint(eng_a)
    fp_b = fp.capture_stochastic_fingerprint(eng_b)
    assert fp_a.matches(fp_b)


def test_advancing_the_rng_changes_the_digest_and_breaks_the_match():
    eng = _full_engine(seed=8666)
    fp_before = fp.capture_stochastic_fingerprint(eng)
    eng._hazard_rng.random()  # simulate a real draw advancing the stream
    fp_after = fp.capture_stochastic_fingerprint(eng)
    assert fp_before.hazard_rng_digest != fp_after.hazard_rng_digest
    assert not fp_before.matches(fp_after)


def test_changing_B_or_N_em_breaks_the_match():
    eng = _full_engine()
    fp_before = fp.capture_stochastic_fingerprint(eng)
    eng.B = 0.99
    fp_after = fp.capture_stochastic_fingerprint(eng)
    assert not fp_before.matches(fp_after)


def test_changing_threshold_or_event_length_factor_breaks_the_match():
    eng = _full_engine()
    fp_before = fp.capture_stochastic_fingerprint(eng)
    eng.stochastic_event_length_factor = 2.5
    fp_after = fp.capture_stochastic_fingerprint(eng)
    assert not fp_before.matches(fp_after)


def test_changing_mpz_population_breaks_the_match():
    eng = _full_engine()
    fp_before = fp.capture_stochastic_fingerprint(eng)
    eng.mpz_state.mobile[0, 0] += 1.0
    fp_after = fp.capture_stochastic_fingerprint(eng)
    assert not fp_before.matches(fp_after)
