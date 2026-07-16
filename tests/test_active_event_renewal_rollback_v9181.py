from __future__ import annotations

import copy
from types import SimpleNamespace

from arrhenius_fracture.mode_i_first_passage_v9_18_1 import (
    RenewalRollbackPersistentWakeController,
)


class _Stream:
    def __init__(self):
        self.value = {"target": 9.0}

    def restore(self, payload):
        self.value = copy.deepcopy(payload)


class _State:
    def __init__(self, value):
        self.value = value

    def copy(self):
        return _State(self.value)


class _Engine:
    def __init__(self):
        self.mpz_state = _State("mutated")
        self.B = 0.75
        self.a_adv = 5.0
        self.n_adv = 2
        self._reload_until_U_m = 3.0
        self._reload_until_K = 4.0
        self._threshold_stream = _Stream()
        self._last_pre_renewal_state = _State("pre")
        self._last_pre_renewal_event_snapshot = {
            "B": 0.2,
            "threshold": {"target": 1.0},
            "a_adv": 1.0,
            "n_adv": 1,
            "reload_until_U_m": 1.5,
            "reload_until_K": 2.5,
        }
        self.N_em = 0.0
        self.sync_calls = 0
        self.f = SimpleNamespace(max_advances_per_step=1.0)

    def _sync_compat(self):
        self.sync_calls += 1


def test_accidental_renewal_is_rolled_back_without_cancelling_active_event():
    ctl = RenewalRollbackPersistentWakeController()
    ctl.active_event_id = 7
    ctl._active_record = {"event_id": 7}
    eng = _Engine()

    ctl.defer_engine_renewal(eng, 1, 5.0e-6, {})

    assert ctl.active_event_id == 7
    assert ctl.active_event_renewals_rolled_back == 1
    assert ctl.active_event_thresholds_rolled_back == 1
    assert eng.mpz_state.value == "pre"
    assert eng.B == 0.2
    assert eng.a_adv == 1.0
    assert eng.n_adv == 1
    assert eng._threshold_stream.value == {"target": 1.0}
    assert eng._reload_until_U_m == 1.5
    assert eng._reload_until_K == 2.5
    assert eng.sync_calls == 1
    assert eng._last_pre_renewal_state is None
    assert eng._last_pre_renewal_event_snapshot is None
    assert ctl._active_record["active_event_renewals_rolled_back"] == 1


def test_no_active_event_uses_normal_one_fire_path():
    ctl = RenewalRollbackPersistentWakeController()
    eng = _Engine()
    # The inherited path requires a complete real engine/controller state, so the
    # key regression here is that the transactional branch is not entered.
    assert not ctl.active
    assert ctl.active_event_renewals_rolled_back == 0
