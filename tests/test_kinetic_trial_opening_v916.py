from __future__ import annotations

import json
from pathlib import Path

from arrhenius_fracture.coupled_event_audit_v916 import audit_case_v916
from arrhenius_fracture.mode_i_first_passage_v9_16 import KineticTrialEventController


class _State:
    def __init__(self):
        self.advance_calls = []

    def advance(self, distance):
        self.advance_calls.append(float(distance))
        return {
            "wake_mobile": 1.0,
            "wake_retained": 2.0,
            "wake_slip": 3.0,
            "source_sites_refreshed": 4.0,
        }


class _Engine:
    def __init__(self):
        self.mpz_state = _State()
        self.a_adv = 0.0
        self.n_adv = 0
        self.synced = 0

    def _sync_compat(self):
        self.synced += 1


class _Elem:
    def __init__(self, length=5.0e-6):
        self.damage = 1.0
        self.clock = 0.0
        self.length = length
        self.metadata = {
            "v916_intact_reference_energy_J_per_m": 2.0,
            "v916_recoverable_energy_J_per_m": 2.0,
            "v916_opening_jump_m": 1.0e-9,
            "v916_normal_traction_Pa": 1.0e8,
        }


def _finish(ctl, lam, W_emit=0.0):
    ctl.prepare_substep()
    ctl.finish_substep(
        {
            "anisotropic_KJ_Pa_sqrt_m": 12.0e6,
            "mpz_K_shield_pre_renewal_Pa_sqrt_m": 1.0e5,
            "dN_emit": 0.1,
            "dN_emit_raw": 0.1,
            "W_emit": W_emit,
        },
        lambda_c_current=lam,
        lambda_c_raw_current=lam,
        K_cleave=12.0e6,
        T=700.0,
    )


def test_kinetic_event_commits_deferred_mpz_renewal(monkeypatch):
    monkeypatch.setenv("ARRHENIUS_EVENT_RELAXATION_TIME_S", "1")
    monkeypatch.setenv("ARRHENIUS_EVENT_RELAXATION_SUBSTEPS", "4")
    monkeypatch.setenv("ARRHENIUS_EVENT_MIN_RATE_RATIO", "1e-8")
    ctl = KineticTrialEventController()
    eng = _Engine()
    ctl.defer_engine_renewal(eng, 1, 5.0e-6)
    event_id = ctl.schedule_event({"lambda_c": 2.0, "W_emit": 0.0})
    elem = _Elem()
    ctl.register_geometry([elem], [{}])

    for i in range(4):
        _finish(ctl, 2.0, W_emit=0.1 * (i + 1))

    assert event_id == 1
    assert not ctl.active
    assert eng.mpz_state.advance_calls == [5.0e-6]
    assert eng.a_adv == 5.0e-6
    assert eng.n_adv == 1
    assert elem.damage == 1.0
    event = ctl.events[0]
    assert event["status"] == "complete_committed"
    assert event["mpz_renewal_deferred"]
    assert event["mpz_renewal_committed"]
    assert event["cohesive_work_J_per_m"] > 0.0


def test_event_arrests_and_resumes_under_reload(monkeypatch):
    monkeypatch.setenv("ARRHENIUS_EVENT_RELAXATION_TIME_S", "1")
    monkeypatch.setenv("ARRHENIUS_EVENT_RELAXATION_SUBSTEPS", "4")
    monkeypatch.setenv("ARRHENIUS_EVENT_MIN_RATE_RATIO", "0.5")
    monkeypatch.setenv("ARRHENIUS_EVENT_RESUME_RATE_RATIO", "0.6")
    monkeypatch.setenv("ARRHENIUS_EVENT_ARREST_SUBSTEPS", "2")
    ctl = KineticTrialEventController()
    eng = _Engine()
    ctl.defer_engine_renewal(eng, 1, 5.0e-6)
    ctl.schedule_event({"lambda_c": 1.0})
    elem = _Elem()
    ctl.register_geometry([elem], [{}])

    _finish(ctl, 0.1)
    _finish(ctl, 0.1)
    assert ctl.active
    assert ctl.paused
    assert elem.damage < 1.0
    assert not ctl.events[0]["mpz_renewal_committed"]

    resumed = ctl.note_reload_probe(
        lambda_c_current=0.7,
        K_cleave=13.0e6,
        KJ=13.0e6,
        T=700.0,
    )
    assert resumed
    assert not ctl.paused
    assert ctl.events[0]["status"] == "kinetic_relaxation_resumed"


def _complete_event():
    sub = [
        {
            "dt_s": 0.25,
            "remote_displacement_increment_m": 0.0,
            "damage_progress": q,
            "damage_increment": 0.25,
            "lambda_c_current_s-1": r,
            "lambda_c_ratio_to_nucleation": r,
            "cohesive_damage_work_increment_J_per_m": 0.1,
            "cohesive_recoverable_energy_J_per_m": 0.2,
            "tip_emission_work_increment_J_per_m": 0.01,
            "Kshield_Pa_sqrt_m": 1.0e5 + 2.0e3 * i,
        }
        for i, (q, r) in enumerate(zip((0.25, 0.5, 0.75, 1.0), (1.0, 0.9, 0.8, 0.7)))
    ]
    return {
        "event_id": 1,
        "status": "complete_committed",
        "relaxation_completed": True,
        "final_damage": 1.0,
        "relaxation_time_s": 1.0,
        "mpz_renewal_deferred": True,
        "mpz_renewal_committed": True,
        "committed_nfire": 1,
        "committed_distance_m": 5.0e-6,
        "dN_emit_relaxation": 0.4,
        "tip_emission_work_J_per_m": 0.04,
        "cohesive_work_J_per_m": 0.4,
        "substeps": sub,
    }


def test_v916_audit_requires_commit_and_kinetic_rows(tmp_path: Path):
    (tmp_path / "kinetic_trial_event_relaxation_v916.json").write_text(json.dumps({
        "events": [_complete_event()],
        "active_event_id_at_exit": 2,
        "active_event_progress_at_exit": 0.0,
    }))
    out = audit_case_v916(tmp_path, 700.0, 5.0)
    assert out["all_trial_commit_requirements_passed"]
    assert out["target_committed"]
    assert out["kinetic_response_observed"]
    assert out["physical_coupling_relevance_observed"]


def test_v916_audit_rejects_topology_without_physical_commit(tmp_path: Path):
    event = _complete_event()
    event["status"] = "arrested_pending_reload"
    event["relaxation_completed"] = False
    event["mpz_renewal_committed"] = False
    event.pop("committed_distance_m")
    (tmp_path / "kinetic_trial_event_relaxation_v916.json").write_text(json.dumps({
        "events": [event],
        "active_event_id_at_exit": 1,
        "active_event_progress_at_exit": 0.5,
        "active_event_paused_at_exit": True,
    }))
    out = audit_case_v916(tmp_path, 700.0, 5.0)
    assert not out["all_trial_commit_requirements_passed"]
    assert not out["target_committed"]
