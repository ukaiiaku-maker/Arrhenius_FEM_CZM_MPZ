from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from arrhenius_fracture.coupled_event_audit_v917 import audit_case_v917
from arrhenius_fracture.mode_i_first_passage_v9_17 import HazardClockTrialEventController


class _State:
    def __init__(self):
        self.cfg = SimpleNamespace(source_refresh_length_m=54.8e-6)
        self.dx = 0.5e-6
        self.site_capacity = np.array([14.0, 14.0])
        self.available_sites = np.array([0.0, 0.0])
        self.advance_calls = []

    def advance(self, distance):
        self.advance_calls.append(float(distance))
        self.available_sites[:] = self.site_capacity
        return {
            "wake_mobile": 1.0,
            "wake_retained": 2.0,
            "wake_slip": 3.0,
            "source_sites_refreshed": 28.0,
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


def _finish(ctl, lam, dN=0.0, retained=0.0):
    ctl.prepare_substep()
    ctl.finish_substep(
        {
            "anisotropic_KJ_Pa_sqrt_m": 12.0e6,
            "mpz_K_shield_pre_renewal_Pa_sqrt_m": retained * 1.0e5,
            "lambda_e": 1.0e8,
            "dN_emit": dN,
            "dN_emit_raw": dN,
            "mpz_mobile_count": dN,
            "mpz_retained_count": retained,
            "mpz_emitted_total": dN,
            "mpz_available_site_fraction": 0.1,
            "W_emit": dN,
        },
        lambda_c_current=lam,
        lambda_c_raw_current=1.0e4,
        K_cleave=12.0e6,
        T=700.0,
    )


def test_absolute_hazard_clock_completes_in_inverse_rate_time(monkeypatch):
    monkeypatch.setenv("ARRHENIUS_EVENT_TARGET_DQ", "0.25")
    ctl = HazardClockTrialEventController()
    ctl.external_dt_s = 8.4
    eng = _Engine()
    ctl.defer_engine_renewal(eng, 1, 5.0e-6)
    ctl.schedule_event({"lambda_c": 0.1, "lambda_c_raw": 1.0e4, "W_emit": 0.0})
    elem = _Elem()
    ctl.register_geometry([elem], [{}])

    for _ in range(4):
        _finish(ctl, 0.1, dN=0.2, retained=0.1)

    assert not ctl.active
    event = ctl.events[0]
    assert event["status"] == "complete_committed"
    assert abs(event["relaxation_time_s"] - 10.0) < 1.0e-9
    assert abs(event["opening_hazard_integral"] - 1.0) < 1.0e-12
    assert event["source_sites_refreshed_during_opening"] > 0.0
    assert eng.mpz_state.advance_calls == [5.0e-6]
    assert eng.n_adv == 1
    # advance() attempted a full legacy refresh, but v9.17 restored the
    # incrementally refreshed availability instead of double refreshing.
    assert np.all(eng.mpz_state.available_sites < eng.mpz_state.site_capacity)


def test_hazard_timescale_resumes_external_loading_without_rate_ratio_threshold(monkeypatch):
    monkeypatch.setenv("ARRHENIUS_EVENT_TARGET_DQ", "0.1")
    ctl = HazardClockTrialEventController()
    ctl.external_dt_s = 8.4
    eng = _Engine()
    ctl.defer_engine_renewal(eng, 1, 5.0e-6)
    ctl.schedule_event({"lambda_c": 1.0e-3, "lambda_c_raw": 1.0, "W_emit": 0.0})
    ctl.register_geometry([_Elem()], [{}])

    _finish(ctl, 1.0e-3)
    assert ctl.active
    assert ctl.paused
    assert 0.0 < ctl.progress < 0.1

    resumed = ctl.note_reload_probe(
        lambda_c_current=0.02,
        K_cleave=13.0e6,
        KJ=13.0e6,
        T=700.0,
    )
    assert resumed
    assert not ctl.paused


def test_v917_audit_checks_absolute_hazard_and_refresh(tmp_path: Path):
    rows = []
    for i, q in enumerate((0.25, 0.5, 0.75, 1.0), 1):
        rows.append({
            "substep": i,
            "dt_s": 2.5,
            "damage_progress": q,
            "cleavage_hazard_increment": 0.25,
            "remote_displacement_increment_m": 0.0,
            "lambda_c_current_s-1": 0.1,
            "lambda_c_raw_current_s-1": 1.0e4,
            "lambda_e_current_s-1": 1.0e8,
            "dN_emit": 0.2,
            "retained_count": 0.1,
            "Kshield_Pa_sqrt_m": i * 2.0e3,
        })
    event = {
        "event_id": 1,
        "status": "complete_committed",
        "progress_law": "dq_dt=lambda_c_absolute",
        "relaxation_completed": True,
        "final_damage": 1.0,
        "opening_hazard_integral": 1.0,
        "relaxation_time_s": 10.0,
        "mpz_renewal_deferred": True,
        "mpz_renewal_committed": True,
        "committed_nfire": 1,
        "committed_distance_m": 5.0e-6,
        "source_sites_refreshed_during_opening": 1.0,
        "dN_emit_relaxation": 0.8,
        "substeps": rows,
    }
    (tmp_path / "absolute_hazard_event_relaxation_v917.json").write_text(
        json.dumps({"events": [event]})
    )
    out = audit_case_v917(tmp_path, 700.0, 5.0)
    assert out["software_commit_gate_passed"]
    assert out["hazard_clock_consistency_passed"]
    assert out["source_refresh_during_opening_observed"]
    assert out["physical_coupling_relevance_observed"]
    assert out["retained_state_at_commit_observed"]
