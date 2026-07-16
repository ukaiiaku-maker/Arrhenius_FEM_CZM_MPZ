from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from arrhenius_fracture.coupled_event_audit_v915 import audit_case_v915
from arrhenius_fracture.mode_i_first_passage_v9_15 import (
    CoupledEventRelaxationController,
    _LoadingProxy,
)


class _Elem:
    def __init__(self, length=2.5e-6):
        self.damage = 1.0
        self.clock = 0.0
        self.length = length
        self.metadata = {}


def test_controller_uses_positive_time_and_zero_remote_increment(monkeypatch):
    monkeypatch.setenv("ARRHENIUS_EVENT_RELAXATION_TIME_S", "1e-6")
    monkeypatch.setenv("ARRHENIUS_EVENT_RELAXATION_SUBSTEPS", "4")
    ctl = CoupledEventRelaxationController()
    event_id = ctl.schedule_event({"anisotropic_KJ_Pa_sqrt_m": 12e6})
    elems = [_Elem(), _Elem()]
    rows = [{}, {}]
    ctl.register_geometry(elems, rows)

    class Loading:
        dt = 8.4
        dU_top = 2e-7

    loading = _LoadingProxy(Loading(), ctl)
    for k in range(4):
        assert loading.dt * 1e-8 == 2.5e-7
        assert loading.dU_top * 1e-8 == 0.0
        ctl.finish_substep({
            "anisotropic_KJ_Pa_sqrt_m": 12e6 - k * 1e5,
            "dN_emit": 0.25,
        })
    assert event_id == 1
    assert not ctl.active
    assert all(e.damage == 1.0 for e in elems)
    rec = ctl.events[0]
    assert rec["relaxation_completed"]
    assert rec["n_relaxation_substeps_completed"] == 4
    assert rec["dN_emit_relaxation"] == 1.0
    assert all(r["physical_event_id"] == 1 for r in rows)


def test_audit_accepts_coupled_analysis_event_and_pending_guard(tmp_path: Path):
    pd.DataFrame({
        "step": [1, 2, 3, 4],
        "n_fire": [1, 0, 1, 0],
        "crack_extension_m": [5e-6, 5e-6, 10e-6, 10e-6],
    }).to_csv(tmp_path / "steps_0700K.csv", index=False)
    d = tmp_path / "czm_0700K"
    d.mkdir()
    (d / "cohesive_elements.csv").write_text("index\n")
    (d / "czm_advance_log.json").write_text(json.dumps([
        {"physical_event_id": 1, "length_m": 5e-6},
        {"physical_event_id": 2, "length_m": 5e-6},
    ]))
    (tmp_path / "field_snapshot_manifest_700K.json").write_text("{}")
    sub = [
        {
            "dt_s": 2.5e-7,
            "remote_displacement_increment_m": 0.0,
            "KJ_Pa_sqrt_m": 12e6,
            "dN_emit": 0.1,
        }
        for _ in range(4)
    ]
    (tmp_path / "coupled_event_relaxation_v915.json").write_text(json.dumps({
        "event_relaxation_substeps": 4,
        "active_event_id_at_exit": 3,
        "geometry_retry_attempts": 1,
        "geometry_retry_successes": 1,
        "events": [
            {
                "event_id": 1,
                "status": "complete",
                "relaxation_completed": True,
                "n_relaxation_substeps_completed": 4,
                "relaxation_time_s": 1e-6,
                "final_damage": 1.0,
                "dN_emit_relaxation": 0.4,
                "substeps": sub,
            },
            {
                "event_id": 2,
                "status": "complete",
                "relaxation_completed": True,
                "n_relaxation_substeps_completed": 4,
                "relaxation_time_s": 1e-6,
                "final_damage": 1.0,
                "dN_emit_relaxation": 0.4,
                "substeps": sub,
            },
            {"event_id": 3, "status": "relaxation_pending", "substeps": []},
        ],
    }))
    out = audit_case_v915(tmp_path, 700.0, 10.0)
    assert out["all_coupled_event_requirements_passed"]
    assert out["n_physical_events_in_analysis"] == 2
    assert out["active_guard_event_at_exit"] == 3


def test_audit_rejects_zero_time_operator_split(tmp_path: Path):
    pd.DataFrame({
        "step": [1],
        "n_fire": [1],
        "crack_extension_m": [5e-6],
    }).to_csv(tmp_path / "steps_0700K.csv", index=False)
    d = tmp_path / "czm_0700K"
    d.mkdir()
    (d / "cohesive_elements.csv").write_text("index\n")
    (d / "czm_advance_log.json").write_text(json.dumps([
        {"physical_event_id": 1, "length_m": 5e-6},
    ]))
    (tmp_path / "field_snapshot_manifest_700K.json").write_text("{}")
    (tmp_path / "coupled_event_relaxation_v915.json").write_text(json.dumps({
        "event_relaxation_substeps": 1,
        "events": [{
            "event_id": 1,
            "status": "complete",
            "relaxation_completed": True,
            "n_relaxation_substeps_completed": 1,
            "relaxation_time_s": 0.0,
            "final_damage": 1.0,
            "substeps": [{"dt_s": 0.0, "remote_displacement_increment_m": 0.0}],
        }],
    }))
    out = audit_case_v915(tmp_path, 700.0, 5.0)
    assert not out["all_coupled_event_requirements_passed"]
