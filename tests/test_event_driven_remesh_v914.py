from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

import arrhenius_fracture.mode_i_first_passage_v9_14 as mode_i_v914
from arrhenius_fracture.mode_i_first_passage_v9_14 import (
    _ControlledScalar,
    _PostEventEquilibriumController,
    _forward_mode_i_plane,
)
from arrhenius_fracture.remesh_audit_v914 import audit_case
import run_mpz_v9_14_event_driven_remesh as runner


def test_forward_plane_is_exact_mode_i():
    plane = _forward_mode_i_plane(45.0)[0]
    assert list(plane["t"]) == [1.0, 0.0]
    assert list(plane["n"]) == [0.0, 1.0]


def test_post_event_controller_zeroes_one_dt_and_du_pair():
    ctl = _PostEventEquilibriumController()
    dt = _ControlledScalar(8.4, ctl, "dt")
    du = _ControlledScalar(2e-7, ctl, "dU")
    assert dt * 0.5 == 4.2
    assert du * 0.5 == 1e-7
    event_id = ctl.schedule_physical_event()
    assert event_id == 1
    # Repeated low-level calls while the event is pending must not create a new
    # physical event or a second correction requirement.
    assert ctl.schedule_physical_event() == 1
    assert ctl.events_scheduled == 1
    assert ctl.duplicate_schedule_calls == 1
    assert dt * 1.0 == 0.0
    assert du * 1.0 == 0.0
    assert not ctl.pending
    assert ctl.corrections_consumed == 1
    assert ctl.corrected_event_ids == [1]
    assert dt * 0.5 == 4.2


def test_geometry_veto_cancels_pending_physical_event():
    ctl = _PostEventEquilibriumController()
    assert ctl.schedule_physical_event() == 1
    assert ctl.cancel_pending_event() == 1
    assert not ctl.pending
    assert ctl.events_cancelled == 1
    assert ctl.cancelled_event_ids == [1]


def test_v914_entry_injects_only_supported_branch_controls(monkeypatch, tmp_path):
    captured = {}

    def fake_main(argv):
        captured["argv"] = list(argv)
        return []

    monkeypatch.setattr(mode_i_v914._base, "main", fake_main)
    mode_i_v914.main(["--out", str(tmp_path)])
    argv = captured["argv"]
    assert "--crystal-aniso" in argv
    assert "--crack-backend" in argv
    assert "adaptive_czm" in argv
    assert "--adaptive-events" in argv
    assert "--no-crystal-branch" not in argv
    assert "--crystal-branch" not in argv


class _Args:
    parameter_root = Path("mpz_v9_11_parameters")
    T_K = 700.0
    target_extension_um = 25.0
    steps = 100
    nx = 18
    ny = 36
    tip_h_fine = 1e-6
    tip_ratio = 1.2
    dU = 2e-7
    dt = 8.4
    n_stagger = 2
    print_every = 25
    adaptive_event_target = 0.15
    da_phys_um = 5.0
    mpz_length_um = 100.0
    mpz_n_bins = 200
    crystal_theta_deg = 45.0
    save_snapshots = 2
    snapshot_cols = 2
    snapshot_by_extension_um = 10.0


def test_v914_command_selects_driver_and_adds_one_guard_increment():
    cmd = runner._build_command_v914(_Args(), "DBTT", Path("runs/x"), True)
    assert cmd[1] == "run_mpz_v9_14_mode_i_rcurve.py"
    i = cmd.index("--target-extension-um")
    assert float(cmd[i + 1]) == 30.0


def _write_diagnostics(tmp_path: Path, rows: list[dict]):
    out = tmp_path / "czm_0700K"
    out.mkdir()
    (out / "czm_advance_log.json").write_text(json.dumps(rows))
    (out / "cohesive_elements.csv").write_text("index,front_id\n")
    (tmp_path / "field_snapshot_manifest_700K.json").write_text("{}")


def test_audit_rejects_missing_same_load_correction(tmp_path):
    pd.DataFrame({
        "step": [1, 2, 3],
        "Uapp_m": [1.0e-4, 1.1e-4, 1.2e-4],
        "n_fire": [0, 1, 0],
        "adaptive_frac": [1.0, 1e-8, 1.0],
        "dt_cur_s": [8.4, 1e-8, 8.4],
        "crack_extension_m": [0.0, 5e-6, 5e-6],
    }).to_csv(tmp_path / "steps_0700K.csv", index=False)
    _write_diagnostics(tmp_path, [
        {"physical_event_id": 1, "physical_subsegment_index": 0, "length_m": 2e-6},
        {"physical_event_id": 1, "physical_subsegment_index": 1, "length_m": 3e-6},
    ])
    (tmp_path / "post_event_equilibrium_audit_v914.json").write_text(json.dumps({
        "events_scheduled": 1,
        "corrections_consumed": 0,
        "corrected_event_ids": [],
        "cancelled_event_ids": [],
        "pending_event_id": 1,
    }))
    out = audit_case(tmp_path, 700.0, analysis_target_extension_um=5.0)
    assert out["requirements_1_to_4_passed"]
    assert out["n_czm_advances"] == 1
    assert out["n_czm_subsegments_in_analysis"] == 2
    assert not out["requirement_5_passed"]


def test_audit_accepts_two_corrected_events_and_one_pending_guard(tmp_path):
    pd.DataFrame({
        "step": [1, 2, 3, 4, 5, 6, 7],
        "Uapp_m": [1.0e-4, 1.1e-4, 1.1e-4, 1.2e-4, 1.2e-4, 1.3e-4, 1.3e-4],
        "n_fire": [0, 1, 0, 1, 0, 1, 0],
        "adaptive_frac": [1.0, 1e-8, 1.0, 1e-8, 1.0, 1e-8, 1.0],
        "dt_cur_s": [8.4, 1e-8, 0.0, 1e-8, 0.0, 1e-8, 0.0],
        "crack_extension_m": [0.0, 5e-6, 5e-6, 10e-6, 10e-6, 15e-6, 15e-6],
    }).to_csv(tmp_path / "steps_0700K.csv", index=False)
    _write_diagnostics(tmp_path, [
        {"physical_event_id": 1, "physical_subsegment_index": 0, "length_m": 2e-6},
        {"physical_event_id": 1, "physical_subsegment_index": 1, "length_m": 3e-6},
        {"physical_event_id": 2, "physical_subsegment_index": 0, "length_m": 5e-6},
        {"physical_event_id": 3, "physical_subsegment_index": 0, "length_m": 5e-6},
    ])
    (tmp_path / "post_event_equilibrium_audit_v914.json").write_text(json.dumps({
        "events_scheduled": 3,
        "corrections_consumed": 2,
        "events_cancelled": 0,
        "scheduled_event_ids": [1, 2, 3],
        "corrected_event_ids": [1, 2],
        "cancelled_event_ids": [],
        "pending_event_id": 3,
        "pending_at_exit": True,
    }))
    out = audit_case(tmp_path, 700.0, analysis_target_extension_um=10.0)
    assert out["n_physical_events_in_analysis"] == 2
    assert out["n_czm_subsegments_in_analysis"] == 3
    assert out["n_guard_fire_rows"] == 1
    assert out["guard_event_ids"] == [3]
    assert out["terminal_pending_guard_event_allowed"]
    assert out["requirement_5_passed"]
    assert out["all_five_requirements_passed"]
