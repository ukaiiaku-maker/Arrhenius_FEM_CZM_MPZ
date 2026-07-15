from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from arrhenius_fracture.material_rcurve_audit_v913 import audit_campaign


def _case(root: Path, cls: str, *, rc: int = 0, complete: bool = True, event_statistics: str = "deterministic", curve=(1.0, 1.1, 1.2), yoffset: float = 0.0):
    case = root / "seed_1" / "tip_only" / cls / "T700_th45"
    case.mkdir(parents=True)
    pd.DataFrame({
        "crack_extension_m": [0.0, 50e-6, 100e-6 if complete else 5e-6],
        "mpz_K_shield_Pa_sqrt_m": [0.0, 1e5, 2e5],
        "mpz_retained_count": [0.0, 1.0, 2.0],
        "mpz_mobile_count": [0.0, 2.0, 1.0],
        "mpz_local_slip_count": [0.0, 1.0, 2.0],
    }).to_csv(case / "steps_0700K.csv", index=False)
    (case / "anisotropic_calibrated_tip_first_passage_summary.json").write_text(json.dumps({
        "control_state": "first_passage" if complete else "right_censored_endpoint",
        "Kc_first_existing_MPa_sqrt_m": 20.0 if complete else None,
    }))
    (case / "rcurve_run_audit.json").write_text(json.dumps({"target_extension_um": 100.0}))
    (case / "v9_13_run_config.json").write_text(json.dumps({
        "event_statistics": event_statistics, "stochastic_emission": False,
    }))
    (case / "v9_13_case_summary.json").write_text(json.dumps({"subprocess_returncode": rc}))
    pd.DataFrame([{
        "n_raw_topology_events": 3, "n_independent_load_events": 3,
        "n_unstable_same_load_cascades": 0, "fraction_topology_events_in_cascades": 0.0,
        "largest_same_load_jump_um": 5.0,
    }]).to_csv(case / "R_curve_cascade_metrics.csv", index=False)
    pd.DataFrame({
        "crack_extension_um": [0.0, 50.0, 100.0],
        "KJ_MPa_sqrt_m": [20.0 * v for v in curve],
    }).to_csv(case / "R_curve_load_events_clustered.csv", index=False)
    pd.DataFrame({
        "crack_extension_after_um": [5.0, 55.0, 105.0],
        "KJ_MPa_sqrt_m": [20.0 * v for v in curve],
    }).to_csv(case / "R_curve_topology_events_raw.csv", index=False)
    pd.DataFrame({
        "x_m": [0.0, 1e-4, 2e-4],
        "y_m": [yoffset, yoffset, yoffset],
    }).to_csv(case / "crack_path_700K.csv", index=False)
    (case / "mpz_state_snapshots_0700K.json").write_text(json.dumps({
        "snapshots": [{"fronts": [{"state": {"emitted_total": 3.0}}]}],
        "final_fronts": [],
    }))
    for name in ("field_snapshots_700K.png", "field_snapshots_tip_zoom_700K.png", "field_snapshot_manifest_700K.json"):
        (case / name).write_text("x")
    if complete:
        (case / ".long_growth_complete").write_text("")
    return case


def test_right_censored_failed_case_cannot_pass(tmp_path):
    _case(tmp_path, "ceramic", rc=1, complete=False)
    _case(tmp_path, "weakT", curve=(1.0, 1.2, 0.9), yoffset=1e-6)
    _case(tmp_path, "DBTT", curve=(1.0, 0.9, 1.3), yoffset=2e-6)
    out = audit_campaign(tmp_path, 1, 700.0)
    assert not out["execution_gate_passed"]
    assert not out["completion_gate_passed"]
    assert not out["material_transfer_gate_passed"]
    assert "ceramic" in out["failed_solver_cases"]
    assert "ceramic" in out["incomplete_cases"]


def test_single_class_does_not_pass_vacuously(tmp_path):
    _case(tmp_path, "ceramic")
    out = audit_campaign(tmp_path, 1, 700.0, classes=["ceramic"])
    assert out["n_pairwise_comparisons"] == 0
    assert not out["pairwise_coverage_gate_passed"]
    assert not out["material_transfer_gate_passed"]


def test_stochastic_realization_is_not_deterministic_gate(tmp_path):
    _case(tmp_path, "ceramic", event_statistics="stochastic")
    _case(tmp_path, "weakT", event_statistics="stochastic", curve=(1.0, 1.2, 0.9), yoffset=1e-6)
    _case(tmp_path, "DBTT", event_statistics="stochastic", curve=(1.0, 0.9, 1.3), yoffset=2e-6)
    out = audit_campaign(tmp_path, 1, 700.0)
    assert not out["deterministic_mean_protocol_gate_passed"]
    assert out["interpretation"] == "stochastic_ensemble_not_deterministic_transfer_gate"


def test_complete_distinct_deterministic_campaign_can_pass(tmp_path):
    _case(tmp_path, "ceramic", curve=(1.0, 1.01, 0.99))
    _case(tmp_path, "weakT", curve=(1.0, 1.25, 0.90), yoffset=1e-6)
    _case(tmp_path, "DBTT", curve=(1.0, 0.90, 1.35), yoffset=2e-6)
    out = audit_campaign(tmp_path, 1, 700.0)
    assert out["execution_gate_passed"]
    assert out["completion_gate_passed"]
    assert out["field_output_gate_passed"]
    assert out["pairwise_coverage_gate_passed"]
    assert out["deterministic_mean_protocol_gate_passed"]
    assert out["material_transfer_gate_passed"]
