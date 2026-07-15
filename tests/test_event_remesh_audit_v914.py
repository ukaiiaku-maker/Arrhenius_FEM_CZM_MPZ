from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from arrhenius_fracture.event_remesh_audit_v914 import audit_campaign


def _case(root: Path, cls: str, curve, yoffset: float = 0.0, complete: bool = True):
    case = root / "seed_1" / "tip_only" / cls / "T700_th45"
    case.mkdir(parents=True)
    pd.DataFrame({
        "crack_extension_m": [0.0, 50e-6, 100e-6 if complete else 5e-6],
        "n_fire": [0.0, 1.0, 1.0],
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
    (case / "v9_14_run_config.json").write_text(json.dumps({
        "schema": "v9.14_event_remesh_material_transfer",
        "event_statistics": "deterministic",
        "stochastic_emission": False,
        "adaptive_event_coordinate": "absolute_integrated_hazard_action",
    }))
    (case / "v9_14_case_summary.json").write_text(json.dumps({
        "schema": "v9.14_case_summary",
        "subprocess_returncode": 0,
    }))
    pd.DataFrame([{
        "n_raw_topology_events": 2,
        "n_independent_load_events": 3,
        "n_unstable_same_load_cascades": 0,
        "fraction_topology_events_in_cascades": 0.0,
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
    for name in (
        "field_snapshots_700K.png",
        "field_snapshots_tip_zoom_700K.png",
        "field_snapshot_manifest_700K.json",
    ):
        (case / name).write_text("x")
    if complete:
        (case / ".long_growth_complete").write_text("")

    czm = case / "czm_0700K"
    czm.mkdir()
    (czm / "event_remesh_audit_v914.json").write_text(json.dumps({
        "n_events": 2,
        "n_failed_event_attempts": 0,
        "all_parent_maps_valid": True,
        "all_patch_targets_satisfied": True,
        "all_preexisting_cohesive_states_unchanged": True,
        "all_events_one_physical_cohesive_event": True,
        "all_post_event_equilibria_completed": True,
        "max_parent_relative_area_conservation_error": 1e-14,
        "max_relative_total_area_error": 1e-15,
    }))
    (case / "event_equilibrium_audit_v914.json").write_text(json.dumps({
        "n_post_event_equilibria": 2,
        "all_same_time": True,
        "all_zero_hazard_increment": True,
        "all_J_recomputed": True,
        "all_MPZ_profiles_recomputed": True,
        "max_relative_boundary_displacement_drift": 0.0,
        "max_relative_rho_area_integral_error": 1e-14,
        "max_relative_ep_area_integral_error": 1e-14,
        "max_relative_total_mesh_area_error": 1e-15,
    }))
    return case


def test_complete_conservative_campaign_passes_numerical_gate(tmp_path):
    _case(tmp_path, "ceramic", (1.0, 1.01, 0.99), 0.0)
    _case(tmp_path, "weakT", (1.0, 1.25, 0.90), 1e-6)
    _case(tmp_path, "DBTT", (1.0, 0.90, 1.35), 2e-6)
    out = audit_campaign(tmp_path, 1, 700.0)
    assert out["numerical_event_remesh_gate_passed"]
    assert not out["failed_numerical_remesh_cases"]
    assert out["material_transfer_gate_passed_v914"]
    assert out["material_audit_v913"]["deterministic_mean_protocol_gate_passed"]
    assert all(case["all_MPZ_profiles_recomputed"] for case in out["cases"])


def test_missing_same_load_equilibrium_fails_numerical_gate(tmp_path):
    _case(tmp_path, "ceramic", (1.0, 1.01, 0.99), 0.0)
    weak = _case(tmp_path, "weakT", (1.0, 1.25, 0.90), 1e-6)
    _case(tmp_path, "DBTT", (1.0, 0.90, 1.35), 2e-6)
    (weak / "event_equilibrium_audit_v914.json").unlink()
    out = audit_campaign(tmp_path, 1, 700.0)
    assert not out["numerical_event_remesh_gate_passed"]
    assert "weakT" in out["failed_numerical_remesh_cases"]
    assert not out["material_transfer_gate_passed_v914"]


def test_missing_post_event_mpz_profile_fails_numerical_gate(tmp_path):
    _case(tmp_path, "ceramic", (1.0, 1.01, 0.99), 0.0)
    weak = _case(tmp_path, "weakT", (1.0, 1.25, 0.90), 1e-6)
    _case(tmp_path, "DBTT", (1.0, 0.90, 1.35), 2e-6)
    path = weak / "event_equilibrium_audit_v914.json"
    payload = json.loads(path.read_text())
    payload["all_MPZ_profiles_recomputed"] = False
    path.write_text(json.dumps(payload))
    out = audit_campaign(tmp_path, 1, 700.0)
    assert not out["numerical_event_remesh_gate_passed"]
    assert "weakT" in out["failed_numerical_remesh_cases"]
