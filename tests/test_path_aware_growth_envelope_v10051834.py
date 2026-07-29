from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pytest

from arrhenius_fracture.config import GeometryConfig, MeshConfig
from arrhenius_fracture import mode_i_first_passage_v9_18_5_2 as v91852
from arrhenius_fracture import mode_i_first_passage_v9_18_5_3 as v91853
from arrhenius_fracture import mode_i_first_passage_v10_0_5_13_2_barrier_only as v1005132
from arrhenius_fracture import mode_i_first_passage_v10_0_5_18_3_4_four_class_path_aware_growth_envelope as entry
from arrhenius_fracture import path_aware_growth_envelope_v10051834 as envelope
from arrhenius_fracture import physical_refinement_mesh_v100510 as physical


def _set_env(monkeypatch, target_um="1000", theta="30", minimum_forward="0.05"):
    monkeypatch.setenv("ARRHENIUS_COMMITTED_TARGET_EXTENSION_UM", target_um)
    monkeypatch.setenv("ARRHENIUS_PHYSICAL_DA_UM", "5")
    monkeypatch.setenv("ARRHENIUS_CORRIDOR_PROCESS_ZONE_UM", "50")
    monkeypatch.setenv("ARRHENIUS_CORRIDOR_GUARD_UM", "10")
    monkeypatch.setenv("ARRHENIUS_CRYSTAL_THETA_DEG", theta)
    monkeypatch.setenv("ARRHENIUS_MIN_GLOBAL_FORWARD", minimum_forward)
    monkeypatch.setenv("ARRHENIUS_MAX_CORRIDOR_H_OVER_LPZ", "0.25")
    monkeypatch.setenv("ARRHENIUS_MIN_INITIAL_TRIANGLE_QUALITY", "0.035")
    monkeypatch.setenv("ARRHENIUS_MAX_TIP_H_OVER_DA", "0.75")
    monkeypatch.setenv("ARRHENIUS_ENVELOPE_SUPPORT_H_OVER_LPZ", "0.245")
    monkeypatch.setenv("ARRHENIUS_ENVELOPE_SUPPORT_CANDIDATES", "4")
    monkeypatch.setenv("ARRHENIUS_ENVELOPE_DIRECTION_SAMPLES", "181")
    monkeypatch.setenv("ARRHENIUS_ENVELOPE_BUFFER_LPZ", "1.0")
    monkeypatch.setenv("ARRHENIUS_ENVELOPE_AUDIT_SAMPLE_LPZ", "0.5")
    monkeypatch.setenv("ARRHENIUS_ENVELOPE_AUDIT_NEAREST_ELEMENTS", "16")


def _build_production_envelope(monkeypatch, target_um=1000.0):
    _set_env(monkeypatch, str(target_um))
    geom = GeometryConfig()
    mesh_cfg = MeshConfig(
        nx=36,
        ny=72,
        tip_h_fine=2.5e-6,
        tip_ratio=1.15,
    )
    physical.configure_physical_refinement_v100510(330.0e-6)
    envelope.path_aware_long_growth_envelope_mesh._original = (
        physical.make_physical_refinement_mesh_v100510
    )
    try:
        mesh = envelope.path_aware_long_growth_envelope_mesh(
            geom,
            mesh_cfg,
            seed=42,
        )
        audit = dict(v91852._STARTUP_AUDIT)
    finally:
        physical.clear_physical_refinement_v100510()
    return geom, mesh, audit


def test_global_forward_envelope_covers_continuous_production_cone():
    directions = envelope._global_forward_directions(0.05, 181)
    assert directions.shape == (181, 2)
    assert np.min(directions[:, 0]) == pytest.approx(0.05)
    assert np.max(directions[:, 0]) == pytest.approx(1.0)
    assert np.max(np.abs(directions[:, 1])) == pytest.approx(
        np.sqrt(1.0 - 0.05**2)
    )


def test_1000um_path_aware_envelope_resolves_extreme_and_switched_paths(
    monkeypatch,
):
    geom, mesh, audit = _build_production_envelope(monkeypatch, 1000.0)

    assert audit["schema"] == envelope.CORRIDOR_SCHEMA
    assert audit["path_aware_growth_envelope"] is True
    assert audit["directional_support_basis"].startswith(
        "complete_global_forward_gate"
    )
    assert audit["minimum_global_forward_cosine"] == pytest.approx(0.05)
    assert audit["production_physical_provider_detected"] is True
    assert audit["path_aware_physical_refinement_active"] is True
    assert audit["swept_physical_refinement_active"] is False
    assert (
        audit["path_aware_physical_refinement_policy"]
        == envelope.PATH_AWARE_PHYSICAL_POLICY
    )
    assert audit["full_requested_reachable_envelope_covered"] is True
    assert audit["minimum_initial_triangle_quality"] >= 0.035
    assert audit["maximum_sampled_hbar_tip_over_L_pz"] <= 0.25
    assert audit["tip_h_over_da_enforced_as_veto"] is False
    assert audit["reachable_envelope_y_max_m"] > 330.0e-6
    assert audit["reachable_envelope_y_min_m"] < -330.0e-6
    assert audit["reachable_envelope_y_max_m"] == pytest.approx(
        (1010.0e-6) * np.sqrt(1.0 - 0.05**2),
        rel=1.0e-6,
    )
    assert mesh.nn < 25000
    assert mesh.ne < 50000
    assert np.all(np.isfinite(mesh.area_e))
    assert np.all(mesh.area_e > 0.0)

    vertices = np.asarray(audit["reachable_envelope_vertices_m"], dtype=float)
    hull = envelope.ConvexHull(vertices)
    resolution = envelope._path_envelope_resolution(
        mesh,
        hull,
        5.0e-6,
        50.0e-6,
    )
    assert resolution["maximum_sampled_hbar_tip_over_L_pz"] <= 0.25

    length = 1010.0e-6
    extreme = np.sqrt(1.0 - 0.05**2)
    start = np.array([geom.a0, 0.0])
    points = np.asarray(
        [
            start + length * np.array([0.05, extreme]),
            start + length * np.array([0.05, -extreme]),
            start + length * np.array([1.0, 0.0]),
            start
            + 0.45 * length * np.array([0.05, extreme])
            + 0.55 * length * np.array([0.05, -extreme]),
        ]
    )
    excess = envelope._signed_support_distance(points, hull)
    assert np.max(excess) <= 1.0e-12


def test_nominal_bcc_trace_audit_uses_theta_without_restricting_support(
    monkeypatch,
):
    _, _, audit = _build_production_envelope(monkeypatch, 400.0)
    traces = audit["nominal_BCC_100_trace_audit"]
    angles = sorted(round(float(row["angle_deg"]), 6) for row in traces)
    assert angles == [-60.0, 30.0]
    assert all(row["globally_admissible"] for row in traces)
    assert audit["direction_sample_count"] == 181


def test_entrypoint_propagates_direction_contract_and_patches_both_slots(
    monkeypatch,
    tmp_path,
):
    observed = {}

    def fake_base(argv):
        observed["target"] = os.environ.get(
            "ARRHENIUS_COMMITTED_TARGET_EXTENSION_UM"
        )
        observed["da"] = os.environ.get("ARRHENIUS_PHYSICAL_DA_UM")
        observed["lpz"] = os.environ.get(
            "ARRHENIUS_CORRIDOR_PROCESS_ZONE_UM"
        )
        observed["theta"] = os.environ.get("ARRHENIUS_CRYSTAL_THETA_DEG")
        observed["forward"] = os.environ.get(
            "ARRHENIUS_MIN_GLOBAL_FORWARD"
        )
        observed["v91853"] = v91853._quality_selected_corridor_mesh
        observed["v1005132"] = (
            v1005132._quality_selected_corridor_mesh_v1005132
        )
        return "ok"

    monkeypatch.setattr(entry._base, "main", fake_base)
    original_v91853 = v91853._quality_selected_corridor_mesh
    original_v1005132 = v1005132._quality_selected_corridor_mesh_v1005132
    result = entry.main(
        [
            "--target-crack-extension-um",
            "1000",
            "--da-phys",
            "5e-6",
            "--mpz-length-um",
            "50",
            "--crystal-theta-deg",
            "30",
            "--min-global-forward",
            "0.05",
            "--out",
            str(tmp_path),
        ]
    )

    assert result == "ok"
    assert observed["target"] == "1000.0"
    assert float(observed["da"]) == pytest.approx(5.0)
    assert observed["lpz"] == "50.0"
    assert observed["theta"] == "30.0"
    assert observed["forward"] == "0.05"
    assert (
        observed["v91853"]
        is envelope.path_aware_long_growth_envelope_mesh
    )
    assert (
        observed["v1005132"]
        is envelope.path_aware_long_growth_envelope_mesh
    )
    assert v91853._quality_selected_corridor_mesh is original_v91853
    assert (
        v1005132._quality_selected_corridor_mesh_v1005132
        is original_v1005132
    )


def test_committed_endpoint_audit_certifies_actual_geometry_events(tmp_path):
    corridor = {
        "reachable_envelope_vertices_m": [
            [0.5e-3, 0.0],
            [0.55e-3, -1.0e-3],
            [1.51e-3, 0.0],
            [0.55e-3, 1.0e-3],
        ],
        "physical_da_um": 5.0,
        "minimum_global_forward_cosine": 0.05,
        "maximum_sampled_hbar_tip_over_L_pz": 0.245,
        "maximum_hbar_tip_over_L_pz_required": 0.25,
    }
    events = [
        {
            "inserted": True,
            "moved_m": 5.0e-6,
            "p0_m": [0.5e-3, 0.0],
            "actual_p1_m": [0.504e-3, 3.0e-6],
        },
        {
            "inserted": True,
            "moved_m": 5.0e-6,
            "p0_m": [0.504e-3, 3.0e-6],
            "actual_p1_m": [0.508e-3, 0.0],
        },
    ]
    (tmp_path / "stochastic_geometry_events_v10_0_5_16.json").write_text(
        json.dumps(events)
    )
    audit = entry._committed_endpoint_audit(tmp_path, corridor)
    assert audit["committed_endpoint_audit_available"] is True
    assert audit["committed_endpoint_count"] == 2
    assert audit["all_committed_endpoints_inside_certified_envelope"] is True
    assert audit["all_committed_segments_pass_global_forward_gate"] is True
    assert audit["committed_endpoint_support_certified"] is True


def test_runner_and_package_contracts():
    root = Path(__file__).resolve().parents[1]
    runner = (
        root
        / "run_v10_0_5_18_3_4_four_class_focused_long_growth.sh"
    ).read_text()
    assert (
        "mode_i_first_passage_v10_0_5_18_3_4_four_class_"
        "path_aware_growth_envelope"
    ) in runner
    assert "ARRHENIUS_ENVELOPE_SUPPORT_H_OVER_LPZ" in runner
    assert "path_aware_growth_envelope=true" in runner
    assert "child_area_ratio_floor_relaxed=false" in runner

    pyproject = (root / "pyproject.toml").read_text()
    assert 'version = "10.0.5.18.3.4"' in pyproject
