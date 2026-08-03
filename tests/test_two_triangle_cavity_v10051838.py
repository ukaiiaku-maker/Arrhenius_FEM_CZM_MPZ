from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from arrhenius_fracture.crack_backend import AdaptiveCZMBackend
from arrhenius_fracture.mesh import make_boundary_data, rebuild_tri_mesh
from arrhenius_fracture import quality_aware_czm_patch_v10051835 as base
from arrhenius_fracture import quality_aware_czm_two_triangle_cavity_v10051838 as cavity
from arrhenius_fracture import (
    mode_i_first_passage_v10_0_5_18_3_8_two_triangle_cavity as entry
)


def _geom():
    return SimpleNamespace(
        Lx=2.0,
        Ly=3.0,
        a0=0.5,
        notch_half_thickness=1.0e-3,
    )


def _backend():
    return AdaptiveCZMBackend(
        geom=_geom(),
        penalty_normal_Pa_per_m=1.0e18,
        penalty_tangent_Pa_per_m=1.0e18,
        max_angle_error_deg=35.0,
        event_damage=1.0,
        min_area_ratio=0.08,
        min_triangle_quality=0.035,
        max_node_move_factor=1.75,
        max_hrefine_subsegments=64,
    )


def _state():
    nodes = np.asarray(
        [
            [0.5, -0.8],
            [0.0, 0.0],
            [1.0, 0.0],
            [0.5, np.sqrt(3.0) / 2.0],
        ],
        dtype=float,
    )
    elems = np.asarray([[0, 2, 1], [1, 2, 3]], dtype=int)
    mesh = rebuild_tri_mesh(nodes, elems, tip_centers=[[0.5, -0.8]])
    return base.State(
        mesh=mesh,
        boundary=make_boundary_data(mesh, _geom()),
        damage=np.zeros(mesh.nn),
        displacement=np.zeros(mesh.ndof),
        parent_to_root=np.arange(mesh.ne, dtype=int),
    )


def _boundary_edges(elems):
    counts = {}
    for conn in np.asarray(elems, dtype=int):
        for i, j in ((conn[0], conn[1]), (conn[1], conn[2]), (conn[2], conn[0])):
            edge = tuple(sorted((int(i), int(j))))
            counts[edge] = counts.get(edge, 0) + 1
    return {edge for edge, count in counts.items() if count == 1}


def test_two_triangle_cavity_removes_near_edge_and_preserves_outer_boundary(monkeypatch):
    monkeypatch.setenv("ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY", "0.035")
    monkeypatch.setenv("ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO", "0.08")
    backend = _backend()
    state = _state()
    p0 = np.asarray([0.5, -0.8])
    target = np.asarray([0.5, 0.01])

    old_boundary = _boundary_edges(state.mesh.elems)
    refined, record = cavity._two_triangle_cavity_candidate(
        backend, state, p0, target, 0, 1
    )
    assert refined is not None, record
    assert record["refinement_kind"] == "exact_endpoint_convex_two_triangle_delaunay"
    assert record["min_triangle_quality"] >= 0.035
    assert record["min_immediate_child_area_ratio"] >= 0.08
    assert record["n_removed_elements"] == 2
    assert record["n_new_cavity_elements"] == 4
    assert record["n_net_new_elements"] == 2
    assert record["exact_endpoint_connected_tip_ids"]
    assert _boundary_edges(refined.mesh.elems) == old_boundary

    target_id = record["exact_endpoint_node_id"]
    assert np.linalg.norm(refined.mesh.nodes[target_id] - target) <= 1.0e-14
    assert len(refined.damage) == refined.mesh.nn
    assert len(refined.displacement) == refined.mesh.ndof
    assert len(refined.parent_to_root) == refined.mesh.ne
    assert set(refined.parent_to_root.tolist()) <= {0, 1}


def test_preinserted_exact_endpoint_commits_with_unmodified_backend(monkeypatch):
    monkeypatch.setenv("ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY", "0.035")
    monkeypatch.setenv("ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO", "0.08")
    backend = _backend()
    state = _state()
    p0 = np.asarray([0.5, -0.8])
    target = np.asarray([0.5, 0.01])
    direction = (target - p0) / np.linalg.norm(target - p0)

    refined, record = cavity._two_triangle_cavity_candidate(
        backend, state, p0, target, 0, 1
    )
    assert refined is not None, record

    result = backend.advance(
        mesh=refined.mesh,
        boundary=refined.boundary,
        damage=refined.damage,
        displacement=refined.displacement,
        p0=p0,
        p1=target,
        direction=direction,
        front_id=0,
    )
    assert result.inserted is True, result.reason
    assert result.moved == pytest.approx(float(np.linalg.norm(target - p0)))
    assert np.linalg.norm(backend.tip_nodes[0][2] - target) <= 1.0e-12
    assert backend.advance_log[-1]["node_move_m"] == pytest.approx(0.0)


def test_fallback_activates_only_after_one_edge_cavity_exhaustion(monkeypatch):
    backend = _backend()
    state = _state()
    p0 = np.asarray([0.5, -0.8])
    target = np.asarray([0.5, 0.01])

    def exhausted(*args, **kwargs):
        return None, {
            "level": 1,
            "accepted": False,
            "reason": "no_target_tip_triangle",
            "target_cavity_failure": {
                "reason": "no_quality_safe_adjacent_target_cavity_bisection"
            },
        }

    monkeypatch.setattr(cavity, "_ORIGINAL_REFINE_ONCE", exhausted)
    refined, record = cavity.refine_once(backend, state, p0, target, 0, 1)
    assert refined is not None, record
    assert record["refinement_kind"] == "exact_endpoint_convex_two_triangle_delaunay"


def test_v10051838_output_and_runner_contract(monkeypatch, tmp_path):
    observed = {}

    def fake_edge(argv):
        observed["install"] = entry._local37.install_ray_handoff
        for old_name in entry._RENAMES:
            payload = {"schema": "v10.0.5.18.3.7_test"}
            if old_name == entry._local37.PRODUCTION_MANIFEST:
                payload["physics_contract"] = {}
            (tmp_path / old_name).write_text(json.dumps(payload))
        return "ok"

    monkeypatch.setattr(entry._edge37, "main", fake_edge)
    assert entry.main(["--out", str(tmp_path)]) == "ok"
    assert observed["install"] is entry.install_two_triangle_cavity

    manifest = json.loads((tmp_path / entry.PRODUCTION_MANIFEST).read_text())
    physics = manifest["physics_contract"]
    assert manifest["point_release"] == "10.0.5.18.3.8"
    assert physics["two_triangle_exact_endpoint_cavity_active"] is True
    assert physics["old_near_diagonal_removed"] is True
    assert physics["equal_length_partition_retry"] is False
    assert physics["triangle_quality_floor_relaxed"] is False

    root = Path(__file__).resolve().parents[1]
    runner = (root / "run_v10_0_5_18_3_8_four_class_two_triangle_cavity.sh").read_text()
    assert "mode_i_first_passage_v10_0_5_18_3_8_two_triangle_cavity" in runner
    assert "ARRHENIUS_QUALITY_AWARE_PARTITIONS" not in runner
    assert 'version = "10.0.5.18.3.8"' in (root / "pyproject.toml").read_text()
