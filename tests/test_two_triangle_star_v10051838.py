from __future__ import annotations

import json
from types import SimpleNamespace

import numpy as np
from scipy.spatial import Delaunay

from arrhenius_fracture.crack_backend import AdaptiveCZMBackend
from arrhenius_fracture.mesh import make_boundary_data, rebuild_tri_mesh
from arrhenius_fracture import quality_aware_czm_patch_v10051835 as base
from arrhenius_fracture import quality_aware_czm_two_triangle_cavity_v10051838 as cavity
from arrhenius_fracture import quality_aware_czm_two_triangle_star_v10051838 as star
from arrhenius_fracture import (
    mode_i_first_passage_v10_0_5_18_3_8_two_triangle_star as entry
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


def test_endpoint_star_replaces_delaunay_sector_omission():
    points = np.asarray(
        [
            [-0.47711663, -1.23376595],
            [0.03239987, -0.67828931],
            [0.64847162, -1.06655914],
            [0.50113899, -0.30504866],
            [0.17297781, -0.60552101],
        ],
        dtype=float,
    )
    unconstrained = np.asarray(Delaunay(points).simplices, dtype=int)
    assert any(4 not in conn for conn in unconstrained)

    constrained = star.EndpointStarTriangulation(points).simplices
    assert constrained.shape == (4, 3)
    assert all(4 in conn for conn in constrained)
    outer = {
        tuple(sorted((int(conn[0]), int(conn[1]))))
        for conn in constrained
    }
    assert len(outer) == 4


def test_constrained_star_cavity_commits_exact_endpoint(monkeypatch):
    monkeypatch.setenv("ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY", "0.035")
    monkeypatch.setenv("ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO", "0.08")
    monkeypatch.setattr(cavity, "Delaunay", star.EndpointStarTriangulation)

    backend = _backend()
    state = _state()
    p0 = np.asarray([0.5, -0.8])
    target = np.asarray([0.5, 0.01])
    direction = (target - p0) / np.linalg.norm(target - p0)

    refined, record = cavity._two_triangle_cavity_candidate(
        backend, state, p0, target, 0, 1
    )
    assert refined is not None, record
    assert record["n_new_cavity_elements"] == 4
    assert record["min_triangle_quality"] >= 0.035
    assert record["min_immediate_child_area_ratio"] >= 0.08

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
    assert np.linalg.norm(backend.tip_nodes[0][2] - target) <= 1.0e-12
    assert backend.advance_log[-1]["node_move_m"] == 0.0


def test_star_entry_annotates_outputs(monkeypatch, tmp_path):
    def fake_base(argv):
        payload = {
            "physics_contract": {},
            "policy": {},
        }
        for name in (
            entry._base.PRODUCTION_MANIFEST,
            entry._base.RESOLUTION_AUDIT,
            entry._base.RETRY_AUDIT,
        ):
            (tmp_path / name).write_text(json.dumps(payload))
        assert cavity.Delaunay is star.EndpointStarTriangulation
        return "ok"

    monkeypatch.setattr(entry._base, "main", fake_base)
    assert entry.main(["--out", str(tmp_path)]) == "ok"

    manifest = json.loads(
        (tmp_path / entry._base.PRODUCTION_MANIFEST).read_text()
    )
    assert manifest["two_triangle_constrained_endpoint_star_active"] is True
    assert manifest["exact_endpoint_connected_to_every_cavity_boundary_edge"] is True
    assert manifest["triangle_quality_floor_relaxed"] is False
    assert manifest["constitutive_physics_changed"] is False
