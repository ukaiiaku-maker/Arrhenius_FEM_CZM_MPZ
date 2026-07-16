from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from arrhenius_fracture.crack_backend import AdaptiveCZMBackend, CrackAdvanceResult
from arrhenius_fracture.mesh import rebuild_tri_mesh, make_boundary_data
from arrhenius_fracture.mode_i_first_passage_v9_18_4 import (
    _mechanically_regularized_insert_target,
    _mechanical_topology_issues,
    _finite_solve_dirichlet,
)


def _geom():
    return SimpleNamespace(Lx=2.0, Ly=2.0, a0=0.0, notch_half_thickness=0.0)


def _mesh():
    nodes = np.array([
        [0.0, 0.0],
        [1.0, 0.0],
        [0.0, 1.0],
        [1.0, 1.0],
        [0.0, -1.0],
        [1.0, -1.0],
    ])
    elems = np.array([
        [0, 1, 2],
        [1, 3, 2],
        [4, 5, 1],
        [4, 1, 0],
    ], dtype=int)
    return rebuild_tri_mesh(nodes, elems, tip_centers=[[0.0, 0.0]])


def test_edge_target_is_regularized_at_same_length(monkeypatch):
    mesh = _mesh()
    backend = AdaptiveCZMBackend(geom=_geom())
    u = np.zeros(mesh.ndof)
    monkeypatch.setenv("ARRHENIUS_EDGE_FRONT_REGULARIZATION_ANGLES_DEG", "0.05 0.1 0.2")
    p0 = np.array([0.0, 0.0])
    target = np.array([0.5, 0.5])
    refined, u1, q, reason, meta, parent_map = _mechanically_regularized_insert_target(
        backend, mesh, u, p0, target, 0
    )
    assert reason == "ok"
    assert refined is not None
    assert meta["target_location_case"] == "same_length_interior_regularization"
    assert abs(np.linalg.norm(q - p0) - np.linalg.norm(target - p0)) < 1e-12
    assert 0.0 < abs(meta["regularization_angle_deg"]) <= 0.2
    assert np.all(refined.area_e > 0.0)
    assert np.all(np.isfinite(refined.dNdx_e))
    assert u1.shape == (refined.ndof,)
    assert parent_map.shape == (refined.ne,)


def test_orphan_node_is_rejected():
    mesh0 = _mesh()
    nodes = np.vstack([mesh0.nodes, [[0.25, 0.25]]])
    mesh = rebuild_tri_mesh(nodes, mesh0.elems, tip_centers=[[0.0, 0.0]])
    bnd = make_boundary_data(mesh, _geom())
    result = CrackAdvanceResult(
        mesh=mesh,
        boundary=bnd,
        damage=np.zeros(mesh.nn),
        displacement=np.zeros(mesh.ndof),
        moved=1.0,
        inserted=True,
    )
    backend = AdaptiveCZMBackend(geom=_geom())
    issues = _mechanical_topology_issues(backend, result)
    assert any(item.startswith("orphan_bulk_nodes:") for item in issues)


def test_unanchored_disconnected_component_is_rejected():
    nodes = np.array([
        [0.0, -1.0], [1.0, -1.0], [0.0, 0.0],
        [0.0, 0.5], [1.0, 0.5], [0.0, 1.0],
    ])
    elems = np.array([[0, 1, 2], [3, 4, 5]], dtype=int)
    mesh = rebuild_tri_mesh(nodes, elems, tip_centers=[[0.0, 0.0]])
    bnd = make_boundary_data(mesh, _geom())
    result = CrackAdvanceResult(
        mesh=mesh,
        boundary=bnd,
        damage=np.zeros(mesh.nn),
        displacement=np.zeros(mesh.ndof),
        moved=1.0,
        inserted=True,
    )
    backend = AdaptiveCZMBackend(geom=_geom())
    issues = _mechanical_topology_issues(backend, result)
    assert any(item.startswith("unanchored_bulk_components:") for item in issues)


def test_nonfinite_mechanics_fails_before_mpz(monkeypatch):
    prior = getattr(_finite_solve_dirichlet, "_original", None)
    _finite_solve_dirichlet._original = lambda *a, **k: (np.array([0.0, np.nan]), np.nan)
    try:
        with pytest.raises(RuntimeError, match="non-finite FEM solution"):
            _finite_solve_dirichlet(None)
    finally:
        if prior is None:
            delattr(_finite_solve_dirichlet, "_original")
        else:
            _finite_solve_dirichlet._original = prior
