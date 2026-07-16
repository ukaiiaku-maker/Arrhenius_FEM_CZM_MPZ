from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from arrhenius_fracture.crack_backend import AdaptiveCZMBackend, CrackAdvanceResult
from arrhenius_fracture.mesh import rebuild_tri_mesh
from arrhenius_fracture.mode_i_first_passage_v9_18_3 import (
    _advance_with_identical_veto_guard,
    _edge_aware_insert_target_in_incident_triangle,
)


def _mesh():
    nodes = np.array([
        [0.0, 0.0],
        [1.0, 0.0],
        [0.0, 1.0],
        [1.0, 1.0],
    ])
    elems = np.array([
        [0, 1, 2],
        [1, 3, 2],
    ], dtype=int)
    return rebuild_tri_mesh(nodes, elems, tip_centers=[[0.0, 0.0]])


def _backend():
    geom = SimpleNamespace(Lx=2.0, Ly=2.0, a0=0.0, notch_half_thickness=0.0)
    return AdaptiveCZMBackend(geom=geom)


def test_target_on_shared_edge_uses_exact_edge_split_without_degenerate_child():
    mesh = _mesh()
    backend = _backend()
    u = np.zeros(mesh.ndof)
    refined, u1, q, reason, meta, parent_map = (
        _edge_aware_insert_target_in_incident_triangle(
            backend, mesh, u, np.array([0.0, 0.0]),
            np.array([0.5, 0.5]), 0,
        )
    )
    assert reason == "ok"
    assert refined is not None
    assert refined.ne == 4
    assert np.all(refined.area_e > 0.0)
    assert np.all(np.isfinite(refined.dNdx_e))
    assert np.allclose(q, [0.5, 0.5])
    assert meta["target_location_case"] == "existing_edge_exact_split"
    assert parent_map.shape == (4,)
    assert u1.shape == (2 * refined.nn,)


def test_target_on_existing_vertex_reuses_vertex_without_refinement():
    mesh = _mesh()
    backend = _backend()
    u = np.zeros(mesh.ndof)
    refined, u1, q, reason, meta, parent_map = (
        _edge_aware_insert_target_in_incident_triangle(
            backend, mesh, u, np.array([0.0, 0.0]),
            np.array([1.0, 0.0]), 0,
        )
    )
    assert reason == "ok"
    assert refined.ne == mesh.ne
    assert refined.nn == mesh.nn
    assert np.allclose(q, [1.0, 0.0])
    assert meta["target_location_case"] == "existing_vertex"
    assert np.array_equal(parent_map, np.arange(mesh.ne))
    assert u1.shape == u.shape


def test_strict_interior_target_keeps_three_child_refinement():
    mesh = _mesh()
    backend = _backend()
    u = np.zeros(mesh.ndof)
    refined, _, q, reason, meta, parent_map = (
        _edge_aware_insert_target_in_incident_triangle(
            backend, mesh, u, np.array([0.0, 0.0]),
            np.array([0.2, 0.2]), 0,
        )
    )
    assert reason == "ok"
    assert refined.ne == 4
    assert np.all(refined.area_e > 0.0)
    assert np.allclose(q, [0.2, 0.2])
    assert meta["target_location_case"] == "strict_triangle_interior"
    assert parent_map.shape == (4,)


def test_identical_veto_guard_fails_fast(monkeypatch):
    mesh = _mesh()
    fake = SimpleNamespace()

    def always_veto(self, *args, **kwargs):
        return CrackAdvanceResult(
            mesh=mesh,
            boundary=None,
            damage=np.zeros(mesh.nn),
            displacement=np.zeros(mesh.ndof),
            moved=0.0,
            inserted=False,
            reason="local_hrefine_error:degenerate elements after topology update: [2126]",
        )

    prior = getattr(_advance_with_identical_veto_guard, "_original", None)
    _advance_with_identical_veto_guard._original = always_veto
    monkeypatch.setenv("ARRHENIUS_MAX_IDENTICAL_GEOMETRY_VETOES", "3")
    kwargs = {
        "mesh": mesh,
        "p0": np.array([0.0, 0.0]),
        "p1": np.array([0.5, 0.0]),
        "front_id": 0,
    }
    try:
        _advance_with_identical_veto_guard(fake, **kwargs)
        _advance_with_identical_veto_guard(fake, **kwargs)
        with pytest.raises(RuntimeError, match="repeated identical geometry veto"):
            _advance_with_identical_veto_guard(fake, **kwargs)
    finally:
        if prior is None:
            delattr(_advance_with_identical_veto_guard, "_original")
        else:
            _advance_with_identical_veto_guard._original = prior
