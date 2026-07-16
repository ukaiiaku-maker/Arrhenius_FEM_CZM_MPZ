from __future__ import annotations

from types import SimpleNamespace

import numpy as np
from scipy import sparse

from arrhenius_fracture.crack_backend import AdaptiveCZMBackend, CrackAdvanceResult
from arrhenius_fracture.mesh import rebuild_tri_mesh
from arrhenius_fracture.mode_i_first_passage_v9_18_5 import (
    _DynamicStepHorizon,
    _RUNTIME,
    _component_anchored_solve,
    _corridor_centers,
    _strict_quality_advance,
)


def test_dynamic_horizon_stops_only_after_controller_request():
    controller = SimpleNamespace(v9185_stop_requested=False)
    horizon = _DynamicStepHorizon(100, controller)
    assert 5 < horizon
    assert int(horizon) == 100
    assert horizon * 2 == 200
    controller.v9185_stop_requested = True
    assert not (5 < horizon)
    assert int(horizon) == 100


def test_corridor_centers_cover_target_and_guard(monkeypatch):
    monkeypatch.setenv("ARRHENIUS_COMMITTED_TARGET_EXTENSION_UM", "100")
    monkeypatch.setenv("ARRHENIUS_CORRIDOR_GUARD_UM", "10")
    monkeypatch.setenv("ARRHENIUS_CORRIDOR_CENTER_SPACING_UM", "25")
    geom = SimpleNamespace(a0=5.0e-4, Lx=2.0e-3)
    mesh_cfg = SimpleNamespace(tip_h_fine=1.0e-6)
    centers = _corridor_centers(geom, mesh_cfg)
    assert np.isclose(centers[0, 0], geom.a0)
    assert centers[-1, 0] >= geom.a0 + 110.0e-6 - 1.0e-15
    assert np.allclose(centers[:, 1], 0.0)
    assert np.max(np.diff(centers[:, 0])) <= 25.0e-6 + 1.0e-15


def _disconnected_mesh():
    nodes = np.array([
        [0.0, -1.0], [1.0, -1.0], [0.0, -0.2],
        [0.0, 0.2], [1.0, 1.0], [0.0, 1.0],
    ])
    elems = np.array([[0, 1, 2], [3, 4, 5]], dtype=int)
    return rebuild_tri_mesh(nodes, elems, tip_centers=[[0.0, 0.0]])


def test_component_solve_adds_minimal_x_anchor_to_upper_body():
    mesh = _disconnected_mesh()
    _RUNTIME["mesh"] = mesh
    _RUNTIME["component_anchor_history"] = []
    bnd = SimpleNamespace(
        top_nodes=np.array([4, 5], dtype=int),
        bot_nodes=np.array([0, 1], dtype=int),
        left_bot=0,
        right_bot=1,
    )
    K = sparse.eye(mesh.ndof, format="csr")
    u = np.zeros(mesh.ndof)
    R = np.zeros(mesh.ndof)
    u_new, force = _component_anchored_solve(K, R, u, bnd, 1.0e-6, -1.0e-6)
    assert np.all(np.isfinite(u_new))
    assert np.isfinite(force)
    assert _RUNTIME["component_anchor_history"]
    anchors = _RUNTIME["component_anchor_history"][-1]
    assert len(anchors) == 1
    assert anchors[0]["node"] in {4, 5}


class _FakeNetwork:
    elements = []


class _FakeBackend:
    min_triangle_quality = 0.035
    min_area_ratio = 0.08
    cohesive_network = _FakeNetwork()
    advance_log = []

    def __init__(self):
        self.rolled_back = False

    def _transaction_snapshot(self):
        return {"snapshot": True}

    def _transaction_rollback(self, snapshot):
        assert snapshot["snapshot"]
        self.rolled_back = True

    @staticmethod
    def _triangle_quality(nodes, elems):
        return AdaptiveCZMBackend._triangle_quality(nodes, elems)


def test_strict_quality_gate_rolls_back_sliver(monkeypatch):
    old = rebuild_tri_mesh(
        np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]]),
        np.array([[0, 1, 2]], dtype=int),
        tip_centers=[[0.0, 0.0]],
    )
    sliver = rebuild_tri_mesh(
        np.array([[0.0, 0.0], [1.0, 0.0], [1.0e-5, 1.0e-5]]),
        np.array([[0, 1, 2]], dtype=int),
        tip_centers=[[0.0, 0.0]],
    )
    damage = np.zeros(old.nn)
    displacement = np.zeros(old.ndof)

    def fake_original(self, *args, **kwargs):
        return CrackAdvanceResult(
            mesh=sliver,
            boundary=kwargs["boundary"],
            damage=damage,
            displacement=displacement,
            moved=1.0,
            inserted=True,
            reason="ok",
            elem_parent_map=np.array([0], dtype=int),
        )

    backend = _FakeBackend()
    prior = getattr(_strict_quality_advance, "_original", None)
    _strict_quality_advance._original = fake_original
    monkeypatch.setenv("ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY", "0.035")
    monkeypatch.setenv("ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO", "0.08")
    monkeypatch.setenv("ARRHENIUS_MAX_TIP_H_OVER_DA", "100")
    try:
        result = _strict_quality_advance(
            backend,
            mesh=old,
            boundary=object(),
            damage=damage,
            displacement=displacement,
            p0=np.array([0.0, 0.0]),
            p1=np.array([1.0, 0.0]),
            direction=np.array([1.0, 0.0]),
            front_id=0,
        )
    finally:
        if prior is None:
            delattr(_strict_quality_advance, "_original")
        else:
            _strict_quality_advance._original = prior
    assert not result.inserted
    assert "v9185_quality_veto" in result.reason
    assert backend.rolled_back
