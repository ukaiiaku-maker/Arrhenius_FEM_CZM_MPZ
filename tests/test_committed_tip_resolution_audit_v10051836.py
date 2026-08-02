from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from arrhenius_fracture import committed_tip_resolution_audit_v10051836 as audit
from arrhenius_fracture.numerical_resilience_v1005183 import (
    RetryingAdaptiveCZMBackendV1005183,
)


class Mesh:
    def __init__(self, nodes, elems, areas=None):
        self.nodes = np.asarray(nodes, dtype=float)
        self.elems = np.asarray(elems, dtype=int)
        self.nn = len(self.nodes)
        self.ne = len(self.elems)
        self.hbar = 1.0
        self.hbar_tip = 0.75
        if areas is None:
            tri = self.nodes[self.elems]
            areas = 0.5 * np.abs(
                (tri[:, 1, 0] - tri[:, 0, 0])
                * (tri[:, 2, 1] - tri[:, 0, 1])
                - (tri[:, 1, 1] - tri[:, 0, 1])
                * (tri[:, 2, 0] - tri[:, 0, 0])
            )
        self.area_e = np.asarray(areas, dtype=float)


class SyntheticBackend(audit.AuditedRetryingAdaptiveCZMBackendV10051836):
    def __init__(self):
        # Avoid constructing the production cohesive backend; the parent
        # ``advance`` method is replaced by the tests below.
        self.tip_nodes = {}

    @staticmethod
    def _triangle_quality(nodes, elems):
        tri = np.asarray(nodes)[np.asarray(elems)]
        l2 = (
            np.sum((tri[:, 1] - tri[:, 0]) ** 2, axis=1)
            + np.sum((tri[:, 2] - tri[:, 1]) ** 2, axis=1)
            + np.sum((tri[:, 0] - tri[:, 2]) ** 2, axis=1)
        )
        a2 = np.abs(
            (tri[:, 1, 0] - tri[:, 0, 0])
            * (tri[:, 2, 1] - tri[:, 0, 1])
            - (tri[:, 1, 1] - tri[:, 0, 1])
            * (tri[:, 2, 0] - tri[:, 0, 0])
        )
        return 2.0 * np.sqrt(3.0) * a2 / np.maximum(l2, 1.0e-300)

    @staticmethod
    def _incident_elements(elems, node_ids):
        return np.where(
            np.any(np.isin(np.asarray(elems), list(node_ids)), axis=1)
        )[0]

    def _tip_geometric_node_ids(self, mesh, point, front_id):
        distance = np.linalg.norm(
            np.asarray(mesh.nodes) - np.asarray(point)[None, :], axis=1
        )
        return np.where(distance <= np.min(distance) + 1.0e-12)[0].tolist()


def _configure():
    audit.configure(
        tip_h_fine_m=2.5e-6,
        tip_ratio=1.15,
        triangle_quality_floor=0.035,
        child_area_ratio_floor=0.08,
        event_minimum_factor=0.5,
        event_maximum_factor=4.0,
    )
    audit.reset_audit()


def test_event_length_atoms_are_reported_as_threshold_clipping():
    _configure()
    lower = audit.event_length_diagnostics(0.1, 5.0e-6)
    upper = audit.event_length_diagnostics(8.0, 5.0e-6)
    interior = audit.event_length_diagnostics(1.25, 5.0e-6)

    assert lower["event_length_clip_state"] == "lower"
    assert lower["realized_event_length_factor"] == pytest.approx(
        0.45946801912497465
    )
    assert lower["realized_event_length_m"] == pytest.approx(
        2.2973400956248734e-6
    )
    assert upper["event_length_clip_state"] == "upper"
    assert upper["realized_event_length_m"] == pytest.approx(
        18.378720764998987e-6
    )
    assert interior["event_length_clip_state"] == "none"


def test_tip_h_fine_contract_does_not_equate_input_spacing_with_hbar_tip():
    _configure()
    contract = audit.tip_h_fine_contract()
    assert contract["tip_h_fine_m"] == pytest.approx(2.5e-6)
    assert contract["tip_h_fine_equals_one_ring_mean_edge"] is False
    assert "first radial-ring increment" in contract["tip_h_fine_semantics"]
    assert "nearest max(4,2_percent)" in contract["initial_hbar_tip_semantics"]


def test_observer_records_pre_event_state_and_returns_base_result_unchanged(
    monkeypatch,
):
    _configure()
    mesh = Mesh(
        [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
        [[0, 1, 2]],
    )
    sentinel = SimpleNamespace(inserted=False, moved=0.0, reason="sentinel")
    calls = []

    def fake_parent(self, **kwargs):
        calls.append(kwargs)
        return sentinel

    monkeypatch.setattr(
        RetryingAdaptiveCZMBackendV1005183,
        "advance",
        fake_parent,
    )
    backend = SyntheticBackend()
    kwargs = {
        "mesh": mesh,
        "boundary": object(),
        "damage": np.zeros(mesh.nn),
        "displacement": np.zeros(2 * mesh.nn),
        "p0": np.array([0.0, 0.0]),
        "p1": np.array([0.2, 0.2]),
        "direction": np.array([1.0, 1.0]) / np.sqrt(2.0),
        "front_id": 0,
    }
    result = backend.advance(**kwargs)

    assert result is sentinel
    assert len(calls) == 1
    assert calls[0] is not kwargs
    payload = audit.audit_payload()
    assert payload["event_count"] == 1
    row = payload["events"][0]
    assert row["committed_tip_point_m"] == [0.0, 0.0]
    assert row["requested_endpoint_m"] == [0.2, 0.2]
    assert row["committed_tip_one_ring"]["available"] is True
    assert row["target_parent"]["available"] is True
    assert row["backend_result"]["inserted"] is False
    np.testing.assert_array_equal(mesh.elems, np.array([[0, 1, 2]]))


def test_joint_margin_classifies_area_only_failure(monkeypatch):
    _configure()
    old = Mesh(
        [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
        [[0, 1, 2]],
        areas=[1.0],
    )
    new = Mesh(
        [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
        [[0, 1, 2]],
        areas=[0.0267],
    )
    result = SimpleNamespace(
        inserted=True,
        moved=0.1,
        reason="ok",
        mesh=new,
        elem_parent_map=np.array([0], dtype=int),
    )
    backend = SyntheticBackend()
    monkeypatch.setattr(
        audit._v9185,
        "_affected_elements",
        lambda old_mesh, new_mesh: np.array([0], dtype=int),
    )
    monkeypatch.setattr(
        backend,
        "_triangle_quality",
        lambda nodes, elems: np.array([0.0752]),
    )

    metrics = audit._candidate_quality(
        backend,
        old,
        result,
        np.array([0.0, 0.0]),
        0,
    )

    assert metrics["joint_margin_classification"] == "area_only_failure"
    assert metrics["quality_margin"] > 1.0
    assert metrics["area_margin"] < 1.0
