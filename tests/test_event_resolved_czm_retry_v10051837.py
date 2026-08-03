from __future__ import annotations

from types import SimpleNamespace
from pathlib import Path

import numpy as np
import pytest

from arrhenius_fracture import event_resolved_czm_retry_v10051837 as event
from arrhenius_fracture import quality_aware_czm_patch_v10051835 as base
from arrhenius_fracture import quality_aware_czm_ray_handoff_v10051835 as handoff
from arrhenius_fracture import quality_aware_czm_retry_v10051835 as local


class Mesh:
    def __init__(self, nodes, elems, hbar_tip=1.0):
        self.nodes = np.asarray(nodes, dtype=float)
        self.elems = np.asarray(elems, dtype=int)
        self.nn = len(self.nodes)
        self.ne = len(self.elems)
        tri = self.nodes[self.elems]
        self.area_e = 0.5 * np.abs(
            (tri[:, 1, 0] - tri[:, 0, 0]) * (tri[:, 2, 1] - tri[:, 0, 1])
            - (tri[:, 1, 1] - tri[:, 0, 1]) * (tri[:, 2, 0] - tri[:, 0, 0])
        )
        self.hbar = float(hbar_tip)
        self.hbar_tip = float(hbar_tip)


class RefineBackend:
    def __init__(self):
        self.geom = SimpleNamespace(Lx=1.0, Ly=1.0)
        self.min_triangle_quality = 0.035
        self.min_area_ratio = 0.08

    @staticmethod
    def _tip_geometric_node_ids(mesh, p0, front_id):
        distance = np.linalg.norm(mesh.nodes - np.asarray(p0)[None, :], axis=1)
        dmin = float(np.min(distance))
        return np.where(distance <= dmin + 1.0e-12)[0].astype(int).tolist()

    @staticmethod
    def _incident_elements(elems, node_ids):
        return np.where(np.any(np.isin(np.asarray(elems), list(node_ids)), axis=1))[0]

    @staticmethod
    def _triangle_quality(nodes, elems):
        tri = np.asarray(nodes)[np.asarray(elems)]
        l2 = (
            np.sum((tri[:, 1] - tri[:, 0]) ** 2, axis=1)
            + np.sum((tri[:, 2] - tri[:, 1]) ** 2, axis=1)
            + np.sum((tri[:, 0] - tri[:, 2]) ** 2, axis=1)
        )
        a2 = np.abs(
            (tri[:, 1, 0] - tri[:, 0, 0]) * (tri[:, 2, 1] - tri[:, 0, 1])
            - (tri[:, 1, 1] - tri[:, 0, 1]) * (tri[:, 2, 0] - tri[:, 0, 0])
        )
        return 2.0 * np.sqrt(3.0) * a2 / np.maximum(l2, 1.0e-300)

    @staticmethod
    def _insert_point_on_edge(mesh, displacement, q, edge_i, edge_j):
        q = np.asarray(q, dtype=float)
        nodes = np.vstack([mesh.nodes, q[None, :]])
        new_id = len(mesh.nodes)
        parents = np.where(
            np.any(mesh.elems == int(edge_i), axis=1)
            & np.any(mesh.elems == int(edge_j), axis=1)
        )[0]
        elems = mesh.elems.copy()
        appended = []
        appended_parent = []
        for parent in parents:
            conn = [int(value) for value in mesh.elems[int(parent)]]
            third = [value for value in conn if value not in (edge_i, edge_j)][0]
            elems[int(parent)] = [third, int(edge_i), new_id]
            appended.append([third, new_id, int(edge_j)])
            appended_parent.append(int(parent))
        elems = np.vstack([elems, np.asarray(appended, dtype=int)])
        refined = Mesh(nodes, elems, hbar_tip=mesh.hbar_tip)
        parent_map = np.concatenate([
            np.arange(mesh.ne, dtype=int),
            np.asarray(appended_parent, dtype=int),
        ])
        u = np.asarray(displacement, dtype=float).reshape(-1, 2)
        uq = 0.5 * (u[int(edge_i)] + u[int(edge_j)])
        u = np.vstack([u, uq[None, :]]).reshape(-1)
        return refined, u, "ok", {"n_new_bulk_elements": len(appended)}, parent_map


def _state(mesh):
    return base.State(
        mesh=mesh,
        boundary=object(),
        damage=np.zeros(mesh.nn),
        displacement=np.zeros(2 * mesh.nn),
        parent_to_root=np.arange(mesh.ne, dtype=int),
    )


def _target_metrics(backend, state, target):
    elem = handoff._target_triangle_any(state.mesh, target)
    assert elem is not None
    return base._predicted_exact_target_metrics(backend, state.mesh, elem, target)


def test_archived_step67_reaches_joint_feasible_endpoint(monkeypatch):
    monkeypatch.setattr(base, "make_boundary_data", lambda mesh, geom: object())
    monkeypatch.setenv("ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY", "0.035")
    monkeypatch.setenv("ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO", "0.08")

    p0 = np.array([0.0005861866667123705, 1.200896179212802e-05])
    target = np.array([0.0005887734388763731, 1.2527425261258408e-05])
    midpoint_1 = np.array([0.0005971874999999998, 1.0369721710082628e-05])
    midpoint_2 = np.array([0.0005987187499999999, 1.3021924509172447e-05])
    vertex_a = 2.0 * midpoint_2 - midpoint_1
    vertex_b = 2.0 * midpoint_1 - vertex_a

    backend = RefineBackend()
    state = _state(Mesh([p0, vertex_a, vertex_b], [[0, 1, 2]], hbar_tip=1.1196217271646168e-05))
    initial = _target_metrics(backend, state, target)
    assert initial["predicted_min_triangle_quality"] == pytest.approx(0.021082302341594004)
    assert initial["predicted_min_child_area_ratio"] == pytest.approx(0.017275666997289947)

    records = []
    final = initial
    for level in range(1, 9):
        state, record = handoff.refine_once(backend, state, p0, target, 0, level)
        assert state is not None, record
        records.append(record)
        assert record["min_triangle_quality"] >= 0.035
        assert record["min_immediate_child_area_ratio"] >= 0.08
        final = _target_metrics(backend, state, target)
        if final["predicted_min_triangle_quality"] >= 0.035 and final["predicted_min_child_area_ratio"] >= 0.08:
            break

    assert len(records) <= 4
    assert final["predicted_min_triangle_quality"] >= 0.035
    assert final["predicted_min_child_area_ratio"] >= 0.08


def test_archived_step71_connectivity_ranked_cavity_reaches_joint_feasible_endpoint(monkeypatch):
    monkeypatch.setattr(base, "make_boundary_data", lambda mesh, geom: object())
    monkeypatch.setenv("ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY", "0.035")
    monkeypatch.setenv("ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO", "0.08")

    p0 = np.array([0.00061963283714399, 1.236159753472365e-06])
    full_target = np.array([0.0006219277841311654, 1.1313273108774453e-06])
    target = 0.5 * (p0 + full_target)
    midpoint_a = np.array([0.000622191418571995, -2.15366766549196e-06])
    midpoint_b = np.array([0.000619128918571995, 3.150737932687678e-06])
    vertex_a = 2.0 * midpoint_a - p0
    vertex_b = 2.0 * midpoint_b - p0

    backend = RefineBackend()
    state = _state(Mesh([p0, vertex_a, vertex_b], [[0, 1, 2]], hbar_tip=1.115236968592306e-05))
    initial = _target_metrics(backend, state, target)
    assert initial["predicted_min_triangle_quality"] == pytest.approx(0.013649178395857288)
    assert initial["predicted_min_child_area_ratio"] == pytest.approx(0.07125645537506213)

    records = []
    final = initial
    for level in range(1, 9):
        state, record = handoff.refine_once(backend, state, p0, target, 0, level)
        assert state is not None, record
        records.append(record)
        final = _target_metrics(backend, state, target)
        if final["predicted_min_triangle_quality"] >= 0.035 and final["predicted_min_child_area_ratio"] >= 0.08:
            break

    assert len(records) <= 4
    assert records[0]["refinement_kind"] == "endpoint_centered_tip_radial_ray_handoff"
    assert any(row["refinement_kind"] == "endpoint_centered_adjacent_target_cavity_midpoint" for row in records[1:])
    assert final["predicted_min_triangle_quality"] >= 0.035
    assert final["predicted_min_child_area_ratio"] >= 0.08


def test_long_event_uses_exact_mesh_crossings_not_equal_partitions(monkeypatch):
    event.reset_audit()
    mesh = SimpleNamespace(ne=1, hbar_tip=1.0, hbar=1.0)
    root = local._State(
        mesh=mesh,
        boundary=object(),
        damage=np.zeros(1),
        displacement=np.zeros(2),
        parent_to_root=np.array([0], dtype=int),
    )

    class Backend:
        def __init__(self):
            self.advance_log = []
            self.tip_nodes = {}
            self.max_hrefine_subsegments = 16
            self.rollbacks = 0
        def _transaction_snapshot(self):
            return {}
        def _transaction_rollback(self, snap):
            self.rollbacks += 1
        def _ray_exit_edge(self, mesh, p0, direction, front_id, remaining):
            x = float(p0[0])
            q = np.array([min(np.floor(x + 1.0), 3.0), 0.0])
            return float(q[0] - x), 0, 1, q, 0, 0.5

    backend = Backend()
    monkeypatch.setattr(event, "_local_endpoint_reachable", lambda self, mesh, p0, p1, front_id: float(p1[0] - p0[0]) <= 1.0 + 1.0e-12)

    calls = []
    def fake_segment(original, self, state, root_kwargs, p0, p1, direction, **kwargs):
        moved = float(np.linalg.norm(np.asarray(p1) - np.asarray(p0)))
        calls.append((kwargs["segment_kind"], moved))
        self.tip_nodes[0] = (0, 1, np.asarray(p1, dtype=float).copy())
        self.advance_log.append({})
        result = SimpleNamespace(
            mesh=state.mesh,
            boundary=state.boundary,
            damage=state.damage,
            displacement=state.displacement,
            elem_parent_map=state.parent_to_root.copy(),
            moved=moved,
            angle_error_deg=0.0,
        )
        return result, 0, None

    monkeypatch.setattr(local, "_segment_retry", fake_segment)
    result, failure, audit = event._exact_ray_transaction(
        lambda *a, **k: None,
        backend,
        {
            "mesh": root.mesh,
            "boundary": root.boundary,
            "damage": root.damage,
            "displacement": root.displacement,
            "front_id": 0,
        },
        np.array([0.0, 0.0]),
        np.array([3.0, 0.0]),
        np.array([1.0, 0.0]),
    )

    assert failure is None
    assert result is not None and result.inserted is True
    assert result.moved == pytest.approx(3.0)
    assert [kind for kind, _ in calls] == [
        "exact_mesh_ray_crossing",
        "exact_mesh_ray_crossing",
        "exact_physical_endpoint",
    ]
    assert audit["equal_length_partition_retry"] is False
    assert audit["segment_count"] == 3
    assert audit["endpoint_error_m"] == pytest.approx(0.0)


def test_runner_and_package_contracts():
    root = Path(__file__).resolve().parents[1]
    runner = (root / "run_v10_0_5_18_3_7_four_class_event_resolved_local_cavity.sh").read_text()
    assert "mode_i_first_passage_v10_0_5_18_3_7_event_resolved_local_cavity" in runner
    assert "ARRHENIUS_EVENT_RESOLVED_PATCH_LEVELS" in runner
    assert "ARRHENIUS_EVENT_RESOLVED_MAX_RAY_SEGMENTS" in runner
    assert "ARRHENIUS_QUALITY_AWARE_PARTITIONS" not in runner
    pyproject = (root / "pyproject.toml").read_text()
    assert 'version = "10.0.5.18.3.7"' in pyproject
