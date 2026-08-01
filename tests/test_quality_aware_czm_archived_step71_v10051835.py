from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from arrhenius_fracture import quality_aware_czm_patch_v10051835 as base
from arrhenius_fracture import quality_aware_czm_ray_handoff_v10051835 as handoff


class Mesh:
    def __init__(self, nodes, elems, hbar_tip=1.0):
        self.nodes = np.asarray(nodes, dtype=float)
        self.elems = np.asarray(elems, dtype=int)
        self.nn = len(self.nodes)
        self.ne = len(self.elems)
        tri = self.nodes[self.elems]
        self.area_e = 0.5 * np.abs(
            (tri[:, 1, 0] - tri[:, 0, 0])
            * (tri[:, 2, 1] - tri[:, 0, 1])
            - (tri[:, 1, 1] - tri[:, 0, 1])
            * (tri[:, 2, 0] - tri[:, 0, 0])
        )
        self.hbar = float(hbar_tip)
        self.hbar_tip = float(hbar_tip)


class Backend:
    def __init__(self):
        self.geom = SimpleNamespace(Lx=1.0, Ly=1.0)
        self.min_triangle_quality = 0.035
        self.min_area_ratio = 0.08

    @staticmethod
    def _tip_geometric_node_ids(mesh, p0, front_id):
        distance = np.linalg.norm(
            mesh.nodes - np.asarray(p0)[None, :], axis=1
        )
        dmin = float(np.min(distance))
        return np.where(distance <= dmin + 1.0e-12)[0].astype(int).tolist()

    @staticmethod
    def _incident_elements(elems, node_ids):
        return np.where(
            np.any(np.isin(np.asarray(elems), list(node_ids)), axis=1)
        )[0]

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
            third = [
                value
                for value in conn
                if value not in (edge_i, edge_j)
            ][0]
            elems[int(parent)] = [third, int(edge_i), new_id]
            appended.append([third, new_id, int(edge_j)])
            appended_parent.append(int(parent))
        elems = np.vstack([elems, np.asarray(appended, dtype=int)])
        refined = Mesh(nodes, elems, hbar_tip=mesh.hbar_tip)
        parent_map = np.concatenate(
            [
                np.arange(mesh.ne, dtype=int),
                np.asarray(appended_parent, dtype=int),
            ]
        )
        u = np.asarray(displacement, dtype=float).reshape(-1, 2)
        uq = 0.5 * (u[int(edge_i)] + u[int(edge_j)])
        u = np.vstack([u, uq[None, :]]).reshape(-1)
        return (
            refined,
            u,
            "ok",
            {"n_new_bulk_elements": len(appended)},
            parent_map,
        )


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
    return base._predicted_exact_target_metrics(
        backend, state.mesh, elem, target
    )


def test_archived_step71_partition_segment_reaches_both_quality_floors(
    monkeypatch,
):
    monkeypatch.setattr(base, "make_boundary_data", lambda mesh, geom: object())
    monkeypatch.setenv("ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY", "0.035")
    monkeypatch.setenv("ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO", "0.08")

    # Reconstructed from the step-71 audit:
    # p0 is event-20 start, target is partition-2 segment 0, and the two
    # midpoint records identify the two original tip-radial vertices.
    p0 = np.array([0.00061963283714399, 1.236159753472365e-06])
    full_target = np.array(
        [0.0006219277841311654, 1.1313273108774453e-06]
    )
    target = 0.5 * (p0 + full_target)
    midpoint_a = np.array(
        [0.000622191418571995, -2.15366766549196e-06]
    )
    midpoint_b = np.array(
        [0.000619128918571995, 3.150737932687678e-06]
    )
    vertex_a = 2.0 * midpoint_a - p0
    vertex_b = 2.0 * midpoint_b - p0

    backend = Backend()
    state = _state(
        Mesh(
            [p0, vertex_a, vertex_b],
            [[0, 1, 2]],
            hbar_tip=1.115236968592306e-05,
        )
    )

    initial = _target_metrics(backend, state, target)
    assert initial["predicted_min_triangle_quality"] == np.testing.assert_allclose(
        initial["predicted_min_triangle_quality"],
        0.013649178395857288,
        rtol=1.0e-10,
        atol=1.0e-14,
    )
    np.testing.assert_allclose(
        initial["predicted_min_child_area_ratio"],
        0.07125645537506213,
        rtol=1.0e-10,
        atol=1.0e-14,
    )

    records = []
    final = initial
    for level in range(1, 9):
        state, record = handoff.refine_once(
            backend,
            state,
            p0,
            target,
            front_id=0,
            level=level,
        )
        assert state is not None, record
        records.append(record)
        assert record["min_triangle_quality"] >= 0.035
        assert record["min_immediate_child_area_ratio"] >= 0.08
        final = _target_metrics(backend, state, target)
        if (
            final["predicted_min_triangle_quality"] >= 0.035
            and final["predicted_min_child_area_ratio"] >= 0.08
        ):
            break

    assert len(records) <= 4
    assert records[0]["refinement_kind"] == (
        "endpoint_centered_tip_radial_ray_handoff"
    )
    assert any(
        row["refinement_kind"]
        == "endpoint_centered_adjacent_target_cavity_midpoint"
        for row in records[1:]
    )
    assert final["predicted_min_triangle_quality"] >= 0.035
    assert final["predicted_min_child_area_ratio"] >= 0.08
