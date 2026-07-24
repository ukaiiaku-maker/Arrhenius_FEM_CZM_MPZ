from __future__ import annotations

import numpy as np

from arrhenius_fracture.adaptive_czm_tip_support_v1005144 import (
    _repair_pair_support,
    _split_segment_topology_supported_v1005144,
    reset_tip_support_audit_v1005144,
    tip_support_audit_v1005144,
)
from arrhenius_fracture.config import GeometryConfig
from arrhenius_fracture.crack_backend import AdaptiveCZMBackend
from arrhenius_fracture.mesh import rebuild_tri_mesh


def incidence(mesh):
    return np.bincount(mesh.elems.ravel(), minlength=mesh.nn)


def test_pair_support_repair_assigns_one_triangle_to_orphan_copy():
    nodes = np.array(
        [
            [0.0, 0.0],
            [0.0, 0.0],
            [1.0, 1.0],
            [-1.0, 1.0],
            [-1.0, -1.0],
            [1.0, -1.0],
        ],
        dtype=float,
    )
    # Node 1 is the intentionally orphaned minus-side copy.
    elems = np.array(
        [
            [0, 2, 3],
            [0, 3, 4],
            [0, 4, 5],
            [0, 5, 2],
        ],
        dtype=int,
    )
    repaired, rows = _repair_pair_support(
        nodes=nodes,
        elems=elems,
        point=np.array([0.0, 0.0]),
        plus_id=0,
        minus_id=1,
        p0=np.array([-1.0, 0.0]),
        p1=np.array([0.0, 0.0]),
        label="leading",
    )
    counts = np.bincount(repaired.ravel(), minlength=len(nodes))
    assert counts[0] > 0
    assert counts[1] > 0
    assert len(rows) == 1
    assert rows[0]["missing_node"] == 1


def test_splitter_reuses_authoritative_trailing_pair_and_supports_new_tip():
    scale = 1.0e-6
    p0 = np.array([0.5, 0.0]) * scale
    p1 = np.array([1.0, 0.0]) * scale
    nodes = np.array(
        [
            p0,
            p0,
            p1,
            np.array([0.75, 0.5]) * scale,
            np.array([0.75, -0.5]) * scale,
        ]
    )
    elems = np.array(
        [
            [0, 2, 3],
            [1, 4, 2],
        ],
        dtype=int,
    )
    mesh = rebuild_tri_mesh(nodes, elems, tip_centers=[p0])
    geom = GeometryConfig(Lx=2.0e-6, Ly=2.0e-6, a0=0.5e-6)
    backend = AdaptiveCZMBackend(geom=geom)
    backend.tip_nodes[0] = (0, 1, p0.copy())
    displacement = np.zeros(2 * mesh.nn)

    reset_tip_support_audit_v1005144()
    new_mesh, new_u, cohesive = _split_segment_topology_supported_v1005144(
        backend,
        mesh,
        displacement,
        p0,
        p1,
        0,
    )
    counts = incidence(new_mesh)
    assert np.all(counts > 0)
    assert cohesive.plus_nodes[0] == 0
    assert cohesive.minus_nodes[0] == 1
    assert counts[cohesive.plus_nodes[1]] > 0
    assert counts[cohesive.minus_nodes[1]] > 0
    assert new_u.shape == (2 * new_mesh.nn,)
    assert backend.tip_nodes[0][0] == cohesive.plus_nodes[1]
    assert backend.tip_nodes[0][1] == cohesive.minus_nodes[1]
    audit = tip_support_audit_v1005144()
    assert audit["authoritative_pair_reuses"] == 1


def test_every_committed_cohesive_endpoint_is_bulk_supported():
    scale = 1.0e-6
    p0 = np.array([0.5, 0.0]) * scale
    p1 = np.array([1.0, 0.0]) * scale
    nodes = np.array(
        [
            p0,
            p0,
            p1,
            np.array([0.75, 0.5]) * scale,
            np.array([0.75, -0.5]) * scale,
        ]
    )
    elems = np.array([[0, 2, 3], [1, 4, 2]], dtype=int)
    mesh = rebuild_tri_mesh(nodes, elems, tip_centers=[p0])
    backend = AdaptiveCZMBackend(
        geom=GeometryConfig(Lx=2.0e-6, Ly=2.0e-6, a0=0.5e-6)
    )
    backend.tip_nodes[0] = (0, 1, p0.copy())
    new_mesh, _, _ = _split_segment_topology_supported_v1005144(
        backend, mesh, np.zeros(2 * mesh.nn), p0, p1, 0
    )
    counts = incidence(new_mesh)
    for element in backend.cohesive_network.elements:
        for node_id in element.nodes4:
            assert 0 <= node_id < new_mesh.nn
            assert counts[node_id] > 0
