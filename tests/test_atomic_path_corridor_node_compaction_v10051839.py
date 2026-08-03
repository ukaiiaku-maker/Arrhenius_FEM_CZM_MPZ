from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from arrhenius_fracture.atomic_path_corridor_adaptive_quality_v10051839 import (
    AdaptiveQualityAtomicPathCorridorCZMBackendV10051839,
)
from arrhenius_fracture.atomic_path_corridor_local_scale_v10051839 import (
    apply_backend_node_remap,
    compact_candidate_nodes,
)
from arrhenius_fracture.cohesive import CohesiveElement


def _backend():
    geom = SimpleNamespace(
        Lx=2.0,
        Ly=2.0,
        a0=0.25,
        notch_half_thickness=1.0e-3,
    )
    return AdaptiveQualityAtomicPathCorridorCZMBackendV10051839(
        geom=geom,
        min_area_ratio=0.08,
        min_triangle_quality=0.035,
    )


def test_unprotected_orphan_is_compacted_without_changing_elements():
    backend = _backend()
    nodes = np.asarray(
        [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [0.25, 0.25]],
        dtype=float,
    )
    elems = np.asarray([[0, 1, 2]], dtype=int)
    u = np.arange(8, dtype=float)
    d = np.asarray([0.0, 0.1, 0.2, 0.3])

    state, record = compact_candidate_nodes(
        backend, nodes, elems, u, d, 1.0e-12
    )
    assert state is not None, record
    assert state["nodes"].shape == (3, 2)
    assert np.array_equal(state["elems"], elems)
    assert state["old_to_new"].tolist() == [0, 1, 2, -1]
    assert record["retired_orphan_node_ids_first20"] == [3]
    assert record["protected_orphan_replacement_count"] == 0


def test_protected_orphan_maps_to_supported_same_side_copy():
    backend = _backend()
    nodes = np.asarray(
        [
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
            [0.0, 0.0],  # protected orphan, same geometry/state as node 0
        ],
        dtype=float,
    )
    elems = np.asarray([[0, 1, 2]], dtype=int)
    u = np.asarray(
        [
            [0.25, -0.10],
            [0.0, 0.0],
            [0.0, 0.0],
            [0.25, -0.10],
        ],
        dtype=float,
    )
    d = np.asarray([0.4, 0.0, 0.0, 0.4])
    backend.tip_nodes[0] = (3, 3, nodes[3].copy())
    backend.cohesive_network.add(
        CohesiveElement(
            plus_nodes=(3, 1),
            minus_nodes=(3, 2),
            normal=np.asarray([0.0, 1.0]),
            tangent=np.asarray([1.0, 0.0]),
            length=1.0,
            front_id=0,
            event_index=0,
        )
    )

    state, record = compact_candidate_nodes(
        backend, nodes, elems, u.reshape(-1), d, 1.0e-12
    )
    assert state is not None, record
    assert state["old_to_new"].tolist() == [0, 1, 2, 0]
    assert record["protected_orphan_replacements"] == {"3": 0}

    apply_backend_node_remap(backend, state["old_to_new"])
    assert backend.tip_nodes[0][0:2] == (0, 0)
    cohesive = backend.cohesive_network.elements[0]
    assert cohesive.plus_nodes == (0, 1)
    assert cohesive.minus_nodes == (0, 2)
    assert cohesive.metadata["v10051839_node_compaction_remapped"] is True


def test_unresolved_protected_orphan_is_rejected_without_mutation():
    backend = _backend()
    nodes = np.asarray(
        [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0], [0.4, 0.4]],
        dtype=float,
    )
    elems = np.asarray([[0, 1, 2]], dtype=int)
    backend.tip_nodes[0] = (3, 3, nodes[3].copy())
    before = backend._transaction_snapshot()

    state, record = compact_candidate_nodes(
        backend,
        nodes,
        elems,
        np.zeros(8),
        np.zeros(4),
        1.0e-12,
    )
    assert state is None
    assert record["reason"] == "protected_orphan_has_no_supported_coincident_copy"
    after = backend._transaction_snapshot()
    assert before["tip_nodes"] == after["tip_nodes"]
    assert before["n_cohesive"] == after["n_cohesive"]
