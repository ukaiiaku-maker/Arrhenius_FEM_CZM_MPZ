"""Adaptive-CZM crack-tip endpoint support repair for v10.0.5.14.4.

Repeated cohesive insertion can leave several coincident topological node IDs at
the active crack tip.  The inherited splitter selected the first two IDs rather
than the authoritative plus/minus pair stored by the backend.  On the 17th
5-micrometre advance this can reassign every incident bulk triangle away from a
previously referenced cohesive endpoint, producing an orphan bulk node and an
unsupported cohesive endpoint.

This point release changes topology bookkeeping only:

* reuse the authoritative active-tip plus/minus pair;
* partition the endpoint star by the new crack tangent;
* if roundoff leaves one side unsupported, move one quality-valid incident
  triangle to that side while preserving support on the donor side;
* validate all node incidence and cohesive endpoint support before returning.

The production quality veto remains active and still fails closed if the repair
cannot construct a valid two-sided crack-tip star.
"""
from __future__ import annotations

import copy
from contextlib import contextmanager
from typing import Any, Iterator

import numpy as np

from .cohesive import CohesiveElement
from .crack_backend import AdaptiveCZMBackend
from .mesh import rebuild_tri_mesh

MODEL_ID = "adaptive_czm_authoritative_tip_pair_support_v10_0_5_14_4"

_AUDIT: dict[str, Any] = {
    "split_calls": 0,
    "authoritative_pair_reuses": 0,
    "new_endpoint_duplicates": 0,
    "support_repairs": 0,
    "support_repair_rows": [],
}


def reset_tip_support_audit_v1005144() -> None:
    _AUDIT.clear()
    _AUDIT.update(
        {
            "split_calls": 0,
            "authoritative_pair_reuses": 0,
            "new_endpoint_duplicates": 0,
            "support_repairs": 0,
            "support_repair_rows": [],
        }
    )


def tip_support_audit_v1005144() -> dict[str, Any]:
    return copy.deepcopy(_AUDIT)


def _node_incidence(elems: np.ndarray, n_nodes: int) -> np.ndarray:
    return np.bincount(np.asarray(elems, dtype=int).ravel(), minlength=int(n_nodes))


def _signed_twice_area(nodes: np.ndarray, conn: np.ndarray) -> float:
    tri = np.asarray(nodes, dtype=float)[np.asarray(conn, dtype=int)]
    return float(
        (tri[1, 0] - tri[0, 0]) * (tri[2, 1] - tri[0, 1])
        - (tri[1, 1] - tri[0, 1]) * (tri[2, 0] - tri[0, 0])
    )


def _repair_pair_support(
    *,
    nodes: np.ndarray,
    elems: np.ndarray,
    point: np.ndarray,
    plus_id: int,
    minus_id: int,
    p0: np.ndarray,
    p1: np.ndarray,
    label: str,
) -> tuple[np.ndarray, list[dict[str, Any]]]:
    """Ensure both topological copies have a nondegenerate bulk incident star."""
    connectivity = np.asarray(elems, dtype=int).copy()
    repairs: list[dict[str, Any]] = []
    tangent = np.asarray(p1, dtype=float) - np.asarray(p0, dtype=float)
    if float(np.linalg.norm(tangent)) <= 1.0e-30:
        raise RuntimeError("cannot support a zero-length cohesive endpoint pair")

    for missing_id, donor_id, desired_sign in (
        (int(plus_id), int(minus_id), 1.0),
        (int(minus_id), int(plus_id), -1.0),
    ):
        incidence = _node_incidence(connectivity, len(nodes))
        if incidence[missing_id] > 0:
            continue
        donor_elements = np.where(np.any(connectivity == donor_id, axis=1))[0]
        if len(donor_elements) < 2:
            raise RuntimeError(
                f"{label} endpoint pair cannot support both sides: "
                f"missing={missing_id}, donor={donor_id}, "
                f"donor_incidence={len(donor_elements)}"
            )
        candidates: list[tuple[float, float, int, np.ndarray]] = []
        for element in donor_elements:
            trial = connectivity[int(element)].copy()
            trial[trial == donor_id] = missing_id
            if len(set(int(value) for value in trial)) != 3:
                continue
            area = _signed_twice_area(nodes, trial)
            old_area = _signed_twice_area(nodes, connectivity[int(element)])
            if abs(area) <= 1.0e-24 or old_area * area <= 0.0:
                continue
            centroid = np.mean(np.asarray(nodes)[trial], axis=0)
            rel = centroid - np.asarray(point, dtype=float)
            cross = float(tangent[0] * rel[1] - tangent[1] * rel[0])
            preferred = desired_sign * cross
            # Prefer the geometrically correct side, then a large area margin.
            candidates.append((preferred, abs(area), int(element), trial))
        if not candidates:
            raise RuntimeError(
                f"{label} endpoint pair has no orientation-preserving support repair: "
                f"missing={missing_id}, donor={donor_id}"
            )
        candidates.sort(key=lambda row: (row[0] >= 0.0, row[0], row[1]), reverse=True)
        preferred, area, element, trial = candidates[0]
        connectivity[element] = trial
        repairs.append(
            {
                "endpoint": label,
                "missing_node": missing_id,
                "donor_node": donor_id,
                "reassigned_element": element,
                "desired_cross_sign": desired_sign,
                "selected_signed_cross_score": preferred,
                "selected_twice_area_m2": area,
            }
        )

    incidence = _node_incidence(connectivity, len(nodes))
    if incidence[int(plus_id)] <= 0 or incidence[int(minus_id)] <= 0:
        raise RuntimeError(
            f"{label} endpoint support repair did not produce a two-sided star"
        )
    return connectivity, repairs


def _split_segment_topology_supported_v1005144(
    self: AdaptiveCZMBackend,
    mesh,
    displacement: np.ndarray,
    p0: np.ndarray,
    p1: np.ndarray,
    front_id: int,
):
    """Corrected copy of the topology splitter with authoritative tip pairing."""
    _AUDIT["split_calls"] += 1
    nodes = np.asarray(mesh.nodes, dtype=float).copy()
    elems = np.asarray(mesh.elems, dtype=int).copy()
    u = np.asarray(displacement, dtype=float).reshape(-1, 2).copy()
    p0 = np.asarray(p0, dtype=float)
    p1 = np.asarray(p1, dtype=float)
    tol = max(1.0e-12, 1.0e-6 * max(mesh.hbar_tip, mesh.hbar, 1.0e-12))

    ids0 = self._coincident_ids(nodes, p0, tol)
    ids1 = self._coincident_ids(nodes, p1, tol)
    if not ids0 or not ids1:
        raise RuntimeError("CZM edge endpoints are not represented by mesh nodes")

    def duplicate(base: int, point: np.ndarray) -> int:
        nonlocal nodes, u
        new_id = int(len(nodes))
        nodes = np.vstack([nodes, np.asarray(point, dtype=float)[None, :]])
        u = np.vstack([u, u[int(base)][None, :]])
        _AUDIT["new_endpoint_duplicates"] += 1
        return new_id

    authoritative = self.tip_nodes.get(int(front_id))
    if authoritative is not None:
        a, b, xy = authoritative
        a = int(a)
        b = int(b)
        if (
            a != b
            and 0 <= a < len(nodes)
            and 0 <= b < len(nodes)
            and float(np.linalg.norm(nodes[a] - p0)) <= tol
            and float(np.linalg.norm(nodes[b] - p0)) <= tol
            and float(np.linalg.norm(np.asarray(xy, dtype=float) - p0)) <= tol
        ):
            p0_plus, p0_minus = a, b
            _AUDIT["authoritative_pair_reuses"] += 1
        elif len(ids0) >= 2:
            incidence0 = _node_incidence(elems, len(nodes))
            ordered = sorted(ids0, key=lambda nid: (-int(incidence0[int(nid)]), int(nid)))
            p0_plus, p0_minus = int(ordered[0]), int(ordered[1])
        else:
            p0_plus = int(ids0[0])
            p0_minus = duplicate(p0_plus, p0)
    elif len(ids0) >= 2:
        incidence0 = _node_incidence(elems, len(nodes))
        ordered = sorted(ids0, key=lambda nid: (-int(incidence0[int(nid)]), int(nid)))
        p0_plus, p0_minus = int(ordered[0]), int(ordered[1])
    else:
        p0_plus = int(ids0[0])
        p0_minus = duplicate(p0_plus, p0)

    if len(ids1) >= 2:
        incidence1 = _node_incidence(elems, len(nodes))
        ordered = sorted(ids1, key=lambda nid: (-int(incidence1[int(nid)]), int(nid)))
        p1_plus, p1_minus = int(ordered[0]), int(ordered[1])
    else:
        p1_plus = int(ids1[0])
        p1_minus = duplicate(p1_plus, p1)

    # Partition the full geometric endpoint stars onto the authoritative pairs.
    for point, plus_id, minus_id in (
        (p0, p0_plus, p0_minus),
        (p1, p1_plus, p1_minus),
    ):
        coincident = set(self._coincident_ids(nodes, point, tol)) | {
            int(plus_id),
            int(minus_id),
        }
        incident = self._incident_elements(elems, coincident)
        side = self._classify_element_sides(
            rebuild_tri_mesh(nodes, elems, tip_centers=[p1], validate=False),
            p0,
            p1,
            incident,
        )
        for local_index, element in enumerate(incident):
            target = int(plus_id if side[local_index] >= 0.0 else minus_id)
            for local_node in range(3):
                if int(elems[int(element), local_node]) in coincident:
                    elems[int(element), local_node] = target

    repairs: list[dict[str, Any]] = []
    elems, rows = _repair_pair_support(
        nodes=nodes,
        elems=elems,
        point=p0,
        plus_id=p0_plus,
        minus_id=p0_minus,
        p0=p0,
        p1=p1,
        label="trailing",
    )
    repairs.extend(rows)
    elems, rows = _repair_pair_support(
        nodes=nodes,
        elems=elems,
        point=p1,
        plus_id=p1_plus,
        minus_id=p1_minus,
        p0=p0,
        p1=p1,
        label="leading",
    )
    repairs.extend(rows)

    new_mesh = rebuild_tri_mesh(nodes, elems, tip_centers=[p1])
    incidence = _node_incidence(new_mesh.elems, new_mesh.nn)
    orphan = np.where(incidence <= 0)[0]
    if orphan.size:
        raise RuntimeError(
            "v10.0.5.14.4 endpoint repair left orphan bulk nodes: "
            f"{orphan[:10].tolist()}"
        )

    tangent = p1 - p0
    length = float(np.linalg.norm(tangent))
    tangent /= max(length, 1.0e-300)
    normal = np.array([-tangent[1], tangent[0]], dtype=float)
    element = CohesiveElement(
        plus_nodes=(p0_plus, p1_plus),
        minus_nodes=(p0_minus, p1_minus),
        normal=normal,
        tangent=tangent,
        length=length,
        damage=self.event_damage,
        front_id=int(front_id),
        event_index=int(self.event_counter),
        barrier_kind="exp_floor",
        metadata={
            "tip_support_model": MODEL_ID,
            "tip_support_repairs": copy.deepcopy(repairs),
            "authoritative_trailing_pair": [p0_plus, p0_minus],
            "leading_pair": [p1_plus, p1_minus],
        },
    )
    self.cohesive_network.add(element)
    self.tip_nodes[int(front_id)] = (p1_plus, p1_minus, p1.copy())
    self.event_counter += 1

    bad_endpoints: list[int] = []
    for cohesive in self.cohesive_network.elements:
        for node_id in cohesive.nodes4:
            if (
                int(node_id) < 0
                or int(node_id) >= new_mesh.nn
                or incidence[int(node_id)] <= 0
            ):
                bad_endpoints.append(int(node_id))
    if bad_endpoints:
        raise RuntimeError(
            "v10.0.5.14.4 endpoint repair left unsupported cohesive endpoints: "
            f"{bad_endpoints[:10]}"
        )

    if repairs:
        _AUDIT["support_repairs"] += len(repairs)
        _AUDIT["support_repair_rows"].extend(copy.deepcopy(repairs))
    return new_mesh, u.reshape(-1), element


@contextmanager
def installed_tip_support_repair_v1005144() -> Iterator[None]:
    """Temporarily install the corrected topology splitter."""
    reset_tip_support_audit_v1005144()
    old = AdaptiveCZMBackend._split_segment_topology
    AdaptiveCZMBackend._split_segment_topology = (
        _split_segment_topology_supported_v1005144
    )
    try:
        yield
    finally:
        AdaptiveCZMBackend._split_segment_topology = old


__all__ = [
    "MODEL_ID",
    "installed_tip_support_repair_v1005144",
    "reset_tip_support_audit_v1005144",
    "tip_support_audit_v1005144",
    "_repair_pair_support",
]
