from __future__ import annotations

import json
from collections import Counter
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.spatial import Delaunay

from arrhenius_fracture.atomic_path_corridor_czm_v10051839 import (
    audit_payload,
    reset_audit,
)
from arrhenius_fracture.atomic_path_corridor_local_scale_v10051839 import (
    CertifiedAtomicPathCorridorCZMBackendV10051839,
    install,
    restore,
)
from arrhenius_fracture.mesh import make_boundary_data, rebuild_tri_mesh


@pytest.fixture(autouse=True)
def _certified_local_scale_stack():
    saved = install()
    try:
        yield
    finally:
        restore(saved)


def _geom():
    return SimpleNamespace(
        Lx=20.0,
        Ly=16.0,
        a0=2.0,
        notch_half_thickness=1.0e-3,
    )


def _mesh():
    # A healthy inherited unstructured mesh with the physical tip inserted as
    # an exact vertex but no pre-existing edge to the requested endpoint.
    x = np.linspace(0.0, 20.0, 9)
    y = np.linspace(-8.0, 8.0, 8)
    xx, yy = np.meshgrid(x, y)
    nodes = np.c_[xx.ravel(), yy.ravel()]
    p0 = np.array([5.0, 0.35])
    nodes = np.vstack([nodes, p0])
    elems = Delaunay(nodes).simplices
    cent = nodes[elems].mean(axis=1)
    keep = (
        (cent[:, 0] >= 0.0)
        & (cent[:, 0] <= 20.0)
        & (cent[:, 1] >= -8.0)
        & (cent[:, 1] <= 8.0)
    )
    return rebuild_tri_mesh(nodes, elems[keep], tip_centers=[p0]), p0


def _backend():
    return CertifiedAtomicPathCorridorCZMBackendV10051839(
        geom=_geom(),
        penalty_normal_Pa_per_m=1.0e18,
        penalty_tangent_Pa_per_m=1.0e18,
        max_angle_error_deg=35.0,
        event_damage=1.0,
        min_area_ratio=0.08,
        min_triangle_quality=0.035,
        max_node_move_factor=1.75,
        max_hrefine_subsegments=512,
    )


def _failure_summary() -> str:
    audit = audit_payload()
    event = audit["events"][-1] if audit.get("events") else {}
    attempts = event.get("corridor_attempts", [])
    histogram = Counter(
        (str(row.get("stage")), str(row.get("reason"))) for row in attempts
    )
    ranked = [
        {"stage": stage, "reason": reason, "count": count}
        for (stage, reason), count in histogram.most_common()
    ]
    tail = attempts[-12:]
    return json.dumps(
        {
            "failure": event.get("failure"),
            "attempt_count": len(attempts),
            "rejection_histogram": ranked,
            "last_attempts": tail,
        },
        indent=2,
        sort_keys=True,
        default=str,
    )


def test_atomic_corridor_commits_exact_endpoint_and_length(monkeypatch):
    monkeypatch.setenv("ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY", "0.035")
    monkeypatch.setenv("ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO", "0.08")
    reset_audit()
    mesh, p0 = _mesh()
    backend = _backend()
    direction = np.array([0.985, 0.1726])
    direction /= np.linalg.norm(direction)
    requested = 5.601322683862894
    p1 = p0 + requested * direction
    damage = np.zeros(mesh.nn)
    displacement = np.zeros(mesh.ndof)

    result = backend.advance(
        mesh=mesh,
        boundary=make_boundary_data(mesh, _geom()),
        damage=damage,
        displacement=displacement,
        p0=p0,
        p1=p1,
        direction=direction,
        front_id=0,
    )

    assert result.inserted is True, _failure_summary()
    assert abs(float(result.moved) - requested) <= 1.0e-10
    assert np.linalg.norm(backend.tip_nodes[0][2] - p1) <= 1.0e-10
    assert result.elem_parent_map is not None
    assert len(result.elem_parent_map) == result.mesh.ne
    assert np.all(result.elem_parent_map >= 0)
    assert np.all(result.elem_parent_map < mesh.ne)

    q = backend._triangle_quality(result.mesh.nodes, result.mesh.elems)
    assert float(np.min(q)) >= 0.035
    incidence = np.bincount(result.mesh.elems.ravel(), minlength=result.mesh.nn)
    assert np.all(incidence > 0)

    certificate = result.v10051839_quality_certificate
    assert certificate["min_local_support_area_ratio"] >= 0.08
    assert certificate["min_triangle_quality"] >= 0.035
    assert certificate["area_ratio_threshold_relaxed"] is False
    assert (
        certificate["old_parent_ratio_is_state_transfer_diagnostic_not_immediate_child_gate"]
        is True
    )

    audit = audit_payload()
    assert audit["event_count"] == 1
    assert audit["success_count"] == 1
    event = audit["events"][0]
    assert event["exact_endpoint_preserved"] is True
    assert event["atomic_transaction"] is True
    assert event["equal_length_physical_partition_retry"] is False
    assert event["topology_subsegment_count"] >= 1
    assert abs(event["length_error_m"]) <= 1.0e-10
    assert abs(event["endpoint_error_m"]) <= 1.0e-10


def test_short_exact_event_does_not_create_residual_ray_fragment(monkeypatch):
    monkeypatch.setenv("ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY", "0.035")
    monkeypatch.setenv("ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO", "0.08")
    reset_audit()
    mesh, p0 = _mesh()
    backend = _backend()
    direction = np.array([0.9906378371905247, 0.13651620975722836])
    requested = 2.2973400956249
    p1 = p0 + requested * direction

    result = backend.advance(
        mesh=mesh,
        boundary=make_boundary_data(mesh, _geom()),
        damage=np.zeros(mesh.nn),
        displacement=np.zeros(mesh.ndof),
        p0=p0,
        p1=p1,
        direction=direction,
        front_id=0,
    )
    assert result.inserted is True, _failure_summary()
    event = audit_payload()["events"][0]
    lengths = [row["length_m"] for row in event["topology_subsegments"]]
    assert abs(sum(lengths) - requested) <= 1.0e-10
    # Numerical support is generated for the complete event at once; there is
    # no inherited-edge residual such as the archived 0.1666 micrometre tail.
    assert min(lengths) > 1.0e-6


def test_failed_corridor_rolls_back_backend_state(monkeypatch):
    monkeypatch.setenv("ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY", "0.999999")
    monkeypatch.setenv("ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO", "0.999999")
    reset_audit()
    mesh, p0 = _mesh()
    backend = _backend()
    p1 = p0 + np.array([5.0, 0.75])
    before = backend._transaction_snapshot()

    result = backend.advance(
        mesh=mesh,
        boundary=make_boundary_data(mesh, _geom()),
        damage=np.zeros(mesh.nn),
        displacement=np.zeros(mesh.ndof),
        p0=p0,
        p1=p1,
        direction=(p1 - p0) / np.linalg.norm(p1 - p0),
        front_id=0,
    )
    assert result.inserted is False
    assert result.mesh is mesh
    after = backend._transaction_snapshot()
    assert after["n_cohesive"] == before["n_cohesive"]
    assert after["event_counter"] == before["event_counter"]
    assert after["n_log"] == before["n_log"]
    event = audit_payload()["events"][0]
    assert event["success"] is False
    assert event["failure"] == "no_feasible_atomic_path_corridor"
