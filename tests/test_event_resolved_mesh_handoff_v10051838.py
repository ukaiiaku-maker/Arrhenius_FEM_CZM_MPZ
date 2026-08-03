from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from arrhenius_fracture import event_resolved_czm_retry_v10051837 as event
from arrhenius_fracture import quality_aware_czm_retry_v10051835 as local


def _state(tag: str = "root") -> local._State:
    mesh = SimpleNamespace(ne=1, nn=3, hbar_tip=1.0e-2, tag=tag)
    return local._State(
        mesh=mesh,
        boundary=SimpleNamespace(tag=tag),
        damage=np.zeros(3),
        displacement=np.zeros(6),
        parent_to_root=np.asarray([0], dtype=int),
    )


class _SnapshotBackend:
    def __init__(self):
        self.rollbacks = []

    def _transaction_snapshot(self):
        return {"tip": "unchanged"}

    def _transaction_rollback(self, snap):
        self.rollbacks.append(snap)


def test_segment_retry_returns_mesh_only_handoff_before_second_endpoint_trial(
    monkeypatch,
):
    local.reset_audit()
    backend = _SnapshotBackend()
    root = _state("root")
    refined = _state("refined")
    calls = []

    def rejected(*args, **kwargs):
        calls.append("endpoint_trial")
        return SimpleNamespace(inserted=False, reason="quality_veto")

    def handoff_refinement(*args, **kwargs):
        return refined, {
            "accepted": True,
            "level": 1,
            "refinement_kind": "endpoint_centered_adjacent_target_cavity_midpoint",
            "exact_ray_march_required_after_refinement": True,
            "predicted_exact_target_pass": False,
        }

    monkeypatch.setattr(local, "_call_on_state", rejected)
    monkeypatch.setattr(local, "_refine_once", handoff_refinement)
    monkeypatch.setattr(local, "_patch_levels", lambda: 4)

    result, levels, failure = local._segment_retry(
        lambda *args, **kwargs: None,
        backend,
        root,
        {"front_id": 0},
        np.asarray([0.0, 0.0]),
        np.asarray([1.0, 0.0]),
        np.asarray([1.0, 0.0]),
    )

    assert isinstance(result, local.MeshOnlyRayHandoff)
    assert result.state is refined
    assert result.record["exact_ray_march_required_after_refinement"] is True
    assert levels == 1
    assert failure is None
    assert calls == ["endpoint_trial"]
    assert len(backend.rollbacks) == 2


class _TransactionBackend:
    max_hrefine_subsegments = 16

    def __init__(self):
        self.tip_nodes = {0: (None, None, np.asarray([0.0, 0.0]))}
        self.advance_log = []
        self.rollbacks = []

    def _transaction_snapshot(self):
        return {"tip": np.asarray(self.tip_nodes[0][2]).copy()}

    def _transaction_rollback(self, snap):
        self.rollbacks.append(snap)

    def _ray_exit_edge(self, mesh, point, direction, front_id, remaining):
        assert mesh.tag == "refined"
        assert point[0] == pytest.approx(0.0)
        return (
            0.4,
            0,
            1,
            np.asarray([0.4, 0.0]),
            0,
            0.5,
        )


def test_outer_transaction_restarts_at_exact_ray_crossing_after_handoff(
    monkeypatch,
):
    local.reset_audit()
    event.reset_audit()
    backend = _TransactionBackend()
    root = _state("root")
    refined = _state("refined")
    calls = []

    monkeypatch.setattr(event, "_local_endpoint_reachable", lambda *args: True)
    monkeypatch.setattr(event, "_max_crossings", lambda self: 16)
    monkeypatch.setattr(local, "_patch_levels", lambda: 4)

    def success_result(state, p0, p1):
        moved = float(np.linalg.norm(np.asarray(p1) - np.asarray(p0)))
        backend.tip_nodes[0] = (None, None, np.asarray(p1, dtype=float).copy())
        return SimpleNamespace(
            mesh=state.mesh,
            boundary=state.boundary,
            damage=state.damage,
            displacement=state.displacement,
            elem_parent_map=np.asarray([0], dtype=int),
            moved=moved,
            inserted=True,
            angle_error_deg=0.0,
            reason="ok",
        )

    def segment_retry(
        original,
        self,
        state,
        root_kwargs,
        p0,
        p1,
        direction,
        **metadata,
    ):
        calls.append(
            {
                "kind": metadata["segment_kind"],
                "p0": np.asarray(p0).copy(),
                "p1": np.asarray(p1).copy(),
            }
        )
        if len(calls) == 1:
            return (
                local.MeshOnlyRayHandoff(
                    refined,
                    {
                        "accepted": True,
                        "exact_ray_march_required_after_refinement": True,
                        "predicted_exact_target_pass": False,
                    },
                ),
                1,
                None,
            )
        return success_result(state, p0, p1), 0, None

    monkeypatch.setattr(local, "_segment_retry", segment_retry)

    kwargs = {
        "mesh": root.mesh,
        "boundary": root.boundary,
        "damage": root.damage,
        "displacement": root.displacement,
        "front_id": 0,
    }
    result, failure, record = event._exact_ray_transaction(
        lambda *args, **kwargs: None,
        backend,
        kwargs,
        np.asarray([0.0, 0.0]),
        np.asarray([1.0, 0.0]),
        np.asarray([1.0, 0.0]),
    )

    assert failure is None
    assert result is not None and result.inserted is True
    assert result.moved == pytest.approx(1.0)
    assert record["success"] is True
    assert record["mesh_only_handoff_count"] == 1
    assert record["segment_count"] == 2
    assert [row["segment_kind"] for row in record["ray_segments"]] == [
        "exact_mesh_ray_crossing",
        "exact_physical_endpoint",
    ]
    assert record["mesh_only_handoffs"][0]["physical_motion_m"] == 0.0
    assert record["mesh_only_handoffs"][0][
        "force_exact_ray_crossing_restart"
    ] is True
    assert calls[0]["kind"] == "exact_physical_endpoint"
    assert calls[1]["kind"] == "exact_mesh_ray_crossing"
    assert calls[1]["p1"][0] == pytest.approx(0.4)
    assert calls[2]["kind"] == "exact_physical_endpoint"
    assert record["length_error_m"] <= 1.0e-12
    assert record["endpoint_error_m"] <= 1.0e-12


def test_audit_declares_mesh_only_handoff_has_no_physical_length():
    payload = event.audit_payload()
    assert payload["mesh_only_refinement_handoff_restart"] is True
    assert payload["mesh_only_handoff_consumes_physical_length"] is False
    assert payload["equal_length_partition_retry"] is False
    assert payload["triangle_quality_floor_relaxed"] is False
    assert payload["child_area_ratio_floor_relaxed"] is False
