from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from arrhenius_fracture import quality_aware_czm_patch_v10051835 as base
from arrhenius_fracture import quality_aware_czm_ray_handoff_v10051835 as handoff
from arrhenius_fracture import quality_aware_czm_retry_v10051835 as retry


def _state(mesh):
    return base.State(
        mesh=mesh,
        boundary=object(),
        damage=np.zeros(3),
        displacement=np.zeros(6),
        parent_to_root=np.array([0], dtype=int),
    )


def test_base_refinement_success_is_preserved(monkeypatch):
    accepted_state = object()
    accepted_record = {"accepted": True, "refinement_kind": "base"}
    monkeypatch.setattr(
        handoff,
        "_ORIGINAL_REFINE_ONCE",
        lambda *args, **kwargs: (accepted_state, accepted_record),
    )

    state, record = handoff.refine_once(
        object(),
        object(),
        np.array([0.0, 0.0]),
        np.array([1.0, 0.0]),
        0,
        1,
    )

    assert state is accepted_state
    assert record is accepted_record


def test_safe_forward_split_can_hand_off_to_exact_ray_march(monkeypatch):
    root_mesh = SimpleNamespace(
        elems=np.array([[0, 1, 2]], dtype=int),
        nodes=np.array(
            [[0.0, 0.0], [4.0, -4.0], [-1.0, 3.0]],
            dtype=float,
        ),
    )
    root = _state(root_mesh)
    forward_mesh = SimpleNamespace(name="forward")
    backward_mesh = SimpleNamespace(name="backward")
    forward_state = _state(forward_mesh)
    backward_state = _state(backward_mesh)

    monkeypatch.setattr(
        handoff,
        "_ORIGINAL_REFINE_ONCE",
        lambda *args, **kwargs: (
            None,
            {
                "accepted": False,
                "reason": "no_quality_safe_tip_radial_midpoint_bisection",
            },
        ),
    )
    monkeypatch.setattr(
        base,
        "_target_tip_triangle",
        lambda self, mesh, p0, target, front_id: (
            0 if mesh is root_mesh else None
        ),
    )
    monkeypatch.setattr(
        base,
        "_tip_node_in_triangle",
        lambda mesh, elem_id, p0: 0,
    )

    def fake_split(self, state, edge_i, edge_j, qfloor, afloor):
        if edge_j == 1:
            return forward_state, {
                "edge_i": 0,
                "edge_j": 1,
                "midpoint_m": [2.0, -2.0],
                "min_triangle_quality": 0.09,
                "min_immediate_child_area_ratio": 0.5,
            }
        return backward_state, {
            "edge_i": 0,
            "edge_j": 2,
            "midpoint_m": [-0.5, 1.5],
            "min_triangle_quality": 0.10,
            "min_immediate_child_area_ratio": 0.5,
        }

    monkeypatch.setattr(base, "_midpoint_split_candidate", fake_split)

    backend = SimpleNamespace(
        min_triangle_quality=0.035,
        min_area_ratio=0.08,
    )
    state, record = handoff.refine_once(
        backend,
        root,
        np.array([0.0, 0.0]),
        np.array([1.0, -0.1]),
        0,
        1,
    )

    assert state is forward_state
    assert record["accepted"] is True
    assert record["refinement_kind"] == (
        "endpoint_centered_tip_radial_ray_handoff"
    )
    assert record["target_remains_in_tip_triangle"] is False
    assert record["exact_ray_march_required_after_refinement"] is True
    assert record["ray_forward_projection_m"] > 0.0
    assert record["ray_alignment_cosine"] > 0.0
    assert record["predicted_exact_target_pass"] is False


def test_non_midpoint_failure_remains_fail_closed(monkeypatch):
    failure = {"accepted": False, "reason": "no_target_tip_triangle"}
    monkeypatch.setattr(
        handoff,
        "_ORIGINAL_REFINE_ONCE",
        lambda *args, **kwargs: (None, failure),
    )

    state, record = handoff.refine_once(
        object(),
        object(),
        np.array([0.0, 0.0]),
        np.array([1.0, 0.0]),
        0,
        1,
    )

    assert state is None
    assert record is failure


def test_install_changes_only_retry_refinement_hook(monkeypatch):
    sentinel = object()
    monkeypatch.setattr(retry, "_refine_once", sentinel)

    handoff.install()

    assert retry._refine_once is handoff.refine_once
    assert base.refine_once is handoff._ORIGINAL_REFINE_ONCE
