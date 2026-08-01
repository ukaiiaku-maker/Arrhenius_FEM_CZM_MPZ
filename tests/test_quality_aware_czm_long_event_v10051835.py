from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from arrhenius_fracture import quality_aware_czm_retry_v10051835 as quality


def _backend():
    return SimpleNamespace(
        advance_log=[],
        tip_nodes={},
        _transaction_snapshot=lambda: {},
        _transaction_rollback=lambda snapshot: None,
    )


def _root_state():
    mesh = SimpleNamespace(ne=1)
    return quality._State(
        mesh=mesh,
        boundary=object(),
        damage=np.zeros(1),
        displacement=np.zeros(2),
        parent_to_root=np.array([0], dtype=int),
    )


def test_distant_endpoint_failure_falls_through_to_partition_retry(monkeypatch):
    quality.reset_audit()
    backend = _backend()
    root = _root_state()
    accepted = SimpleNamespace(inserted=True)
    partition_calls = []

    monkeypatch.setattr(
        quality,
        "_segment_retry",
        lambda *args, **kwargs: (None, 0, "no_target_tip_triangle"),
    )

    def fake_partition(
        original,
        self,
        state,
        root_kwargs,
        p0,
        direction,
        requested,
        partitions,
    ):
        partition_calls.append(partitions)
        return accepted, None, {
            "max_patch_levels": 0,
            "patched_segments": 0,
            "failed_segment_index": None,
        }

    monkeypatch.setattr(quality, "_partition_retry", fake_partition)
    quality.quality_aware_strict_advance_v10051835._original = (
        lambda *args, **kwargs: None
    )

    result = quality.quality_aware_strict_advance_v10051835(
        backend,
        mesh=root.mesh,
        boundary=root.boundary,
        damage=root.damage,
        displacement=root.displacement,
        p0=np.array([0.0, 0.0]),
        p1=np.array([17.65, 0.0]),
        direction=np.array([1.0, 0.0]),
        front_id=0,
    )

    assert result is accepted
    assert partition_calls == [2]
    success = quality.audit_payload()["retry_successes"][-1]
    assert success["retry_kind"] == "partition"
    assert success["partitions"] == 2


def test_each_partition_segment_receives_local_patch_retry(monkeypatch):
    quality.reset_audit()
    backend = _backend()
    root = _root_state()
    calls = []

    def fake_segment(
        original,
        self,
        state,
        root_kwargs,
        p0,
        p1,
        direction,
        *,
        partitions,
        segment_index,
        segment_count,
    ):
        moved = float(np.linalg.norm(np.asarray(p1) - np.asarray(p0)))
        calls.append(
            {
                "partitions": partitions,
                "segment_index": segment_index,
                "segment_count": segment_count,
                "moved": moved,
            }
        )
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
        patch_levels = 2 if segment_index == 1 else 0
        return result, patch_levels, None

    monkeypatch.setattr(quality, "_segment_retry", fake_segment)

    result, failure, summary = quality._partition_retry(
        lambda *args, **kwargs: None,
        backend,
        root,
        {"front_id": 0},
        np.array([0.0, 0.0]),
        np.array([1.0, 0.0]),
        requested=17.65,
        partitions=2,
    )

    assert failure is None
    assert result is not None
    assert result.inserted is True
    assert result.moved == 17.65
    assert len(calls) == 2
    assert [row["segment_index"] for row in calls] == [0, 1]
    assert all(row["segment_count"] == 2 for row in calls)
    assert summary["patched_segments"] == 1
    assert summary["max_patch_levels"] == 2
    assert all(
        row["v10051835_retry_kind"] == "partition_with_segment_patch"
        for row in backend.advance_log
    )


def test_failed_raw_trial_subsegments_are_not_reported_as_committed(monkeypatch):
    quality.reset_audit()
    backend = _backend()

    def fake_legacy(self, **kwargs):
        from arrhenius_fracture import mode_i_first_passage_v9_18_5_6 as v91856

        v91856._AUDIT["accepted_events"].append({"trial": True})
        v91856._AUDIT["resolution_warnings"].append({"trial": True})
        v91856._AUDIT["quality_vetoes"].append({"accepted": False})
        return SimpleNamespace(inserted=False, reason="quality_veto")

    quality.configure_legacy_quality_wrapper(fake_legacy)
    result = quality._checked_call(
        lambda *args, **kwargs: None,
        backend,
        {},
    )

    assert result.inserted is False
    payload = quality.audit_payload()
    assert payload["candidate_vetoes"] == [{"accepted": False}]

    from arrhenius_fracture import mode_i_first_passage_v9_18_5_6 as v91856

    assert v91856._AUDIT["accepted_events"] == []
    assert v91856._AUDIT["resolution_warnings"] == []
