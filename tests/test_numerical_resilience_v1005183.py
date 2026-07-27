from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from arrhenius_fracture.crack_backend import (
    AdaptiveCZMBackend,
    CrackAdvanceResult,
)
from arrhenius_fracture.numerical_resilience_v1005183 import (
    CZM_RETRY_SCHEMA,
    PROBE_RETRY_SCHEMA,
    RetryingAdaptiveCZMBackendV1005183,
    make_robust_process_zone_traction_probe,
)


def test_tensor_probe_retries_with_expanded_physical_annulus():
    calls = []

    def original(
        mesh,
        sigma_gp,
        d,
        tip,
        crack_direction,
        radius_m,
        annulus_half_width=0.45,
        sector_half_angle_deg=40.0,
        damage_cutoff=0.85,
        min_elements=4,
    ):
        calls.append((radius_m, damage_cutoff, min_elements))
        if radius_m < 1.2:
            return {
                "reliable": False,
                "n_elements": 1,
                "probe_radius_m": radius_m,
            }
        return {
            "reliable": True,
            "n_elements": 5,
            "probe_radius_m": radius_m,
            "sigma_nn_Pa": 2.0,
            "tau_tn_Pa": 1.0,
            "stress_tensor": np.array([[2.0, 1.0], [1.0, 3.0]]),
        }

    robust = make_robust_process_zone_traction_probe(original)
    result = robust(
        object(),
        np.zeros((3, 1)),
        np.zeros(1),
        np.zeros(2),
        np.array([1.0, 0.0]),
        1.0,
    )

    assert result["reliable"] is True
    assert result["probe_retry_active"] is True
    assert result["probe_retry_level"] == 1
    assert result["probe_retry_schema"] == PROBE_RETRY_SCHEMA
    assert calls[0][0] == pytest.approx(1.0)
    assert calls[1][0] == pytest.approx(1.25)


def test_tensor_probe_returns_best_failed_attempt_without_fabricating_stress():
    def original(
        mesh,
        sigma_gp,
        d,
        tip,
        crack_direction,
        radius_m,
        annulus_half_width=0.45,
        sector_half_angle_deg=40.0,
        damage_cutoff=0.85,
        min_elements=4,
    ):
        return {
            "reliable": False,
            "n_elements": int(round(radius_m * 10.0)),
            "probe_radius_m": radius_m,
        }

    robust = make_robust_process_zone_traction_probe(original)
    result = robust(
        object(),
        np.zeros((3, 1)),
        np.zeros(1),
        np.zeros(2),
        np.array([1.0, 0.0]),
        1.0,
    )

    assert result["reliable"] is False
    assert result["probe_retry_exhausted"] is True
    assert result["probe_retry_schema"] == PROBE_RETRY_SCHEMA
    assert "stress_tensor" not in result


def test_atomic_czm_partition_retry_preserves_total_event_length(monkeypatch):
    original_advance = AdaptiveCZMBackend.advance

    def fake_advance(self, **kwargs):
        mesh = kwargs["mesh"]
        boundary = kwargs["boundary"]
        damage = kwargs["damage"]
        displacement = kwargs["displacement"]
        p0 = np.asarray(kwargs["p0"], dtype=float)
        p1 = np.asarray(kwargs["p1"], dtype=float)
        front_id = int(kwargs.get("front_id", 0))
        length = float(np.linalg.norm(p1 - p0))
        if length > 0.6:
            return CrackAdvanceResult(
                mesh,
                boundary,
                damage,
                displacement,
                0.0,
                False,
                reason="single_segment_quality_failure",
            )
        self.tip_nodes[front_id] = (0, 0, p1.copy())
        self.advance_log.append(
            {
                "front_id": front_id,
                "event_index": len(self.advance_log),
                "x0": float(p0[0]),
                "y0": float(p0[1]),
                "x1": float(p1[0]),
                "y1": float(p1[1]),
                "length_m": length,
                "reason": "ok",
            }
        )
        return CrackAdvanceResult(
            mesh,
            boundary,
            damage,
            displacement,
            length,
            True,
            angle_error_deg=0.0,
            selected_edge_length=length,
            reason="ok",
            elem_parent_map=None,
        )

    monkeypatch.setattr(AdaptiveCZMBackend, "advance", fake_advance)

    backend = RetryingAdaptiveCZMBackendV1005183(
        geom=SimpleNamespace(Lx=10.0, Ly=10.0)
    )
    mesh = SimpleNamespace(ne=1)
    boundary = object()
    damage = np.zeros(1)
    displacement = np.zeros(2)

    result = backend.advance(
        mesh=mesh,
        boundary=boundary,
        damage=damage,
        displacement=displacement,
        p0=np.array([0.0, 0.0]),
        p1=np.array([1.0, 0.0]),
        direction=np.array([1.0, 0.0]),
        front_id=0,
    )

    assert result.inserted is True
    assert result.moved == pytest.approx(1.0)
    assert result.reason == "partition_retry_2"
    assert backend.partition_retry_attempts == 1
    assert backend.partition_retry_successes == 1
    assert backend.partition_retry_failures == 0
    assert len(backend.advance_log) == 2
    assert sum(row["length_m"] for row in backend.advance_log) == pytest.approx(1.0)
    for row in backend.advance_log:
        assert row["physical_event_partition_retry"] is True
        assert row["physical_event_partition_count"] == 2
        assert row["physical_event_retry_schema"] == CZM_RETRY_SCHEMA

    monkeypatch.setattr(AdaptiveCZMBackend, "advance", original_advance)


def test_failed_partition_retry_rolls_back_backend_bookkeeping(monkeypatch):
    original_advance = AdaptiveCZMBackend.advance

    def fake_advance(self, **kwargs):
        mesh = kwargs["mesh"]
        boundary = kwargs["boundary"]
        damage = kwargs["damage"]
        displacement = kwargs["displacement"]
        return CrackAdvanceResult(
            mesh,
            boundary,
            damage,
            displacement,
            0.0,
            False,
            reason="always_fail",
        )

    monkeypatch.setattr(AdaptiveCZMBackend, "advance", fake_advance)

    backend = RetryingAdaptiveCZMBackendV1005183(
        geom=SimpleNamespace(Lx=10.0, Ly=10.0)
    )
    mesh = SimpleNamespace(ne=1)
    result = backend.advance(
        mesh=mesh,
        boundary=object(),
        damage=np.zeros(1),
        displacement=np.zeros(2),
        p0=np.array([0.0, 0.0]),
        p1=np.array([1.0, 0.0]),
        direction=np.array([1.0, 0.0]),
        front_id=0,
    )

    assert result.inserted is False
    assert result.moved == 0.0
    assert "partition_retry_exhausted" in result.reason
    assert backend.partition_retry_attempts == 3
    assert backend.partition_retry_successes == 0
    assert backend.partition_retry_failures == 3
    assert backend.advance_log == []
    assert backend.tip_nodes == {}

    monkeypatch.setattr(AdaptiveCZMBackend, "advance", original_advance)
