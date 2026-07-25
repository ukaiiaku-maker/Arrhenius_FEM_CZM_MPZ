from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from arrhenius_fracture import crack_backend as cb
from arrhenius_fracture import adaptive_czm_quality_subdivision_v100518 as uniform
from arrhenius_fracture import adaptive_czm_quality_partition_v100519 as adaptive


class _Mesh:
    def __init__(self):
        self.ne = 1
        self.nn = 3
        self.hbar_tip = 1.0
        self.nodes = np.array([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
        self.elems = np.array([[0, 1, 2]], dtype=int)
        self.area_e = np.array([0.5])


class _Backend:
    def __init__(self):
        self.advance_log = []
        self.cohesive_network = SimpleNamespace(elements=[])
        self._v91856_consecutive_geometry_vetoes = 0
        self._v91856_last_veto_reason = None

    def _transaction_snapshot(self):
        return {"n_log": len(self.advance_log)}

    def _transaction_rollback(self, snap):
        del self.advance_log[snap["n_log"] :]


def _row(issues):
    return {
        "active_tip_h_mean_m": 0.25,
        "min_triangle_quality": 0.1,
        "min_child_area_ratio": 0.1,
        "child_area_ratio_reference": "test",
        "active_tip_h_over_da": 0.5,
        "accepted": not issues,
        "issues": list(issues),
    }


def _original(self, **kwargs):
    moved = float(np.linalg.norm(np.asarray(kwargs["p1"]) - np.asarray(kwargs["p0"])))
    self.advance_log.append(
        {
            "x0": float(kwargs["p0"][0]),
            "x1": float(kwargs["p1"][0]),
            "length_m": moved,
        }
    )
    return cb.CrackAdvanceResult(
        mesh=kwargs["mesh"],
        boundary=kwargs["boundary"],
        damage=kwargs["damage"],
        displacement=kwargs["displacement"],
        moved=moved,
        inserted=True,
        reason="ok",
        elem_parent_map=np.arange(kwargs["mesh"].ne),
    )


def test_adaptive_partition_recovers_when_uniform_counts_are_insufficient(monkeypatch):
    backend = _Backend()
    mesh = _Mesh()

    def assess(self, old_mesh, result, kwargs, log_start):
        moved = float(np.linalg.norm(np.asarray(kwargs["p1"]) - np.asarray(kwargs["p0"])))
        issues = ["triangle_quality=test"] if moved > 0.30 else []
        return issues, _row(issues), False

    monkeypatch.setattr(uniform, "_assess", assess)
    monkeypatch.setattr(uniform, "_subdivision_counts", lambda: (2,))
    monkeypatch.setenv("ARRHENIUS_QUALITY_PARTITION_MAX_DEPTH", "4")
    uniform._quality_subdividing_advance_v100518._original = _original

    with adaptive.installed_adaptive_quality_partition_v100519():
        result = uniform._quality_subdividing_advance_v100518(
            backend,
            mesh=mesh,
            boundary=object(),
            damage=np.zeros(mesh.nn),
            displacement=np.zeros(2 * mesh.nn),
            p0=np.array([0.0, 0.0]),
            p1=np.array([1.0, 0.0]),
            direction=np.array([1.0, 0.0]),
            front_id=0,
        )

    assert result.inserted is True
    assert np.isclose(result.moved, 1.0)
    assert result.reason == "ok_adaptive_quality_partition_4_leaves_depth_2"
    assert len(backend.advance_log) == 4
    assert np.allclose(
        [row["length_m"] for row in backend.advance_log],
        [0.25, 0.25, 0.25, 0.25],
    )
    assert backend.advance_log[-1]["v100519_adaptive_quality_partition_recovery"] is True
    assert backend.advance_log[-1]["v100519_partition_leaf_count"] == 4


def test_adaptive_partition_rolls_back_failed_physical_event_once(monkeypatch):
    backend = _Backend()
    mesh = _Mesh()

    monkeypatch.setattr(
        uniform,
        "_assess",
        lambda *args, **kwargs: (["bad"], _row(["bad"]), False),
    )
    monkeypatch.setattr(uniform, "_subdivision_counts", lambda: (2,))
    monkeypatch.setenv("ARRHENIUS_QUALITY_PARTITION_MAX_DEPTH", "2")

    calls = {"physical_vetoes": 0}

    def record_once(self, kwargs, result):
        calls["physical_vetoes"] += 1
        self._v91856_consecutive_geometry_vetoes += 1
        return result

    monkeypatch.setattr(uniform._v91856, "_record_or_raise", record_once)
    uniform._quality_subdividing_advance_v100518._original = _original

    with adaptive.installed_adaptive_quality_partition_v100519():
        result = uniform._quality_subdividing_advance_v100518(
            backend,
            mesh=mesh,
            boundary=object(),
            damage=np.zeros(mesh.nn),
            displacement=np.zeros(2 * mesh.nn),
            p0=np.array([0.0, 0.0]),
            p1=np.array([1.0, 0.0]),
            direction=np.array([1.0, 0.0]),
            front_id=0,
        )

    assert result.inserted is False
    assert backend.advance_log == []
    assert calls["physical_vetoes"] == 1
    assert backend._v91856_consecutive_geometry_vetoes == 1


def test_v100519_wrapper_keeps_quality_floors_and_exact_event_contract():
    text = open(
        "arrhenius_fracture/mode_i_first_passage_v10_0_5_19_four_class_standalone.py"
    ).read()
    assert "triangle_quality_floor_relaxed\": False" in text
    assert "child_area_ratio_floor_relaxed\": False" in text
    assert "adaptive_czm_exact_event_endpoint_preserved\": True" in text
    assert "adaptive_czm_exact_event_length_preserved\": True" in text
    assert "physical_event_atomic\": True" in text
