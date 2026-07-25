from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from arrhenius_fracture import crack_backend as cb
from arrhenius_fracture import adaptive_czm_quality_subdivision_v100518 as repair


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


def test_atomic_subdivision_recovers_one_physical_event(monkeypatch):
    backend = _Backend()
    mesh = _Mesh()

    def original(self, **kwargs):
        moved = float(
            np.linalg.norm(np.asarray(kwargs["p1"]) - np.asarray(kwargs["p0"]))
        )
        self.advance_log.append(
            {"x1": float(kwargs["p1"][0]), "y1": float(kwargs["p1"][1])}
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

    def assess(self, old_mesh, result, kwargs, log_start):
        moved = float(
            np.linalg.norm(np.asarray(kwargs["p1"]) - np.asarray(kwargs["p0"]))
        )
        issues = ["child_area_ratio=test"] if moved > 0.6 else []
        return issues, _row(issues), False

    monkeypatch.setattr(repair, "_assess", assess)
    monkeypatch.setattr(repair, "_subdivision_counts", lambda: (2,))
    repair._quality_subdividing_advance_v100518._original = original

    result = repair._quality_subdividing_advance_v100518(
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
    assert result.reason == "ok_atomic_quality_subdivision_2"
    assert len(backend.advance_log) == 2
    assert backend.advance_log[-1]["v100518_subdivision_count"] == 2


def test_recursive_quality_veto_is_not_counted_as_physical_attempt(monkeypatch):
    backend = _Backend()
    mesh = _Mesh()

    def original(self, **kwargs):
        self.advance_log.append({})
        return cb.CrackAdvanceResult(
            mesh=kwargs["mesh"],
            boundary=kwargs["boundary"],
            damage=kwargs["damage"],
            displacement=kwargs["displacement"],
            moved=0.5,
            inserted=True,
            reason="ok",
            elem_parent_map=np.arange(kwargs["mesh"].ne),
        )

    monkeypatch.setattr(
        repair,
        "_assess",
        lambda *a, **k: (["bad"], _row(["bad"]), False),
    )
    monkeypatch.setattr(
        repair._v91856,
        "_record_or_raise",
        lambda *a, **k: (_ for _ in ()).throw(
            AssertionError("nested veto counted")
        ),
    )
    repair._quality_subdividing_advance_v100518._original = original

    result = repair._quality_subdividing_advance_v100518(
        backend,
        mesh=mesh,
        boundary=object(),
        damage=np.zeros(mesh.nn),
        displacement=np.zeros(2 * mesh.nn),
        p0=np.array([0.0, 0.0]),
        p1=np.array([0.5, 0.0]),
        direction=np.array([1.0, 0.0]),
        front_id=0,
        _subdepth=3,
    )

    assert result.inserted is False
    assert backend._v91856_consecutive_geometry_vetoes == 0


def test_v100518_smoke_keeps_exact_failed_run_seeds_and_no_pf_dependency():
    text = open(
        "run_v10_0_5_18_frozen_four_class_300_1000K_20um_smoke.sh"
    ).read()
    assert "mode_i_first_passage_v10_0_5_18_four_class_standalone" in text
    assert (
        "SEED_NAMESPACE=${SEED_NAMESPACE:-v10.0.5.17-frozen-four-class}"
        in text
    )
    assert "PFROOT" not in text
    assert "triangle_quality_floor_relaxed=false" in text
    assert "child_area_ratio_floor_relaxed=false" in text
