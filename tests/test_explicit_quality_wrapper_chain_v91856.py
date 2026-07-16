from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from arrhenius_fracture import crack_backend as cb
from arrhenius_fracture import mode_i_first_passage_v9_18_5_6 as m


class _FakeBackend:
    min_triangle_quality = 0.035
    min_area_ratio = 0.08

    def __init__(self, quality=0.5):
        self.quality = float(quality)
        self.cohesive_network = SimpleNamespace(elements=[])
        self.advance_log = [{}]
        self.tip_nodes = {}
        self.rollback_count = 0

    def _transaction_snapshot(self):
        return {"token": 1}

    def _transaction_rollback(self, snap):
        assert snap == {"token": 1}
        self.rollback_count += 1

    def _triangle_quality(self, nodes, elems):
        return np.full(len(elems), self.quality, dtype=float)


def _mesh():
    nodes = np.array([[0.0, 0.0], [0.1, 0.0], [0.0, 1.0]], dtype=float)
    elems = np.array([[0, 1, 2]], dtype=int)
    area = np.array([0.05], dtype=float)
    return SimpleNamespace(
        nodes=nodes,
        elems=elems,
        area_e=area,
        nn=3,
        ne=1,
        hbar=0.7,
        hbar_tip=0.7,
    )


def _kwargs(mesh):
    return {
        "mesh": mesh,
        "boundary": object(),
        "damage": np.zeros(mesh.nn),
        "displacement": np.zeros(2 * mesh.nn),
        "p0": np.array([0.0, 0.0]),
        "p1": np.array([0.1, 0.0]),
        "front_id": 0,
    }


def test_resolution_excess_is_warning_but_all_quality_checks_execute(monkeypatch):
    monkeypatch.setenv("ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY", "0.035")
    monkeypatch.setenv("ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO", "0.08")
    monkeypatch.setenv("ARRHENIUS_MAX_TIP_H_OVER_DA", "0.75")
    mesh = _mesh()
    backend = _FakeBackend(quality=0.5)
    result = cb.CrackAdvanceResult(
        mesh=mesh,
        boundary=object(),
        damage=np.zeros(mesh.nn),
        displacement=np.zeros(2 * mesh.nn),
        moved=0.1,
        inserted=True,
        elem_parent_map=np.array([0], dtype=int),
    )
    m._AUDIT["accepted_events"] = []
    m._AUDIT["resolution_warnings"] = []
    m._AUDIT["quality_vetoes"] = []
    m._strict_quality_advance_v91856._original = lambda self, *a, **k: result

    out = m._strict_quality_advance_v91856(backend, **_kwargs(mesh))

    assert out.inserted is True
    assert backend.rollback_count == 0
    assert len(m._AUDIT["accepted_events"]) == 1
    row = m._AUDIT["accepted_events"][0]
    assert row["min_triangle_quality"] == pytest.approx(0.5)
    assert row["min_child_area_ratio"] == pytest.approx(1.0)
    assert row["active_tip_h_over_da"] > 0.75
    assert row["resolution_warning"] is True
    assert len(m._AUDIT["resolution_warnings"]) == 1
    assert backend.advance_log[-1]["v91856_quality_gate_passed"] is True


def test_triangle_quality_floor_remains_a_hard_veto(monkeypatch):
    monkeypatch.setenv("ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY", "0.035")
    monkeypatch.setenv("ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO", "0.08")
    monkeypatch.setenv("ARRHENIUS_MAX_IDENTICAL_GEOMETRY_VETOES", "12")
    mesh = _mesh()
    backend = _FakeBackend(quality=0.01)
    result = cb.CrackAdvanceResult(
        mesh=mesh,
        boundary=object(),
        damage=np.zeros(mesh.nn),
        displacement=np.zeros(2 * mesh.nn),
        moved=0.1,
        inserted=True,
        elem_parent_map=np.array([0], dtype=int),
    )
    m._AUDIT["accepted_events"] = []
    m._AUDIT["quality_vetoes"] = []
    m._strict_quality_advance_v91856._original = lambda self, *a, **k: result

    out = m._strict_quality_advance_v91856(backend, **_kwargs(mesh))

    assert out.inserted is False
    assert out.reason.startswith("v91856_quality_veto:triangle_quality=")
    assert backend.rollback_count == 1
    assert len(m._AUDIT["accepted_events"]) == 0
    assert len(m._AUDIT["quality_vetoes"]) == 1


def test_consecutive_veto_guard_raises_without_signature_matching(monkeypatch):
    monkeypatch.setenv("ARRHENIUS_MAX_IDENTICAL_GEOMETRY_VETOES", "3")
    backend = _FakeBackend()
    kwargs = _kwargs(_mesh())
    r1 = cb.CrackAdvanceResult(
        mesh=kwargs["mesh"], boundary=kwargs["boundary"], damage=kwargs["damage"],
        displacement=kwargs["displacement"], moved=0.0, inserted=False,
        reason="first_reason",
    )
    r2 = cb.CrackAdvanceResult(
        mesh=kwargs["mesh"], boundary=kwargs["boundary"], damage=kwargs["damage"],
        displacement=kwargs["displacement"], moved=0.0, inserted=False,
        reason="slightly_different_reason",
    )
    m._record_or_raise(backend, kwargs, r1)
    m._record_or_raise(backend, kwargs, r2)
    with pytest.raises(RuntimeError, match="count=3/3"):
        m._record_or_raise(backend, kwargs, r1)


def test_main_installs_directly_into_v9185_slot():
    source = open(m.__file__, encoding="utf-8").read()
    assert "_v9185._strict_quality_advance = _strict_quality_advance_v91856" in source
    assert "return _v91853.main(user_args)" in source
    assert "return _v91854.main(user_args)" not in source
