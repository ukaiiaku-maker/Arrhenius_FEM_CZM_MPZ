from __future__ import annotations

import os
from types import SimpleNamespace

import pytest

from arrhenius_fracture import mode_i_first_passage_v9_18_5_4 as v54
from arrhenius_fracture import mode_i_first_passage_v9_18_5_5 as v55


class DummyBackend:
    pass


def result(inserted: bool, reason: str = "ok"):
    return SimpleNamespace(inserted=inserted, reason=reason)


def test_consecutive_guard_does_not_depend_on_exact_reason(monkeypatch):
    monkeypatch.setenv("ARRHENIUS_MAX_IDENTICAL_GEOMETRY_VETOES", "3")
    backend = DummyBackend()
    kwargs = {"front_id": 0, "p0": [0.5e-3, 0.0], "p1": [0.505e-3, 0.0]}

    v55._consecutive_veto_guard(backend, kwargs, result(False, "reason_a"))
    v55._consecutive_veto_guard(backend, kwargs, result(False, "reason_b"))
    with pytest.raises(RuntimeError, match="count=3/3"):
        v55._consecutive_veto_guard(backend, kwargs, result(False, "reason_c"))


def test_success_resets_consecutive_guard(monkeypatch):
    monkeypatch.setenv("ARRHENIUS_MAX_IDENTICAL_GEOMETRY_VETOES", "2")
    backend = DummyBackend()
    kwargs = {"front_id": 0, "p0": [0.0, 0.0], "p1": [1.0, 0.0]}
    v55._consecutive_veto_guard(backend, kwargs, result(False, "first"))
    v55._consecutive_veto_guard(backend, kwargs, result(True, "ok"))
    v55._consecutive_veto_guard(backend, kwargs, result(False, "new_first"))
    assert backend._v91855_consecutive_geometry_vetoes == 1


def test_runtime_resolution_threshold_is_audit_only_and_restored(monkeypatch):
    monkeypatch.setenv("ARRHENIUS_MAX_TIP_H_OVER_DA", "0.75")
    old_rows = list(v54._AUDIT.get("accepted_events", []))
    v54._AUDIT["accepted_events"] = []
    v55._AUDIT["resolution_warnings"] = []

    def inherited(self, *args, **kwargs):
        assert os.environ["ARRHENIUS_MAX_TIP_H_OVER_DA"] == "inf"
        v54._AUDIT["accepted_events"].append({"active_tip_h_over_da": 1.31})
        return result(True, "ok")

    v55._quality_with_resolution_audit_only._original = inherited
    try:
        out = v55._quality_with_resolution_audit_only(
            DummyBackend(),
            front_id=0,
            p0=[0.5e-3, 0.0],
            p1=[0.505e-3, 0.0],
        )
        assert out.inserted
        assert os.environ["ARRHENIUS_MAX_TIP_H_OVER_DA"] == "0.75"
        row = v54._AUDIT["accepted_events"][-1]
        assert row["v91855_resolution_threshold_enforced_as_veto"] is False
        assert row["v91855_resolution_warning"] is True
        assert len(v55._AUDIT["resolution_warnings"]) == 1
    finally:
        v54._AUDIT["accepted_events"] = old_rows
