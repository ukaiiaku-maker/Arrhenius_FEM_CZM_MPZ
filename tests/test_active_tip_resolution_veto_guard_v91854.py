from types import SimpleNamespace

import numpy as np
import pytest

from arrhenius_fracture import mode_i_first_passage_v9_18_5_4 as m


def test_active_tip_one_ring_ignores_broad_stored_hbar_tip():
    # The first triangle is the active endpoint one-ring; the second is a remote
    # coarse element.  The stored value reproduces the rejected run's broad
    # nearest-2%-style hbar_tip, but the active one-ring is well below da.
    nodes = np.array([
        [0.0, 0.0], [1.0, 0.0], [0.0, 1.0],
        [10.0, 0.0], [20.0, 0.0], [10.0, 10.0],
    ])
    elems = np.array([[0, 1, 2], [3, 4, 5]], dtype=int)
    mesh = SimpleNamespace(
        nodes=nodes,
        elems=elems,
        hbar=8.0,
        hbar_tip=5.48628,
    )
    audit = m._active_tip_one_ring_resolution(mesh, np.array([0.0, 0.0]))
    assert audit["legacy_stored_hbar_tip_m"] == pytest.approx(5.48628)
    assert audit["active_tip_h_mean_m"] == pytest.approx((1.0 + 1.0 + np.sqrt(2.0)) / 3.0)
    assert audit["active_tip_h_mean_m"] / 5.0 < 0.75
    assert audit["legacy_stored_hbar_tip_m"] / 5.0 > 1.0


def test_identical_veto_guard_raises_at_configured_limit(monkeypatch):
    monkeypatch.setenv("ARRHENIUS_MAX_IDENTICAL_GEOMETRY_VETOES", "3")
    backend = SimpleNamespace()
    kwargs = {
        "front_id": 0,
        "p0": np.array([0.515e-3, 0.0]),
        "p1": np.array([0.520e-3, 0.0]),
    }
    veto = SimpleNamespace(inserted=False, reason="v91854_quality_veto:active_tip_h_over_da")
    assert m._record_veto_or_raise(backend, kwargs, veto) is veto
    assert m._record_veto_or_raise(backend, kwargs, veto) is veto
    with pytest.raises(RuntimeError, match="count=3/3"):
        m._record_veto_or_raise(backend, kwargs, veto)


def test_success_resets_identical_veto_counter(monkeypatch):
    monkeypatch.setenv("ARRHENIUS_MAX_IDENTICAL_GEOMETRY_VETOES", "3")
    backend = SimpleNamespace()
    kwargs = {"front_id": 0, "p0": np.zeros(2), "p1": np.ones(2)}
    veto = SimpleNamespace(inserted=False, reason="same")
    accepted = SimpleNamespace(inserted=True, reason="ok")
    m._record_veto_or_raise(backend, kwargs, veto)
    assert backend._v91854_identical_veto_count == 1
    m._record_veto_or_raise(backend, kwargs, accepted)
    assert backend._v91854_identical_veto_count == 0
    assert backend._v91854_last_veto_signature is None
