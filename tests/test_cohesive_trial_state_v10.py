from __future__ import annotations

import copy
from types import SimpleNamespace

import numpy as np
import pytest

from arrhenius_fracture.cohesive import (
    CohesiveElement,
    CohesiveNetwork,
    cohesive_contribution,
)
from arrhenius_fracture.cohesive_trial_state import KineticCZMTransactionSnapshot


def make_element(damage=0.0, clock=0.0):
    return CohesiveElement(
        plus_nodes=(0, 1),
        minus_nodes=(2, 3),
        normal=np.array([0.0, 1.0]),
        tangent=np.array([1.0, 0.0]),
        length=1.0,
        damage=damage,
        clock=clock,
        status="trial" if damage < 1.0 else "committed",
    )


def test_clock_linear_damage_and_monotonicity():
    elem = make_element()
    elem.set_clock_damage(0.25, "clock_linear")
    assert elem.damage == pytest.approx(0.25)
    assert elem.clock == pytest.approx(0.25)
    with pytest.raises(ValueError):
        elem.set_clock_damage(0.20, "clock_linear")


def test_abrupt_mode_reproduces_old_behavior():
    elem = make_element()
    elem.set_clock_damage(0.9, "abrupt")
    assert elem.damage == 0.0
    elem.set_clock_damage(1.0, "abrupt")
    assert elem.damage == 1.0


def test_damage_one_removes_tension_and_shear_but_keeps_contact():
    network = CohesiveNetwork(elements=[make_element(damage=1.0, clock=1.0)])
    u = np.zeros(8)
    u[0] = u[2] = 1e-6
    u[1] = u[3] = 1e-6
    _, R = cohesive_contribution(network, u, 8)
    assert np.linalg.norm(R) == pytest.approx(0.0)

    u = np.zeros(8)
    u[1] = u[3] = -1e-6
    _, R = cohesive_contribution(network, u, 8)
    assert np.linalg.norm(R) > 0.0


def test_transaction_restores_backend_and_front_state():
    elem = make_element(damage=0.2, clock=0.2)
    backend = SimpleNamespace(
        cohesive_network=CohesiveNetwork(elements=[elem]),
        tip_nodes={0: (1, 2, np.array([0.0, 0.0]))},
        event_counter=4,
        advance_log=[{"event": 1}],
    )
    backend._transaction_snapshot = lambda: {
        "n_cohesive": len(backend.cohesive_network.elements),
        "tip_nodes": copy.deepcopy(backend.tip_nodes),
        "event_counter": backend.event_counter,
        "n_log": len(backend.advance_log),
    }
    engine = SimpleNamespace(value=3)
    engine.snapshot_kinetic_state = lambda: {"value": engine.value}
    engine.restore_kinetic_state = lambda state: setattr(engine, "value", state["value"])
    snap = KineticCZMTransactionSnapshot.capture(
        backend=backend,
        mesh={"nodes": [1]},
        boundary={"b": 2},
        displacement=np.array([1.0, 2.0]),
        damage=np.array([0.0]),
        front_engine=engine,
        bulk_history={"rho": np.array([4.0])},
    )
    backend.cohesive_network.elements.append(make_element())
    backend.tip_nodes[0] = (8, 9, np.array([1.0, 1.0]))
    backend.event_counter = 9
    backend.advance_log.append({"event": 2})
    engine.value = 10
    snap.restore_backend_and_front(backend, engine)
    assert len(backend.cohesive_network.elements) == 1
    assert backend.event_counter == 4
    assert len(backend.advance_log) == 1
    assert engine.value == 3
