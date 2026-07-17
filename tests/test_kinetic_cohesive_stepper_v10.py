from __future__ import annotations

import copy
from types import SimpleNamespace

from arrhenius_fracture.kinetic_cohesive_stepper import (
    KineticCohesiveStepper,
    KineticCohesiveStepperConfig,
)


class FakeEngine:
    def __init__(self, rate=1.0):
        self.B = 0.0
        self.rate = rate

    def snapshot_kinetic_state(self):
        return {"B": self.B}

    def restore_kinetic_state(self, state):
        self.B = state["B"]

    def integrate_kinetics(self, Kopen, Kcleave, T, dt, system_weights=None):
        need = max(1.0 - self.B, 0.0)
        action = self.rate * dt
        fired = action >= need and need > 0.0
        if fired:
            consumed = need / self.rate
            self.B = 0.0
        else:
            consumed = dt
            self.B += action
        return {
            "fired": fired,
            "n_fire": 1 if fired else 0,
            "dt_consumed_s": consumed,
            "dt_unused_s": max(dt - consumed, 0.0),
        }


class FakeBackend:
    def __init__(self):
        tx = SimpleNamespace(metadata={"last_progress": 0.0})
        self.trial = SimpleNamespace(progress=0.0, event_index=1, transaction=tx)
        self.active_trials = {0: self.trial}
        self.commits = 0

    def active_trial(self, front_id):
        return self.active_trials.get(front_id)

    def update_trial_segment(self, front_id, B, coupling=None):
        self.trial.progress = B
        self.trial.transaction.metadata["last_progress"] = B
        return []

    def active_trial_diagnostics(self, front_id):
        return {"trial_cohesive_damage": self.trial.progress}

    def commit_trial_segment(self, front_id):
        self.commits += 1
        del self.active_trials[front_id]
        return []

    def rollback_trial_segment(self, front_id, front_engine=None):
        del self.active_trials[front_id]
        return {}


def patch_trial_snapshot(monkeypatch):
    monkeypatch.setattr(
        KineticCohesiveStepper,
        "_trial_step_snapshot",
        staticmethod(lambda backend, front_id: {
            "progress": backend.trial.progress,
            "metadata": copy.deepcopy(backend.trial.transaction.metadata),
        }),
    )
    monkeypatch.setattr(
        KineticCohesiveStepper,
        "_restore_trial_step",
        staticmethod(lambda backend, front_id, snap: (
            setattr(backend.trial, "progress", snap["progress"]),
            setattr(backend.trial.transaction, "metadata", copy.deepcopy(snap["metadata"])),
        )),
    )


def mechanics():
    return {"K_open_Pa_sqrt_m": 1.0, "K_cleave_input_Pa_sqrt_m": 1.0}


def test_large_clock_linear_damage_increment_is_rejected(monkeypatch):
    patch_trial_snapshot(monkeypatch)
    backend = FakeBackend()
    engine = FakeEngine(rate=1.0)
    stepper = KineticCohesiveStepper(
        KineticCohesiveStepperConfig(maximum_damage_change=0.05)
    )
    result = stepper.advance(
        backend=backend,
        front_engine=engine,
        front_id=0,
        T_K=700.0,
        dt_s=0.2,
        solve_mechanics=mechanics,
    )
    assert not result.accepted
    assert engine.B == 0.0
    assert backend.trial.progress == 0.0
    assert result.recommended_dt_s < 0.2


def test_one_checkpoint_and_unused_time_are_returned(monkeypatch):
    patch_trial_snapshot(monkeypatch)
    backend = FakeBackend()
    engine = FakeEngine(rate=2.0)
    stepper = KineticCohesiveStepper(
        KineticCohesiveStepperConfig(
            opening_coupling="abrupt",
            maximum_damage_change=1.0,
        )
    )
    result = stepper.advance(
        backend=backend,
        front_engine=engine,
        front_id=0,
        T_K=700.0,
        dt_s=1.0,
        solve_mechanics=mechanics,
    )
    assert result.accepted
    assert result.committed
    assert backend.commits == 1
    assert result.dt_consumed_s == 0.5
    assert result.dt_unused_s == 0.5
