"""Staggered FEM/kinetic coupling for one active trial cohesive segment.

The stepper performs the specified predictor--corrector without introducing a
cohesive failure criterion.  Mechanics supplies directional K; the front engine
integrates the Arrhenius state; the trial interface receives ``damage=B`` (or
the abrupt compatibility mapping).  One outer call can commit at most one
cohesive checkpoint and returns unused physical time.
"""
from __future__ import annotations

import copy
from dataclasses import asdict, dataclass
import math
from typing import Any, Callable, Mapping

import numpy as np

from .cohesive_trial_state import (
    KineticTrialAdaptiveCZMBackend,
    normalize_opening_coupling,
)

MechanicsSolve = Callable[[], Mapping[str, Any]]
SnapshotCallback = Callable[[], Any]
RestoreCallback = Callable[[Any], None]
FullRollbackCallback = Callable[[Mapping[str, Any]], None]


@dataclass
class KineticCohesiveStepperConfig:
    opening_coupling: str = "clock_linear"
    maximum_damage_change: float = 0.05
    correction_skip_threshold: float = 1.0e-4
    rejection_safety_factor: float = 0.8

    def validate(self) -> "KineticCohesiveStepperConfig":
        self.opening_coupling = normalize_opening_coupling(self.opening_coupling)
        if not (0.0 < float(self.maximum_damage_change) <= 1.0):
            raise ValueError("maximum_damage_change must lie in (0, 1]")
        if not (0.0 <= float(self.correction_skip_threshold) <= 1.0):
            raise ValueError("correction_skip_threshold must lie in [0, 1]")
        if not (0.0 < float(self.rejection_safety_factor) < 1.0):
            raise ValueError("rejection_safety_factor must lie in (0, 1)")
        return self


@dataclass
class KineticCohesiveStepResult:
    accepted: bool
    committed: bool
    front_id: int
    trial_event_id: int
    progress_before: float
    progress_after: float
    damage_before: float
    damage_after: float
    dt_requested_s: float
    dt_consumed_s: float
    dt_unused_s: float
    recommended_dt_s: float | None
    mechanics_predictor: dict[str, Any]
    mechanics_corrector: dict[str, Any] | None
    kinetics: dict[str, Any]
    trial_diagnostics: dict[str, Any]
    rejection_reason: str | None = None

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


class KineticCohesiveStepper:
    """Predictor--corrector driver for one active trial cohesive event."""

    def __init__(self, config: KineticCohesiveStepperConfig | None = None):
        self.config = (config or KineticCohesiveStepperConfig()).validate()
        self.accepted_steps = 0
        self.rejected_steps = 0
        self.committed_events = 0

    @staticmethod
    def _drives(mechanics: Mapping[str, Any]) -> tuple[float, float, np.ndarray | None]:
        def first(*names: str) -> float:
            for name in names:
                if name in mechanics and mechanics[name] is not None:
                    value = float(mechanics[name])
                    if math.isfinite(value):
                        return value
            raise KeyError(f"mechanics result is missing all drive keys {names}")

        K_open = first(
            "K_open_Pa_sqrt_m",
            "KJ_Pa_sqrt_m",
            "anisotropic_KJ_Pa_sqrt_m",
        )
        K_cleave = first(
            "K_cleave_input_Pa_sqrt_m",
            "K_cleave_Pa_sqrt_m",
            "KJ_Pa_sqrt_m",
            "anisotropic_KJ_Pa_sqrt_m",
        )
        weights = mechanics.get("slip_system_weights")
        if weights is not None:
            weights = np.asarray(weights, dtype=float)
        return K_open, K_cleave, weights

    @staticmethod
    def _trial_step_snapshot(
        backend: KineticTrialAdaptiveCZMBackend,
        front_id: int,
    ) -> dict[str, Any]:
        trial = backend.active_trial(front_id)
        if trial is None:
            raise RuntimeError(f"front {front_id} has no active trial interface")
        return {
            "elements": {
                i: copy.deepcopy(backend.cohesive_network.elements[i])
                for i in trial.element_indices
            },
            "logs": {i: copy.deepcopy(backend.advance_log[i]) for i in trial.log_indices},
            "metadata": copy.deepcopy(trial.transaction.metadata),
        }

    @staticmethod
    def _restore_trial_step(
        backend: KineticTrialAdaptiveCZMBackend,
        front_id: int,
        snapshot: Mapping[str, Any],
    ) -> None:
        trial = backend.active_trial(front_id)
        if trial is None:
            raise RuntimeError(f"front {front_id} lost its active trial interface")
        for i, elem in snapshot["elements"].items():
            backend.cohesive_network.elements[int(i)] = copy.deepcopy(elem)
        for i, row in snapshot["logs"].items():
            backend.advance_log[int(i)] = copy.deepcopy(row)
        trial.transaction.metadata = copy.deepcopy(snapshot["metadata"])

    @staticmethod
    def _progress_after(
        engine: Any,
        result: Mapping[str, Any],
        progress_before: float,
    ) -> float:
        if bool(result.get("fired", False)):
            return 1.0
        q = float(getattr(engine, "B", progress_before))
        return float(np.clip(q, progress_before, 1.0))

    @staticmethod
    def _damage_for_progress(progress: float, coupling: str) -> float:
        if coupling == "abrupt":
            return 1.0 if progress >= 1.0 - 1.0e-12 else 0.0
        return float(np.clip(progress, 0.0, 1.0))

    def _recommended_dt(self, dt_s: float, increment: float) -> float:
        ratio = float(self.config.maximum_damage_change) / max(float(increment), 1.0e-300)
        return max(
            float(dt_s) * ratio * float(self.config.rejection_safety_factor),
            np.finfo(float).tiny,
        )

    def advance(
        self,
        *,
        backend: KineticTrialAdaptiveCZMBackend,
        front_engine: Any,
        front_id: int,
        T_K: float,
        dt_s: float,
        solve_mechanics: MechanicsSolve,
        external_snapshot: SnapshotCallback | None = None,
        external_restore: RestoreCallback | None = None,
        on_full_rollback: FullRollbackCallback | None = None,
    ) -> KineticCohesiveStepResult:
        fid = int(front_id)
        trial = backend.active_trial(fid)
        if trial is None:
            raise RuntimeError(f"front {fid} has no active trial interface")
        dt = max(float(dt_s), 0.0)
        q0 = float(trial.progress)
        damage0 = self._damage_for_progress(q0, self.config.opening_coupling)
        trial_step0 = self._trial_step_snapshot(backend, fid)
        engine0 = front_engine.snapshot_kinetic_state()
        external0 = external_snapshot() if external_snapshot is not None else None

        try:
            predictor = dict(solve_mechanics())
            Kopen0, Kcleave0, weights0 = self._drives(predictor)
            provisional = front_engine.integrate_kinetics(
                Kopen0,
                Kcleave0,
                float(T_K),
                dt,
                system_weights=weights0,
            )
            q_provisional = self._progress_after(front_engine, provisional, q0)
            provisional_increment = max(q_provisional - q0, 0.0)

            corrector: dict[str, Any] | None = None
            final = provisional
            q_final = q_provisional
            if provisional_increment > float(self.config.correction_skip_threshold):
                self._restore_trial_step(backend, fid, trial_step0)
                front_engine.restore_kinetic_state(copy.deepcopy(engine0))
                if external_restore is not None:
                    external_restore(copy.deepcopy(external0))

                q_mid = 0.5 * (q0 + q_provisional)
                backend.update_trial_segment(
                    fid, q_mid, coupling=self.config.opening_coupling
                )
                corrector = dict(solve_mechanics())
                Kopen_mid, Kcleave_mid, weights_mid = self._drives(corrector)

                front_engine.restore_kinetic_state(copy.deepcopy(engine0))
                final = front_engine.integrate_kinetics(
                    Kopen_mid,
                    Kcleave_mid,
                    float(T_K),
                    dt,
                    system_weights=weights_mid,
                )
                q_final = self._progress_after(front_engine, final, q0)

            increment = max(q_final - q0, 0.0)
            damage_final = self._damage_for_progress(
                q_final, self.config.opening_coupling
            )
            damage_increment = max(damage_final - damage0, 0.0)
            too_large = (
                self.config.opening_coupling == "clock_linear"
                and damage_increment > float(self.config.maximum_damage_change) + 1.0e-14
            )
            if too_large:
                self._restore_trial_step(backend, fid, trial_step0)
                front_engine.restore_kinetic_state(copy.deepcopy(engine0))
                if external_restore is not None:
                    external_restore(copy.deepcopy(external0))
                self.rejected_steps += 1
                recommended = self._recommended_dt(dt, damage_increment)
                return KineticCohesiveStepResult(
                    accepted=False,
                    committed=False,
                    front_id=fid,
                    trial_event_id=int(trial.event_index),
                    progress_before=q0,
                    progress_after=q0,
                    damage_before=damage0,
                    damage_after=damage0,
                    dt_requested_s=dt,
                    dt_consumed_s=0.0,
                    dt_unused_s=dt,
                    recommended_dt_s=recommended,
                    mechanics_predictor=predictor,
                    mechanics_corrector=corrector,
                    kinetics=dict(final),
                    trial_diagnostics=backend.active_trial_diagnostics(fid),
                    rejection_reason=(
                        "trial_damage_increment_exceeds_limit:"
                        f"{damage_increment:.6e}>{self.config.maximum_damage_change:.6e}"
                    ),
                )

            backend.update_trial_segment(
                fid, q_final, coupling=self.config.opening_coupling
            )
            committed = bool(final.get("fired", False))
            diagnostics = backend.active_trial_diagnostics(fid)
            if committed:
                if int(final.get("n_fire", 0)) != 1:
                    raise RuntimeError(
                        "kinetic cohesive step must commit exactly one checkpoint"
                    )
                diagnostics = backend.active_trial_diagnostics(fid)
                backend.commit_trial_segment(fid)
                self.committed_events += 1

            self.accepted_steps += 1
            return KineticCohesiveStepResult(
                accepted=True,
                committed=committed,
                front_id=fid,
                trial_event_id=int(trial.event_index),
                progress_before=q0,
                progress_after=q_final,
                damage_before=damage0,
                damage_after=damage_final,
                dt_requested_s=dt,
                dt_consumed_s=float(final.get("dt_consumed_s", dt)),
                dt_unused_s=float(final.get("dt_unused_s", 0.0)),
                recommended_dt_s=None,
                mechanics_predictor=predictor,
                mechanics_corrector=corrector,
                kinetics=dict(final),
                trial_diagnostics=diagnostics,
            )
        except Exception:
            if fid in backend.active_trials:
                restored = backend.rollback_trial_segment(
                    fid, front_engine=front_engine
                )
                if on_full_rollback is not None:
                    on_full_rollback(restored)
            if external_restore is not None:
                external_restore(copy.deepcopy(external0))
            raise

    def audit_payload(self) -> dict[str, Any]:
        return {
            "schema": "kinetic_cohesive_stepper_v10_0",
            "config": asdict(self.config),
            "accepted_steps": int(self.accepted_steps),
            "rejected_steps": int(self.rejected_steps),
            "committed_events": int(self.committed_events),
            "one_topology_event_per_geometry_solve": True,
            "unused_physical_time_returned": True,
            "cohesive_failure_criterion_added": False,
        }


__all__ = [
    "KineticCohesiveStepperConfig",
    "KineticCohesiveStepResult",
    "KineticCohesiveStepper",
]
