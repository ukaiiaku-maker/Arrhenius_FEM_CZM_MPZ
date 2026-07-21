"""Adapter for coupling the v9.12 emergent-GND state to a 2-D crack solver.

The adapter intentionally does not own the FEM/CZM mechanical solve. The 2-D
solver supplies the current signed-J-derived K, temperature, physical time
increment, and accepted crack advance. The adapter returns the state-generated
resolved GND stress and active-zone shielding to be used in the same mechanical
fixed point as the production solver.
"""
from __future__ import annotations

import copy
from dataclasses import asdict
from typing import Any

import numpy as np

from .emergent_gnd_dbtt_v912 import (
    CandidateParameters,
    CommonPhysics,
    EmergentGNDState,
)


class EmergentGND2DAdapter:
    """Stateful v9.12 constitutive adapter for one active crack front."""

    state_model = "emergent_gnd_peierls_taylor_v912"

    def __init__(
        self,
        candidate: CandidateParameters,
        physics: CommonPhysics,
        *,
        gnd_stress_kernel_Pa_per_m2: np.ndarray | None = None,
        shielding_kernel_MPa_sqrt_m_per_m2: np.ndarray | None = None,
    ) -> None:
        self.candidate = candidate
        self.physics = physics
        self.state = EmergentGNDState(candidate, physics)
        self.state.set_mechanical_kernels(
            gnd_stress_kernel_Pa_per_m2=gnd_stress_kernel_Pa_per_m2,
            shielding_kernel_MPa_sqrt_m_per_m2=(
                shielding_kernel_MPa_sqrt_m_per_m2
            ),
        )
        self.accepted_steps = 0
        self.accepted_crack_events = 0

    def constitutive_feedback(self) -> dict[str, Any]:
        tau = self.state.tau_gnd_Pa()
        return {
            "state_model": self.state_model,
            "K_shield_MPa_sqrt_m": self.state.K_shield_MPa_sqrt_m(),
            "tau_gnd_Pa": tau.copy(),
            "signed_gnd_m2": self.state.signed_gnd_m2().copy(),
            "forest_density_m2": self.state.forest_density_m2().copy(),
            "source_available_fraction": self.state.source_available_fraction(),
            "explicit_N_sat_active": False,
            "independent_backstress_law_active": False,
            "constitutive_K_shield_cap_active": False,
        }

    def evolve_trial(
        self,
        *,
        dt_s: float,
        T_K: float,
        K_applied_MPa_sqrt_m: float,
    ) -> dict[str, Any]:
        """Evaluate one trial increment without committing the state."""
        trial = copy.deepcopy(self.state)
        kinetics = trial.advance_time(dt_s, K_applied_MPa_sqrt_m, T_K)
        tau = trial.tau_gnd_Pa()
        return {
            "trial_state": trial,
            "kinetics": kinetics,
            "K_shield_MPa_sqrt_m": trial.K_shield_MPa_sqrt_m(),
            "tau_gnd_Pa": tau,
            "source_available_fraction": trial.source_available_fraction(),
        }

    def commit_trial(
        self,
        trial_payload: dict[str, Any],
        *,
        accepted_crack_advance_m: float = 0.0,
    ) -> dict[str, Any]:
        trial = trial_payload.get("trial_state")
        if not isinstance(trial, EmergentGNDState):
            raise TypeError("trial_payload does not contain an EmergentGNDState")
        self.state = trial
        self.accepted_steps += 1
        if accepted_crack_advance_m > 0.0:
            self.state.translate_tip(accepted_crack_advance_m)
            self.accepted_crack_events += 1
        return self.constitutive_feedback()

    def evolve_and_commit(
        self,
        *,
        dt_s: float,
        T_K: float,
        K_applied_MPa_sqrt_m: float,
        accepted_crack_advance_m: float = 0.0,
    ) -> dict[str, Any]:
        trial = self.evolve_trial(
            dt_s=dt_s,
            T_K=T_K,
            K_applied_MPa_sqrt_m=K_applied_MPa_sqrt_m,
        )
        return self.commit_trial(
            trial,
            accepted_crack_advance_m=accepted_crack_advance_m,
        )

    def audit(self) -> dict[str, Any]:
        feedback = self.constitutive_feedback()
        feedback.update(
            {
                "candidate_id": self.candidate.candidate_id,
                "accepted_steps": self.accepted_steps,
                "accepted_crack_events": self.accepted_crack_events,
                "extension_m": self.state.extension_m,
                "time_s": self.state.time_s,
                "physics": asdict(self.physics),
            }
        )
        return feedback


__all__ = ["EmergentGND2DAdapter"]
