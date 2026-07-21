"""Stiff-safe integration policy for the v9.12 emergent-GND state.

The homogeneous reactions use bounded/exact updates.  In the spatial model,
Peierls transport and forest encounter/storage are both proportional to the
same velocity and can be many orders of magnitude faster than the crack-growth
protocol.  Advancing those two operators sequentially makes the retained state
depend on the arbitrary split order and feedback timestep.

This implementation therefore advances mobile transport, mobile-to-retained
encounter storage, Taylor release, and first-order recovery in one coupled,
positivity-preserving backward-Euler finite-volume solve.  Emission is applied
at the midpoint of a symmetric split, and opposite-sign annihilation uses the
exact frozen-coefficient pair-reaction solution.  The timestep resolves only
constitutive feedback, not every Peierls jump or cell transit.

No candidate-specific cap, backstress multiplier, shielding coefficient, or
saturation law is introduced.
"""
from __future__ import annotations

import math
import os
from typing import Mapping

import numpy as np
from scipy.linalg import solve_banded

from .emergent_gnd_state_v912 import EmergentGNDState as _BaseState


class EmergentGNDState(_BaseState):
    """Emergent-GND state with coupled stiff spatial integration."""

    max_feedback_substep_s = float(
        os.environ.get("MPZ_V912_MAX_FEEDBACK_SUBSTEP_S", "0.1")
    )
    coupled_operator_substeps = int(
        os.environ.get("MPZ_V912_COUPLED_OPERATOR_SUBSTEPS", "2")
    )

    def integration_metadata(self) -> dict[str, float | int | str]:
        return {
            "spatial_integrator": "coupled_mobile_retained_backward_euler_v1",
            "max_feedback_substep_s": float(self.max_feedback_substep_s),
            "coupled_operator_substeps": int(
                max(self.coupled_operator_substeps, 1)
            ),
        }

    def _substep(
        self,
        rates: Mapping[str, np.ndarray],
        remaining_s: float,
    ) -> float:
        del rates
        remaining = max(float(remaining_s), 0.0)
        if remaining <= 0.0:
            return 0.0
        if remaining <= self.c.min_substep_s:
            return remaining
        dt = min(
            remaining,
            max(self.max_feedback_substep_s, self.c.min_substep_s),
        )
        if not np.isfinite(dt) or dt <= 0.0:
            raise RuntimeError(
                "invalid emergent-GND feedback substep: "
                f"dt={dt!r}, remaining={remaining!r}"
            )
        return dt

    @staticmethod
    def _transport_diagonals(
        velocity: np.ndarray,
        dx: float,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return lower/diagonal/upper coefficients of conservative upwind A."""
        v = np.asarray(velocity, dtype=float)
        n = v.size
        face_v = np.empty(n + 1, dtype=float)
        if n > 1:
            face_v[1:-1] = 0.5 * (v[:-1] + v[1:])
        face_v[0] = v[0]
        face_v[-1] = v[-1]
        diagonal = (
            -np.maximum(face_v[1:], 0.0)
            + np.minimum(face_v[:-1], 0.0)
        ) / dx
        lower = np.maximum(face_v[1:n], 0.0) / dx
        upper = -np.minimum(face_v[1:n], 0.0) / dx
        return lower, diagonal, upper

    @classmethod
    def _coupled_banded_matrix(
        cls,
        velocity: np.ndarray,
        encounter_s: np.ndarray,
        taylor_s: np.ndarray,
        recovery_s: float,
        dx: float,
        dt: float,
    ) -> np.ndarray:
        """Build (I-dt*B) for interleaved [mobile_j, retained_j] unknowns."""
        velocity = np.asarray(velocity, dtype=float)
        encounter = np.maximum(np.asarray(encounter_s, dtype=float), 0.0)
        taylor = np.maximum(np.asarray(taylor_s, dtype=float), 0.0)
        n = velocity.size
        if encounter.shape != (n,) or taylor.shape != (n,):
            raise ValueError("coupled rate arrays must match spatial velocity")

        lower_A, diagonal_A, upper_A = cls._transport_diagonals(
            velocity, dx
        )
        size = 2 * n
        banded = np.zeros((5, size), dtype=float)  # lower=upper=2

        for j in range(n):
            mobile = 2 * j
            retained = mobile + 1

            banded[2, mobile] = 1.0 - dt * (
                diagonal_A[j] - encounter[j]
            )
            banded[1, retained] = -dt * taylor[j]
            banded[3, mobile] = -dt * encounter[j]
            banded[2, retained] = 1.0 + dt * (
                taylor[j] + max(float(recovery_s), 0.0)
            )

            if j > 0:
                column = mobile - 2
                banded[4, column] = -dt * lower_A[j - 1]
            if j < n - 1:
                column = mobile + 2
                banded[0, column] = -dt * upper_A[j]
        return banded

    def _coupled_mobile_retained(
        self,
        rates: Mapping[str, np.ndarray],
        dt: float,
    ) -> None:
        """Advance transport/storage/release/recovery as one stiff operator."""
        if dt <= 0.0:
            return
        n_substeps = max(int(self.coupled_operator_substeps), 1)
        h = float(dt) / float(n_substeps)
        recovery = float(rates["recovery_rate_s"])
        velocity_base = np.asarray(rates["velocity_m_s"], dtype=float)
        encounter = np.asarray(rates["encounter_s"], dtype=float)
        taylor = np.asarray(rates["taylor_completion_s"], dtype=float)

        for system in range(self.c.n_systems):
            for q in range(2):
                sign = -1.0 if q == 0 else 1.0
                velocity = sign * velocity_base[system]
                banded = self._coupled_banded_matrix(
                    velocity,
                    encounter[system],
                    taylor[system],
                    recovery,
                    self.dx,
                    h,
                )
                state = np.empty(2 * self.c.n_bins, dtype=float)
                state[0::2] = np.maximum(self.mobile_m2[system, q], 0.0)
                state[1::2] = np.maximum(self.retained_m2[system, q], 0.0)
                for _ in range(n_substeps):
                    state = solve_banded(
                        (2, 2),
                        banded,
                        state,
                        overwrite_ab=False,
                        overwrite_b=False,
                        check_finite=False,
                    )
                    state = np.maximum(state, 0.0)
                self.mobile_m2[system, q] = state[0::2]
                self.retained_m2[system, q] = state[1::2]

        if not (
            np.all(np.isfinite(self.mobile_m2))
            and np.all(np.isfinite(self.retained_m2))
        ):
            raise RuntimeError(
                "coupled emergent-GND operator produced nonfinite state"
            )

    def _emit_exact(
        self,
        rates: Mapping[str, np.ndarray],
        dt: float,
    ) -> float:
        emitted_per_m = 0.0
        emission = np.asarray(rates["emission_rate_s"], dtype=float)
        for system, sign in enumerate(self.c.emission_signs):
            q = 1 if sign > 0 else 0
            fraction = 1.0 - np.exp(
                -np.minimum(emission[system, q] * dt, 700.0)
            )
            emitted = self.source_available_m2[system] * fraction
            self.source_available_m2[system] -= emitted
            self.mobile_m2[system, q] += emitted
            emitted_per_m += float(np.sum(emitted) * self.dx)
        return emitted_per_m

    def _annihilate_exact(
        self,
        rates: Mapping[str, np.ndarray],
        dt: float,
    ) -> float:
        """Exact A+B annihilation for frozen velocity in each system/bin."""
        if dt <= 0.0:
            return 0.0
        capture = self.c.annihilation_capture_radius_b * self.c.b_m
        velocity = np.abs(np.asarray(rates["velocity_m_s"], dtype=float))
        annihilated_per_m = 0.0

        for system in range(self.c.n_systems):
            a = np.maximum(self.retained_m2[system, 0], 0.0)
            b = np.maximum(self.retained_m2[system, 1], 0.0)
            lo = np.minimum(a, b)
            hi = np.maximum(a, b)
            difference = hi - lo
            coefficient = 4.0 * capture * velocity[system]
            new_lo = lo.copy()

            active = (lo > 0.0) & (coefficient > 0.0)
            near_equal = active & (
                difference <= 1.0e-12 * np.maximum(hi, 1.0)
            )
            new_lo[near_equal] = lo[near_equal] / (
                1.0
                + coefficient[near_equal]
                * lo[near_equal]
                * float(dt)
            )

            unequal = active & ~near_equal
            if np.any(unequal):
                exponent = np.minimum(
                    coefficient[unequal]
                    * difference[unequal]
                    * float(dt),
                    700.0,
                )
                ratio = (
                    lo[unequal]
                    / np.maximum(hi[unequal], 1.0e-300)
                    * np.exp(-exponent)
                )
                new_lo[unequal] = (
                    difference[unequal]
                    * ratio
                    / np.maximum(1.0 - ratio, 1.0e-300)
                )

            removed = np.clip(lo - new_lo, 0.0, lo)
            self.retained_m2[system, 0] = a - removed
            self.retained_m2[system, 1] = b - removed
            annihilated_per_m += float(2.0 * np.sum(removed) * self.dx)
        return annihilated_per_m

    def _advance_spatial_step(
        self,
        dt: float,
        K_MPa_sqrt_m: float,
        T_K: float,
    ) -> dict[str, float]:
        """Symmetric split around a coupled stiff transport/storage operator."""
        totals = {"emitted_per_m": 0.0, "annihilated_per_m": 0.0}
        rates_start = self.local_rates(K_MPa_sqrt_m, T_K)
        totals["annihilated_per_m"] += self._annihilate_exact(
            rates_start, 0.5 * dt
        )
        self._coupled_mobile_retained(rates_start, 0.5 * dt)

        rates_mid = self.local_rates(K_MPa_sqrt_m, T_K)
        totals["emitted_per_m"] += self._emit_exact(rates_mid, dt)
        self._coupled_mobile_retained(rates_mid, 0.5 * dt)

        rates_end = self.local_rates(K_MPa_sqrt_m, T_K)
        totals["annihilated_per_m"] += self._annihilate_exact(
            rates_end, 0.5 * dt
        )
        return totals

    def advance_time(
        self,
        duration_s: float,
        K_MPa_sqrt_m: float,
        T_K: float,
    ) -> dict[str, float]:
        # Retain the exact bounded homogeneous implementation for the 0-D gate.
        if self.c.n_bins == 1:
            return super().advance_time(duration_s, K_MPa_sqrt_m, T_K)

        remaining = max(float(duration_s), 0.0)
        totals = {"emitted_per_m": 0.0, "annihilated_per_m": 0.0}
        steps = 0
        while remaining > 0.0:
            steps += 1
            if steps > 10_000_000:
                raise RuntimeError(
                    "coupled emergent-GND integration exceeded max feedback steps"
                )
            dt = self._substep({}, remaining)
            increment = self._advance_spatial_step(
                dt, K_MPa_sqrt_m, T_K
            )
            totals["emitted_per_m"] += increment["emitted_per_m"]
            totals["annihilated_per_m"] += increment["annihilated_per_m"]
            self.time_s += dt
            remaining -= dt
        return totals


__all__ = ["EmergentGNDState"]
