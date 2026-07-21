"""Stiff-safe integration policy for the v9.12 emergent-GND state.

The base state uses exact bounded updates for source depletion, mobile/retained
exchange, and first-order recovery.  Those rates therefore must not directly
set an explicit Euler timestep.

The spatial mobile-density equation is discretized by the same conservative
first-order upwind finite-volume operator as the base implementation, but the
transport operator is advanced implicitly.  This removes the explicit CFL
restriction while preserving nonnegative densities, zero inflow, and free
outflow at the MPZ boundaries.  A small fixed number of backward-Euler
substeps limits large-step temporal diffusion without resolving every Peierls
jump or grid-cell transit.

The policy is common numerical physics and is not candidate-specific.
"""
from __future__ import annotations

import os
from typing import Mapping

import numpy as np
from scipy.linalg import solve_banded

from .emergent_gnd_state_v912 import EmergentGNDState as _BaseState


class EmergentGNDState(_BaseState):
    """Base emergent-GND state with stiff-safe reaction/transport splitting."""

    max_feedback_substep_s = float(
        os.environ.get("MPZ_V912_MAX_FEEDBACK_SUBSTEP_S", "0.1")
    )
    implicit_transport_substeps = int(
        os.environ.get("MPZ_V912_IMPLICIT_TRANSPORT_SUBSTEPS", "4")
    )

    def _substep(
        self,
        rates: Mapping[str, np.ndarray],
        remaining_s: float,
    ) -> float:
        del rates  # Transport is implicit; only feedback re-evaluation limits dt.
        remaining = max(float(remaining_s), 0.0)
        if remaining <= 0.0:
            return 0.0

        # A final remainder below min_substep_s is a floating-point closure of
        # the requested physical interval, not an unresolved transport scale.
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

    @classmethod
    def _advect(
        cls,
        field: np.ndarray,
        velocity: np.ndarray,
        dx: float,
        dt: float,
    ) -> np.ndarray:
        """Advance conservative upwind transport without an explicit CFL limit.

        The face flux convention matches the base donor-cell implementation:
        zero external inflow and free outflow at both MPZ boundaries.  For a
        frozen velocity field, backward Euler applied to the upwind generator
        is positivity preserving.  Repeating a few equal implicit substeps
        improves the approximation when the Courant number is large.
        """
        field = np.asarray(field, dtype=float)
        velocity = np.asarray(velocity, dtype=float)
        if velocity.shape != field.shape:
            raise ValueError("velocity field must match mobile density shape")
        if not np.isfinite(dx) or dx <= 0.0:
            raise ValueError("dx must be finite and positive")
        if not np.isfinite(dt) or dt < 0.0:
            raise ValueError("dt must be finite and nonnegative")
        if dt == 0.0 or field.shape[-1] <= 1:
            return np.maximum(field.copy(), 0.0)

        n_substeps = max(int(cls.implicit_transport_substeps), 1)
        h = float(dt) / float(n_substeps)
        n = field.shape[-1]
        out = np.empty_like(field)

        for system in range(field.shape[0]):
            for q in range(field.shape[1]):
                f = np.maximum(field[system, q], 0.0)
                v = velocity[system, q]

                # Cell-face velocities.  Boundary fluxes use zero inflow and
                # the local boundary-cell velocity for outward motion.
                face_v = np.empty(n + 1, dtype=float)
                face_v[1:-1] = 0.5 * (v[:-1] + v[1:])
                face_v[0] = v[0]
                face_v[-1] = v[-1]

                # Upwind semi-discrete generator A, represented by its three
                # diagonals.  The implicit system is (I - h A) f_new = f_old.
                diagonal_A = (
                    -np.maximum(face_v[1:], 0.0)
                    + np.minimum(face_v[:-1], 0.0)
                ) / dx
                lower_A = np.maximum(face_v[1:n], 0.0) / dx
                upper_A = -np.minimum(face_v[1:n], 0.0) / dx

                banded = np.zeros((3, n), dtype=float)
                banded[0, 1:] = -h * upper_A
                banded[1, :] = 1.0 - h * diagonal_A
                banded[2, :-1] = -h * lower_A

                transported = f.copy()
                for _ in range(n_substeps):
                    transported = solve_banded(
                        (1, 1),
                        banded,
                        transported,
                        overwrite_ab=False,
                        overwrite_b=False,
                        check_finite=False,
                    )
                    transported = np.maximum(transported, 0.0)
                out[system, q] = transported

        if not np.all(np.isfinite(out)):
            raise RuntimeError("implicit emergent-GND transport produced nonfinite state")
        return out


__all__ = ["EmergentGNDState"]
