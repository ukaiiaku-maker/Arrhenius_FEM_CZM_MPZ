"""Stiff-safe integration policy for the v9.12 emergent-GND state.

The base state uses exact bounded updates for source depletion, mobile/retained
exchange, and first-order recovery.  Those rates therefore must not directly
set an explicit Euler timestep.  Doing so makes a near-prefactor Arrhenius rate
request trillions of unnecessary substeps even though the corresponding state
update is already exponential and bounded.

This subclass limits only:

* constitutive feedback re-evaluation over a common physical time interval; and
* the explicit 1-D transport CFL condition.

The policy is common numerical physics and is not candidate-specific.
"""
from __future__ import annotations

import os
from typing import Mapping

import numpy as np

from .emergent_gnd_state_v912 import EmergentGNDState as _BaseState


class EmergentGNDState(_BaseState):
    """Base emergent-GND state with stiff-safe operator splitting."""

    max_feedback_substep_s = float(
        os.environ.get("MPZ_V912_MAX_FEEDBACK_SUBSTEP_S", "0.1")
    )
    transport_cfl = float(os.environ.get("MPZ_V912_TRANSPORT_CFL", "0.25"))

    def _substep(
        self,
        rates: Mapping[str, np.ndarray],
        remaining_s: float,
    ) -> float:
        remaining = max(float(remaining_s), 0.0)
        if remaining <= 0.0:
            return 0.0

        # Emission, Taylor exchange, and recovery use exact bounded exponential
        # updates in advance_time.  Re-evaluate their state feedback at a common
        # physical interval instead of resolving every elementary attempt.
        dt = min(remaining, max(self.max_feedback_substep_s, self.c.min_substep_s))

        # Only the spatial donor-cell transport operator is explicit.  The
        # homogeneous 0-D gate has no transport CFL restriction.
        if self.c.n_bins > 1:
            vmax = float(np.max(np.abs(rates["velocity_m_s"])))
            if vmax > 0.0:
                dt = min(dt, self.transport_cfl * self.dx / vmax)

        if not np.isfinite(dt) or dt <= 0.0:
            raise RuntimeError(
                "invalid emergent-GND substep: "
                f"dt={dt!r}, remaining={remaining!r}"
            )
        if dt < self.c.min_substep_s:
            raise RuntimeError(
                "emergent-GND transport CFL is below min_substep_s: "
                f"required_dt={dt:.6e}, min_substep_s={self.c.min_substep_s:.6e}. "
                "Classify this candidate as unresolved/stiff or refine the MPZ grid policy."
            )
        return dt


__all__ = ["EmergentGNDState"]
