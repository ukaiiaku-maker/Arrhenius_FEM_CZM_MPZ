"""Count-aware timestep guard for the v9.10.4.3 reduced front.

The v9.10.4.1 population-aware guard disabled rate limits only after a
population became exactly negligible.  Continuous crack-advance refresh can
leave a tiny available-source trace, and very large microscopic rates then
force hundreds of thousands of steps even though exhausting that trace changes
an insignificant fraction of the finite source capacity.

This guard limits the absolute source or active-dislocation count that may
change in one step, measured relative to each system's finite source capacity.
The constitutive rates and exact exponential updates are unchanged.
"""
from __future__ import annotations

import math

import numpy as np

from .reduced_campaign_front_v9104 import ReducedCampaignFront


def _dt_for_count_change(
    population: np.ndarray,
    rate_s: np.ndarray,
    allowed_count: np.ndarray,
) -> float:
    """Return a timestep that limits exponential depletion by absolute count.

    If the entire current population is smaller than the allowed count, exact
    depletion of that population in one step is acceptable and imposes no
    timestep restriction.
    """
    pop, rate, allowed = np.broadcast_arrays(
        np.maximum(np.asarray(population, dtype=float), 0.0),
        np.maximum(np.asarray(rate_s, dtype=float), 0.0),
        np.maximum(np.asarray(allowed_count, dtype=float), 0.0),
    )
    active = (rate > 0.0) & (pop > allowed) & (allowed > 0.0)
    if not np.any(active):
        return math.inf
    fraction = np.clip(allowed[active] / pop[active], 0.0, 1.0 - 1.0e-15)
    limits = -np.log1p(-fraction) / rate[active]
    finite = limits[np.isfinite(limits) & (limits > 0.0)]
    return float(np.min(finite)) if finite.size else math.inf


def _count_aware_choose_dt(self: ReducedCampaignFront, rates):
    h = self.s.max_dK_substep_MPa_sqrt_m / max(self.s.Kdot_MPa_sqrt_m_s, 1.0e-30)

    lam_c = max(float(rates["lambda_c_s"]), 0.0)
    if lam_c > 0.0:
        h = min(h, self.s.max_action_substep / lam_c)
        h = min(h, max(1.0 - self.B, self.s.event_tolerance) / lam_c)

    if self.plasticity_active:
        capacity = np.maximum(np.asarray(self.capacity, dtype=float), 0.0)

        emit_allowed = max(float(self.s.max_emit_fraction_substep), 0.0) * capacity
        emit_limit = _dt_for_count_change(
            self.available,
            np.asarray(rates["lambda_e_s"], dtype=float),
            emit_allowed,
        )
        if np.isfinite(emit_limit):
            h = min(h, emit_limit)

        active_count = np.maximum(self.mobile, 0.0) + np.maximum(self.retained, 0.0)
        exchange_allowed = max(float(self.s.max_exchange_fraction_substep), 0.0) * capacity
        recovery_rate = max(float(self.p.get("retained_recovery_rate_s", 0.0)), 0.0)
        kinetic_rate = (
            np.asarray(rates["encounter_rate_s"], dtype=float)
            + np.asarray(rates["taylor_rate_s"], dtype=float)
            + np.asarray(rates["velocity_m_s"], dtype=float)
            / max(self.s.L_pz_m, 1.0e-30)
            + recovery_rate
        )
        exchange_limit = _dt_for_count_change(
            active_count,
            kinetic_rate,
            exchange_allowed,
        )
        if np.isfinite(exchange_limit):
            h = min(h, exchange_limit)

    return max(float(h), self.s.min_dt_s)


ReducedCampaignFront._choose_dt = _count_aware_choose_dt


__all__ = ["_count_aware_choose_dt", "_dt_for_count_change"]
