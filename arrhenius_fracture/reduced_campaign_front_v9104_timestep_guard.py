"""Population-aware timestep guard for the v9.10.4 reduced front.

Microscopic Peierls/Taylor rates can remain enormous after all mobile/retained
content has vanished, and emission rates can remain enormous after the finite
source inventory is depleted. Such inactive rates must not freeze the global
K/time integration. This focused guard installs the population-aware v9.10.4
step selector without changing any constitutive rate.
"""
from __future__ import annotations

import math

import numpy as np

from .reduced_campaign_front_v9104 import ReducedCampaignFront


def _population_aware_choose_dt(self: ReducedCampaignFront, rates):
    h = self.s.max_dK_substep_MPa_sqrt_m / max(self.s.Kdot_MPa_sqrt_m_s, 1.0e-30)
    lam_c = max(float(rates["lambda_c_s"]), 0.0)
    if lam_c > 0.0:
        h = min(h, self.s.max_action_substep / lam_c)
        h = min(h, max(1.0 - self.B, self.s.event_tolerance) / lam_c)
    if self.plasticity_active:
        if float(np.sum(np.maximum(self.available, 0.0))) > 1.0e-14:
            lam_e = float(np.max(np.asarray(rates["lambda_e_s"], dtype=float)))
            if lam_e > 0.0 and self.s.max_emit_fraction_substep < 1.0:
                h = min(
                    h,
                    -math.log(max(1.0 - self.s.max_emit_fraction_substep, 1.0e-12)) / lam_e,
                )
        active_count = float(
            np.sum(np.maximum(self.mobile, 0.0) + np.maximum(self.retained, 0.0))
        )
        if active_count > 1.0e-14:
            kinetic = float(
                np.max(
                    np.asarray(rates["encounter_rate_s"], dtype=float)
                    + np.asarray(rates["taylor_rate_s"], dtype=float)
                    + np.asarray(rates["velocity_m_s"], dtype=float)
                    / max(self.s.L_pz_m, 1.0e-30)
                )
            )
            if kinetic > 0.0:
                h = min(h, self.s.max_exchange_fraction_substep / kinetic)
    return max(float(h), self.s.min_dt_s)


ReducedCampaignFront._choose_dt = _population_aware_choose_dt

__all__ = ["_population_aware_choose_dt"]
