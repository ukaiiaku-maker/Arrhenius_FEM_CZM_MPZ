"""v9.11 MPZ state: v9.10.2 kinetics coupled to non-double-counted 2-D profiles."""
from __future__ import annotations

import numpy as np

from .moving_process_zone_v9102 import MovingProcessZoneState as _V9102State
from .process_zone_2d_v911 import ProcessZone2DProfile


class MovingProcessZoneState(_V9102State):
    """Independent-shape MPZ with optional 2-D forest/stress-shape inputs.

    The 2-D scalar density augments the local forest density used by Peierls/Taylor
    transport. It is not interpreted as signed shielding. The direct K shield
    remains the unresolved retained-line integral inherited from v9.10.2.
    """

    state_model = "moving_pz_v911_independent_shapes_2d_profile"

    def __init__(self, cfg):
        super().__init__(cfg)
        self._profile_2d: ProcessZone2DProfile | None = None

    def set_2d_profile(self, profile: ProcessZone2DProfile | None) -> None:
        self._profile_2d = profile

    @staticmethod
    def _resample(values, n):
        q = np.asarray(values, float).reshape(-1)
        if q.size == n:
            return q.copy()
        old = (np.arange(q.size, dtype=float) + 0.5) / max(q.size, 1)
        new = (np.arange(n, dtype=float) + 0.5) / n
        return np.interp(new, old, q)

    def local_forest_density_m2(self) -> np.ndarray:
        local = np.asarray(super().local_forest_density_m2(), float)
        if self._profile_2d is None:
            return local
        floor = float(self.cfg.pt_forest_density_floor_m2)
        bulk = self._resample(self._profile_2d.forest_density_m2, self.n_bins)
        retained_excess = np.maximum(local - floor, 0.0)
        return np.maximum(bulk, floor) + retained_excess

    def local_stress_profile_Pa(self, tip_stress_Pa: float) -> np.ndarray:
        if self._profile_2d is None or not self._profile_2d.reliable:
            return super().local_stress_profile_Pa(tip_stress_Pa)
        shape = self._resample(self._profile_2d.stress_shape, self.n_bins)
        return max(float(tip_stress_Pa), 0.0) * np.maximum(shape, 0.0)

    def diagnostics(self, G_shear, nu, b, r0, c_blunt):
        out = super().diagnostics(G_shear, nu, b, r0, c_blunt)
        if self._profile_2d is None:
            out.update({
                "mpz_2d_profile_active": 0.0,
                "bulk_scalar_rho_used_for_signed_shielding": 0.0,
            })
        else:
            out.update(self._profile_2d.diagnostics())
            out["mpz_2d_profile_active"] = 1.0
            out["bulk_scalar_rho_used_for_signed_shielding"] = 0.0
        out.update({
            "mpz_shielding_source": "unresolved_retained_line_integral",
            "bulk_plastic_shielding_source": "FEM_stress_redistribution_already_in_J",
            "bulk_K_shield_subtracted_again": 0.0,
            "explicit_GND_backstress_active": 0.0,
        })
        return out


__all__ = ["MovingProcessZoneState"]
