"""2-D validation front engine for the v9.10.2/v9.10.3 MPZ parameterizations."""
from __future__ import annotations

from .mpz_front_engine import MovingProcessZoneFrontEngine as _BaseEngine
from .moving_process_zone_v911 import MovingProcessZoneState


class MovingProcessZone2DFrontEngine(_BaseEngine):
    # Keep the compatibility token expected by the mixed-mode v8 wrapper.
    state_model = "moving_pz"
    state_model_detail = "moving_pz_v911_independent_shapes_2d_profile"

    def reset(self):
        self.mpz_state = MovingProcessZoneState(self.mpz_config)
        self.N_em = 0.0
        self.B = 0.0
        self.a_adv = 0.0
        self.n_adv = 0
        self.W_emit = 0.0
        self.t = 0.0
        self.K_prev = None
        self._lambda_c_prev = None
        self._K_cleave_prev = None
        self._last_pre_renewal_state = None

    def set_2d_process_zone_profile(self, profile) -> None:
        self.mpz_state.set_2d_profile(profile)

    def step_drives(self, K_cleave, K_emit, T, dt, metadata=None):
        out = super().step_drives(K_cleave, K_emit, T, dt, metadata=metadata)
        out.update({
            "front_state_model_detail": self.state_model_detail,
            "independent_shape_all_four_active": 1.0,
            "bulk_pt_model_v9102_active": 1.0,
            "explicit_GND_backstress_active": 0.0,
            "source_inventory_model": "finite_site_v9102",
        })
        return out


__all__ = ["MovingProcessZone2DFrontEngine"]
