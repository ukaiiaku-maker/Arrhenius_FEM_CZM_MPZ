"""2-D v9.11 MPZ front engine with stochastic Arrhenius first passage."""
from __future__ import annotations

import copy

import numpy as np

from .mpz_front_engine import MovingProcessZoneFrontEngine as _BaseEngine
from .moving_process_zone_v911 import MovingProcessZoneState
from .stochastic_kinetics_v911 import (
    HazardThresholdStream,
    normalize_event_statistics,
)


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
        self._last_pre_renewal_event_snapshot = None

        mode = normalize_event_statistics(
            getattr(self.mpz_config, "event_statistics", "deterministic")
        )
        seed = int(getattr(self.mpz_config, "stochastic_seed", 1))
        stream = int(getattr(self.mpz_config, "stochastic_cleavage_stream", 91101))
        self.event_statistics = mode
        self.stochastic_seed = seed
        self._threshold_stream = HazardThresholdStream(mode, seed=seed, stream=stream)

    @property
    def B_target(self) -> float:
        return float(self._threshold_stream.target)

    def set_2d_process_zone_profile(self, profile) -> None:
        self.mpz_state.set_2d_profile(profile)

    def _renew(self, dt):
        """Consume deterministic or Exp(1) renewal thresholds.

        ``B`` remains accumulated integrated hazard. In stochastic mode the event
        surface is ``B_target ~ Exp(1)``; after a firing the used threshold is
        subtracted and a fresh threshold is drawn. Residual action is retained,
        which is required when adaptive topology permits only one edge insertion
        per equilibrium solve.
        """
        self._sync_compat()
        Npre = float(self.N_em)
        Ksh_pre = self.K_shield()
        rpre = self.r_eff()
        mobile_pre = float(self.mpz_state.mobile_count)
        retained_pre = float(self.mpz_state.retained_count)
        site_fraction_pre = float(self.mpz_state.available_site_fraction)
        local_slip_pre = float(self.mpz_state.local_slip_count())
        emitted_total_pre = float(self.mpz_state.emitted_total)
        self._last_pre_renewal_state = self.mpz_state.copy()
        self._last_pre_renewal_event_snapshot = {
            "B": float(self.B),
            "threshold": self._threshold_stream.snapshot(),
            "a_adv": float(self.a_adv),
            "n_adv": int(self.n_adv),
        }

        if not np.isfinite(self.B):
            self.B = 0.0
        target_before = float(self.B_target)
        max_fire = float(getattr(self.f, "max_advances_per_step", float("inf")))
        limit = 10_000_000 if not np.isfinite(max_fire) else max(int(max_fire), 0)
        self.B, nfire, crossed = self._threshold_stream.consume(self.B, max_events=limit)
        fired = nfire >= 1

        wake = {
            "wake_mobile": 0.0,
            "wake_retained": 0.0,
            "wake_slip": 0.0,
            "source_sites_refreshed": 0.0,
        }
        if fired:
            distance = self.f.da * nfire
            wake = self.mpz_state.advance(distance)
            self.a_adv += distance
            self.n_adv += nfire
        self._sync_compat()

        ready_after_limit = bool(self.B + 1.0e-15 >= self.B_target)
        return {
            "fired": bool(fired),
            "n_fire": int(nfire),
            "n_fire_available": int(nfire + (1 if ready_after_limit else 0)),
            "v_crack": self.f.da * nfire / dt if dt > 0 else 0.0,
            "N_em_pre_renewal": Npre,
            "N_em_retained": float(self.N_em),
            "N_em_shed_to_wake": float(wake["wake_retained"]),
            "sigma_back_pre_renewal": float(
                Ksh_pre / np.sqrt(2.0 * np.pi * max(rpre, 1.0e-30))
            ),
            "r_eff_pre_renewal": float(rpre),
            "mpz_K_shield_pre_renewal_Pa_sqrt_m": float(Ksh_pre),
            "mpz_mobile_pre_renewal": mobile_pre,
            "mpz_retained_pre_renewal": retained_pre,
            "mpz_available_site_fraction_pre_renewal": site_fraction_pre,
            "mpz_local_slip_pre_renewal": local_slip_pre,
            "mpz_emitted_total_pre_renewal": emitted_total_pre,
            "dG_emb_pre_renewal_eV": 0.0,
            "mpz_wake_mobile_block": float(wake["wake_mobile"]),
            "mpz_wake_retained_block": float(wake["wake_retained"]),
            "mpz_wake_slip_block": float(wake["wake_slip"]),
            "mpz_source_sites_refreshed_on_advance": float(
                wake["source_sites_refreshed"]
            ),
            "event_statistics": self.event_statistics,
            "stochastic_first_passage_active": float(
                self.event_statistics == "stochastic"
            ),
            "B_target_pre_renewal": target_before,
            "B_thresholds_crossed": ";".join(f"{x:.17g}" for x in crossed),
            "B_target": float(self.B_target),
            "B_fraction_of_target": float(
                self.B / max(self.B_target, 1.0e-300)
            ),
            "stochastic_event_index": float(self._threshold_stream.event_index),
        }

    def restore_geometry_veto(self, n_restore: int) -> None:
        """Rollback MPZ, hazard action, threshold and RNG after a rejected edge."""
        snap = self._last_pre_renewal_event_snapshot
        if snap is None:
            return super().restore_geometry_veto(n_restore)
        if self._last_pre_renewal_state is not None:
            self.mpz_state = self._last_pre_renewal_state.copy()
        self.B = float(snap["B"])
        self._threshold_stream.restore(snap["threshold"])
        self.a_adv = float(snap["a_adv"])
        self.n_adv = int(snap["n_adv"])
        self._sync_compat()
        self._last_pre_renewal_event_snapshot = None

    def clone_split(self, daughter_fraction=0.5):
        child = super().clone_split(daughter_fraction)
        if self.event_statistics == "stochastic":
            # Exponential residual lives are memoryless. Once one physical tip
            # becomes two, give the daughters independent future event streams.
            child._threshold_stream = self._threshold_stream.fork()
            self.B = 0.0
            child.B = 0.0
            self._last_pre_renewal_event_snapshot = None
            child._last_pre_renewal_event_snapshot = None
        return child

    def step_drives(self, K_cleave, K_emit, T, dt, metadata=None):
        out = super().step_drives(K_cleave, K_emit, T, dt, metadata=metadata)
        out.update({
            "front_state_model_detail": self.state_model_detail,
            "independent_shape_all_four_active": 1.0,
            "bulk_pt_model_v9102_active": 1.0,
            "explicit_GND_backstress_active": 0.0,
            "source_inventory_model": "finite_site_v9102",
            "event_statistics": self.event_statistics,
            "stochastic_first_passage_active": float(
                self.event_statistics == "stochastic"
            ),
            "B_target": float(self.B_target),
            "B_fraction_of_target": float(
                self.B / max(self.B_target, 1.0e-300)
            ),
            "stochastic_event_index": float(self._threshold_stream.event_index),
        })
        return out

    def export_process_zone_state(self):
        payload = super().export_process_zone_state()
        payload["first_passage_v911"] = {
            "event_statistics": self.event_statistics,
            "stochastic_seed": self.stochastic_seed,
            "B": float(self.B),
            "threshold_stream": copy.deepcopy(self._threshold_stream.state_dict()),
        }
        return payload


__all__ = ["MovingProcessZone2DFrontEngine"]
