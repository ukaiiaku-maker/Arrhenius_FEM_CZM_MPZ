"""2-D v9.11 MPZ front engine with stochastic Arrhenius first passage."""
from __future__ import annotations

import copy
import os

import numpy as np

from .bulk_remesh_transfer_v911 import install_bulk_remesh_transfer_patch
from .mpz_front_engine import MovingProcessZoneFrontEngine as _BaseEngine
from .moving_process_zone_v911 import MovingProcessZoneState
from .stochastic_kinetics_v911 import (
    HazardThresholdStream,
    normalize_event_statistics,
)


install_bulk_remesh_transfer_patch()


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
            getattr(
                self.mpz_config,
                "event_statistics",
                os.environ.get("ARRHENIUS_EVENT_STATISTICS", "deterministic"),
            )
        )
        seed = int(
            getattr(
                self.mpz_config,
                "stochastic_seed",
                os.environ.get("ARRHENIUS_STOCHASTIC_SEED", "1"),
            )
        )
        stream = int(getattr(self.mpz_config, "stochastic_cleavage_stream", 91101))
        self.event_statistics = mode
        self.stochastic_seed = seed
        self._threshold_stream = HazardThresholdStream(mode, seed=seed, stream=stream)

        control = str(
            getattr(
                self.mpz_config,
                "propagation_control",
                os.environ.get("ARRHENIUS_PROPAGATION_CONTROL", "raw"),
            )
        ).strip().lower().replace("-", "_")
        if control in {"reload", "event", "event_to_event", "event_reload"}:
            control = "event_reload"
        if control not in {"raw", "event_reload"}:
            raise ValueError(
                f"unknown propagation control {control!r}; expected raw or event_reload"
            )
        self.propagation_control = control
        self.reload_relative_U = max(
            float(os.environ.get("ARRHENIUS_RELOAD_RELATIVE_U", "1e-4")), 0.0
        )
        self.reload_absolute_U_m = max(
            float(os.environ.get("ARRHENIUS_RELOAD_ABSOLUTE_U_M", "1e-12")), 0.0
        )
        self.reload_relative_K = max(
            float(os.environ.get("ARRHENIUS_RELOAD_RELATIVE_K", "1e-4")), 0.0
        )
        self.reload_absolute_K = max(
            float(os.environ.get("ARRHENIUS_RELOAD_ABSOLUTE_K_PA_SQRT_M", "1e3")), 0.0
        )
        self._reload_until_U_m = None
        self._reload_until_K = None
        self._reload_gate_count = 0

    @property
    def B_target(self) -> float:
        return float(self._threshold_stream.target)

    def set_2d_process_zone_profile(self, profile) -> None:
        self.mpz_state.set_2d_profile(profile)

    def _current_remote_U(self) -> float | None:
        context = getattr(self, "_mm", None)
        loading = getattr(context, "last_loading", {}) if context is not None else {}
        try:
            value = float(loading.get("U_total_m"))
        except (TypeError, ValueError):
            return None
        return value if np.isfinite(value) else None

    def _reload_gate_active(self, K_cleave: float) -> bool:
        if self.propagation_control != "event_reload":
            return False
        if self._reload_until_U_m is None and self._reload_until_K is None:
            return False
        U = self._current_remote_U()
        if U is not None and self._reload_until_U_m is not None:
            ready = U + 1.0e-15 >= float(self._reload_until_U_m)
        else:
            ready = float(K_cleave) + 1.0e-9 >= float(self._reload_until_K or 0.0)
        if ready:
            self._reload_until_U_m = None
            self._reload_until_K = None
            return False
        return True

    def _set_reload_floor(self, K_cleave: float) -> None:
        if self.propagation_control != "event_reload":
            return
        U = self._current_remote_U()
        if U is not None:
            self._reload_until_U_m = max(
                U * (1.0 + self.reload_relative_U),
                U + self.reload_absolute_U_m,
            )
        self._reload_until_K = max(
            float(K_cleave) * (1.0 + self.reload_relative_K),
            float(K_cleave) + self.reload_absolute_K,
        )

    def predict_clock_increment_drives(self, K_cleave, K_emit, T, dt):
        if self._reload_gate_active(float(K_cleave)):
            return 0.0
        dB = float(
            super().predict_clock_increment_drives(K_cleave, K_emit, T, dt)
        )
        remaining = max(float(self.B_target) - float(self.B), 1.0e-12)
        return max(dB / remaining, 0.0)

    def _renew(self, dt):
        """Consume deterministic or Exp(1) renewal thresholds.

        ``B`` remains accumulated integrated hazard. In stochastic mode the event
        surface is ``B_target ~ Exp(1)``; after a firing the used threshold is
        subtracted and a fresh threshold is drawn. The residual action is
        retained, which is required when adaptive topology permits only one edge
        insertion per equilibrium solve.
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
            "reload_until_U_m": self._reload_until_U_m,
            "reload_until_K": self._reload_until_K,
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
        self._reload_until_U_m = snap.get("reload_until_U_m")
        self._reload_until_K = snap.get("reload_until_K")
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
        gated = self._reload_gate_active(float(K_cleave))
        if gated:
            # Continue emission/Peierls/Taylor evolution at the actual load while
            # holding only the next cleavage renewal until a measurable reload.
            # This is an explicit continuation protocol, not a material barrier.
            saved_B = float(self.B)
            saved_nu0 = float(self.f.nu0_c)
            self.B = 0.0
            self.f.nu0_c = 0.0
            try:
                out = super().step_drives(
                    K_cleave, K_emit, T, dt, metadata=metadata
                )
            finally:
                self.f.nu0_c = saved_nu0
            self.B = saved_B
            out.update({
                "fired": False,
                "n_fire": 0,
                "n_fire_available": 0,
                "v_crack": 0.0,
                "B": float(self.B),
                "lambda_c": 0.0,
                "lambda_c_raw": 0.0,
            })
            self._reload_gate_count += 1
        else:
            out = super().step_drives(
                K_cleave, K_emit, T, dt, metadata=metadata
            )
            if bool(out.get("fired", False)):
                self._set_reload_floor(float(K_cleave))

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
            "propagation_control": self.propagation_control,
            "event_reload_gate_active": float(gated),
            "event_reload_gate_count": float(self._reload_gate_count),
            "event_reload_until_U_m": self._reload_until_U_m,
            "event_reload_until_K_Pa_sqrt_m": self._reload_until_K,
        })
        return out

    def export_process_zone_state(self):
        payload = super().export_process_zone_state()
        payload["first_passage_v911"] = {
            "event_statistics": self.event_statistics,
            "stochastic_seed": self.stochastic_seed,
            "B": float(self.B),
            "threshold_stream": copy.deepcopy(self._threshold_stream.state_dict()),
            "propagation_control": self.propagation_control,
            "reload_until_U_m": self._reload_until_U_m,
            "reload_until_K_Pa_sqrt_m": self._reload_until_K,
            "reload_gate_count": self._reload_gate_count,
        }
        return payload


__all__ = ["MovingProcessZone2DFrontEngine"]
