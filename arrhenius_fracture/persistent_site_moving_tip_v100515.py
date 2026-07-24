"""PF-equivalent continuous moving-tip coupling for FEM/CZM v10.0.5.15.

The outer FEM solution remains fixed over one accepted mechanics interval.  The
front-local constitutive state is integrated with the same Strang sequence used
by the PF kinetic moving-tip cell:

    half plastic update
    -> recompute cleavage rate from the evolved shielding/blunting state
    -> fractional crack-tip/MPZ translation da = da_checkpoint*dB
    -> half plastic update.

A cohesive topology event is requested only when the cumulative fractional
progress reaches one checkpoint.  The MPZ has already translated by exactly the
same physical checkpoint distance, so the cohesive crack cannot run ahead of a
rapidly developing plastic zone.
"""
from __future__ import annotations

import copy
import math
from typing import Any

import numpy as np

from .config import EV_TO_J
from .persistent_site_front_engine_v100514 import (
    PersistentSiteMovingProcessZoneFrontEngineV100514,
)
from .persistent_site_pf_update_v100515 import (
    PF_REFERENCE_COMMIT,
    PF_UPDATE_MAP,
    evolve_pf_v10222,
)

MODEL_ID = "FEM_CZM_persistent_site_PF_moving_tip_v10_0_5_15"
COUPLING_SCHEME = "PF_v10_1_kinetic_tip_Strang_fractional_moving_frame"


class PersistentSitePFMovingTipFrontEngineV100515(
    PersistentSiteMovingProcessZoneFrontEngineV100514
):
    """Persistent signed MPZ with PF update order and continuous tip motion."""

    pf_update_map_active = True
    kinetic_tip_cell_active = True
    max_action_substep = 0.02
    max_translation_substep_m = 1.0e-7
    min_substep_s = 1.0e-15
    max_internal_steps = 20000

    def reset(self):
        super().reset()
        self.micro_advance_total_m = 0.0
        self.checkpoint_advance_total_m = 0.0
        self.kinetic_packet_count_mean_total = 0.0
        self.kinetic_packet_variance_total_m2 = 0.0
        self._checkpoint_origin_snapshot: dict[str, Any] | None = None
        self._checkpoint_fired_last_step = False
        self._geometry_veto_snapshot: dict[str, Any] | None = None
        self.kinetic_prediction_calls = 0
        self.kinetic_internal_substeps_total = 0

    @staticmethod
    def _sum_numeric(target: dict[str, float], source: dict[str, Any]) -> None:
        for key, value in source.items():
            if isinstance(value, (bool, np.bool_)):
                continue
            if isinstance(value, (int, float, np.integer, np.floating)):
                target[key] = target.get(key, 0.0) + float(value)

    def _capture_state(self) -> dict[str, Any]:
        return {
            "mpz_state": self.mpz_state.copy(),
            "N_em": float(self.N_em),
            "B": float(self.B),
            "a_adv": float(self.a_adv),
            "n_adv": int(self.n_adv),
            "W_emit": float(self.W_emit),
            "t": float(self.t),
            "K_prev": self.K_prev,
            "lambda_c_prev": self._lambda_c_prev,
            "K_cleave_prev": self._K_cleave_prev,
            "micro_advance_total_m": float(self.micro_advance_total_m),
            "checkpoint_advance_total_m": float(self.checkpoint_advance_total_m),
            "packet_count_mean_total": float(self.kinetic_packet_count_mean_total),
            "packet_variance_total_m2": float(self.kinetic_packet_variance_total_m2),
            "accepted_constitutive_steps": int(self.accepted_constitutive_steps),
            "rejected_constitutive_trials": int(self.rejected_constitutive_trials),
        }

    def _restore_state(self, snapshot: dict[str, Any]) -> None:
        self.mpz_state = snapshot["mpz_state"].copy()
        self.N_em = float(snapshot["N_em"])
        self.B = float(snapshot["B"])
        self.a_adv = float(snapshot["a_adv"])
        self.n_adv = int(snapshot["n_adv"])
        self.W_emit = float(snapshot["W_emit"])
        self.t = float(snapshot["t"])
        self.K_prev = snapshot["K_prev"]
        self._lambda_c_prev = snapshot["lambda_c_prev"]
        self._K_cleave_prev = snapshot["K_cleave_prev"]
        self.micro_advance_total_m = float(snapshot["micro_advance_total_m"])
        self.checkpoint_advance_total_m = float(
            snapshot["checkpoint_advance_total_m"]
        )
        self.kinetic_packet_count_mean_total = float(
            snapshot["packet_count_mean_total"]
        )
        self.kinetic_packet_variance_total_m2 = float(
            snapshot["packet_variance_total_m2"]
        )
        self.accepted_constitutive_steps = int(
            snapshot["accepted_constitutive_steps"]
        )
        self.rejected_constitutive_trials = int(
            snapshot["rejected_constitutive_trials"]
        )
        self._sync_compat()

    def _opening_stress(self, K_emit: float) -> tuple[float, float, float, bool]:
        K_shield = float(self.K_shield())
        radius = float(self.r_eff())
        uncapped = max(float(K_emit) - K_shield, 0.0) / math.sqrt(
            2.0 * math.pi * max(radius, 1.0e-30)
        )
        capped = uncapped
        cap_active = bool(self.f.sigma_cap > 0.0 and uncapped > self.f.sigma_cap)
        if self.f.sigma_cap > 0.0:
            capped = min(capped, float(self.f.sigma_cap))
        return capped, uncapped, K_shield, cap_active

    def _plastic_half_step(
        self,
        *,
        dt_s: float,
        T_K: float,
        K_emit: float,
        drive_factors: np.ndarray,
        tau_signed_Pa: np.ndarray,
    ) -> tuple[dict[str, Any], float]:
        dt = max(float(dt_s), 0.0)
        if dt <= 0.0:
            return {
                "dN_emit": 0.0,
                "dN_trapped": 0.0,
                "dN_released": 0.0,
                "dN_escaped": 0.0,
                "dN_recovered": 0.0,
                "transport_integrator": PF_UPDATE_MAP,
            }, 0.0
        opening, _, _, _ = self._opening_stress(K_emit)
        result = evolve_pf_v10222(
            self.mpz_state,
            dt_s=dt,
            T_K=float(T_K),
            opening_stress_Pa=opening,
            drive_factors=drive_factors,
            tau_signed_Pa=tau_signed_Pa,
            emission_rate_function=self._emission_rate_per_site,
        )
        self.W_emit += (
            opening
            * float(self.b)
            * float(self.f.L_pz)
            * max(float(result.get("dN_emit", 0.0)), 0.0)
        )
        self._sync_compat()
        return result, opening

    def _substep_limit(self, dt_remaining: float, lambda_c_s: float) -> float:
        remaining = max(float(dt_remaining), 0.0)
        h = remaining
        lam = max(float(lambda_c_s), 0.0)
        if lam > 0.0 and math.isfinite(lam):
            h = min(h, float(self.max_action_substep) / lam)
            translation_rate = max(float(self.f.da) * lam, 1.0e-300)
            h = min(h, float(self.max_translation_substep_m) / translation_rate)
            action_remaining = max(1.0 - float(self.B), 0.0)
            if action_remaining > 0.0:
                h = min(h, action_remaining / lam)
        return max(min(h, remaining), min(float(self.min_substep_s), remaining))

    def _integrate_coupled(
        self,
        *,
        K_cleave: float,
        K_emit: float,
        T_K: float,
        dt_s: float,
        drive_factors: np.ndarray,
        tau_signed_Pa: np.ndarray,
    ) -> dict[str, Any]:
        dt_requested = max(float(dt_s), 0.0)
        remaining = dt_requested
        consumed = 0.0
        dB_total = 0.0
        da_total = 0.0
        packet_mean = 0.0
        packet_variance = 0.0
        plastic_totals: dict[str, float] = {}
        advance_totals: dict[str, float] = {}
        fired = False
        microsteps = 0
        last_lambda = 0.0
        last_lambda_raw = 0.0
        last_Gc = 0.0
        last_sigma = self.sigma_tip(K_cleave)
        last_opening = self._opening_stress(K_emit)[0]

        while remaining > 0.0:
            microsteps += 1
            if microsteps > int(self.max_internal_steps):
                raise RuntimeError(
                    "PF moving-tip cell exceeded max_internal_steps; reduce the "
                    "outer timestep without changing the constitutive update map"
                )

            sigma0 = self.sigma_tip(K_cleave)
            lambda0, _, _ = self.lambda_cleave(sigma0, T_K)
            lambda0 = max(float(lambda0), 0.0) if math.isfinite(lambda0) else 0.0
            h = self._substep_limit(remaining, lambda0)
            microstep_start = self._capture_state()

            first, _ = self._plastic_half_step(
                dt_s=0.5 * h,
                T_K=T_K,
                K_emit=K_emit,
                drive_factors=drive_factors,
                tau_signed_Pa=tau_signed_Pa,
            )
            sigma_mid = self.sigma_tip(K_cleave)
            lambda_mid, raw_mid, Gc_mid = self.lambda_cleave(sigma_mid, T_K)
            lambda_mid = (
                max(float(lambda_mid), 0.0) if math.isfinite(lambda_mid) else 0.0
            )
            action_remaining = max(1.0 - float(self.B), 0.0)

            # Match the PF transactional event localization: repeat the first half
            # once at the shortened event interval when shielding changed lambda.
            if lambda_mid > 0.0 and lambda_mid * h > action_remaining + 1.0e-12:
                self._restore_state(microstep_start)
                h = max(action_remaining / lambda_mid, float(self.min_substep_s))
                h = min(h, remaining)
                first, _ = self._plastic_half_step(
                    dt_s=0.5 * h,
                    T_K=T_K,
                    K_emit=K_emit,
                    drive_factors=drive_factors,
                    tau_signed_Pa=tau_signed_Pa,
                )
                sigma_mid = self.sigma_tip(K_cleave)
                lambda_mid, raw_mid, Gc_mid = self.lambda_cleave(sigma_mid, T_K)
                lambda_mid = (
                    max(float(lambda_mid), 0.0)
                    if math.isfinite(lambda_mid)
                    else 0.0
                )

            action_remaining = max(1.0 - float(self.B), 0.0)
            dB = min(lambda_mid * h, action_remaining)
            da = float(self.f.da) * dB
            advance = self.mpz_state.advance(da) if da > 0.0 else {}

            second, last_opening = self._plastic_half_step(
                dt_s=0.5 * h,
                T_K=T_K,
                K_emit=K_emit,
                drive_factors=drive_factors,
                tau_signed_Pa=tau_signed_Pa,
            )
            self._sum_numeric(plastic_totals, first)
            self._sum_numeric(plastic_totals, second)
            self._sum_numeric(advance_totals, advance)

            packet_rate = (
                float(self.f.da)
                / max(float(self.b), 1.0e-30)
                * lambda_mid
            )
            packet_n = packet_rate * h
            packet_var = float(self.b) ** 2 * packet_n
            self.B += dB
            self.micro_advance_total_m += da
            self.kinetic_packet_count_mean_total += packet_n
            self.kinetic_packet_variance_total_m2 += packet_var
            dB_total += dB
            da_total += da
            packet_mean += packet_n
            packet_variance += packet_var
            consumed += h
            remaining = max(remaining - h, 0.0)
            self.t += h
            last_lambda = lambda_mid
            last_lambda_raw = float(raw_mid)
            last_Gc = float(Gc_mid)
            last_sigma = self.sigma_tip(K_cleave)

            if self.B >= 1.0 - 1.0e-10:
                # The physical MPZ has already moved by the complete checkpoint
                # distance over this and earlier fractional increments.
                self.B = max(self.B - 1.0, 0.0)
                self.a_adv += float(self.f.da)
                self.checkpoint_advance_total_m += float(self.f.da)
                self.n_adv += 1
                fired = True
                break
            if h <= 0.0:
                break

        self.kinetic_internal_substeps_total += microsteps
        return {
            "fired": fired,
            "n_fire": 1 if fired else 0,
            "v_crack": da_total / consumed if consumed > 0.0 else 0.0,
            "dB": dB_total,
            "da": da_total,
            "dt_consumed": consumed,
            "dt_unused": max(dt_requested - consumed, 0.0),
            "packet_mean": packet_mean,
            "packet_variance_m2": packet_variance,
            "lambda_c": last_lambda,
            "lambda_c_raw": last_lambda_raw,
            "Gc_J": last_Gc,
            "sigma_tip": last_sigma,
            "sigma_emit_tip": last_opening,
            "plastic": plastic_totals,
            "advance": advance_totals,
            "microsteps": microsteps,
        }

    def predict_clock_increment_drives(self, K_cleave, K_emit, T, dt):
        """Non-mutating PF coupled predictor used by the outer event controller."""
        factors, tau_signed, _ = self._two_channel_drive()
        trial = copy.copy(self)
        trial.mpz_state = self.mpz_state.copy()
        trial._checkpoint_origin_snapshot = None
        trial._geometry_veto_snapshot = None
        trial._checkpoint_fired_last_step = False
        predicted = trial._integrate_coupled(
            K_cleave=float(K_cleave),
            K_emit=float(K_emit),
            T_K=float(T),
            dt_s=max(float(dt), 0.0),
            drive_factors=factors,
            tau_signed_Pa=tau_signed,
        )
        self.kinetic_prediction_calls += 1
        return float(max(predicted["dB"], 0.0))

    def step_drives(
        self,
        K_cleave,
        K_emit,
        T,
        dt,
        metadata: dict[str, Any] | None = None,
    ):
        factors, tau_signed, drive_metadata = self._two_channel_drive()
        step_snapshot = self._capture_state()
        if self._checkpoint_fired_last_step:
            # No geometry-veto callback occurred, so the previous cohesive
            # checkpoint was accepted.  This is the origin of the next renewal.
            self._checkpoint_origin_snapshot = step_snapshot
            self._checkpoint_fired_last_step = False
        elif self._checkpoint_origin_snapshot is None:
            self._checkpoint_origin_snapshot = step_snapshot
        self._geometry_veto_snapshot = None

        N_pre = float(self.N_em)
        radius_pre = float(self.r_eff())
        Kshield_pre = float(self.K_shield())
        try:
            coupled = self._integrate_coupled(
                K_cleave=float(K_cleave),
                K_emit=float(K_emit),
                T_K=float(T),
                dt_s=max(float(dt), 0.0),
                drive_factors=factors,
                tau_signed_Pa=tau_signed,
            )
        except BaseException:
            rejected = int(self.rejected_constitutive_trials) + 1
            self._restore_state(step_snapshot)
            self.rejected_constitutive_trials = rejected
            raise

        self.accepted_constitutive_steps += 1
        self._lambda_c_prev = float(coupled["lambda_c"])
        self._K_cleave_prev = float(K_cleave)
        self.K_prev = float(K_cleave)
        self._sync_compat()
        if coupled["fired"]:
            self._geometry_veto_snapshot = copy.deepcopy(
                self._checkpoint_origin_snapshot
            )
            self._checkpoint_fired_last_step = True

        plastic = coupled["plastic"]
        advance = coupled["advance"]
        state_diag = self.mpz_state.diagnostics(
            self.G, self.nu, self.b, self.f.r0, self.f.c_blunt
        )
        aggregate = np.asarray(
            self.mpz_state.last_emission.get(
                "aggregate_hazard_initial_by_system_s", np.zeros(2)
            ),
            dtype=float,
        )
        opening, opening_uncapped, _, cap_active = self._opening_stress(K_emit)
        sig_cleave = self.sigma_tip(K_cleave)
        G_emit = self.eb.G_barrier(
            np.asarray([max(opening, 0.0)]), float(T), self.b
        )[0]
        out = {
            "fired": bool(coupled["fired"]),
            "n_fire": int(coupled["n_fire"]),
            "v_crack": float(coupled["v_crack"]),
            "N_em_pre_renewal": N_pre,
            "N_em_retained": float(self.N_em),
            "N_em_shed_to_wake": float(
                advance.get("wake_mobile", 0.0)
                + advance.get("wake_retained", 0.0)
            ),
            "sigma_back_pre_renewal": float(
                Kshield_pre / math.sqrt(2.0 * math.pi * max(radius_pre, 1.0e-30))
            ),
            "r_eff_pre_renewal": radius_pre,
            "mpz_K_shield_pre_renewal_Pa_sqrt_m": Kshield_pre,
            "sigma_tip": float(sig_cleave),
            "sigma_emit_tip": float(opening),
            "sigma_back": float(self.sigma_back()),
            "lambda_e": float(np.sum(aggregate)),
            "lambda_c": float(coupled["lambda_c"]),
            "lambda_c_raw": float(coupled["lambda_c_raw"]),
            "B": float(self.B),
            "N_em": float(self.N_em),
            "r_eff": float(self.r_eff()),
            "dG_emb_eV": 0.0,
            "G_cleave_eff_eV": float(coupled["Gc_J"] / EV_TO_J),
            **self.cleavage_diagnostics(sig_cleave, T),
            "G_emit_eV": float(G_emit / EV_TO_J),
            "W_emit": float(self.W_emit),
            "sigma_tip_uncapped": float(opening_uncapped),
            "sigma_cap_active": bool(cap_active),
            "dN_emit_raw": float(plastic.get("dN_emit", 0.0)),
            "dN_cap_active": False,
            "N_sat_factor": 1.0,
            "N_sat_active": False,
            "front_state_model_code": 3.0,
            "persistent_site_source_active": True,
            "persistent_source_inventory_active": False,
            "source_depletion_active": False,
            "source_refresh_active": False,
            "available_site_fraction": 1.0,
            "source_sites_refreshed": 0.0,
            "trial_commit_state_active": True,
            "accepted_constitutive_steps": self.accepted_constitutive_steps,
            "rejected_constitutive_trials": self.rejected_constitutive_trials,
            "two_channel_drive_factors": factors.tolist(),
            "two_channel_tau_signed_Pa": tau_signed.tolist(),
            "two_channel_drive_reliable": True,
            "two_channel_names": drive_metadata.get("two_channel_names"),
            "pf_update_map_active": True,
            "pf_reference_commit": PF_REFERENCE_COMMIT,
            "kinetic_tip_cell_active": True,
            "kinetic_coupling_scheme": COUPLING_SCHEME,
            "kinetic_micro_advance_step_m": float(coupled["da"]),
            "kinetic_micro_advance_total_m": float(self.micro_advance_total_m),
            "kinetic_checkpoint_progress_m": float(self.B * self.f.da),
            "kinetic_checkpoint_committed_total_m": float(
                self.checkpoint_advance_total_m
            ),
            "kinetic_dt_consumed_s": float(coupled["dt_consumed"]),
            "kinetic_dt_unused_s": float(coupled["dt_unused"]),
            "kinetic_internal_substeps": int(coupled["microsteps"]),
            "kinetic_prediction_calls": int(self.kinetic_prediction_calls),
            "microstructure_advance_precedes_cohesive_checkpoint": True,
            "mpz_advance_repeated_at_checkpoint": False,
            **plastic,
            **advance,
            **state_diag,
        }
        if metadata:
            out.update(metadata)
        return out

    def restore_geometry_veto(self, n_restore: int) -> None:
        """Fail closed if the cohesive backend rejects a continuous-tip checkpoint."""
        if self._geometry_veto_snapshot is not None:
            rejected = int(self.rejected_constitutive_trials) + 1
            self._restore_state(self._geometry_veto_snapshot)
            self.rejected_constitutive_trials = rejected
        raise RuntimeError(
            "adaptive-CZM rejected a PF continuous moving-tip checkpoint; the "
            "renewal-origin constitutive state was restored and the run is stopped "
            "rather than leaving the MPZ ahead of the cohesive crack"
        )

    @classmethod
    def audit_payload(cls) -> dict[str, Any]:
        payload = super().audit_payload()
        payload.update(
            {
                "model_id": MODEL_ID,
                "PF_reference_commit": PF_REFERENCE_COMMIT,
                "PF_update_map": PF_UPDATE_MAP,
                "transport_operator_order": (
                    "emit_exchange_zero_recovery_scalar_advection"
                ),
                "kinetic_tip_cell": True,
                "coupling_scheme": COUPLING_SCHEME,
                "continuous_fractional_MPZ_translation": True,
                "microstructure_advance_precedes_cohesive_checkpoint": True,
                "checkpoint_MPZ_advance_repeated": False,
                "nonmutating_coupled_clock_prediction": True,
                "geometry_veto_policy": "restore_renewal_origin_and_fail_closed",
            }
        )
        return payload


__all__ = [
    "MODEL_ID",
    "COUPLING_SCHEME",
    "PersistentSitePFMovingTipFrontEngineV100515",
]
