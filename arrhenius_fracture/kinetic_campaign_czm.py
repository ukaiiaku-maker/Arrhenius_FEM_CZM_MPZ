"""PF-equivalent moving-tip kinetic state for the FEM/CZM solver.

The new ``kinetic_campaign_czm`` state keeps the existing v9.11 spatial
Peierls--Taylor transport machinery but replaces the front-local closure with
three separated stress channels, the bounded PF campaign source budget, and
continuous cleavage-action translation inside one coarse CZM checkpoint.

This module contains no cohesive critical traction/opening/Gc law. The only
fracture clock is the Arrhenius cleavage action ``B``.
"""
from __future__ import annotations

import copy
from dataclasses import asdict, dataclass
import math
from typing import Any, Mapping

import numpy as np

from .config import KB, EV_TO_J
from .mpz_front_engine_v911 import MovingProcessZone2DFrontEngine
from .moving_process_zone_v911 import MovingProcessZoneState as _TransportState
from .pf_equivalent_material_manifest import PFEquivalentMaterialManifest

MODEL_ID = "FEM_CZM_kinetic_campaign_czm_pf_v10_1_7_1"
STATE_MODEL = "kinetic_campaign_czm"


@dataclass
class KineticCampaignCZMConfig:
    """Numerical coupling controls; none are material fit parameters."""

    max_action_substep: float = 0.02
    max_translation_substep_m: float = 1.0e-7
    min_substep_s: float = 1.0e-15
    max_internal_steps: int = 20000
    coupling_scheme: str = "strang"
    wake_shielding: bool = False
    active_shielding: bool = True
    signed_active_shielding: bool = True
    mobile_shield_fraction: float = 1.0
    backstress_scale: float = 1.0
    source_refresh_scale: float = 1.0

    def validate(self) -> "KineticCampaignCZMConfig":
        if not (0.0 < float(self.max_action_substep) <= 1.0):
            raise ValueError("max_action_substep must lie in (0, 1]")
        if float(self.max_translation_substep_m) <= 0.0:
            raise ValueError("max_translation_substep_m must be positive")
        if float(self.min_substep_s) <= 0.0:
            raise ValueError("min_substep_s must be positive")
        if int(self.max_internal_steps) < 1:
            raise ValueError("max_internal_steps must be positive")
        if str(self.coupling_scheme) != "strang":
            raise ValueError("only coupling_scheme='strang' is supported")
        if float(self.backstress_scale) < 0.0:
            raise ValueError("backstress_scale must be nonnegative")
        if float(self.source_refresh_scale) <= 0.0:
            raise ValueError("source_refresh_scale must be positive")
        return self


def apply_pf_manifest_to_mpz_config(
    cfg: Any,
    manifest: PFEquivalentMaterialManifest,
    kinetic_cfg: KineticCampaignCZMConfig,
) -> Any:
    """Install the PF campaign material fields without adding fitted defaults."""

    emit0 = max(float(manifest.emission.G00_eV), 1.0e-30)
    cfg.n_systems = 2
    cfg.source_sites_per_system = float(manifest.source_sites_per_system)
    cfg.source_recovery_rate_s = 0.0
    cfg.source_refresh_length_m = float(manifest.source_refresh_length_m)
    cfg.source_bin_count = max(2, int(round(0.02 * int(cfg.n_bins))))
    cfg.shielding_orientation_factors = tuple(manifest.orientation_factors)
    cfg.mobile_shield_fraction = float(kinetic_cfg.mobile_shield_fraction)
    cfg.retained_recovery_nu0_s = float(manifest.retained_recovery_rate_s)
    cfg.retained_recovery_barrier_eV = 0.0
    cfg.retained_recovery_activation_volume_b3 = 0.0
    cfg.mobile_recovery_rate_s = 0.0
    cfg.pair_annihilation_rate_per_count_s = 0.0

    cfg.pt_emit_G00_eV = float(manifest.emission.G00_eV)
    cfg.pt_emit_gT_eV_per_K = float(manifest.emission.gT_eV_per_K)
    cfg.pt_emit_sigc0_Pa = float(manifest.emission.sigc0_Pa)
    cfg.pt_emit_sT_Pa_per_K = float(manifest.emission.sT_Pa_per_K)
    cfg.pt_emit_Tref_K = float(manifest.emission.Tref_K)
    cfg.pt_emit_exp_a = float(manifest.emission.alpha)
    cfg.pt_emit_exp_n = float(manifest.emission.exponent)
    cfg.pt_emit_floor_frac = float(manifest.emission.floor_fraction)
    cfg.pt_emit_floor_min_eV = float(manifest.emission.floor_min_eV)
    cfg.pt_emit_floor_max_frac = float(manifest.emission.floor_max_fraction)

    cfg.pt_peierls_energy_ratio = float(manifest.peierls.H0_eV) / emit0
    cfg.pt_peierls_entropy_ratio = float(manifest.peierls.activation_entropy_kB)
    cfg.pt_peierls_exp_a = float(manifest.peierls.alpha)
    cfg.pt_peierls_exp_n = float(manifest.peierls.exponent)
    cfg.pt_peierls_stress_ratio = 1.0
    cfg.pt_peierls_stress_fraction = 1.0 / math.sqrt(3.0)
    cfg.pt_peierls_nu0_s = float(manifest.peierls.attempt_frequency_s)

    cfg.pt_taylor_energy_ratio = float(manifest.taylor.H0_eV) / emit0
    cfg.pt_taylor_entropy_ratio = float(manifest.taylor.activation_entropy_kB)
    cfg.pt_taylor_exp_a = float(manifest.taylor.alpha)
    cfg.pt_taylor_exp_n = float(manifest.taylor.exponent)
    cfg.pt_taylor_stress_ratio = 1.0
    cfg.pt_taylor_stress_fraction = 1.0 / math.sqrt(3.0)
    cfg.pt_taylor_nu0_s = float(manifest.taylor.attempt_frequency_s)
    cfg.pt_taylor_corr_rho_c = float(manifest.taylor_corr_rho_c_m2)
    cfg.pt_taylor_renewal_time_s = 1.0
    cfg.pt_taylor_m_exponent = 1.0
    cfg.pt_taylor_m_scale = float(manifest.taylor_corr_scale)
    cfg.pt_taylor_m_cap = math.inf
    cfg.pt_encounter_efficiency = float(manifest.encounter_efficiency)
    cfg.pt_forest_density_floor_m2 = 5.0e12
    cfg.pt_mobile_fraction = 0.01
    cfg.pt_mobile_saturation_density_m2 = math.inf
    cfg.pt_mobile_density_floor_m2 = 0.0
    cfg.pt_jump_fraction = 1.0
    cfg.pt_jump_length_min_m = 0.0
    cfg.pt_taylor_phi_max = math.inf
    return cfg


class CampaignKineticMPZState(_TransportState):
    """v9.11 spatial transport with the exact PF bounded source law."""

    state_model = STATE_MODEL

    def __init__(
        self,
        cfg: Any,
        manifest: PFEquivalentMaterialManifest,
        *,
        b: float,
        G_Pa: float,
        kinetic_cfg: KineticCampaignCZMConfig,
    ) -> None:
        apply_pf_manifest_to_mpz_config(cfg, manifest, kinetic_cfg)
        super().__init__(cfg)
        self.manifest = manifest
        self._campaign_b = abs(float(b))
        self._campaign_G_Pa = max(float(G_Pa), 0.0)
        self._kinetic_cfg = copy.deepcopy(kinetic_cfg).validate()
        self.site_capacity = np.full(
            self.n_systems, max(float(manifest.source_sites_per_system), 0.0)
        )
        self.available_sites = self.site_capacity.copy()

        self.cumulative_emitted = 0.0
        self.cumulative_refreshed = 0.0
        self.cumulative_trapped = 0.0
        self.cumulative_released = 0.0
        self.cumulative_recovered = 0.0
        self.cumulative_escaped = 0.0
        self.mobile_residence_count_s = 0.0
        self.retained_residence_count_s = 0.0
        self.active_residence_count_s = 0.0
        self.last_source_refresh_fraction = 0.0
        self.last_source_refreshed = 0.0
        self.last_emitted_per_system = np.zeros(self.n_systems)
        self.last_emission_rate_per_system_s = np.zeros(self.n_systems)
        self.last_sigma_back_per_system_Pa = np.zeros(self.n_systems)
        self.last_sigma_emit_per_system_Pa = np.zeros(self.n_systems)
        self.last_local_density_per_system_m2 = np.zeros(self.n_systems)

    def local_backstress_density_m2(self) -> np.ndarray:
        length = max(float(self.cfg.blunting_length_m), float(self.dx), 1.0e-12)
        weights = np.exp(-np.asarray(self.x, dtype=float) / length)
        norm = max(float(np.sum(weights)), 1.0e-30)
        count = np.maximum(self.mobile, 0.0) + np.maximum(self.retained, 0.0)
        near_count = np.sum(count * weights[None, :], axis=1) / norm
        width = max(float(self.cfg.blunting_length_m), float(self.dx), 1.0e-12)
        return np.maximum(near_count / max(float(self.dx) * width, 1.0e-30), 0.0)

    def taylor_backstress_Pa(self) -> np.ndarray:
        rho = self.local_backstress_density_m2()
        tau = (
            float(self._kinetic_cfg.backstress_scale)
            * self._campaign_G_Pa
            * self._campaign_b
            * np.sqrt(np.maximum(rho, 0.0))
        )
        resolved = max(abs(float(self.cfg.pt_taylor_stress_fraction)), 1.0e-6)
        sigma = tau / resolved
        self.last_local_density_per_system_m2 = rho.copy()
        self.last_sigma_back_per_system_Pa = sigma.copy()
        return sigma

    @staticmethod
    def _weights(system_weights: np.ndarray | None, n_systems: int) -> np.ndarray:
        if system_weights is None:
            return np.ones(n_systems, dtype=float)
        raw = np.maximum(np.asarray(system_weights, dtype=float).reshape(-1), 0.0)
        if raw.size < n_systems:
            raw = np.pad(raw, (0, n_systems - raw.size), mode="edge")
        raw = raw[:n_systems]
        return raw / np.max(raw) if np.max(raw) > 0.0 else np.zeros_like(raw)

    def emit_exact(
        self,
        dt_s: float,
        sigma_opening_Pa: float,
        T_K: float,
        system_weights: np.ndarray | None = None,
    ) -> dict[str, Any]:
        dt = max(float(dt_s), 0.0)
        weights = self._weights(system_weights, self.n_systems)
        sigma_back = self.taylor_backstress_Pa()
        sigma_emit = np.maximum(weights * max(float(sigma_opening_Pa), 0.0) - sigma_back, 0.0)
        rates = np.maximum(self.manifest.emission.rate(sigma_emit, T_K), 0.0)
        available0 = np.maximum(self.available_sites.copy(), 0.0)
        probability = 1.0 - np.exp(-np.minimum(rates * dt, 700.0))
        emitted_system = np.minimum(available0 * probability * weights, available0)
        self.available_sites = np.clip(
            available0 - emitted_system, 0.0, self.site_capacity
        )
        nsrc = max(min(int(self.cfg.source_bin_count), self.n_bins), 1)
        self.mobile[:, :nsrc] += emitted_system[:, None] / nsrc
        self.accumulated_slip[:, :nsrc] += emitted_system[:, None] / nsrc
        emitted = float(np.sum(emitted_system))
        self.emitted_total += emitted
        self.cumulative_emitted += emitted
        self.last_emitted_per_system = emitted_system.copy()
        self.last_emission_rate_per_system_s = rates.copy()
        self.last_sigma_emit_per_system_Pa = sigma_emit.copy()
        return {
            "dN_emit": emitted,
            "dN_emit_per_system": emitted_system.tolist(),
            "lambda_emit_per_system_s-1": rates.tolist(),
            "sigma_emission_backstress_per_system_Pa": sigma_back.tolist(),
            "sigma_emission_effective_per_system_Pa": sigma_emit.tolist(),
            "emission_backstress_density_per_system_m2": (
                self.last_local_density_per_system_m2.tolist()
            ),
        }

    def evolve_campaign(
        self,
        dt_s: float,
        T_K: float,
        sigma_opening_Pa: float,
        b: float,
        system_weights: np.ndarray | None = None,
    ) -> dict[str, Any]:
        dt = max(float(dt_s), 0.0)
        emitted = self.emit_exact(dt, sigma_opening_Pa, T_K, system_weights)
        transport = super().evolve(
            dt,
            T_K,
            max(float(sigma_opening_Pa), 0.0),
            b,
            emission_hazard_integral=0.0,
            system_weights=system_weights,
        )
        transport["dN_emit"] = float(emitted["dN_emit"])
        transport.update({k: v for k, v in emitted.items() if k != "dN_emit"})
        trapped = max(float(transport.get("dN_trapped", 0.0)), 0.0)
        released = max(
            float(transport.get("dN_released", transport.get("dN_detrapped", 0.0))),
            0.0,
        )
        recovered = max(float(transport.get("dN_recovered", 0.0)), 0.0)
        escaped = max(float(transport.get("dN_escaped", 0.0)), 0.0)
        self.cumulative_trapped += trapped
        self.cumulative_released += released
        self.cumulative_recovered += recovered
        self.cumulative_escaped += escaped
        self.mobile_residence_count_s += self.mobile_count * dt
        self.retained_residence_count_s += self.retained_count * dt
        self.active_residence_count_s += self.active_count * dt
        transport["dN_released"] = released
        return transport

    def advance_campaign(self, distance_m: float) -> dict[str, float]:
        before = self.available_sites.copy()
        result = super().advance(distance_m)
        self.available_sites = before.copy()
        d = max(float(distance_m), 0.0)
        Lref = max(
            float(self.manifest.source_refresh_length_m)
            * float(self._kinetic_cfg.source_refresh_scale),
            float(self.dx),
            1.0e-12,
        )
        fraction = 1.0 - math.exp(-min(d / Lref, 700.0))
        increment = (self.site_capacity - self.available_sites) * fraction
        self.available_sites = np.clip(
            self.available_sites + increment, 0.0, self.site_capacity
        )
        refreshed = float(np.sum(increment))
        self.cumulative_refreshed += refreshed
        self.last_source_refresh_fraction = fraction
        self.last_source_refreshed = refreshed
        result["source_sites_refreshed"] = refreshed
        result["source_refresh_fraction"] = fraction
        result["source_refresh_length_m"] = Lref
        return {k: float(v) for k, v in result.items()}

    def diagnostics_campaign(self) -> dict[str, Any]:
        mobile = max(float(self.mobile_count), 0.0)
        retained = max(float(self.retained_count), 0.0)
        active = mobile + retained
        return {
            "source_budget_total": float(np.sum(self.site_capacity)),
            "source_budget_remaining": float(np.sum(self.available_sites)),
            "source_budget_consumed": float(np.sum(self.site_capacity - self.available_sites)),
            "source_sites_refreshed_step": float(self.last_source_refreshed),
            "source_refresh_fraction": float(self.last_source_refresh_fraction),
            "mobile_count": mobile,
            "retained_count": retained,
            "active_count": active,
            "retained_fraction": retained / active if active > 0.0 else 0.0,
            "cumulative_emitted": float(self.cumulative_emitted),
            "cumulative_refreshed": float(self.cumulative_refreshed),
            "cumulative_trapped": float(self.cumulative_trapped),
            "cumulative_released": float(self.cumulative_released),
            "cumulative_recovered": float(self.cumulative_recovered),
            "cumulative_escaped": float(self.cumulative_escaped),
            "mobile_residence_count_s": float(self.mobile_residence_count_s),
            "retained_residence_count_s": float(self.retained_residence_count_s),
            "active_residence_count_s": float(self.active_residence_count_s),
            "emission_backstress_density_m2": float(
                np.mean(self.last_local_density_per_system_m2)
            ),
        }


class CampaignCalibratedCZMFrontEngine(MovingProcessZone2DFrontEngine):
    """PF-equivalent front engine with separated stress channels and moving B."""

    state_model = STATE_MODEL
    state_model_detail = "pf_v10_1_7_1_campaign_calibrated_continuous_tip"

    def __init__(
        self,
        fcfg: Any,
        cleave_barrier: Any,
        emit_barrier: Any,
        G_shear: float,
        nu: float,
        b: float,
        mpz_config: Any,
        manifest: PFEquivalentMaterialManifest,
        kinetic_config: KineticCampaignCZMConfig | None = None,
    ) -> None:
        self.manifest = manifest
        self.kinetic_config = copy.deepcopy(
            kinetic_config or KineticCampaignCZMConfig()
        ).validate()
        apply_pf_manifest_to_mpz_config(mpz_config, manifest, self.kinetic_config)
        super().__init__(
            fcfg, cleave_barrier, emit_barrier, G_shear, nu, b, mpz_config
        )
        self.mpz_state = CampaignKineticMPZState(
            mpz_config,
            manifest,
            b=b,
            G_Pa=G_shear,
            kinetic_cfg=self.kinetic_config,
        )
        self.f.c_blunt = float(manifest.c_blunt)
        self.f.sigma_cap = 0.0
        self.f.dN_cap = math.inf
        self.f.N_sat = math.inf
        self.f.recover_k = 0.0
        self.f.k_shield = 0.0
        self.f.chi_shield = 0.0
        self.f.v_emb_b3 = 0.0
        self.micro_advance_total_m = 0.0
        self.checkpoint_advance_total_m = 0.0
        self._last_pre_checkpoint_snapshot: dict[str, Any] | None = None
        self._last_channels: dict[str, Any] = {}
        self._sync_compat()

    def snapshot_kinetic_state(self) -> dict[str, Any]:
        return {
            "mpz_state": self.mpz_state.copy(),
            "B": float(self.B),
            "a_adv": float(self.a_adv),
            "n_adv": int(self.n_adv),
            "W_emit": float(self.W_emit),
            "t": float(self.t),
            "K_prev": self.K_prev,
            "lambda_prev": self._lambda_c_prev,
            "K_cleave_prev": self._K_cleave_prev,
            "micro_advance_total_m": float(self.micro_advance_total_m),
            "checkpoint_advance_total_m": float(self.checkpoint_advance_total_m),
            "last_channels": copy.deepcopy(self._last_channels),
        }

    def restore_kinetic_state(self, payload: Mapping[str, Any]) -> None:
        self.mpz_state = copy.deepcopy(payload["mpz_state"])
        self.B = float(payload["B"])
        self.a_adv = float(payload["a_adv"])
        self.n_adv = int(payload["n_adv"])
        self.W_emit = float(payload["W_emit"])
        self.t = float(payload["t"])
        self.K_prev = payload.get("K_prev")
        self._lambda_c_prev = payload.get("lambda_prev")
        self._K_cleave_prev = payload.get("K_cleave_prev")
        self.micro_advance_total_m = float(payload["micro_advance_total_m"])
        self.checkpoint_advance_total_m = float(payload["checkpoint_advance_total_m"])
        self._last_channels = copy.deepcopy(payload.get("last_channels", {}))
        self._sync_compat()

    def _active_shielding_raw(self) -> float:
        state = self.mpz_state
        core = max(float(state.cfg.shielding_core_m), 0.25 * abs(self.b), 1.0e-12)
        kernel = (
            self.G
            * self.b
            / max(1.0 - self.nu, 1.0e-6)
            / np.sqrt(2.0 * np.pi * np.maximum(state.x, core))
        )
        signed = state.retained + float(state.cfg.mobile_shield_fraction) * state.mobile
        return float(
            np.sum(state.orientation_factors[:, None] * signed * kernel[None, :])
        )

    def active_K_shielding(self) -> float:
        if not self.kinetic_config.active_shielding:
            return 0.0
        raw = self._active_shielding_raw()
        if not self.kinetic_config.signed_active_shielding:
            raw = max(raw, 0.0)
        cap = max(float(self.manifest.max_K_shield_MPa_sqrt_m), 0.0) * 1.0e6
        return float(np.clip(raw, -cap, cap)) if cap > 0.0 else float(raw)

    def K_shield(self) -> float:
        return self.active_K_shielding()

    def r_eff(self) -> float:
        return self.mpz_state.blunted_radius(
            self.f.r0, float(self.manifest.c_blunt), self.b
        )

    def sigma_opening_tip(self, K_open: float) -> float:
        return float(
            max(float(K_open), 0.0)
            / math.sqrt(2.0 * math.pi * max(float(self.r_eff()), 1.0e-30))
        )

    def sigma_cleavage_tip(self, K_cleave: float) -> float:
        K_eff = max(float(K_cleave) - self.active_K_shielding(), 0.0)
        return float(
            K_eff / math.sqrt(2.0 * math.pi * max(float(self.r_eff()), 1.0e-30))
        )

    def stress_channels(
        self,
        K_open: float,
        K_cleave: float,
        system_weights: np.ndarray | None = None,
    ) -> dict[str, Any]:
        opening = self.sigma_opening_tip(K_open)
        cleavage = self.sigma_cleavage_tip(K_cleave)
        weights = CampaignKineticMPZState._weights(
            system_weights, self.mpz_state.n_systems
        )
        back = self.mpz_state.taylor_backstress_Pa()
        emission = np.maximum(weights * opening - back, 0.0)
        out = {
            "K_open_Pa_sqrt_m": float(K_open),
            "K_cleave_input_Pa_sqrt_m": float(K_cleave),
            "K_shield_raw_Pa_sqrt_m": float(self._active_shielding_raw()),
            "K_shield_effective_Pa_sqrt_m": float(self.active_K_shielding()),
            "sigma_opening_tip_Pa": opening,
            "sigma_cleave_eff_Pa": cleavage,
            "sigma_emission_backstress_Pa": float(np.mean(back)),
            "sigma_emission_effective_Pa": float(np.mean(emission)),
            "sigma_emission_backstress_per_system_Pa": back.tolist(),
            "sigma_emission_effective_per_system_Pa": emission.tolist(),
            "slip_system_directional_weights": weights.tolist(),
        }
        self._last_channels = copy.deepcopy(out)
        return out

    def _cleavage_rate(self, sigma_cleave: float, T_K: float) -> tuple[float, float, float]:
        Gstar = float(self.cb.G_barrier(np.array([max(sigma_cleave, 0.0)]), T_K, self.b)[0])
        raw = self.f.nu0_c * math.exp(
            float(np.clip(-Gstar / max(KB * T_K, 1.0e-30), -700.0, 0.0))
        )
        m = max(float(self.f.m_hits), 1.0)
        if m > 1.0 + 1.0e-12:
            from scipy.special import gammainc

            tau = max(float(self.f.tau_c), 1.0e-30)
            effective = float(gammainc(m, min(raw * tau, 1.0e12)) / tau)
        else:
            effective = float(raw)
        return effective, float(raw), Gstar

    @staticmethod
    def _sum_numeric(target: dict[str, float], source: Mapping[str, Any]) -> None:
        for key, value in source.items():
            if isinstance(value, (bool, np.bool_)):
                continue
            if isinstance(value, (int, float, np.integer, np.floating)):
                target[key] = target.get(key, 0.0) + float(value)

    def _plastic_half_step(
        self,
        dt_s: float,
        T_K: float,
        sigma_opening: float,
        system_weights: np.ndarray | None,
    ) -> dict[str, Any]:
        out = self.mpz_state.evolve_campaign(
            dt_s,
            T_K,
            sigma_opening,
            self.b,
            system_weights,
        )
        self.W_emit += (
            max(float(sigma_opening), 0.0)
            * self.b
            * self.f.L_pz
            * max(float(out.get("dN_emit", 0.0)), 0.0)
        )
        self._sync_compat()
        return out

    def _substep_limit(self, remaining: float, lam: float) -> float:
        h = max(float(remaining), 0.0)
        if lam > 0.0 and math.isfinite(lam):
            h = min(h, float(self.kinetic_config.max_action_substep) / lam)
            velocity = max(float(self.f.da) * lam, 1.0e-300)
            h = min(
                h,
                float(self.kinetic_config.max_translation_substep_m) / velocity,
            )
            remaining_action = max(1.0 - float(self.B), 0.0)
            if remaining_action > 0.0:
                h = min(h, remaining_action / lam)
        minimum = min(float(self.kinetic_config.min_substep_s), remaining)
        return max(min(h, remaining), minimum)

    def integrate_kinetics(
        self,
        K_open: float,
        K_cleave: float,
        T_K: float,
        dt_s: float,
        *,
        system_weights: np.ndarray | None = None,
    ) -> dict[str, Any]:
        requested = max(float(dt_s), 0.0)
        remaining = requested
        consumed = 0.0
        dB_total = 0.0
        da_total = 0.0
        totals: dict[str, float] = {}
        advance_totals: dict[str, float] = {}
        fired = False
        internal = 0
        last_lam = 0.0
        last_raw = 0.0
        last_G = 0.0
        last_channels = self.stress_channels(K_open, K_cleave, system_weights)

        while remaining > 0.0:
            internal += 1
            if internal > int(self.kinetic_config.max_internal_steps):
                raise RuntimeError(
                    "kinetic_campaign_czm exceeded max_internal_steps; reduce the outer dt"
                )
            channels0 = self.stress_channels(K_open, K_cleave, system_weights)
            lam0, raw0, G0 = self._cleavage_rate(
                float(channels0["sigma_cleave_eff_Pa"]), T_K
            )
            h = self._substep_limit(remaining, lam0)
            snapshot = self.snapshot_kinetic_state()
            first = self._plastic_half_step(
                0.5 * h,
                T_K,
                float(channels0["sigma_opening_tip_Pa"]),
                system_weights,
            )
            channels_mid = self.stress_channels(K_open, K_cleave, system_weights)
            lam_mid, raw_mid, G_mid = self._cleavage_rate(
                float(channels_mid["sigma_cleave_eff_Pa"]), T_K
            )
            remaining_action = max(1.0 - float(self.B), 0.0)
            if lam_mid > 0.0 and lam_mid * h > remaining_action + 1.0e-12:
                self.restore_kinetic_state(snapshot)
                h = min(max(remaining_action / lam_mid, self.kinetic_config.min_substep_s), remaining)
                channels0 = self.stress_channels(K_open, K_cleave, system_weights)
                first = self._plastic_half_step(
                    0.5 * h,
                    T_K,
                    float(channels0["sigma_opening_tip_Pa"]),
                    system_weights,
                )
                channels_mid = self.stress_channels(K_open, K_cleave, system_weights)
                lam_mid, raw_mid, G_mid = self._cleavage_rate(
                    float(channels_mid["sigma_cleave_eff_Pa"]), T_K
                )

            dB = min(lam_mid * h, max(1.0 - float(self.B), 0.0))
            da = float(self.f.da) * dB
            advance = self.mpz_state.advance_campaign(da) if da > 0.0 else {}
            channels1 = self.stress_channels(K_open, K_cleave, system_weights)
            second = self._plastic_half_step(
                0.5 * h,
                T_K,
                float(channels1["sigma_opening_tip_Pa"]),
                system_weights,
            )

            self._sum_numeric(totals, first)
            self._sum_numeric(totals, second)
            self._sum_numeric(advance_totals, advance)
            self.B += dB
            self.micro_advance_total_m += da
            dB_total += dB
            da_total += da
            consumed += h
            remaining = max(remaining - h, 0.0)
            self.t += h
            last_lam, last_raw, last_G = lam_mid, raw_mid, G_mid
            last_channels = channels1

            if self.B >= 1.0 - 1.0e-10:
                self._last_pre_checkpoint_snapshot = snapshot
                self.B = max(self.B - 1.0, 0.0)
                self.a_adv += float(self.f.da)
                self.checkpoint_advance_total_m += float(self.f.da)
                self.n_adv += 1
                fired = True
                break
            if h <= 0.0:
                break

        return {
            "fired": fired,
            "n_fire": 1 if fired else 0,
            "dB": dB_total,
            "micro_advance_step_m": da_total,
            "micro_advance_total_m": float(self.micro_advance_total_m),
            "checkpoint_committed_total_m": float(self.checkpoint_advance_total_m),
            "dt_consumed_s": consumed,
            "dt_unused_s": max(requested - consumed, 0.0),
            "internal_substeps": internal,
            "lambda_c_effective_s-1": last_lam,
            "lambda_c_raw_s-1": last_raw,
            "G_cleave_eff_eV": last_G / EV_TO_J,
            "plastic": totals,
            "advance": advance_totals,
            "channels": last_channels,
        }

    def predict_clock_increment_drives(self, K_cleave, K_emit, T, dt):
        channels = self.stress_channels(float(K_emit), float(K_cleave), None)
        lam, _, _ = self._cleavage_rate(
            float(channels["sigma_cleave_eff_Pa"]), float(T)
        )
        return max(float(lam) * max(float(dt), 0.0), 0.0)

    def step_drives(self, K_cleave, K_emit, T, dt, metadata=None):
        md = dict(metadata or {})
        K_open = float(md.get("K_open_Pa_sqrt_m", K_emit))
        weights = md.get("slip_system_weights")
        if weights is not None:
            weights = np.asarray(weights, dtype=float)
        pre = self.snapshot_kinetic_state()
        result = self.integrate_kinetics(
            K_open,
            float(K_cleave),
            float(T),
            float(dt),
            system_weights=weights,
        )
        if result["fired"]:
            self._last_pre_checkpoint_snapshot = pre
        channels = dict(result["channels"])
        plastic = dict(result["plastic"])
        advance = dict(result["advance"])
        state_diag = self.mpz_state.diagnostics_campaign()
        generic_diag = self.mpz_state.diagnostics(
            self.G, self.nu, self.b, self.f.r0, float(self.manifest.c_blunt)
        )
        out = {
            "fired": bool(result["fired"]),
            "n_fire": int(result["n_fire"]),
            "n_fire_available": int(result["n_fire"]),
            "v_crack": (
                float(result["micro_advance_step_m"]) / float(result["dt_consumed_s"])
                if float(result["dt_consumed_s"]) > 0.0
                else 0.0
            ),
            "B": float(self.B),
            "cleavage_clock_B": float(self.B),
            "N_em": float(self.N_em),
            "r_eff": float(self.r_eff()),
            "r_eff_m": float(self.r_eff()),
            "local_slip_count": float(self.mpz_state.local_slip_count()),
            "lambda_c": float(result["lambda_c_effective_s-1"]),
            "lambda_c_raw": float(result["lambda_c_raw_s-1"]),
            "lambda_c_effective_s-1": float(result["lambda_c_effective_s-1"]),
            "lambda_c_raw_s-1": float(result["lambda_c_raw_s-1"]),
            "G_cleave_eff_eV": float(result["G_cleave_eff_eV"]),
            "W_emit": float(self.W_emit),
            "kinetic_campaign_czm_active": True,
            "front_state_model": STATE_MODEL,
            "material_parameter_source": self.manifest.parameter_source,
            "material_class": self.manifest.name,
            "candidate_id": self.manifest.candidate_id,
            "micro_advance_step_m": float(result["micro_advance_step_m"]),
            "micro_advance_total_m": float(result["micro_advance_total_m"]),
            "checkpoint_committed_total_m": float(result["checkpoint_committed_total_m"]),
            "dt_consumed_s": float(result["dt_consumed_s"]),
            "dt_unused_s": float(result["dt_unused_s"]),
            "internal_substeps": int(result["internal_substeps"]),
            "wake_shielding_mechanically_active": False,
            "stored_energy_cleavage_active": False,
            "temporal_source_recycling_active": False,
            "per_step_emission_cap_active": False,
            "temperature_dependent_source_capacity_active": False,
            **channels,
            **plastic,
            **advance,
            **state_diag,
            **generic_diag,
            **md,
        }
        return out

    def step(self, K, T, dt):
        return self.step_drives(K, K, T, dt)

    def restore_geometry_veto(self, n_restore: int) -> None:
        if self._last_pre_checkpoint_snapshot is not None:
            self.restore_kinetic_state(self._last_pre_checkpoint_snapshot)
            self._last_pre_checkpoint_snapshot = None
            return
        super().restore_geometry_veto(n_restore)

    def audit_payload(self) -> dict[str, Any]:
        return {
            "model": MODEL_ID,
            "front_state_model": STATE_MODEL,
            "kinetic_config": asdict(self.kinetic_config),
            "material": self.manifest.as_dict(),
            "stress_channels_separated": True,
            "opening_stress_unshielded": True,
            "cleavage_uses_active_elastic_shielding_only": True,
            "emission_uses_local_taylor_backstress_only": True,
            "continuous_mpz_translation": True,
            "source_refresh_from_advance_only": True,
            "wake_shielding_active": False,
            "stored_energy_cleavage_active": False,
        }


class DevelopedStateDiagnosticCZMFrontEngine(CampaignCalibratedCZMFrontEngine):
    """Semantic alias: diagnostics are integrated by CampaignKineticMPZState."""

    developed_state_diagnostics_active = True


__all__ = [
    "MODEL_ID",
    "STATE_MODEL",
    "KineticCampaignCZMConfig",
    "apply_pf_manifest_to_mpz_config",
    "CampaignKineticMPZState",
    "CampaignCalibratedCZMFrontEngine",
    "DevelopedStateDiagnosticCZMFrontEngine",
]
