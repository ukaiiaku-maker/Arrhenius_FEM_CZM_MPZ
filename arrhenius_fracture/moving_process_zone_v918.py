"""Persistent signed plastic wake for the v9.18 moving process zone.

The active MPZ remains ahead of the current crack tip.  State crossed by a
physical crack increment is conservatively remapped into a second finite-volume
field behind the tip instead of being reduced to uncoupled cumulative scalars.

Only the signed retained/mobile line field contributes to the Mode-I-equivalent
wake shielding.  Wake slip is conserved and reported, but is not interpreted as
bridging, transformation toughening, or an extra blunting law.
"""
from __future__ import annotations

import copy
import math
import os
from typing import Any

import numpy as np

from .moving_process_zone_v911 import MovingProcessZoneState as _V911State


class MovingProcessZoneState(_V911State):
    """v9.11 MPZ plus a conservative spatial state behind the crack tip."""

    state_model = "moving_pz_v918_persistent_signed_wake"

    def __init__(self, cfg):
        super().__init__(cfg)
        default_length_um = self.length_m * 1.0e6
        wake_length_um = float(
            os.environ.get("ARRHENIUS_WAKE_LENGTH_UM", f"{default_length_um:.16g}")
        )
        self.wake_length_m = max(wake_length_um * 1.0e-6, self.dx)
        requested_bins = int(os.environ.get("ARRHENIUS_WAKE_N_BINS", "0") or 0)
        if requested_bins <= 0:
            requested_bins = max(int(round(self.wake_length_m / self.dx)), 1)
        self.wake_n_bins = requested_bins
        self.wake_dx = self.wake_length_m / self.wake_n_bins
        self.wake_x = (np.arange(self.wake_n_bins, dtype=float) + 0.5) * self.wake_dx

        shape = (self.n_systems, self.wake_n_bins)
        self.wake_mobile = np.zeros(shape, dtype=float)
        self.wake_retained = np.zeros(shape, dtype=float)
        self.wake_slip = np.zeros(shape, dtype=float)

        self.wake_discarded_mobile_total = 0.0
        self.wake_discarded_retained_total = 0.0
        self.wake_discarded_slip_total = 0.0
        self.wake_recovered_total = 0.0
        self.wake_trapped_total = 0.0
        self.wake_released_total = 0.0
        self.wake_mobile_transport_loss_total = 0.0

        self.wake_shield_projection = float(
            os.environ.get("ARRHENIUS_WAKE_SHIELD_PROJECTION", "1")
        )
        self.wake_shielding_enabled = (
            os.environ.get("ARRHENIUS_WAKE_SHIELDING", "1").strip().lower()
            not in {"0", "false", "off", "no"}
        )

    @property
    def wake_mobile_count(self) -> float:
        return float(np.sum(self.wake_mobile))

    @property
    def wake_retained_count(self) -> float:
        return float(np.sum(self.wake_retained))

    @property
    def wake_slip_count(self) -> float:
        return float(np.sum(self.wake_slip))

    def _shielding_raw(
        self,
        retained: np.ndarray,
        mobile: np.ndarray,
        x: np.ndarray,
        G_shear: float,
        nu: float,
        b: float,
    ) -> float:
        core = max(float(self.cfg.shielding_core_m), 0.25 * abs(float(b)), 1.0e-12)
        kernel = (
            float(G_shear)
            * float(b)
            / max(1.0 - float(nu), 1.0e-6)
            / np.sqrt(2.0 * np.pi * np.maximum(np.asarray(x, float), core))
        )
        signed = np.asarray(retained, float) + float(
            self.cfg.mobile_shield_fraction
        ) * np.asarray(mobile, float)
        return float(
            np.sum(self.orientation_factors[:, None] * signed * kernel[None, :])
        )

    def active_K_shielding(self, G_shear: float, nu: float, b: float) -> float:
        return float(max(self._shielding_raw(
            self.retained, self.mobile, self.x, G_shear, nu, b
        ), 0.0))

    def wake_K_shielding(self, G_shear: float, nu: float, b: float) -> float:
        if not self.wake_shielding_enabled:
            return 0.0
        raw = self._shielding_raw(
            self.wake_retained,
            self.wake_mobile,
            self.wake_x,
            G_shear,
            nu,
            b,
        )
        return float(max(self.wake_shield_projection * raw, 0.0))

    def shielding_K(self, G_shear: float, nu: float, b: float) -> float:
        return float(
            self.active_K_shielding(G_shear, nu, b)
            + self.wake_K_shielding(G_shear, nu, b)
        )

    @staticmethod
    def _deposit_interval(
        out: np.ndarray,
        left: float,
        right: float,
        mass: np.ndarray,
        source_width: float,
        target_dx: float,
        target_length: float,
    ) -> float:
        left = max(float(left), 0.0)
        right = min(float(right), float(target_length))
        if right <= left or source_width <= 0.0:
            return 0.0
        j0 = max(int(math.floor(left / target_dx)), 0)
        j1 = min(
            int(math.floor((right - 1.0e-15 * target_dx) / target_dx)),
            out.shape[1] - 1,
        )
        deposited = 0.0
        for j in range(j0, j1 + 1):
            overlap = max(
                min(right, (j + 1) * target_dx) - max(left, j * target_dx),
                0.0,
            )
            if overlap <= 0.0:
                continue
            fraction = overlap / source_width
            out[:, j] += mass * fraction
            deposited += float(np.sum(mass) * fraction)
        return deposited

    def _advance_active_field(
        self, field: np.ndarray, distance_m: float
    ) -> tuple[np.ndarray, np.ndarray, float, float]:
        d = max(float(distance_m), 0.0)
        active = np.zeros_like(field)
        crossed = np.zeros((self.n_systems, self.wake_n_bins), dtype=float)
        deposited_active = 0.0
        deposited_wake = 0.0
        total = float(np.sum(field))
        for i in range(self.n_bins):
            mass = np.asarray(field[:, i], float)
            left = i * self.dx - d
            right = (i + 1) * self.dx - d
            if right > 0.0:
                deposited_active += self._deposit_interval(
                    active,
                    max(left, 0.0),
                    right,
                    mass,
                    self.dx,
                    self.dx,
                    self.length_m,
                )
            if left < 0.0:
                # x<0 becomes a positive distance behind the new tip.
                y_left = max(-right, 0.0)
                y_right = -left
                deposited_wake += self._deposit_interval(
                    crossed,
                    y_left,
                    y_right,
                    mass,
                    self.dx,
                    self.wake_dx,
                    self.wake_length_m,
                )
        discarded = max(total - deposited_active - deposited_wake, 0.0)
        return active, crossed, discarded, deposited_wake

    def _shift_wake_field(
        self, field: np.ndarray, distance_m: float
    ) -> tuple[np.ndarray, float]:
        d = max(float(distance_m), 0.0)
        shifted = np.zeros_like(field)
        deposited = 0.0
        total = float(np.sum(field))
        for i in range(self.wake_n_bins):
            mass = np.asarray(field[:, i], float)
            left = i * self.wake_dx + d
            right = (i + 1) * self.wake_dx + d
            deposited += self._deposit_interval(
                shifted,
                left,
                right,
                mass,
                self.wake_dx,
                self.wake_dx,
                self.wake_length_m,
            )
        return shifted, max(total - deposited, 0.0)

    def advance(self, distance_m: float) -> dict[str, float]:
        d = max(float(distance_m), 0.0)
        old_wake_mobile, lost_old_wm = self._shift_wake_field(self.wake_mobile, d)
        old_wake_retained, lost_old_wr = self._shift_wake_field(self.wake_retained, d)
        old_wake_slip, lost_old_ws = self._shift_wake_field(self.wake_slip, d)

        self.mobile, crossed_m, lost_new_m, crossed_m_count = self._advance_active_field(
            self.mobile, d
        )
        self.retained, crossed_r, lost_new_r, crossed_r_count = self._advance_active_field(
            self.retained, d
        )
        self.accumulated_slip, crossed_s, lost_new_s, crossed_s_count = (
            self._advance_active_field(self.accumulated_slip, d)
        )

        self.wake_mobile = np.maximum(old_wake_mobile + crossed_m, 0.0)
        self.wake_retained = np.maximum(old_wake_retained + crossed_r, 0.0)
        self.wake_slip = np.maximum(old_wake_slip + crossed_s, 0.0)

        discarded_m = lost_old_wm + lost_new_m
        discarded_r = lost_old_wr + lost_new_r
        discarded_s = lost_old_ws + lost_new_s
        self.wake_discarded_mobile_total += discarded_m
        self.wake_discarded_retained_total += discarded_r
        self.wake_discarded_slip_total += discarded_s

        self.wake_mobile_total += crossed_m_count
        self.wake_retained_total += crossed_r_count
        self.wake_slip_total += crossed_s_count
        self.advance_total_m += d

        Lrefresh = max(float(self.cfg.source_refresh_length_m), self.dx)
        fresh = min(d / Lrefresh, 1.0)
        source_refreshed = (self.site_capacity - self.available_sites) * fresh
        self.available_sites += source_refreshed

        return {
            "wake_mobile": float(crossed_m_count),
            "wake_retained": float(crossed_r_count),
            "wake_slip": float(crossed_s_count),
            "wake_mobile_discarded": float(discarded_m),
            "wake_retained_discarded": float(discarded_r),
            "wake_slip_discarded": float(discarded_s),
            "wake_mobile_postcommit": self.wake_mobile_count,
            "wake_retained_postcommit": self.wake_retained_count,
            "wake_slip_postcommit": self.wake_slip_count,
            "active_mobile_postcommit": self.mobile_count,
            "active_retained_postcommit": self.retained_count,
            "active_slip_postcommit": float(np.sum(self.accumulated_slip)),
            "source_sites_refreshed": float(np.sum(source_refreshed)),
        }

    def _wake_forest_density(self) -> np.ndarray:
        width = max(float(self.cfg.blunting_length_m), self.wake_dx, 1.0e-12)
        count = np.sum(np.maximum(self.wake_retained, 0.0), axis=0)
        return np.maximum(
            float(self.cfg.pt_forest_density_floor_m2)
            + count / max(self.wake_dx * width, 1.0e-30),
            1.0,
        )

    def _evolve_wake(self, dt_s: float, T_K: float, b: float) -> dict[str, float]:
        dt = max(float(dt_s), 0.0)
        if dt <= 0.0 or (self.wake_mobile_count + self.wake_retained_count) <= 0.0:
            return {
                "wake_dN_trapped": 0.0,
                "wake_dN_released": 0.0,
                "wake_dN_recovered": 0.0,
                "wake_dN_mobile_transport_loss": 0.0,
            }

        rho = self._wake_forest_density()
        model = self._pt_model()
        rates = model.rates(np.zeros(self.wake_n_bins), rho, T_K, b)
        peierls = np.asarray(rates["peierls_rate_s"], float).reshape(-1)
        taylor = np.asarray(rates["taylor_completion_rate_s"], float).reshape(-1)
        jump = np.asarray(rates["jump_length_m"], float).reshape(-1)
        eta = float(getattr(self.cfg, "pt_encounter_efficiency", 1.0))
        encounter = self.encounter_rate_s(peierls, jump, rho, eta)

        self.wake_mobile, self.wake_retained, trapped, released = (
            self._exchange_mobile_retained(
                self.wake_mobile,
                self.wake_retained,
                encounter,
                taylor,
                dt,
            )
        )

        kr = max(float(getattr(self.cfg, "retained_recovery_nu0_s", 0.0)), 0.0)
        km = max(float(getattr(self.cfg, "mobile_recovery_rate_s", 0.0)), 0.0)
        fr = 1.0 - math.exp(-min(kr * dt, 700.0))
        fm = 1.0 - math.exp(-min(km * dt, 700.0))
        rec_r = self.wake_retained * fr
        rec_m = self.wake_mobile * fm
        self.wake_retained -= rec_r
        self.wake_mobile -= rec_m
        recovered = float(np.sum(rec_r) + np.sum(rec_m))

        kpair = max(float(self.cfg.pair_annihilation_rate_per_count_s), 0.0)
        annihilated = 0.0
        if kpair > 0.0 and self.n_systems >= 2:
            for i in range(0, self.n_systems - 1, 2):
                pair = np.minimum(self.wake_retained[i], self.wake_retained[i + 1])
                frac = 1.0 - np.exp(-np.minimum(kpair * pair * dt, 700.0))
                removed = pair * frac
                self.wake_retained[i] -= removed
                self.wake_retained[i + 1] -= removed
                annihilated += 2.0 * float(np.sum(removed))

        mobile_by_bin = np.sum(np.maximum(self.wake_mobile, 0.0), axis=0)
        if float(np.sum(mobile_by_bin)) > 0.0:
            velocity = float(np.sum(jump * peierls * mobile_by_bin) / np.sum(mobile_by_bin))
        else:
            velocity = 0.0
        self.wake_mobile, mobile_lost = self._shift_wake_field(
            self.wake_mobile, max(velocity, 0.0) * dt
        )

        self.wake_mobile = np.maximum(self.wake_mobile, 0.0)
        self.wake_retained = np.maximum(self.wake_retained, 0.0)
        self.wake_trapped_total += float(trapped)
        self.wake_released_total += float(released)
        self.wake_recovered_total += recovered + annihilated
        self.wake_mobile_transport_loss_total += float(mobile_lost)
        self.wake_discarded_mobile_total += float(mobile_lost)
        return {
            "wake_dN_trapped": float(trapped),
            "wake_dN_released": float(released),
            "wake_dN_recovered": float(recovered + annihilated),
            "wake_dN_mobile_transport_loss": float(mobile_lost),
        }

    def evolve(
        self,
        dt_s: float,
        T_K: float,
        stress_Pa: float,
        b: float,
        emission_hazard_integral: float = 0.0,
        system_weights=None,
    ) -> dict[str, float]:
        out = super().evolve(
            dt_s,
            T_K,
            stress_Pa,
            b,
            emission_hazard_integral=emission_hazard_integral,
            system_weights=system_weights,
        )
        out.update(self._evolve_wake(dt_s, T_K, b))
        return out

    def split(self, daughter_fraction: float):
        frac = float(np.clip(daughter_fraction, 0.0, 1.0))
        child = super().split(frac)
        for name in ("wake_mobile", "wake_retained", "wake_slip"):
            original = np.asarray(getattr(self, name), float).copy()
            setattr(child, name, original * frac)
            setattr(self, name, original * (1.0 - frac))
        for name in (
            "wake_discarded_mobile_total",
            "wake_discarded_retained_total",
            "wake_discarded_slip_total",
            "wake_recovered_total",
            "wake_trapped_total",
            "wake_released_total",
            "wake_mobile_transport_loss_total",
        ):
            value = float(getattr(self, name))
            setattr(child, name, value * frac)
            setattr(self, name, value * (1.0 - frac))
        return child

    def state_dict(self) -> dict[str, Any]:
        out = super().state_dict()
        out["persistent_wake_v918"] = {
            "wake_length_m": float(self.wake_length_m),
            "wake_n_bins": int(self.wake_n_bins),
            "wake_mobile": self.wake_mobile.tolist(),
            "wake_retained": self.wake_retained.tolist(),
            "wake_slip": self.wake_slip.tolist(),
            "wake_shield_projection": float(self.wake_shield_projection),
            "wake_shielding_enabled": bool(self.wake_shielding_enabled),
            "wake_discarded_mobile_total": float(self.wake_discarded_mobile_total),
            "wake_discarded_retained_total": float(self.wake_discarded_retained_total),
            "wake_discarded_slip_total": float(self.wake_discarded_slip_total),
            "wake_recovered_total": float(self.wake_recovered_total),
            "wake_trapped_total": float(self.wake_trapped_total),
            "wake_released_total": float(self.wake_released_total),
            "wake_mobile_transport_loss_total": float(self.wake_mobile_transport_loss_total),
        }
        return out

    @classmethod
    def from_state_dict(cls, payload):
        obj = super().from_state_dict(payload)
        wake = dict(payload.get("persistent_wake_v918", {}))
        if wake:
            obj.wake_length_m = float(wake.get("wake_length_m", obj.wake_length_m))
            obj.wake_n_bins = int(wake.get("wake_n_bins", obj.wake_n_bins))
            obj.wake_dx = obj.wake_length_m / max(obj.wake_n_bins, 1)
            obj.wake_x = (np.arange(obj.wake_n_bins, dtype=float) + 0.5) * obj.wake_dx
            shape = (obj.n_systems, obj.wake_n_bins)
            for name in ("wake_mobile", "wake_retained", "wake_slip"):
                arr = np.asarray(wake.get(name, np.zeros(shape)), float)
                if arr.shape != shape:
                    raise ValueError(f"{name} shape {arr.shape} != {shape}")
                setattr(obj, name, np.maximum(arr, 0.0))
            for name in (
                "wake_discarded_mobile_total",
                "wake_discarded_retained_total",
                "wake_discarded_slip_total",
                "wake_recovered_total",
                "wake_trapped_total",
                "wake_released_total",
                "wake_mobile_transport_loss_total",
            ):
                setattr(obj, name, float(wake.get(name, getattr(obj, name))))
            obj.wake_shield_projection = float(
                wake.get("wake_shield_projection", obj.wake_shield_projection)
            )
            obj.wake_shielding_enabled = bool(
                wake.get("wake_shielding_enabled", obj.wake_shielding_enabled)
            )
        return obj

    def diagnostics(self, G_shear, nu, b, r0, c_blunt):
        out = super().diagnostics(G_shear, nu, b, r0, c_blunt)
        active_K = self.active_K_shielding(G_shear, nu, b)
        wake_K = self.wake_K_shielding(G_shear, nu, b)
        out.update({
            "mpz_state_model_v918": self.state_model,
            "mpz_active_K_shield_Pa_sqrt_m": float(active_K),
            "mpz_wake_K_shield_Pa_sqrt_m": float(wake_K),
            "mpz_total_K_shield_Pa_sqrt_m": float(active_K + wake_K),
            "mpz_wake_mobile_count": self.wake_mobile_count,
            "mpz_wake_retained_count": self.wake_retained_count,
            "mpz_wake_slip_count": self.wake_slip_count,
            "mpz_wake_length_m": float(self.wake_length_m),
            "mpz_wake_n_bins": int(self.wake_n_bins),
            "mpz_wake_shielding_enabled": float(self.wake_shielding_enabled),
            "mpz_wake_shield_projection": float(self.wake_shield_projection),
            "mpz_wake_discarded_mobile_total": float(self.wake_discarded_mobile_total),
            "mpz_wake_discarded_retained_total": float(self.wake_discarded_retained_total),
            "mpz_wake_discarded_slip_total": float(self.wake_discarded_slip_total),
            "mpz_wake_recovered_total": float(self.wake_recovered_total),
            "mpz_wake_trapped_total": float(self.wake_trapped_total),
            "mpz_wake_released_total": float(self.wake_released_total),
        })
        return out


__all__ = ["MovingProcessZoneState"]
