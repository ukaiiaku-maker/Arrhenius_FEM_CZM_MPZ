"""Moving one-dimensional crack-tip process-zone state.

This module is the v9 replacement for the legacy scalar ``N_em`` closure.
It is deliberately independent of the FEM/crack-geometry implementation: every
sharp front owns one moving process-zone state, so anisotropy, multifront
branching, crack coalescence, adaptive CZM insertion, cyclic mechanics and
restart logic remain in the production driver.

The active coordinate ``x`` is measured forward from the current crack tip.
When the crack advances, all state is conservatively translated toward the
wake; material that falls behind the new tip is removed from the active zone
and virgin material enters at the far boundary.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
import copy
import math
from typing import Any, Iterable

import numpy as np

from .config import KB, EV_TO_J


@dataclass
class MovingProcessZoneConfig:
    """Material/state parameters for the moving process-zone closure.

    These are loading-protocol independent.  Monotonic, cyclic and dwell
    drivers use the same object and differ only in the imposed K(t).
    """

    length_m: float = 2.0e-6
    n_bins: int = 40
    n_systems: int = 2

    # Finite source-site inventory.  A site emits at most once until recovery.
    source_sites_per_system: float = 200.0
    source_recovery_rate_s: float = 0.0
    source_refresh_length_m: float = 2.5e-7
    source_bin_count: int = 2

    # Signed elastic shielding kernel.  Factors are crystallographic/geometric
    # projections, not fitted back-stress multipliers.  Symmetric systems can
    # both shield Mode I and therefore need not have opposite factors.
    shielding_orientation_factors: tuple[float, ...] = (1.0, 1.0)
    mobile_shield_fraction: float = 0.0
    shielding_core_m: float = 2.5e-10

    # Thermally activated transport and retention.  Stress bias uses tau*V.
    glide_nu0_s: float = 1.0e11
    glide_barrier_eV: float = 0.80
    glide_activation_volume_b3: float = 8.0
    glide_step_m: float = 2.5e-10
    glide_stress_fraction: float = 0.45

    # Emission-derived production transport.  Peierls glide and Taylor release
    # are scaled EXP-floor descendants of the active crack-tip emission barrier.
    # The fixed glide/detrap barriers below remain only as explicit legacy
    # ablations when use_emission_derived_pt is false.
    use_emission_derived_pt: bool = True
    pt_emit_G00_eV: float = 1.94
    pt_emit_gT_eV_per_K: float = 0.003934
    pt_emit_sigc0_Pa: float = 2.298e9
    pt_emit_sT_Pa_per_K: float = -6.564e5
    pt_emit_Tref_K: float = 481.33
    pt_emit_exp_a: float = 0.0845685
    pt_emit_exp_n: float = 1.0
    pt_emit_floor_frac: float = 0.02
    pt_emit_floor_min_eV: float = 1.0e-4
    pt_emit_floor_max_frac: float = 0.95
    pt_peierls_energy_ratio: float = 0.005
    pt_peierls_entropy_ratio: float = 0.005
    pt_peierls_stress_ratio: float = 1.0
    pt_peierls_nu0_s: float = 1.0e12
    pt_taylor_energy_ratio: float = 0.02
    pt_taylor_entropy_ratio: float = 0.02
    pt_taylor_stress_ratio: float = 1.0
    pt_taylor_nu0_s: float = 1.0e11
    pt_taylor_corr_rho_c: float = 1.0e14
    pt_taylor_renewal_time_s: float = 1.0e-9
    pt_taylor_m_exponent: float = 1.0
    pt_taylor_m_scale: float = 1.0
    pt_taylor_m_cap: float = float('inf')
    pt_forest_density_floor_m2: float = 5.0e12
    pt_mobile_fraction: float = 0.01
    pt_mobile_saturation_density_m2: float = 1.0e14
    pt_mobile_density_floor_m2: float = 1.0e6
    pt_jump_fraction: float = 1.0
    pt_jump_length_min_m: float = 2.5e-10
    pt_peierls_stress_fraction: float = 0.5773502691896258
    pt_taylor_stress_fraction: float = 0.5773502691896258
    pt_taylor_phi_max: float = 20.0

    trap_nu0_s: float = 1.0e9
    trap_barrier_eV: float = 0.65
    trap_activation_volume_b3: float = 1.0

    detrap_nu0_s: float = 1.0e10
    detrap_barrier_eV: float = 1.20
    detrap_activation_volume_b3: float = 1.0

    retained_recovery_nu0_s: float = 1.0e9
    retained_recovery_barrier_eV: float = 1.50
    retained_recovery_activation_volume_b3: float = 0.0

    mobile_recovery_rate_s: float = 0.0
    pair_annihilation_rate_per_count_s: float = 0.0

    # Blunting is calculated from accumulated local slip, not retained count.
    blunting_length_m: float = 5.0e-7
    blunting_slip_fraction: float = 1.0

    # Numerical integration only.  This is a CFL control, not a physical cap.
    max_transport_cfl: float = 0.35
    max_transport_substeps: int = 2000

    def normalized_factors(self) -> np.ndarray:
        vals = np.asarray(self.shielding_orientation_factors, dtype=float).reshape(-1)
        if vals.size == 0:
            vals = np.ones(1, dtype=float)
        if vals.size < self.n_systems:
            vals = np.pad(vals, (0, self.n_systems - vals.size), mode="edge")
        return vals[: self.n_systems].copy()

    def as_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["shielding_orientation_factors"] = list(self.shielding_orientation_factors)
        return d


class MovingProcessZoneState:
    """Conservative front-local process-zone inventory."""

    def __init__(self, cfg: MovingProcessZoneConfig):
        self.cfg = copy.deepcopy(cfg)
        self.n_bins = max(int(cfg.n_bins), 4)
        self.n_systems = max(int(cfg.n_systems), 1)
        self.length_m = max(float(cfg.length_m), 1.0e-12)
        self.dx = self.length_m / self.n_bins
        self.x = (np.arange(self.n_bins, dtype=float) + 0.5) * self.dx
        self.orientation_factors = cfg.normalized_factors()

        cap = max(float(cfg.source_sites_per_system), 0.0)
        self.site_capacity = np.full(self.n_systems, cap, dtype=float)
        self.available_sites = self.site_capacity.copy()

        shp = (self.n_systems, self.n_bins)
        self.mobile = np.zeros(shp, dtype=float)
        self.retained = np.zeros(shp, dtype=float)
        self.accumulated_slip = np.zeros(shp, dtype=float)

        self.emitted_total = 0.0
        self.escaped_total = 0.0
        self.recovered_total = 0.0
        self.wake_mobile_total = 0.0
        self.wake_retained_total = 0.0
        self.wake_slip_total = 0.0
        self.time_s = 0.0
        self.advance_total_m = 0.0

    def copy(self) -> "MovingProcessZoneState":
        return copy.deepcopy(self)

    # ------------------------------------------------------------------
    # Conservation, branching and serialization
    # ------------------------------------------------------------------
    def split(self, daughter_fraction: float) -> "MovingProcessZoneState":
        """Conservatively split one tip state into parent and daughter states."""
        frac = float(np.clip(daughter_fraction, 0.0, 1.0))
        child = self.copy()
        for name in ("site_capacity", "available_sites", "mobile", "retained",
                     "accumulated_slip"):
            arr = getattr(self, name)
            setattr(child, name, arr * frac)
            setattr(self, name, arr * (1.0 - frac))
        for name in ("emitted_total", "escaped_total", "recovered_total",
                     "wake_mobile_total", "wake_retained_total", "wake_slip_total"):
            val = float(getattr(self, name))
            setattr(child, name, val * frac)
            setattr(self, name, val * (1.0 - frac))
        child.advance_total_m = 0.0
        return child

    def state_dict(self) -> dict[str, Any]:
        return {
            "schema": "moving_process_zone_v1",
            "config": self.cfg.as_dict(),
            "site_capacity": self.site_capacity.tolist(),
            "available_sites": self.available_sites.tolist(),
            "mobile": self.mobile.tolist(),
            "retained": self.retained.tolist(),
            "accumulated_slip": self.accumulated_slip.tolist(),
            "emitted_total": float(self.emitted_total),
            "escaped_total": float(self.escaped_total),
            "recovered_total": float(self.recovered_total),
            "wake_mobile_total": float(self.wake_mobile_total),
            "wake_retained_total": float(self.wake_retained_total),
            "wake_slip_total": float(self.wake_slip_total),
            "time_s": float(self.time_s),
            "advance_total_m": float(self.advance_total_m),
        }

    @classmethod
    def from_state_dict(cls, payload: dict[str, Any]) -> "MovingProcessZoneState":
        cfgd = dict(payload.get("config", {}))
        if "shielding_orientation_factors" in cfgd:
            cfgd["shielding_orientation_factors"] = tuple(cfgd["shielding_orientation_factors"])
        obj = cls(MovingProcessZoneConfig(**cfgd))
        for name in ("site_capacity", "available_sites", "mobile", "retained",
                     "accumulated_slip"):
            if name in payload:
                setattr(obj, name, np.asarray(payload[name], dtype=float))
        for name in ("emitted_total", "escaped_total", "recovered_total",
                     "wake_mobile_total", "wake_retained_total", "wake_slip_total",
                     "time_s", "advance_total_m"):
            if name in payload:
                setattr(obj, name, float(payload[name]))
        obj._validate_shapes()
        return obj

    def _validate_shapes(self) -> None:
        shp = (self.n_systems, self.n_bins)
        for name in ("mobile", "retained", "accumulated_slip"):
            arr = np.asarray(getattr(self, name), dtype=float)
            if arr.shape != shp:
                raise ValueError(f"{name} shape {arr.shape} != {shp}")
            setattr(self, name, np.maximum(arr, 0.0))
        for name in ("site_capacity", "available_sites"):
            arr = np.asarray(getattr(self, name), dtype=float).reshape(-1)
            if arr.size != self.n_systems:
                raise ValueError(f"{name} length {arr.size} != {self.n_systems}")
            setattr(self, name, np.maximum(arr, 0.0))
        self.available_sites = np.minimum(self.available_sites, self.site_capacity)

    # ------------------------------------------------------------------
    # Derived physical state
    # ------------------------------------------------------------------
    @property
    def retained_count(self) -> float:
        return float(np.sum(self.retained))

    @property
    def mobile_count(self) -> float:
        return float(np.sum(self.mobile))

    @property
    def active_count(self) -> float:
        return self.retained_count + self.mobile_count

    @property
    def available_site_fraction(self) -> float:
        denom = float(np.sum(self.site_capacity))
        return float(np.sum(self.available_sites) / denom) if denom > 0 else 0.0

    def shielding_K(self, G_shear: float, nu: float, b: float) -> float:
        """Mode-I-equivalent K shielding from the signed defect distribution."""
        core = max(float(self.cfg.shielding_core_m), 0.25 * abs(float(b)), 1.0e-12)
        kernel = (float(G_shear) * float(b) /
                  max(1.0 - float(nu), 1.0e-6) /
                  np.sqrt(2.0 * np.pi * np.maximum(self.x, core)))
        signed = self.retained + float(self.cfg.mobile_shield_fraction) * self.mobile
        K = np.sum(self.orientation_factors[:, None] * signed * kernel[None, :])
        return float(max(K, 0.0))

    def local_slip_count(self) -> float:
        Lb = max(float(self.cfg.blunting_length_m), self.dx)
        w = np.exp(-self.x / Lb)
        return float(np.sum(self.accumulated_slip * w[None, :]))

    def blunted_radius(self, r0: float, c_blunt: float, b: float) -> float:
        slip = self.local_slip_count() * max(float(self.cfg.blunting_slip_fraction), 0.0)
        return float(max(float(r0) + max(float(c_blunt), 0.0) * abs(float(b)) * slip,
                         float(r0)))

    # ------------------------------------------------------------------
    # Kinetics
    # ------------------------------------------------------------------
    @staticmethod
    def _activated_rate(nu0: float, barrier_eV: float, activation_volume_b3: float,
                        stress_Pa: float, T_K: float, b: float) -> float:
        H = max(float(barrier_eV), 0.0) * EV_TO_J
        V = max(float(activation_volume_b3), 0.0) * abs(float(b)) ** 3
        G = max(H - max(float(stress_Pa), 0.0) * V, 0.0)
        x = -G / max(KB * float(T_K), 1.0e-30)
        return float(max(float(nu0), 0.0) * math.exp(float(np.clip(x, -700.0, 0.0))))

    def _source_commit_from_hazard(self, hazard_integral: float,
                                   system_weights: np.ndarray | None = None) -> np.ndarray:
        H = max(float(hazard_integral), 0.0)
        p = 1.0 - math.exp(-min(H, 700.0))
        if system_weights is None:
            system_weights = np.ones(self.n_systems, dtype=float)
        w = np.maximum(np.asarray(system_weights, dtype=float).reshape(-1), 0.0)
        if w.size < self.n_systems:
            w = np.pad(w, (0, self.n_systems - w.size), mode="edge")
        w = w[: self.n_systems]
        if np.sum(w) <= 0.0:
            w[:] = 1.0
        # All source systems experience the same per-site hazard but directional
        # activity can weight how much of the site inventory is active.
        active_fraction = w / np.max(w)
        emitted = self.available_sites * active_fraction * p
        self.available_sites = np.maximum(self.available_sites - emitted, 0.0)
        nsrc = max(min(int(self.cfg.source_bin_count), self.n_bins), 1)
        self.mobile[:, :nsrc] += emitted[:, None] / nsrc
        self.accumulated_slip[:, :nsrc] += emitted[:, None] / nsrc
        self.emitted_total += float(np.sum(emitted))
        return emitted

    def _recover_source_sites(self, dt_s: float) -> float:
        k = max(float(self.cfg.source_recovery_rate_s), 0.0)
        if k <= 0.0 or dt_s <= 0.0:
            return 0.0
        frac = 1.0 - math.exp(-min(k * dt_s, 700.0))
        inc = (self.site_capacity - self.available_sites) * frac
        self.available_sites += inc
        return float(np.sum(inc))

    def evolve(self, dt_s: float, T_K: float, stress_Pa: float, b: float,
               emission_hazard_integral: float = 0.0,
               system_weights: np.ndarray | None = None) -> dict[str, float]:
        """Advance source, transport, trapping, recovery and annihilation.

        ``emission_hazard_integral`` is the integral of the per-site emission
        hazard over the accepted loading interval.  It may come from a monotonic
        step, a dwell, or phase quadrature over many fatigue cycles.
        """
        dt_s = max(float(dt_s), 0.0)
        emitted = self._source_commit_from_hazard(emission_hazard_integral, system_weights)
        source_recovered = self._recover_source_sites(dt_s)

        tau = max(float(stress_Pa), 0.0)
        pt_diag: dict[str, float] = {}
        if bool(getattr(self.cfg, 'use_emission_derived_pt', True)):
            from .emission_derived_plasticity import (
                CorrelatedTaylorConfig,
                EmissionDerivedPeierlsTaylorConfig,
                EmissionDerivedPeierlsTaylorModel,
                ExpFloorSurface,
                MechanismScale,
            )
            pt_cfg = EmissionDerivedPeierlsTaylorConfig(
                parent=ExpFloorSurface(
                    G00_eV=self.cfg.pt_emit_G00_eV,
                    gT_eV_per_K=self.cfg.pt_emit_gT_eV_per_K,
                    sigc0_Pa=self.cfg.pt_emit_sigc0_Pa,
                    sT_Pa_per_K=self.cfg.pt_emit_sT_Pa_per_K,
                    Tref_K=self.cfg.pt_emit_Tref_K,
                    a=self.cfg.pt_emit_exp_a,
                    n=self.cfg.pt_emit_exp_n,
                    floor_fraction=self.cfg.pt_emit_floor_frac,
                    floor_min_eV=self.cfg.pt_emit_floor_min_eV,
                    floor_max_fraction=self.cfg.pt_emit_floor_max_frac,
                ),
                peierls=MechanismScale(
                    self.cfg.pt_peierls_energy_ratio,
                    self.cfg.pt_peierls_entropy_ratio,
                    self.cfg.pt_peierls_stress_ratio,
                    self.cfg.pt_peierls_nu0_s,
                ),
                taylor=MechanismScale(
                    self.cfg.pt_taylor_energy_ratio,
                    self.cfg.pt_taylor_entropy_ratio,
                    self.cfg.pt_taylor_stress_ratio,
                    self.cfg.pt_taylor_nu0_s,
                ),
                correlated_taylor=CorrelatedTaylorConfig(
                    self.cfg.pt_taylor_corr_rho_c,
                    self.cfg.pt_taylor_renewal_time_s,
                    self.cfg.pt_taylor_m_exponent,
                    self.cfg.pt_taylor_m_scale,
                    self.cfg.pt_taylor_m_cap,
                ),
                peierls_stress_fraction=self.cfg.pt_peierls_stress_fraction,
                taylor_stress_fraction=self.cfg.pt_taylor_stress_fraction,
                taylor_phi_max=self.cfg.pt_taylor_phi_max,
                mobile_fraction_low_density=self.cfg.pt_mobile_fraction,
                mobile_saturation_density_m2=self.cfg.pt_mobile_saturation_density_m2,
                mobile_density_floor_m2=self.cfg.pt_mobile_density_floor_m2,
                jump_fraction_of_forest_spacing=self.cfg.pt_jump_fraction,
                jump_length_min_m=self.cfg.pt_jump_length_min_m,
            )
            pt_model = EmissionDerivedPeierlsTaylorModel(pt_cfg)
            rho_forest = max(
                float(self.cfg.pt_forest_density_floor_m2)
                + self.retained_count / max(self.length_m ** 2, 1.0e-30),
                1.0,
            )
            rates = pt_model.rates(tau, rho_forest, T_K, b)
            glide_rate = float(np.asarray(rates['peierls_rate_s']))
            detrap_rate = float(np.asarray(rates['taylor_completion_rate_s']))
            series_rate = float(np.asarray(rates['series_rate_s']))
            pt_diag = {
                'peierls_rate_s': glide_rate,
                'taylor_single_hit_rate_s': float(np.asarray(rates['taylor_single_hit_rate_s'])),
                'taylor_completion_rate_s': detrap_rate,
                'peierls_taylor_series_rate_s': series_rate,
                'taylor_m_eff': float(np.asarray(rates['taylor_m_eff'])),
                'G_peierls_eV': float(np.asarray(rates['G_peierls_eV'])),
                'G_taylor_eV': float(np.asarray(rates['G_taylor_eV'])),
                'rho_forest_m2': rho_forest,
            }
        else:
            glide_rate = self._activated_rate(
                self.cfg.glide_nu0_s, self.cfg.glide_barrier_eV,
                self.cfg.glide_activation_volume_b3,
                self.cfg.glide_stress_fraction * tau, T_K, b)
            detrap_rate = self._activated_rate(
                self.cfg.detrap_nu0_s, self.cfg.detrap_barrier_eV,
                self.cfg.detrap_activation_volume_b3, tau, T_K, b)
            series_rate = 1.0 / (
                1.0 / max(glide_rate, 1.0e-300)
                + 1.0 / max(detrap_rate, 1.0e-300)
            )
        velocity = max(float(self.cfg.glide_step_m), 0.0) * glide_rate
        trap_rate = self._activated_rate(
            self.cfg.trap_nu0_s, self.cfg.trap_barrier_eV,
            self.cfg.trap_activation_volume_b3, tau, T_K, b)
        retained_recovery_rate = self._activated_rate(
            self.cfg.retained_recovery_nu0_s,
            self.cfg.retained_recovery_barrier_eV,
            self.cfg.retained_recovery_activation_volume_b3,
            tau, T_K, b)

        # Exact first-order reaction fractions over the full accepted interval.
        # Advection is remapped conservatively over the full travel distance, so
        # no event-count or velocity cap is introduced by the time discretization.
        ftrap = 1.0 - math.exp(-min(trap_rate * dt_s, 700.0)) if dt_s > 0 else 0.0
        fdetrap = 1.0 - math.exp(-min(detrap_rate * dt_s, 700.0)) if dt_s > 0 else 0.0
        frec_r = 1.0 - math.exp(-min(retained_recovery_rate * dt_s, 700.0)) if dt_s > 0 else 0.0
        frec_m = 1.0 - math.exp(-min(max(self.cfg.mobile_recovery_rate_s, 0.0) * dt_s, 700.0)) if dt_s > 0 else 0.0

        dm_to_r = self.mobile * ftrap
        dr_to_m = self.retained * fdetrap
        dr_rec = np.maximum(self.retained - dr_to_m, 0.0) * frec_r
        dm_rec = np.maximum(self.mobile - dm_to_r, 0.0) * frec_m
        self.mobile += dr_to_m - dm_to_r - dm_rec
        self.retained += dm_to_r - dr_to_m - dr_rec
        trapped = float(np.sum(dm_to_r))
        detrapped = float(np.sum(dr_to_m))
        recovered = float(np.sum(dr_rec) + np.sum(dm_rec))

        annihilated = 0.0
        kpair = max(float(self.cfg.pair_annihilation_rate_per_count_s), 0.0)
        if kpair > 0.0 and self.n_systems >= 2 and dt_s > 0.0:
            for i in range(0, self.n_systems - 1, 2):
                pair = np.minimum(self.retained[i], self.retained[i + 1])
                frac = 1.0 - np.exp(-np.minimum(kpair * pair * dt_s, 700.0))
                d = pair * frac
                self.retained[i] -= d
                self.retained[i + 1] -= d
                annihilated += 2.0 * float(np.sum(d))

        self.mobile, escaped = self._advect_forward_field(self.mobile, velocity * dt_s)
        nsub = 1
        self.mobile = np.maximum(self.mobile, 0.0)
        self.retained = np.maximum(self.retained, 0.0)
        self.escaped_total += escaped
        self.recovered_total += recovered + annihilated
        self.time_s += dt_s
        return {
            "dN_emit": float(np.sum(emitted)),
            "dN_source_recovered": source_recovered,
            "dN_trapped": trapped,
            "dN_detrapped": detrapped,
            "dN_escaped": escaped,
            "dN_recovered": recovered,
            "dN_annihilated": annihilated,
            "glide_rate_s": glide_rate,
            "glide_velocity_m_s": velocity,
            "trap_rate_s": trap_rate,
            "detrap_rate_s": detrap_rate,
            "retained_recovery_rate_s": retained_recovery_rate,
            "peierls_taylor_series_rate_s": series_rate,
            "transport_substeps": float(nsub),
            "available_site_fraction": self.available_site_fraction,
            **pt_diag,
        }

    def _advect_forward_field(self, field: np.ndarray, distance_m: float) -> tuple[np.ndarray, float]:
        """Conservatively advect a piecewise-constant field away from the tip."""
        d = max(float(distance_m), 0.0)
        if d <= 0.0:
            return field.copy(), 0.0
        out = np.zeros_like(field)
        lost = 0.0
        for i in range(self.n_bins):
            left = i * self.dx + d
            right = (i + 1) * self.dx + d
            mass = field[:, i]
            if left >= self.length_m:
                lost += float(np.sum(mass))
                continue
            inside_left = max(left, 0.0)
            inside_right = min(right, self.length_m)
            inside_len = max(inside_right - inside_left, 0.0)
            if inside_len < self.dx:
                lost += float(np.sum(mass) * (1.0 - inside_len / self.dx))
            if inside_len <= 0.0:
                continue
            j0 = max(int(math.floor(inside_left / self.dx)), 0)
            j1 = min(int(math.floor((inside_right - 1.0e-15 * self.dx) / self.dx)), self.n_bins - 1)
            for j in range(j0, j1 + 1):
                ol = max(inside_left, j * self.dx)
                orr = min(inside_right, (j + 1) * self.dx)
                frac = max(orr - ol, 0.0) / self.dx
                if frac > 0.0:
                    out[:, j] += mass * frac
        return out, lost

    # ------------------------------------------------------------------
    # Moving-frame crack advance
    # ------------------------------------------------------------------
    def _translate_field(self, field: np.ndarray, distance_m: float) -> tuple[np.ndarray, float]:
        d = max(float(distance_m), 0.0)
        if d <= 0.0:
            return field.copy(), 0.0
        out = np.zeros_like(field)
        lost = 0.0
        # Piecewise-constant finite-volume remap. Old cell [i dx,(i+1)dx]
        # becomes [i dx-d,(i+1)dx-d] in the new tip frame.
        for i in range(self.n_bins):
            left = i * self.dx - d
            right = (i + 1) * self.dx - d
            mass = field[:, i]
            if right <= 0.0 or left >= self.length_m:
                lost += float(np.sum(mass))
                continue
            inside_left = max(left, 0.0)
            inside_right = min(right, self.length_m)
            inside_len = max(inside_right - inside_left, 0.0)
            if inside_len < self.dx:
                lost += float(np.sum(mass) * (1.0 - inside_len / self.dx))
            if inside_len <= 0.0:
                continue
            j0 = max(int(math.floor(inside_left / self.dx)), 0)
            j1 = min(int(math.floor((inside_right - 1.0e-15 * self.dx) / self.dx)), self.n_bins - 1)
            for j in range(j0, j1 + 1):
                ol = max(inside_left, j * self.dx)
                orr = min(inside_right, (j + 1) * self.dx)
                frac = max(orr - ol, 0.0) / self.dx
                if frac > 0.0:
                    out[:, j] += mass * frac
        return out, lost

    def advance(self, distance_m: float) -> dict[str, float]:
        d = max(float(distance_m), 0.0)
        old_mobile = self.mobile
        old_retained = self.retained
        old_slip = self.accumulated_slip
        self.mobile, lost_m = self._translate_field(old_mobile, d)
        self.retained, lost_r = self._translate_field(old_retained, d)
        self.accumulated_slip, lost_s = self._translate_field(old_slip, d)
        self.wake_mobile_total += lost_m
        self.wake_retained_total += lost_r
        self.wake_slip_total += lost_s
        self.advance_total_m += d

        # A moving source zone samples virgin material.  Replenishment is tied to
        # physical advance, not numerical step count.
        Lrefresh = max(float(self.cfg.source_refresh_length_m), self.dx)
        fresh = min(d / Lrefresh, 1.0)
        source_refreshed = (self.site_capacity - self.available_sites) * fresh
        self.available_sites += source_refreshed
        return {
            "wake_mobile": lost_m,
            "wake_retained": lost_r,
            "wake_slip": lost_s,
            "source_sites_refreshed": float(np.sum(source_refreshed)),
        }

    def diagnostics(self, G_shear: float, nu: float, b: float,
                    r0: float, c_blunt: float) -> dict[str, float]:
        return {
            "mpz_mobile_count": self.mobile_count,
            "mpz_retained_count": self.retained_count,
            "mpz_active_count": self.active_count,
            "mpz_available_site_fraction": self.available_site_fraction,
            "mpz_K_shield_Pa_sqrt_m": self.shielding_K(G_shear, nu, b),
            "mpz_local_slip_count": self.local_slip_count(),
            "mpz_r_eff_m": self.blunted_radius(r0, c_blunt, b),
            "mpz_emitted_total": float(self.emitted_total),
            "mpz_escaped_total": float(self.escaped_total),
            "mpz_recovered_total": float(self.recovered_total),
            "mpz_wake_retained_total": float(self.wake_retained_total),
            "mpz_advance_total_m": float(self.advance_total_m),
        }


def parse_float_tuple(value: str | Iterable[float] | None, default=(1.0, 1.0)) -> tuple[float, ...]:
    if value is None:
        return tuple(float(x) for x in default)
    if isinstance(value, str):
        vals = [x for x in value.replace(",", " ").split() if x]
        return tuple(float(x) for x in vals)
    return tuple(float(x) for x in value)


def build_mpz_config_from_namespace(args: Any, *, default_length_m: float | None = None) -> MovingProcessZoneConfig:
    def get(name, default):
        value = getattr(args, name, default)
        return default if value is None else value
    L = float(get("mpz_length_m", default_length_m if default_length_m is not None else 2.0e-6))
    return MovingProcessZoneConfig(
        length_m=L,
        n_bins=int(get("mpz_n_bins", 40)),
        n_systems=int(get("mpz_n_systems", 2)),
        source_sites_per_system=float(get("mpz_source_sites_per_system", 200.0)),
        source_recovery_rate_s=float(get("mpz_source_recovery_rate_s", 0.0)),
        source_refresh_length_m=float(get("mpz_source_refresh_length_m", 2.5e-7)),
        source_bin_count=int(get("mpz_source_bin_count", 2)),
        shielding_orientation_factors=parse_float_tuple(get("mpz_shielding_factors", None)),
        mobile_shield_fraction=float(get("mpz_mobile_shield_fraction", 0.0)),
        shielding_core_m=float(get("mpz_shielding_core_m", 2.5e-10)),
        glide_nu0_s=float(get("mpz_glide_nu0_s", 1.0e11)),
        glide_barrier_eV=float(get("mpz_glide_barrier_eV", 0.80)),
        glide_activation_volume_b3=float(get("mpz_glide_activation_volume_b3", 8.0)),
        glide_step_m=float(get("mpz_glide_step_m", 2.5e-10)),
        glide_stress_fraction=float(get("mpz_glide_stress_fraction", 0.45)),
        use_emission_derived_pt=bool(get("mpz_use_emission_derived_pt", True)),
        pt_emit_G00_eV=float(get("emit_G00_eV", 1.94)),
        pt_emit_gT_eV_per_K=float(get("emit_gT_eV_per_K", 0.003934)),
        pt_emit_sigc0_Pa=float(get("emit_sigc0_GPa", 2.298)) * 1.0e9,
        pt_emit_sT_Pa_per_K=float(get("emit_sT_GPa_per_K", -6.564e-4)) * 1.0e9,
        pt_emit_Tref_K=float(get("emit_Tref_K", 481.33)),
        pt_emit_exp_a=float(get("emit_exp_a", 0.0845685)),
        pt_emit_exp_n=float(get("emit_exp_n", 1.0)),
        pt_emit_floor_frac=float(get("emit_floor_frac", 0.02)),
        pt_emit_floor_min_eV=float(get("emit_floor_min_eV", 1.0e-4)),
        pt_emit_floor_max_frac=float(get("emit_floor_max_frac", 0.95)),
        pt_peierls_energy_ratio=float(get("pt_peierls_energy_ratio", get("peierls_energy_scale", 0.005))),
        pt_peierls_entropy_ratio=float(get("pt_peierls_entropy_ratio", get("peierls_entropy_scale", 0.005))),
        pt_peierls_stress_ratio=float(get("pt_peierls_stress_ratio", get("peierls_stress_scale", 1.0))),
        pt_peierls_nu0_s=float(get("pt_peierls_nu0_s", get("nu0_peierls", 1.0e12))),
        pt_taylor_energy_ratio=float(get("pt_taylor_energy_ratio", get("taylor_energy_scale", 0.02))),
        pt_taylor_entropy_ratio=float(get("pt_taylor_entropy_ratio", get("taylor_entropy_scale", 0.02))),
        pt_taylor_stress_ratio=float(get("pt_taylor_stress_ratio", get("taylor_stress_scale", 1.0))),
        pt_taylor_nu0_s=float(get("pt_taylor_nu0_s", get("nu0_taylor", 1.0e11))),
        pt_taylor_corr_rho_c=float(get("pt_taylor_corr_rho_c", 1.0e14)),
        pt_taylor_renewal_time_s=float(get("pt_taylor_renewal_time_s", 1.0e-9)),
        pt_taylor_m_exponent=float(get("pt_taylor_m_exponent", 1.0)),
        pt_taylor_m_scale=float(get("pt_taylor_m_scale", 1.0)),
        pt_taylor_m_cap=float(get("pt_taylor_m_cap", float('inf'))),
        pt_forest_density_floor_m2=float(get("pt_forest_density_floor_m2", 5.0e12)),
        pt_mobile_fraction=float(get("pt_mobile_fraction", 0.01)),
        pt_mobile_saturation_density_m2=float(get("pt_mobile_saturation_density_m2", 1.0e14)),
        pt_mobile_density_floor_m2=float(get("pt_mobile_density_floor_m2", 1.0e6)),
        pt_jump_fraction=float(get("pt_jump_fraction", 1.0)),
        pt_jump_length_min_m=float(get("pt_jump_length_min_m", 2.5e-10)),
        pt_peierls_stress_fraction=float(get("pt_peierls_stress_fraction", 0.5773502691896258)),
        pt_taylor_stress_fraction=float(get("pt_taylor_stress_fraction", 0.5773502691896258)),
        pt_taylor_phi_max=float(get("pt_taylor_phi_max", 20.0)),
        trap_nu0_s=float(get("mpz_trap_nu0_s", 1.0e9)),
        trap_barrier_eV=float(get("mpz_trap_barrier_eV", 0.65)),
        trap_activation_volume_b3=float(get("mpz_trap_activation_volume_b3", 1.0)),
        detrap_nu0_s=float(get("mpz_detrap_nu0_s", 1.0e10)),
        detrap_barrier_eV=float(get("mpz_detrap_barrier_eV", 1.20)),
        detrap_activation_volume_b3=float(get("mpz_detrap_activation_volume_b3", 1.0)),
        retained_recovery_nu0_s=float(get("mpz_retained_recovery_nu0_s", 1.0e9)),
        retained_recovery_barrier_eV=float(get("mpz_retained_recovery_barrier_eV", 1.50)),
        retained_recovery_activation_volume_b3=float(get("mpz_retained_recovery_activation_volume_b3", 0.0)),
        mobile_recovery_rate_s=float(get("mpz_mobile_recovery_rate_s", 0.0)),
        pair_annihilation_rate_per_count_s=float(get("mpz_pair_annihilation_rate_per_count_s", 0.0)),
        blunting_length_m=float(get("mpz_blunting_length_m", 5.0e-7)),
        blunting_slip_fraction=float(get("mpz_blunting_slip_fraction", 1.0)),
        max_transport_cfl=float(get("mpz_max_transport_cfl", 0.35)),
        max_transport_substeps=int(get("mpz_max_transport_substeps", 2000)),
    )
