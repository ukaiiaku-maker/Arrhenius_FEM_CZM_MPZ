"""Emission-derived Peierls--Taylor kinetics for bulk and process-zone plasticity.

The production plasticity branch uses one fitted crack-tip emission EXP-floor
surface as the parent free-energy landscape. Peierls glide and Taylor obstacle
depinning are mechanism-specific scaled descendants of that same surface.
They are sequential kinetic bottlenecks, not additive quasi-static stresses.

The Taylor branch uses a correlated multi-hit completion probability. This
replaces the independent-site prefactor that grows too rapidly with forest
density and can create a nonphysical high-density stress downturn.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
import math
from typing import Any

import numpy as np
from scipy.special import gammainc

from .config import KB, EV_TO_J


@dataclass(frozen=True)
class ExpFloorSurface:
    """Direct bounded activation free-energy surface."""

    G00_eV: float = 1.94
    gT_eV_per_K: float = 0.003934
    sigc0_Pa: float = 2.298e9
    sT_Pa_per_K: float = -6.564e5
    Tref_K: float = 481.33
    a: float = 0.0845685
    n: float = 1.0
    floor_fraction: float = 0.02
    floor_min_eV: float = 1.0e-4
    floor_max_fraction: float = 0.95

    @classmethod
    def from_fracture_barrier(cls, barrier: Any) -> "ExpFloorSurface":
        return cls(
            G00_eV=float(getattr(barrier, "ef_G00_eV", 1.94)),
            gT_eV_per_K=float(getattr(barrier, "ef_gT_eV_per_K", 0.003934)),
            sigc0_Pa=float(getattr(barrier, "ef_sigc0_Pa", 2.298e9)),
            sT_Pa_per_K=float(getattr(barrier, "ef_sT_Pa_per_K", -6.564e5)),
            Tref_K=float(getattr(barrier, "ef_Tref_K", 481.33)),
            a=float(getattr(barrier, "ef_a", 0.0845685)),
            n=float(getattr(barrier, "ef_n", 1.0)),
            floor_fraction=float(getattr(barrier, "ef_floor_frac", 0.02)),
            floor_min_eV=float(getattr(barrier, "ef_floor_min_eV", 1.0e-4)),
            floor_max_fraction=float(getattr(barrier, "ef_floor_max_frac", 0.95)),
        )

    @classmethod
    def from_namespace(cls, obj: Any, prefix: str = "emit_") -> "ExpFloorSurface":
        get = lambda name, default: getattr(obj, f"{prefix}{name}", default)
        return cls(
            G00_eV=float(get("G00_eV", 1.94)),
            gT_eV_per_K=float(get("gT_eV_per_K", 0.003934)),
            sigc0_Pa=float(get("sigc0_GPa", 2.298)) * 1.0e9,
            sT_Pa_per_K=float(get("sT_GPa_per_K", -6.564e-4)) * 1.0e9,
            Tref_K=float(get("Tref_K", 481.33)),
            a=float(get("exp_a", 0.0845685)),
            n=float(get("exp_n", 1.0)),
            floor_fraction=float(get("floor_frac", 0.02)),
            floor_min_eV=float(get("floor_min_eV", 1.0e-4)),
            floor_max_fraction=float(get("floor_max_frac", 0.95)),
        )

    def as_dict(self) -> dict[str, float]:
        return asdict(self)


@dataclass(frozen=True)
class MechanismScale:
    """Scale one mechanism from the parent emission EXP-floor surface."""

    energy_ratio: float
    entropy_ratio: float
    stress_ratio: float = 1.0
    rate_prefactor_s: float = 1.0e12


@dataclass(frozen=True)
class CorrelatedTaylorConfig:
    """Density-dependent cooperative Taylor completion parameters."""

    rho_c_m2: float = 1.0e14
    renewal_time_s: float = 1.0e-9
    m_exponent: float = 1.0
    m_scale: float = 1.0
    m_cap: float = float("inf")

    def hit_order(self, rho_forest_m2: np.ndarray | float) -> np.ndarray:
        """Continuous cooperative hit order.

        One hit is recovered below the correlation-density crossover. Above
        the crossover the number of obstacles in a correlation length grows as
        sqrt(rho/rho_c); the exponent permits a controlled collective crossover.
        A finite m_cap represents the obstacle count in one correlation domain;
        it is not a cap on total dislocation density.
        """
        rho = np.maximum(np.asarray(rho_forest_m2, dtype=float), 0.0)
        rho_c = max(float(self.rho_c_m2), 1.0)
        obstacle_ratio = np.sqrt(rho / rho_c)
        excess = np.maximum(obstacle_ratio - 1.0, 0.0)
        m = 1.0 + max(float(self.m_scale), 0.0) * np.power(
            excess, max(float(self.m_exponent), 1.0e-9)
        )
        cap = float(self.m_cap)
        if np.isfinite(cap):
            m = np.minimum(m, max(cap, 1.0))
        return np.maximum(m, 1.0)


@dataclass(frozen=True)
class EmissionDerivedPeierlsTaylorConfig:
    """Shared production kinetics for MPZ transport and bulk FEM plasticity."""

    parent: ExpFloorSurface = ExpFloorSurface()
    peierls: MechanismScale = MechanismScale(0.005, 0.005, 1.0, 1.0e12)
    taylor: MechanismScale = MechanismScale(0.02, 0.02, 1.0, 1.0e11)
    correlated_taylor: CorrelatedTaylorConfig = CorrelatedTaylorConfig()

    peierls_stress_fraction: float = 1.0 / math.sqrt(3.0)
    taylor_stress_fraction: float = 1.0 / math.sqrt(3.0)
    taylor_phi_max: float = 20.0
    mobile_fraction_low_density: float = 0.01
    mobile_saturation_density_m2: float = 1.0e14
    mobile_density_floor_m2: float = 1.0e6
    jump_fraction_of_forest_spacing: float = 1.0
    jump_length_min_m: float = 2.5e-10
    equivalent_strain_factor: float = 1.0 / math.sqrt(3.0)
    rate_cap_s: float = float("inf")

    def as_dict(self) -> dict[str, Any]:
        return {
            "parent": self.parent.as_dict(),
            "peierls": asdict(self.peierls),
            "taylor": asdict(self.taylor),
            "correlated_taylor": asdict(self.correlated_taylor),
            **{
                k: v
                for k, v in asdict(self).items()
                if k not in {"parent", "peierls", "taylor", "correlated_taylor"}
            },
        }


class EmissionDerivedPeierlsTaylorModel:
    """Evaluate emission-derived Peierls/Taylor kinetics."""

    def __init__(self, cfg: EmissionDerivedPeierlsTaylorConfig):
        self.cfg = cfg

    @staticmethod
    def _scaled_surface_values(
        parent: ExpFloorSurface,
        scale: MechanismScale,
        sigma_Pa: np.ndarray | float,
        T_K: float,
    ) -> np.ndarray:
        sigma = np.maximum(np.asarray(sigma_Pa, dtype=float), 0.0)
        dT = float(T_K) - float(parent.Tref_K)
        G0 = max(
            float(scale.energy_ratio) * float(parent.G00_eV)
            + float(scale.entropy_ratio) * float(parent.gT_eV_per_K) * dT,
            1.0e-12,
        )
        sigc_parent = float(parent.sigc0_Pa) + float(parent.sT_Pa_per_K) * dT
        sigc = max(float(scale.stress_ratio) * sigc_parent, 1.0)
        raw_floor = max(
            float(parent.floor_min_eV) * max(float(scale.energy_ratio), 1.0e-12),
            float(parent.floor_fraction) * G0,
        )
        floor = min(float(parent.floor_max_fraction) * G0, raw_floor)
        x = sigma / sigc
        return np.maximum(
            floor
            + (G0 - floor)
            * np.exp(
                -max(float(parent.a), 0.0)
                * np.power(x, max(float(parent.n), 1.0e-9))
            ),
            0.0,
        )

    def barrier_eV(
        self, mechanism: str, sigma_Pa: np.ndarray | float, T_K: float
    ) -> np.ndarray:
        key = str(mechanism).lower()
        if key in {"peierls", "p"}:
            scale = self.cfg.peierls
        elif key in {"taylor", "t"}:
            scale = self.cfg.taylor
        elif key in {"emission", "emit", "e"}:
            scale = MechanismScale(1.0, 1.0, 1.0, 1.0)
        else:
            raise ValueError(f"unknown mechanism: {mechanism}")
        return self._scaled_surface_values(
            self.cfg.parent, scale, sigma_Pa, T_K
        )

    @staticmethod
    def _arrhenius_rate(
        barrier_eV: np.ndarray, T_K: float, nu0_s: float
    ) -> np.ndarray:
        exponent = (
            -np.asarray(barrier_eV, dtype=float)
            * EV_TO_J
            / max(KB * float(T_K), 1.0e-30)
        )
        return max(float(nu0_s), 0.0) * np.exp(
            np.clip(exponent, -700.0, 0.0)
        )

    def peierls_rate(
        self, sigma_eq_Pa: np.ndarray | float, T_K: float
    ) -> np.ndarray:
        tau = max(
            float(self.cfg.peierls_stress_fraction), 0.0
        ) * np.maximum(np.asarray(sigma_eq_Pa, dtype=float), 0.0)
        G = self.barrier_eV("peierls", tau, T_K)
        return self._arrhenius_rate(
            G, T_K, self.cfg.peierls.rate_prefactor_s
        )

    def taylor_local_stress(
        self,
        sigma_eq_Pa: np.ndarray | float,
        rho_forest_m2: np.ndarray | float,
        b_m: float,
    ) -> np.ndarray:
        rho = np.maximum(np.asarray(rho_forest_m2, dtype=float), 1.0)
        b = max(abs(float(b_m)), 1.0e-30)
        spacing = 1.0 / (2.0 * np.sqrt(rho))
        phi = np.minimum(
            spacing / b, max(float(self.cfg.taylor_phi_max), 1.0)
        )
        tau = max(
            float(self.cfg.taylor_stress_fraction), 0.0
        ) * np.maximum(np.asarray(sigma_eq_Pa, dtype=float), 0.0)
        return tau * phi

    def taylor_single_hit_rate(
        self,
        sigma_eq_Pa: np.ndarray | float,
        rho_forest_m2: np.ndarray | float,
        T_K: float,
        b_m: float,
    ) -> np.ndarray:
        sigma_local = self.taylor_local_stress(
            sigma_eq_Pa, rho_forest_m2, b_m
        )
        G = self.barrier_eV("taylor", sigma_local, T_K)
        return self._arrhenius_rate(
            G, T_K, self.cfg.taylor.rate_prefactor_s
        )

    def taylor_completion_rate(
        self,
        sigma_eq_Pa: np.ndarray | float,
        rho_forest_m2: np.ndarray | float,
        T_K: float,
        b_m: float,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        h1 = self.taylor_single_hit_rate(
            sigma_eq_Pa, rho_forest_m2, T_K, b_m
        )
        m = self.cfg.correlated_taylor.hit_order(rho_forest_m2)
        tc = max(
            float(self.cfg.correlated_taylor.renewal_time_s), 1.0e-30
        )
        exposure = np.minimum(np.maximum(h1 * tc, 0.0), 1.0e12)
        rate = gammainc(m, exposure) / tc
        return np.maximum(rate, 0.0), m, h1

    @staticmethod
    def series_rate(
        rate_a: np.ndarray | float, rate_b: np.ndarray | float
    ) -> np.ndarray:
        a = np.maximum(np.asarray(rate_a, dtype=float), 0.0)
        b = np.maximum(np.asarray(rate_b, dtype=float), 0.0)
        return 1.0 / (
            1.0 / np.maximum(a, 1.0e-300)
            + 1.0 / np.maximum(b, 1.0e-300)
        )

    def mobile_density(
        self, rho_forest_m2: np.ndarray | float
    ) -> np.ndarray:
        rho = np.maximum(np.asarray(rho_forest_m2, dtype=float), 0.0)
        sat = max(float(self.cfg.mobile_saturation_density_m2), 1.0)
        mobile = (
            max(float(self.cfg.mobile_fraction_low_density), 0.0)
            * rho
            / (1.0 + rho / sat)
        )
        return np.maximum(
            mobile, max(float(self.cfg.mobile_density_floor_m2), 0.0)
        )

    def jump_length(
        self, rho_forest_m2: np.ndarray | float
    ) -> np.ndarray:
        rho = np.maximum(np.asarray(rho_forest_m2, dtype=float), 1.0)
        spacing = 1.0 / (2.0 * np.sqrt(rho))
        return np.maximum(
            max(
                float(self.cfg.jump_fraction_of_forest_spacing), 0.0
            )
            * spacing,
            max(float(self.cfg.jump_length_min_m), 0.0),
        )

    def rates(
        self,
        sigma_eq_Pa: np.ndarray | float,
        rho_forest_m2: np.ndarray | float,
        T_K: float,
        b_m: float,
    ) -> dict[str, np.ndarray]:
        p = self.peierls_rate(sigma_eq_Pa, T_K)
        t, m, t1 = self.taylor_completion_rate(
            sigma_eq_Pa, rho_forest_m2, T_K, b_m
        )
        series = self.series_rate(p, t)
        rho_m = self.mobile_density(rho_forest_m2)
        jump = self.jump_length(rho_forest_m2)
        eq_rate = (
            max(float(self.cfg.equivalent_strain_factor), 0.0)
            * rho_m
            * abs(float(b_m))
            * jump
            * series
        )
        cap = float(self.cfg.rate_cap_s)
        if np.isfinite(cap):
            eq_rate = np.minimum(eq_rate, max(cap, 0.0))
        return {
            "peierls_rate_s": np.asarray(p, dtype=float),
            "taylor_single_hit_rate_s": np.asarray(t1, dtype=float),
            "taylor_completion_rate_s": np.asarray(t, dtype=float),
            "series_rate_s": np.asarray(series, dtype=float),
            "equivalent_plastic_rate_s": np.asarray(eq_rate, dtype=float),
            "taylor_m_eff": np.asarray(m, dtype=float),
            "rho_mobile_m2": np.asarray(rho_m, dtype=float),
            "jump_length_m": np.asarray(jump, dtype=float),
            "taylor_local_stress_Pa": np.asarray(
                self.taylor_local_stress(
                    sigma_eq_Pa, rho_forest_m2, b_m
                ),
                dtype=float,
            ),
            "G_peierls_eV": np.asarray(
                self.barrier_eV(
                    "peierls",
                    max(
                        float(self.cfg.peierls_stress_fraction), 0.0
                    )
                    * np.maximum(
                        np.asarray(sigma_eq_Pa, dtype=float), 0.0
                    ),
                    T_K,
                ),
                dtype=float,
            ),
            "G_taylor_eV": np.asarray(
                self.barrier_eV(
                    "taylor",
                    self.taylor_local_stress(
                        sigma_eq_Pa, rho_forest_m2, b_m
                    ),
                    T_K,
                ),
                dtype=float,
            ),
        }

    def flow_stress(
        self,
        rho_forest_m2: np.ndarray | float,
        T_K: float,
        target_equivalent_rate_s: float,
        b_m: float,
        sigma_max_Pa: float = 40.0e9,
        iterations: int = 80,
    ) -> np.ndarray:
        rho = np.asarray(rho_forest_m2, dtype=float)
        target = max(float(target_equivalent_rate_s), 0.0)
        lo = np.zeros_like(rho)
        hi = np.full_like(rho, max(float(sigma_max_Pa), 1.0))
        rate_hi = self.rates(
            hi, rho, T_K, b_m
        )["equivalent_plastic_rate_s"]
        unresolved = rate_hi < target
        for _ in range(8):
            if not np.any(unresolved):
                break
            hi = np.where(unresolved, hi * 2.0, hi)
            rate_hi = self.rates(
                hi, rho, T_K, b_m
            )["equivalent_plastic_rate_s"]
            unresolved = rate_hi < target
        for _ in range(max(int(iterations), 1)):
            mid = 0.5 * (lo + hi)
            rate = self.rates(
                mid, rho, T_K, b_m
            )["equivalent_plastic_rate_s"]
            below = rate < target
            lo = np.where(below, mid, lo)
            hi = np.where(below, hi, mid)
        out = hi
        return np.where(unresolved, np.nan, out)


def config_from_dislocation_config(
    disl_cfg: Any,
) -> EmissionDerivedPeierlsTaylorConfig:
    """Build the shared model from the FEM dislocation configuration."""
    parent = ExpFloorSurface(
        G00_eV=float(getattr(disl_cfg, "pt_emit_G00_eV", 1.94)),
        gT_eV_per_K=float(
            getattr(disl_cfg, "pt_emit_gT_eV_per_K", 0.003934)
        ),
        sigc0_Pa=float(
            getattr(disl_cfg, "pt_emit_sigc0_Pa", 2.298e9)
        ),
        sT_Pa_per_K=float(
            getattr(disl_cfg, "pt_emit_sT_Pa_per_K", -6.564e5)
        ),
        Tref_K=float(getattr(disl_cfg, "pt_emit_Tref_K", 481.33)),
        a=float(getattr(disl_cfg, "pt_emit_exp_a", 0.0845685)),
        n=float(getattr(disl_cfg, "pt_emit_exp_n", 1.0)),
        floor_fraction=float(
            getattr(disl_cfg, "pt_emit_floor_frac", 0.02)
        ),
        floor_min_eV=float(
            getattr(disl_cfg, "pt_emit_floor_min_eV", 1.0e-4)
        ),
        floor_max_fraction=float(
            getattr(disl_cfg, "pt_emit_floor_max_frac", 0.95)
        ),
    )
    return EmissionDerivedPeierlsTaylorConfig(
        parent=parent,
        peierls=MechanismScale(
            energy_ratio=float(
                getattr(disl_cfg, "pt_peierls_energy_ratio", 0.005)
            ),
            entropy_ratio=float(
                getattr(disl_cfg, "pt_peierls_entropy_ratio", 0.005)
            ),
            stress_ratio=float(
                getattr(disl_cfg, "pt_peierls_stress_ratio", 1.0)
            ),
            rate_prefactor_s=float(
                getattr(disl_cfg, "pt_peierls_nu0_s", 1.0e12)
            ),
        ),
        taylor=MechanismScale(
            energy_ratio=float(
                getattr(disl_cfg, "pt_taylor_energy_ratio", 0.02)
            ),
            entropy_ratio=float(
                getattr(disl_cfg, "pt_taylor_entropy_ratio", 0.02)
            ),
            stress_ratio=float(
                getattr(disl_cfg, "pt_taylor_stress_ratio", 1.0)
            ),
            rate_prefactor_s=float(
                getattr(disl_cfg, "pt_taylor_nu0_s", 1.0e11)
            ),
        ),
        correlated_taylor=CorrelatedTaylorConfig(
            rho_c_m2=float(
                getattr(disl_cfg, "pt_taylor_corr_rho_c", 1.0e14)
            ),
            renewal_time_s=float(
                getattr(
                    disl_cfg, "pt_taylor_renewal_time_s", 1.0e-9
                )
            ),
            m_exponent=float(
                getattr(disl_cfg, "pt_taylor_m_exponent", 1.0)
            ),
            m_scale=float(
                getattr(disl_cfg, "pt_taylor_m_scale", 1.0)
            ),
            m_cap=float(
                getattr(disl_cfg, "pt_taylor_m_cap", float("inf"))
            ),
        ),
        peierls_stress_fraction=float(
            getattr(
                disl_cfg,
                "pt_peierls_stress_fraction",
                1.0 / math.sqrt(3.0),
            )
        ),
        taylor_stress_fraction=float(
            getattr(
                disl_cfg,
                "pt_taylor_stress_fraction",
                1.0 / math.sqrt(3.0),
            )
        ),
        taylor_phi_max=float(
            getattr(disl_cfg, "pt_taylor_phi_max", 20.0)
        ),
        mobile_fraction_low_density=float(
            getattr(disl_cfg, "pt_mobile_fraction", 0.01)
        ),
        mobile_saturation_density_m2=float(
            getattr(
                disl_cfg, "pt_mobile_saturation_density_m2", 1.0e14
            )
        ),
        mobile_density_floor_m2=float(
            getattr(disl_cfg, "pt_mobile_density_floor_m2", 1.0e6)
        ),
        jump_fraction_of_forest_spacing=float(
            getattr(disl_cfg, "pt_jump_fraction", 1.0)
        ),
        jump_length_min_m=float(
            getattr(disl_cfg, "pt_jump_length_min_m", 2.5e-10)
        ),
        equivalent_strain_factor=float(
            getattr(
                disl_cfg,
                "pt_equivalent_strain_factor",
                1.0 / math.sqrt(3.0),
            )
        ),
        rate_cap_s=float(
            getattr(disl_cfg, "dot_ep_max", float("inf"))
        ),
    )
