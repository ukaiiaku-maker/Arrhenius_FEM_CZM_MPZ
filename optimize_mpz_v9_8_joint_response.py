#!/usr/bin/env python3
"""Joint multi-fidelity optimization of fracture and plasticity parameters.

This v9.8 stage optimizes cleavage, emission, Peierls, Taylor, and reduced
process-zone parameters together against one requested response class.  It is
an analytical fidelity-0 optimizer intended to locate physically admissible
basins before the leading candidates are promoted to the spatial moving-PZ and
2-D FEM/CZM solvers.

Global search uses SciPy differential evolution.  Each basin is initialized
from a diverse first-passage atlas row and is then refined locally with bounded
Nelder--Mead.  No mechanism-dominance assumption is imposed.  The only direct
Peierls/Taylor ordering constraint is

    G_P(sigma,T) <= G_T(sigma,T)

on the same resolved-stress grid.  Both mechanisms retain uncapped EXP-floor
kinetics and exact forward-minus-reverse detailed balance.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.optimize import differential_evolution, minimize
from scipy.special import gammainc
from scipy.stats import qmc

from arrhenius_fracture.config import EV_TO_J, KB
from arrhenius_fracture.emission_derived_plasticity import (
    CorrelatedTaylorConfig,
    EmissionDerivedPeierlsTaylorConfig,
    ExpFloorSurface,
)
from arrhenius_fracture.emission_derived_plasticity_v97 import (
    EmissionDerivedPeierlsTaylorModel,
    IndependentEntropyMechanismScale,
)


PARAMETER_NAMES = (
    "cleave_G00_eV",
    "cleave_gT_eV_per_K",
    "cleave_sigc0_GPa",
    "emit_G00_eV",
    "emit_gT_eV_per_K",
    "emit_sigc0_GPa",
    "peierls_H0_eV",
    "delta_H_PT_eV",
    "peierls_activation_entropy_kB",
    "taylor_activation_entropy_kB",
    "log10_taylor_corr_rho_c_m2",
    "log10_taylor_corr_scale",
    "log10_mobile_fraction",
    "log10_source_sites_per_system",
    "log10_recovery_rate_s",
    "log10_source_refresh_length_um",
    "c_blunt",
)

DEFAULT_BOUNDS = {
    "cleave_G00_eV": (0.5, 5.0),
    "cleave_gT_eV_per_K": (-0.003, 0.006),
    "cleave_sigc0_GPa": (1.0, 8.0),
    "emit_G00_eV": (0.5, 5.0),
    "emit_gT_eV_per_K": (-0.003, 0.012),
    "emit_sigc0_GPa": (0.3, 8.0),
    "peierls_H0_eV": (0.02, 8.0),
    "delta_H_PT_eV": (0.0, 12.0),
    "peierls_activation_entropy_kB": (-80.0, 40.0),
    "taylor_activation_entropy_kB": (-80.0, 40.0),
    "log10_taylor_corr_rho_c_m2": (10.0, 17.0),
    "log10_taylor_corr_scale": (-1.0, 1.0),
    "log10_mobile_fraction": (-5.0, -0.3),
    "log10_source_sites_per_system": (0.0, 4.0),
    "log10_recovery_rate_s": (-10.0, 1.0),
    "log10_source_refresh_length_um": (-1.0, 3.0),
    "c_blunt": (0.0, 3.0),
}

SHAPE_COLUMNS = (
    "cleave_sT_GPa_per_K",
    "cleave_exp_a",
    "cleave_exp_n",
    "cleave_floor_frac",
    "emit_sT_GPa_per_K",
    "emit_exp_a",
    "emit_exp_n",
    "emit_floor_frac",
)


def parse_floats(text: str) -> list[float]:
    return [float(x) for x in str(text).replace(",", " ").split() if x]


def bounds_array() -> list[tuple[float, float]]:
    return [DEFAULT_BOUNDS[name] for name in PARAMETER_NAMES]


def vector_to_parameters(x: np.ndarray) -> dict[str, float]:
    values = {name: float(value) for name, value in zip(PARAMETER_NAMES, x)}
    values.update(
        {
            "taylor_H0_eV": values["peierls_H0_eV"]
            + values["delta_H_PT_eV"],
            "taylor_corr_rho_c_m2": 10.0
            ** values["log10_taylor_corr_rho_c_m2"],
            "taylor_corr_scale": 10.0
            ** values["log10_taylor_corr_scale"],
            "mobile_fraction": 10.0 ** values["log10_mobile_fraction"],
            "source_sites_per_system": 10.0
            ** values["log10_source_sites_per_system"],
            "recovery_rate_s": 10.0 ** values["log10_recovery_rate_s"],
            "source_refresh_length_um": 10.0
            ** values["log10_source_refresh_length_um"],
        }
    )
    return values


def huber(value: float, delta: float = 2.0) -> float:
    a = abs(float(value))
    if a <= delta:
        return 0.5 * a * a
    return delta * (a - 0.5 * delta)


def exp_floor_barrier_eV(
    sigma_Pa: np.ndarray | float,
    T_K: np.ndarray | float,
    G00_eV: float,
    gT_eV_per_K: float,
    sigc0_GPa: float,
    sT_GPa_per_K: float,
    a: float,
    n: float,
    floor_fraction: float,
    Tref_K: float,
    floor_min_eV: float = 1.0e-4,
    floor_max_fraction: float = 0.95,
) -> np.ndarray:
    sigma = np.maximum(np.asarray(sigma_Pa, dtype=float), 0.0)
    T = np.asarray(T_K, dtype=float)
    G0 = np.maximum(G00_eV + gT_eV_per_K * (T - Tref_K), 1.0e-12)
    sigc = np.maximum(
        (sigc0_GPa + sT_GPa_per_K * (T - Tref_K)) * 1.0e9,
        1.0,
    )
    floor = np.minimum(
        floor_max_fraction * G0,
        np.maximum(floor_min_eV, floor_fraction * G0),
    )
    return floor + (G0 - floor) * np.exp(
        -max(float(a), 0.0)
        * np.power(sigma / sigc, max(float(n), 1.0e-9))
    )


def cleavage_effective_rate(
    G_eV: np.ndarray,
    T_K: np.ndarray,
    nu0_s: float,
    m_hits: float,
    correlation_time_s: float,
) -> np.ndarray:
    raw = float(nu0_s) * np.exp(
        np.clip(-np.asarray(G_eV) * EV_TO_J / (KB * T_K), -700.0, 0.0)
    )
    if m_hits <= 1.0 + 1.0e-12:
        return raw
    return gammainc(
        float(m_hits),
        np.minimum(raw * float(correlation_time_s), 1.0e12),
    ) / float(correlation_time_s)


def first_passage_proxy(
    p: dict[str, float],
    shape: dict[str, float],
    temperatures: np.ndarray,
    *,
    Kdot: float,
    dK: float,
    Kmax: float,
    r_pz_m: float,
    nu0_c: float,
    nu0_e: float,
    cleavage_hits: float,
    cleavage_correlation_time_s: float,
    Tref_K: float,
) -> dict[str, np.ndarray]:
    T = np.asarray(temperatures, dtype=float)
    B = np.zeros_like(T)
    H = np.zeros_like(T)
    Kc = np.full_like(T, np.nan)
    Hc = np.full_like(T, np.nan)

    def barriers(K: float) -> tuple[np.ndarray, np.ndarray]:
        sigma = K * 1.0e6 / math.sqrt(2.0 * math.pi * r_pz_m)
        Gc = exp_floor_barrier_eV(
            sigma,
            T,
            p["cleave_G00_eV"],
            p["cleave_gT_eV_per_K"],
            p["cleave_sigc0_GPa"],
            shape["cleave_sT_GPa_per_K"],
            shape["cleave_exp_a"],
            shape["cleave_exp_n"],
            shape["cleave_floor_frac"],
            Tref_K,
        )
        Ge = exp_floor_barrier_eV(
            sigma,
            T,
            p["emit_G00_eV"],
            p["emit_gT_eV_per_K"],
            p["emit_sigc0_GPa"],
            shape["emit_sT_GPa_per_K"],
            shape["emit_exp_a"],
            shape["emit_exp_n"],
            shape["emit_floor_frac"],
            Tref_K,
        )
        return Gc, Ge

    Gc0, Ge0 = barriers(0.0)
    lc_prev = cleavage_effective_rate(
        Gc0, T, nu0_c, cleavage_hits, cleavage_correlation_time_s
    )
    le_prev = nu0_e * np.exp(np.clip(-Ge0 * EV_TO_J / (KB * T), -700.0, 0.0))
    Kprev = 0.0
    for istep in range(1, int(math.ceil(Kmax / dK)) + 1):
        K = min(istep * dK, Kmax)
        step = K - Kprev
        Gc, Ge = barriers(K)
        lc = cleavage_effective_rate(
            Gc, T, nu0_c, cleavage_hits, cleavage_correlation_time_s
        )
        le = nu0_e * np.exp(np.clip(-Ge * EV_TO_J / (KB * T), -700.0, 0.0))
        dt = step / Kdot
        dB = 0.5 * (lc_prev + lc) * dt
        dH = 0.5 * (le_prev + le) * dt
        Bnew = B + dB
        Hnew = H + dH
        active = ~np.isfinite(Kc)
        crossed = active & (Bnew >= 1.0)
        if np.any(crossed):
            frac = np.clip(
                (1.0 - B) / np.maximum(Bnew - B, 1.0e-300), 0.0, 1.0
            )
            Kcross = Kprev + frac * step
            Hcross = H + frac * dH
            Kc[crossed] = Kcross[crossed]
            Hc[crossed] = Hcross[crossed]
        B[active] = Bnew[active]
        H[active] = Hnew[active]
        lc_prev = lc
        le_prev = le
        Kprev = K
        if np.all(np.isfinite(Kc)):
            break
    return {"K_intrinsic": Kc, "H_emit_at_Kc": Hc}


def build_pt_model(
    p: dict[str, float],
    shape: dict[str, float],
    *,
    Tref_K: float,
) -> EmissionDerivedPeierlsTaylorModel:
    parent = ExpFloorSurface(
        G00_eV=p["emit_G00_eV"],
        gT_eV_per_K=p["emit_gT_eV_per_K"],
        sigc0_Pa=p["emit_sigc0_GPa"] * 1.0e9,
        sT_Pa_per_K=shape["emit_sT_GPa_per_K"] * 1.0e9,
        Tref_K=Tref_K,
        a=shape["emit_exp_a"],
        n=shape["emit_exp_n"],
        floor_fraction=shape["emit_floor_frac"],
        floor_min_eV=1.0e-4,
        floor_max_fraction=0.95,
    )
    e0 = max(p["emit_G00_eV"], 1.0e-12)
    return EmissionDerivedPeierlsTaylorModel(
        EmissionDerivedPeierlsTaylorConfig(
            parent=parent,
            peierls=IndependentEntropyMechanismScale(
                p["peierls_H0_eV"] / e0,
                p["peierls_activation_entropy_kB"],
                1.0,
                1.0e12,
            ),
            taylor=IndependentEntropyMechanismScale(
                p["taylor_H0_eV"] / e0,
                p["taylor_activation_entropy_kB"],
                1.0,
                1.0e11,
            ),
            correlated_taylor=CorrelatedTaylorConfig(
                rho_c_m2=p["taylor_corr_rho_c_m2"],
                renewal_time_s=1.0,
                m_exponent=1.0,
                m_scale=p["taylor_corr_scale"],
                m_cap=float("inf"),
            ),
            mobile_fraction_low_density=p["mobile_fraction"],
            mobile_saturation_density_m2=float("inf"),
            mobile_density_floor_m2=0.0,
            jump_fraction_of_forest_spacing=1.0,
            jump_length_min_m=0.0,
            taylor_phi_max=float("inf"),
            rate_cap_s=float("inf"),
        )
    )


def barrier_order_margin_eV(
    model: EmissionDerivedPeierlsTaylorModel,
    temperatures: np.ndarray,
    stress_grid_Pa: np.ndarray,
) -> float:
    margin = float("inf")
    for T in np.asarray(temperatures, dtype=float):
        Gp = model._scaled_surface_values(
            model.cfg.parent, model.cfg.peierls, stress_grid_Pa, T
        )
        Gt = model._scaled_surface_values(
            model.cfg.parent, model.cfg.taylor, stress_grid_Pa, T
        )
        margin = min(margin, float(np.min(Gt - Gp)))
    return margin


def load_targets(
    path: Path, target_class: str, temperatures: list[float]
) -> pd.DataFrame:
    df = pd.read_csv(path)
    group = df[df.target_class.astype(str) == target_class].copy()
    if group.empty:
        raise ValueError(f"target class {target_class!r} not found in {path}")
    group = group.sort_values("T_K")
    numeric = [
        "K_init_target",
        "K_init_scale",
        "K_plateau_target",
        "K_plateau_scale",
        "early_rise_per_100um_target",
        "early_rise_scale",
        "plateau_rise_per_100um_target",
        "plateau_rise_scale",
        "delta_KR_min",
        "delta_KR_max",
        "weight",
    ]
    rows = []
    source_T = group.T_K.to_numpy(float)
    for T in temperatures:
        rec: dict[str, float | str] = {"target_class": target_class, "T_K": T}
        for col in numeric:
            rec[col] = float(np.interp(T, source_T, group[col].to_numpy(float)))
        rows.append(rec)
    return pd.DataFrame(rows)


def canonical_rows(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path).copy()
    df["candidate_id"] = "canonical_" + df.target_class.astype(str)
    df["candidate_source"] = "prior_first_passage_reference"
    return df


def select_diverse_seeds(
    atlas: pd.DataFrame,
    canonical: pd.DataFrame,
    target_class: str,
    count: int,
    pool_size: int,
) -> pd.DataFrame:
    score_col = {
        "ceramic": "ceramic_score",
        "weakT": "weakT_score",
        "DBTT": "DBTT_precursor_score",
    }[target_class]
    atlas = atlas.copy()
    atlas["candidate_source"] = "refined_analytic_atlas"
    if "candidate_id" not in atlas:
        atlas["candidate_id"] = [f"atlas_{i:06d}" for i in range(len(atlas))]
    required = [
        "cleave_G00_eV",
        "cleave_gT_eV_per_K",
        "cleave_sigc0_GPa",
        "emit_G00_eV",
        "emit_gT_eV_per_K",
        "emit_sigc0_GPa",
        *SHAPE_COLUMNS,
    ]
    atlas = atlas.dropna(subset=[c for c in required if c in atlas]).copy()
    if score_col in atlas:
        atlas = atlas.sort_values(score_col).head(max(pool_size, count))
    elif "region" in atlas:
        preferred = atlas.region.astype(str).str.contains(target_class, case=False)
        atlas = pd.concat([atlas[preferred], atlas[~preferred]], ignore_index=True)
        atlas = atlas.head(max(pool_size, count))
    else:
        atlas = atlas.head(max(pool_size, count))

    can = canonical[canonical.target_class.astype(str) == target_class].copy()
    pool = pd.concat([can, atlas], ignore_index=True, sort=False)
    pool = pool.drop_duplicates("candidate_id", keep="first").reset_index(drop=True)
    cols = [
        "cleave_G00_eV",
        "cleave_gT_eV_per_K",
        "cleave_sigc0_GPa",
        "emit_G00_eV",
        "emit_gT_eV_per_K",
        "emit_sigc0_GPa",
    ]
    X = pool[cols].to_numpy(float)
    lo = np.nanmin(X, axis=0)
    span = np.maximum(np.nanmax(X, axis=0) - lo, 1.0e-12)
    Z = (X - lo) / span
    selected = [0]
    while len(selected) < min(count, len(pool)):
        remaining = [i for i in range(len(pool)) if i not in selected]
        distances = [
            min(float(np.linalg.norm(Z[i] - Z[j])) for j in selected)
            for i in remaining
        ]
        selected.append(remaining[int(np.argmax(distances))])
    return pool.iloc[selected].reset_index(drop=True)


def shape_from_seed(row: pd.Series) -> dict[str, float]:
    defaults = {
        "cleave_sT_GPa_per_K": 0.0,
        "cleave_exp_a": 0.2,
        "cleave_exp_n": 1.0,
        "cleave_floor_frac": 0.02,
        "emit_sT_GPa_per_K": 0.0,
        "emit_exp_a": 0.2,
        "emit_exp_n": 1.0,
        "emit_floor_frac": 0.02,
    }
    return {name: float(row.get(name, default)) for name, default in defaults.items()}


def seed_vector(row: pd.Series) -> np.ndarray:
    emit0 = max(float(row.emit_G00_eV), 0.5)
    values = {
        "cleave_G00_eV": float(row.cleave_G00_eV),
        "cleave_gT_eV_per_K": float(row.cleave_gT_eV_per_K),
        "cleave_sigc0_GPa": float(row.cleave_sigc0_GPa),
        "emit_G00_eV": emit0,
        "emit_gT_eV_per_K": float(row.emit_gT_eV_per_K),
        "emit_sigc0_GPa": float(row.emit_sigc0_GPa),
        "peierls_H0_eV": min(max(0.5 * emit0, 0.02), 8.0),
        "delta_H_PT_eV": min(max(0.75 * emit0, 0.0), 12.0),
        "peierls_activation_entropy_kB": -20.0,
        "taylor_activation_entropy_kB": -20.0,
        "log10_taylor_corr_rho_c_m2": 14.0,
        "log10_taylor_corr_scale": 0.0,
        "log10_mobile_fraction": -2.0,
        "log10_source_sites_per_system": math.log10(
            max(float(row.get("mpz_source_sites_per_system", 200.0)), 1.0)
        ),
        "log10_recovery_rate_s": -5.0,
        "log10_source_refresh_length_um": math.log10(
            max(float(row.get("mpz_source_refresh_length_m", 2.5e-7)) * 1.0e6, 0.1)
        ),
        "c_blunt": float(row.get("c_blunt", 1.0)),
    }
    x = np.array([values[n] for n in PARAMETER_NAMES], dtype=float)
    bounds = np.asarray(bounds_array(), dtype=float)
    return np.clip(x, bounds[:, 0], bounds[:, 1])


def initial_population(
    x0: np.ndarray,
    popsize: int,
    seed: int,
) -> np.ndarray:
    dim = len(PARAMETER_NAMES)
    n = max(5, int(popsize) * dim)
    sampler = qmc.Sobol(dim, scramble=True, seed=seed)
    m = int(math.ceil(math.log2(n)))
    u = sampler.random_base2(m)[:n]
    bounds = np.asarray(bounds_array(), dtype=float)
    global_pop = qmc.scale(u, bounds[:, 0], bounds[:, 1])
    rng = np.random.default_rng(seed + 17)
    local_count = max(4, n // 2)
    scale = 0.12 * (bounds[:, 1] - bounds[:, 0])
    local = x0 + rng.normal(size=(local_count, dim)) * scale
    local = np.clip(local, bounds[:, 0], bounds[:, 1])
    global_pop[:local_count] = local
    global_pop[0] = x0
    return global_pop


@dataclass
class ObjectiveSettings:
    target_class: str
    temperatures: np.ndarray
    targets: pd.DataFrame
    shape: dict[str, float]
    Kdot: float = 0.005
    dK: float = 0.5
    Kmax: float = 65.0
    r_pz_m: float = 1.0e-6
    b_m: float = 2.74e-10
    G_shear_Pa: float = 160.0e9
    poisson: float = 0.28
    rho0_m2: float = 5.0e12
    nu0_c: float = 1.0e12
    nu0_e: float = 1.0e11
    cleavage_hits: float = 3.0
    cleavage_correlation_time_s: float = 1.0e-6
    Tref_K: float = 481.33
    n_systems: int = 2
    stress_order_max_GPa: float = 30.0
    stress_order_points: int = 31


class JointObjective:
    def __init__(self, settings: ObjectiveSettings):
        self.s = settings
        self.bounds = np.asarray(bounds_array(), dtype=float)

    def evaluate(self, x: np.ndarray, details: bool = False) -> dict[str, Any]:
        x = np.asarray(x, dtype=float)
        below = np.maximum(self.bounds[:, 0] - x, 0.0)
        above = np.maximum(x - self.bounds[:, 1], 0.0)
        bound_penalty = 1.0e5 * float(np.sum(below * below + above * above))
        if bound_penalty > 0.0:
            return {"objective": 1.0e6 + bound_penalty, "bound_penalty": bound_penalty}
        p = vector_to_parameters(x)
        model = build_pt_model(p, self.s.shape, Tref_K=self.s.Tref_K)
        stress_grid = np.linspace(
            0.0, self.s.stress_order_max_GPa * 1.0e9,
            self.s.stress_order_points,
        )
        order_margin = barrier_order_margin_eV(
            model, self.s.temperatures, stress_grid
        )
        raw_barriers = [
            model.raw_zero_stress_barrier_eV(mech, T)
            for mech in ("peierls", "taylor")
            for T in self.s.temperatures
        ]
        min_raw = float(np.min(raw_barriers))
        physical_penalty = 0.0
        if order_margin < -1.0e-9:
            physical_penalty += 1.0e4 + 1.0e5 * abs(order_margin)
        if min_raw <= 0.0:
            physical_penalty += 1.0e4 + 1.0e5 * abs(min_raw)
        if physical_penalty > 0.0:
            return {
                "objective": physical_penalty,
                "physical_penalty": physical_penalty,
                "barrier_order_margin_eV": order_margin,
                "min_raw_barrier_eV": min_raw,
            }

        fp = first_passage_proxy(
            p,
            self.s.shape,
            self.s.temperatures,
            Kdot=self.s.Kdot,
            dK=self.s.dK,
            Kmax=self.s.Kmax,
            r_pz_m=self.s.r_pz_m,
            nu0_c=self.s.nu0_c,
            nu0_e=self.s.nu0_e,
            cleavage_hits=self.s.cleavage_hits,
            cleavage_correlation_time_s=self.s.cleavage_correlation_time_s,
            Tref_K=self.s.Tref_K,
        )
        K0 = fp["K_intrinsic"]
        Hc = fp["H_emit_at_Kc"]
        if not np.all(np.isfinite(K0)) or not np.all(np.isfinite(Hc)):
            return {
                "objective": 5.0e4 + 100.0 * int(np.sum(~np.isfinite(K0))),
                "unresolved_first_passage": True,
                "barrier_order_margin_eV": order_margin,
                "min_raw_barrier_eV": min_raw,
            }

        n_sites = self.s.n_systems * p["source_sites_per_system"]
        K_per_line = (
            self.s.G_shear_Pa
            * self.s.b_m
            / max(1.0 - self.s.poisson, 1.0e-6)
            / math.sqrt(2.0 * math.pi * max(0.5 * self.s.r_pz_m, self.s.b_m))
            / 1.0e6
        )
        rows = []
        for j, T in enumerate(self.s.temperatures):
            sigma_tip = K0[j] * 1.0e6 / math.sqrt(
                2.0 * math.pi * self.s.r_pz_m
            )
            G_emit = float(model.barrier_eV("emission", sigma_tip, T))
            emit_rate = float(model._arrhenius_rate(G_emit, T, self.s.nu0_e))
            N = 0.0
            escape_rate = 0.0
            retention = 0.0
            rates: dict[str, Any] = {}
            for _ in range(24):
                rho = self.s.rho0_m2 + N / max(
                    math.pi * self.s.r_pz_m**2, 1.0e-30
                )
                rates = model.rates(sigma_tip, rho, T, self.s.b_m)
                escape_rate = float(np.asarray(rates["series_rate_s"]))
                denom = emit_rate + escape_rate + p["recovery_rate_s"]
                retention = emit_rate / denom if denom > 0.0 else 0.0
                Nnew = n_sites * retention
                if abs(Nnew - N) <= 1.0e-8 * max(Nnew, 1.0):
                    N = Nnew
                    break
                N = 0.5 * (N + Nnew)
            initial_emitted = n_sites * (1.0 - math.exp(-min(float(Hc[j]), 700.0)))
            N_init = min(initial_emitted * retention, N)
            N_dev = max(N, N_init)

            def resistance_increment(count: float) -> tuple[float, float]:
                shield = K_per_line * max(count, 0.0)
                r_eff = self.s.r_pz_m + p["c_blunt"] * self.s.b_m * max(count, 0.0)
                blunt = K0[j] * (math.sqrt(r_eff / self.s.r_pz_m) - 1.0)
                return shield, blunt

            shield_i, blunt_i = resistance_increment(N_init)
            shield_d, blunt_d = resistance_increment(N_dev)
            K_init = K0[j] + shield_i + blunt_i
            K_plateau = K0[j] + shield_d + blunt_d
            delta_KR = max(K_plateau - K_init, 0.0)
            L = max(p["source_refresh_length_um"], 1.0e-9)
            early = delta_KR * (math.exp(-20.0 / L) - math.exp(-120.0 / L))
            plateau_rise = delta_KR * (
                math.exp(-300.0 / L) - math.exp(-400.0 / L)
            )
            rows.append(
                {
                    "T_K": float(T),
                    "K_intrinsic": float(K0[j]),
                    "K_init_proxy": K_init,
                    "K_plateau_proxy": K_plateau,
                    "delta_KR_proxy": delta_KR,
                    "early_rise_per_100um_proxy": early,
                    "plateau_rise_per_100um_proxy": plateau_rise,
                    "N_initial_proxy": N_init,
                    "N_developed_proxy": N_dev,
                    "retained_fraction_proxy": retention,
                    "emission_rate_s": emit_rate,
                    "escape_rate_s": escape_rate,
                    "peierls_rate_s": float(np.asarray(rates["peierls_rate_s"])),
                    "taylor_completion_rate_s": float(
                        np.asarray(rates["taylor_completion_rate_s"])
                    ),
                    "barrier_bottleneck_log10_P_over_T": math.log10(
                        max(float(np.asarray(rates["peierls_rate_s"])), 1.0e-300)
                        / max(
                            float(np.asarray(rates["taylor_completion_rate_s"])),
                            1.0e-300,
                        )
                    ),
                }
            )

        predictions = pd.DataFrame(rows)
        merged = self.s.targets.merge(predictions, on="T_K", how="left")
        components = {
            "K_init_loss": 0.0,
            "K_plateau_loss": 0.0,
            "early_rise_loss": 0.0,
            "plateau_rise_loss": 0.0,
            "delta_window_loss": 0.0,
        }
        for _, row in merged.iterrows():
            w = float(row.weight)
            components["K_init_loss"] += w * huber(
                (row.K_init_proxy - row.K_init_target)
                / max(row.K_init_scale, 1.0e-9)
            )
            components["K_plateau_loss"] += w * huber(
                (row.K_plateau_proxy - row.K_plateau_target)
                / max(row.K_plateau_scale, 1.0e-9)
            )
            components["early_rise_loss"] += w * huber(
                (row.early_rise_per_100um_proxy - row.early_rise_per_100um_target)
                / max(row.early_rise_scale, 1.0e-9)
            )
            components["plateau_rise_loss"] += w * huber(
                (row.plateau_rise_per_100um_proxy - row.plateau_rise_per_100um_target)
                / max(row.plateau_rise_scale, 1.0e-9)
            )
            if row.delta_KR_proxy < row.delta_KR_min:
                components["delta_window_loss"] += w * huber(
                    row.delta_KR_min - row.delta_KR_proxy
                )
            elif row.delta_KR_proxy > row.delta_KR_max:
                components["delta_window_loss"] += w * huber(
                    row.delta_KR_proxy - row.delta_KR_max
                )
        objective = float(sum(components.values()))
        result: dict[str, Any] = {
            "objective": objective,
            **components,
            "physical_penalty": physical_penalty,
            "barrier_order_margin_eV": order_margin,
            "min_raw_barrier_eV": min_raw,
            "parameters": p,
        }
        if details:
            result["temperature_detail"] = merged.to_dict(orient="records")
        return result

    def __call__(self, x: np.ndarray) -> float:
        return float(self.evaluate(x, details=False)["objective"])


def diverse_shortlist(df: pd.DataFrame, count: int) -> pd.DataFrame:
    if len(df) <= count:
        return df.copy().reset_index(drop=True)
    parameter_cols = [c for c in PARAMETER_NAMES if c in df]
    pool = df.sort_values("objective").head(max(count * 8, count)).copy()
    X = pool[parameter_cols].to_numpy(float)
    bounds = np.asarray(bounds_array(), dtype=float)
    Z = (X - bounds[:, 0]) / np.maximum(bounds[:, 1] - bounds[:, 0], 1.0e-12)
    selected = [0]
    while len(selected) < min(count, len(pool)):
        remaining = [i for i in range(len(pool)) if i not in selected]
        utility = []
        for i in remaining:
            distance = min(np.linalg.norm(Z[i] - Z[j]) for j in selected)
            score_penalty = 0.05 * (
                pool.iloc[i].objective - pool.iloc[0].objective
            )
            utility.append(distance - score_penalty)
        selected.append(remaining[int(np.argmax(utility))])
    return pool.iloc[selected].sort_values("objective").reset_index(drop=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--target-class", choices=["ceramic", "weakT", "DBTT"], required=True)
    ap.add_argument(
        "--prepared-atlas",
        type=Path,
        default=Path(
            "runs/mpz_v9_4_peierls_taylor_search_v1/pt_search_input_joined.csv"
        ),
    )
    ap.add_argument(
        "--atlas",
        type=Path,
        default=Path(
            "runs/mpz_v9_2_analytic_first_passage_atlas/"
            "analytic_first_passage_atlas_shortlist_refined.csv"
        ),
    )
    ap.add_argument(
        "--canonical",
        type=Path,
        default=Path("mpz_v9_6_canonical_first_passage_references.csv"),
    )
    ap.add_argument(
        "--targets", type=Path, default=Path("mpz_three_class_design_targets.csv")
    )
    ap.add_argument("--temperatures", default="300 700 900 1200")
    ap.add_argument("--seed-count", type=int, default=6)
    ap.add_argument("--seed-pool-size", type=int, default=120)
    ap.add_argument("--de-maxiter", type=int, default=40)
    ap.add_argument("--de-popsize", type=int, default=8)
    ap.add_argument("--de-tol", type=float, default=1.0e-3)
    ap.add_argument("--local-maxiter", type=int, default=300)
    ap.add_argument("--local-xatol", type=float, default=1.0e-4)
    ap.add_argument("--local-fatol", type=float, default=1.0e-4)
    ap.add_argument("--seed", type=int, default=98017)
    ap.add_argument("--Kdot", type=float, default=0.005)
    ap.add_argument("--dK", type=float, default=0.5)
    ap.add_argument("--Kmax", type=float, default=65.0)
    ap.add_argument("--r-pz-m", type=float, default=1.0e-6)
    ap.add_argument("--shortlist-count", type=int, default=20)
    ap.add_argument("--resume", action="store_true", default=True)
    ap.add_argument(
        "--out", type=Path, default=Path("runs/mpz_v9_8_joint_response_optimization")
    )
    a = ap.parse_args()

    source = a.prepared_atlas if a.prepared_atlas.exists() else a.atlas
    if not source.exists():
        raise SystemExit(f"atlas not found: {source}")
    atlas = pd.read_csv(source)
    canonical = canonical_rows(a.canonical)
    seeds = select_diverse_seeds(
        atlas, canonical, a.target_class, a.seed_count, a.seed_pool_size
    )
    temperatures = parse_floats(a.temperatures)
    targets = load_targets(a.targets, a.target_class, temperatures)
    out = (a.out / a.target_class).resolve()
    checkpoints = out / "checkpoints"
    checkpoints.mkdir(parents=True, exist_ok=True)
    seeds.to_csv(out / "selected_atlas_seeds.csv", index=False)

    basin_rows: list[dict[str, Any]] = []
    detail_rows: list[dict[str, Any]] = []
    history_rows: list[dict[str, Any]] = []
    for basin_index, seed_row in seeds.iterrows():
        checkpoint = checkpoints / f"basin_{basin_index:03d}.json"
        if a.resume and checkpoint.exists():
            saved = json.loads(checkpoint.read_text())
            basin_rows.append(saved["summary"])
            detail_rows.extend(saved.get("temperature_detail", []))
            history_rows.extend(saved.get("history", []))
            print(f"resumed basin {basin_index}/{len(seeds)-1}", flush=True)
            continue
        shape = shape_from_seed(seed_row)
        settings = ObjectiveSettings(
            target_class=a.target_class,
            temperatures=np.asarray(temperatures, dtype=float),
            targets=targets,
            shape=shape,
            Kdot=a.Kdot,
            dK=a.dK,
            Kmax=a.Kmax,
            r_pz_m=a.r_pz_m,
        )
        objective = JointObjective(settings)
        x0 = seed_vector(seed_row)
        init = initial_population(x0, a.de_popsize, a.seed + 1009 * basin_index)
        history: list[dict[str, Any]] = []

        def callback(xk: np.ndarray, convergence: float) -> bool:
            score = objective(xk)
            rec = {
                "basin_index": int(basin_index),
                "generation": len(history),
                "objective": float(score),
                "convergence": float(convergence),
            }
            history.append(rec)
            checkpoint.write_text(
                json.dumps(
                    {
                        "status": "RUNNING",
                        "basin_index": int(basin_index),
                        "seed_candidate_id": str(seed_row.candidate_id),
                        "history": history,
                    },
                    indent=2,
                )
            )
            print(
                f"class={a.target_class} basin={basin_index} "
                f"generation={len(history)} objective={score:.6g}",
                flush=True,
            )
            return False

        de = differential_evolution(
            objective,
            bounds_array(),
            maxiter=a.de_maxiter,
            popsize=a.de_popsize,
            tol=a.de_tol,
            seed=a.seed + basin_index,
            init=init,
            polish=False,
            updating="immediate",
            workers=1,
            callback=callback,
        )
        local = minimize(
            objective,
            de.x,
            method="Nelder-Mead",
            bounds=bounds_array(),
            options={
                "maxiter": a.local_maxiter,
                "xatol": a.local_xatol,
                "fatol": a.local_fatol,
                "adaptive": True,
            },
        )
        best_x = local.x if local.fun <= de.fun else de.x
        best = objective.evaluate(best_x, details=True)
        p = best.pop("parameters")
        temp_detail = best.pop("temperature_detail")
        summary = {
            "target_class": a.target_class,
            "basin_index": int(basin_index),
            "seed_candidate_id": str(seed_row.candidate_id),
            "seed_candidate_source": str(seed_row.get("candidate_source", "")),
            "seed_region": str(seed_row.get("region", "")),
            "objective": float(best["objective"]),
            "de_objective": float(de.fun),
            "local_objective": float(local.fun),
            "de_success": bool(de.success),
            "local_success": bool(local.success),
            **{name: float(best_x[i]) for i, name in enumerate(PARAMETER_NAMES)},
            **{k: float(v) for k, v in p.items() if k not in PARAMETER_NAMES},
            **{k: v for k, v in best.items() if k != "objective"},
            **{f"shape_{k}": float(v) for k, v in shape.items()},
            "status": "ANALYTICAL_FIDELITY_0_REQUIRES_REDUCED_MPZ_PROMOTION",
        }
        for row in temp_detail:
            row.update(
                {
                    "target_class": a.target_class,
                    "basin_index": int(basin_index),
                    "seed_candidate_id": str(seed_row.candidate_id),
                    "objective": float(best["objective"]),
                }
            )
        payload = {
            "status": "COMPLETE",
            "summary": summary,
            "temperature_detail": temp_detail,
            "history": history,
        }
        checkpoint.write_text(json.dumps(payload, indent=2, allow_nan=True))
        basin_rows.append(summary)
        detail_rows.extend(temp_detail)
        history_rows.extend(history)

    basin_df = pd.DataFrame(basin_rows).sort_values("objective").reset_index(drop=True)
    shortlist = diverse_shortlist(basin_df, a.shortlist_count)
    detail_df = pd.DataFrame(detail_rows)
    history_df = pd.DataFrame(history_rows)
    basin_df.to_csv(out / "joint_response_basin_results.csv", index=False)
    shortlist.to_csv(out / "joint_response_shortlist.csv", index=False)
    detail_df.to_csv(out / "joint_response_temperature_detail.csv", index=False)
    history_df.to_csv(out / "joint_response_generation_history.csv", index=False)
    shortlist.to_csv(out / "joint_response_promotion_manifest.csv", index=False)

    config = vars(a).copy()
    for key, value in list(config.items()):
        if isinstance(value, Path):
            config[key] = str(value)
    config.update(
        {
            "parameter_names": PARAMETER_NAMES,
            "bounds": DEFAULT_BOUNDS,
            "physics_constraints": [
                "G_P(sigma,T) <= G_T(sigma,T) on common resolved-stress grid",
                "positive raw zero-stress Peierls and Taylor barriers",
                "exact detailed balance inherited from v9.7",
                "no constitutive caps or algebraic saturation functions",
            ],
            "mechanism_dominance_assumptions": [],
            "fidelity": "ANALYTICAL_0",
        }
    )
    (out / "joint_response_optimization_config.json").write_text(
        json.dumps(config, indent=2)
    )
    report = {
        "target_class": a.target_class,
        "n_basins": int(len(basin_df)),
        "best_objective": float(basin_df.iloc[0].objective),
        "best_seed_candidate_id": str(basin_df.iloc[0].seed_candidate_id),
        "shortlist_count": int(len(shortlist)),
        "output": str(out),
        "status": "ANALYTICAL_FIDELITY_0_COMPLETE",
    }
    (out / "joint_response_optimization_summary.json").write_text(
        json.dumps(report, indent=2)
    )
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
