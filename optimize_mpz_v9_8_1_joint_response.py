#!/usr/bin/env python3
"""v9.8.1 robustness wrapper for the joint-response optimizer.

The first v9.8 DBTT preflight exposed a data-joining bug: canonical rows acquire
atlas-only columns filled with NaN, so ``Series.get(name, default)`` returns NaN
rather than the intended default. Those NaNs entered the seed vector through
source inventory, source-refresh length, and blunting, then contaminated the
entire differential-evolution population and objective.

This wrapper preserves the v9.8 physics and objective while enforcing finite
seed defaults, finite initial populations, and finite penalty values for every
optimizer evaluation.
"""
from __future__ import annotations

import math
from typing import Any

import numpy as np
import pandas as pd

import optimize_mpz_v9_8_joint_response as core

NONFINITE_PENALTY = 1.0e9


def finite_row_value(row: pd.Series, name: str, default: float) -> float:
    value: Any = row.get(name, default)
    try:
        value = float(value)
    except (TypeError, ValueError):
        return float(default)
    return value if np.isfinite(value) else float(default)


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
    return {
        name: finite_row_value(row, name, default)
        for name, default in defaults.items()
    }


def seed_vector(row: pd.Series) -> np.ndarray:
    emit0 = max(finite_row_value(row, "emit_G00_eV", 1.5), 0.5)
    source_sites = max(
        finite_row_value(row, "mpz_source_sites_per_system", 200.0), 1.0
    )
    refresh_m = max(
        finite_row_value(row, "mpz_source_refresh_length_m", 2.5e-7),
        1.0e-7,
    )
    values = {
        "cleave_G00_eV": finite_row_value(row, "cleave_G00_eV", 2.0),
        "cleave_gT_eV_per_K": finite_row_value(row, "cleave_gT_eV_per_K", 0.0),
        "cleave_sigc0_GPa": finite_row_value(row, "cleave_sigc0_GPa", 4.0),
        "emit_G00_eV": emit0,
        "emit_gT_eV_per_K": finite_row_value(row, "emit_gT_eV_per_K", 0.0),
        "emit_sigc0_GPa": finite_row_value(row, "emit_sigc0_GPa", 2.5),
        "peierls_H0_eV": min(max(0.5 * emit0, 0.02), 8.0),
        "delta_H_PT_eV": min(max(0.75 * emit0, 0.0), 12.0),
        "peierls_activation_entropy_kB": -20.0,
        "taylor_activation_entropy_kB": -20.0,
        "log10_taylor_corr_rho_c_m2": 14.0,
        "log10_taylor_corr_scale": 0.0,
        "log10_mobile_fraction": -2.0,
        "log10_source_sites_per_system": math.log10(source_sites),
        "log10_recovery_rate_s": -5.0,
        "log10_source_refresh_length_um": math.log10(refresh_m * 1.0e6),
        "c_blunt": finite_row_value(row, "c_blunt", 1.0),
    }
    x = np.array([values[name] for name in core.PARAMETER_NAMES], dtype=float)
    bounds = np.asarray(core.bounds_array(), dtype=float)
    x = np.clip(x, bounds[:, 0], bounds[:, 1])
    if not np.all(np.isfinite(x)):
        bad = [
            core.PARAMETER_NAMES[i]
            for i in np.flatnonzero(~np.isfinite(x))
        ]
        raise ValueError(f"non-finite seed parameters after defaults: {bad}")
    return x


_original_initial_population = core.initial_population


def initial_population(x0: np.ndarray, popsize: int, seed: int) -> np.ndarray:
    x0 = np.asarray(x0, dtype=float)
    if not np.all(np.isfinite(x0)):
        raise ValueError("initial seed vector contains NaN or infinity")
    population = _original_initial_population(x0, popsize, seed)
    if not np.all(np.isfinite(population)):
        raise ValueError("generated differential-evolution population is non-finite")
    return population


class JointObjective(core.JointObjective):
    def evaluate(self, x: np.ndarray, details: bool = False) -> dict[str, Any]:
        x = np.asarray(x, dtype=float)
        if not np.all(np.isfinite(x)):
            return {
                "objective": NONFINITE_PENALTY,
                "nonfinite_parameter_vector": True,
            }
        try:
            result = super().evaluate(x, details=details)
        except (FloatingPointError, OverflowError, ValueError, ZeroDivisionError) as exc:
            return {
                "objective": NONFINITE_PENALTY,
                "evaluation_exception": f"{type(exc).__name__}: {exc}",
            }
        objective = float(result.get("objective", NONFINITE_PENALTY))
        if not np.isfinite(objective):
            return {
                "objective": NONFINITE_PENALTY,
                "nonfinite_objective_replaced": True,
            }
        result["objective"] = objective
        return result

    def __call__(self, x: np.ndarray) -> float:
        value = float(self.evaluate(x, details=False)["objective"])
        return value if np.isfinite(value) else NONFINITE_PENALTY


core.shape_from_seed = shape_from_seed
core.seed_vector = seed_vector
core.initial_population = initial_population
core.JointObjective = JointObjective


if __name__ == "__main__":
    core.main()
