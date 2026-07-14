#!/usr/bin/env python3
"""v9.4 detailed-balance wrapper for the v9.3 PT parameter search.

The v9.3 search architecture is retained, but the sampled entropy range and
acceptance logic are corrected for signed forward-minus-reverse kinetics.
"""
from __future__ import annotations

import math

import numpy as np
import pandas as pd
from scipy.stats import qmc

import search_mpz_peierls_taylor_parameters as _legacy


def _log_scale(u: np.ndarray, lo: float, hi: float) -> np.ndarray:
    return 10.0 ** (
        math.log10(lo) + u * (math.log10(hi) - math.log10(lo))
    )


def sample_transport_parameters(n: int, seed: int) -> pd.DataFrame:
    """Sample the post-emission closure without invalid low-T barriers."""
    sampler = qmc.Sobol(10, scramble=True, seed=seed)
    if n > 0 and (n & (n - 1)) == 0:
        u = sampler.random_base2(int(math.log2(n)))
    else:
        u = sampler.random(n)
    sampled = pd.DataFrame({
        "pt_peierls_energy_ratio": _log_scale(u[:, 0], 0.0025, 0.010),
        "pt_taylor_energy_ratio": _log_scale(u[:, 1], 0.010, 0.040),
        "pt_entropy_multiplier": _log_scale(u[:, 2], 0.25, 8.0),
        "pt_taylor_corr_rho_c": _log_scale(u[:, 3], 1.0e9, 1.0e14),
        "pt_taylor_renewal_time_s": _log_scale(u[:, 4], 1.0e-18, 1.0e-8),
        "pt_taylor_m_exponent": 0.5 + 1.5 * u[:, 5],
        "pt_taylor_m_scale": _log_scale(u[:, 6], 0.10, 10.0),
        "pt_taylor_m_cap": 6.0 + 42.0 * u[:, 7],
        "pt_mobile_saturation_density_m2": _log_scale(
            u[:, 8], 1.0e12, 1.0e16
        ),
        "pt_mobile_fraction": _log_scale(u[:, 9], 1.0e-4, 5.0e-2),
    })
    anchors = []
    for entropy_mult in (0.5, 1.0, 2.0, 4.0):
        for m_cap in (15.0, 22.0, 30.0):
            anchors.append({
                "pt_peierls_energy_ratio": 0.005,
                "pt_taylor_energy_ratio": 0.020,
                "pt_entropy_multiplier": entropy_mult,
                "pt_taylor_corr_rho_c": 1.0e11,
                "pt_taylor_renewal_time_s": 1.0e-10,
                "pt_taylor_m_exponent": 1.0,
                "pt_taylor_m_scale": 1.0,
                "pt_taylor_m_cap": m_cap,
                "pt_mobile_saturation_density_m2": 1.0e14,
                "pt_mobile_fraction": 0.01,
            })
    return pd.concat([pd.DataFrame(anchors), sampled], ignore_index=True)


def evaluate_one(
    model,
    rho,
    temperatures,
    strain_rates,
    b,
    min_reference_stress_GPa,
    max_reference_stress_GPa,
    max_stress_GPa,
    slope_tol_GPa_decade,
    drop_tol_fraction,
    zero_stress_threshold_GPa,
):
    out = _legacy.evaluate_one(
        model,
        rho,
        temperatures,
        strain_rates,
        b,
        min_reference_stress_GPa,
        max_reference_stress_GPa,
        max_stress_GPa,
        slope_tol_GPa_decade,
        drop_tol_fraction,
        zero_stress_threshold_GPa,
    )
    raw = [
        model.raw_zero_stress_barrier_eV(mechanism, T)
        for mechanism in ("peierls", "taylor")
        for T in temperatures
    ]
    min_raw = float(np.min(raw))
    barrier_valid = bool(min_raw > 1.0e-8)
    zero_rate = 0.0
    for T in temperatures:
        rates = model.rates(0.0, rho, T, b)
        zero_rate = max(
            zero_rate,
            float(np.max(np.abs(rates["equivalent_plastic_rate_s"]))),
        )
    detailed_balance_valid = bool(zero_rate <= 1.0e-20)

    out["barrier_valid"] = barrier_valid
    out["detailed_balance_valid"] = detailed_balance_valid
    out["min_raw_scaled_G0_eV"] = min_raw
    out["zero_stress_rate_max_s"] = zero_rate
    out["accepted"] = bool(
        out["accepted"] and barrier_valid and detailed_balance_valid
    )
    out["strict_strength_window"] = bool(
        out["strict_strength_window"] and out["accepted"]
    )
    if not barrier_valid:
        out["pt_screen_status"] = "invalid_scaled_zero_stress_barrier"
    elif not detailed_balance_valid:
        out["pt_screen_status"] = "zero_stress_ratchet"
    elif out["strict_strength_window"]:
        out["pt_screen_status"] = "strict_strength_window"
    elif out["accepted"]:
        out["pt_screen_status"] = "monotonic_topology_only"
    return out


def main() -> None:
    _legacy.sample_transport_parameters = sample_transport_parameters
    _legacy.evaluate_one = evaluate_one
    _legacy.main()


if __name__ == "__main__":
    main()
