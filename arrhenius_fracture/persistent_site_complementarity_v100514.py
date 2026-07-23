"""Numerically robust PF v10.2.22 mechanical-blocking complementarity."""
from __future__ import annotations

import math
from typing import Callable


def solve_backstress_limited_activations(
    *,
    multiplicity: float,
    dt_s: float,
    drive_stress_Pa: float,
    rho_initial_m2: float,
    rho_increment_per_activation_m2: float,
    backstress_prefactor_Pa_sqrt_m2: float,
    rate_function: Callable[[float], float],
    tolerance: float = 1.0e-10,
    max_iterations: int = 96,
) -> tuple[float, float, bool]:
    """Return activation increment, mechanical block increment, and block flag.

    The admissible interval endpoint is evaluated at a tolerance-scaled interior
    point. A one-ULP retreat in density can still round through ``sqrt`` to the
    exact blocking stress and incorrectly hide a complementarity solution.
    """
    M = max(float(multiplicity), 0.0)
    dt = max(float(dt_s), 0.0)
    drive = max(float(drive_stress_Pa), 0.0)
    rho0 = max(float(rho_initial_m2), 0.0)
    rho_per = max(float(rho_increment_per_activation_m2), 0.0)
    kback = max(float(backstress_prefactor_Pa_sqrt_m2), 0.0)
    tol = max(float(tolerance), 1.0e-15)
    if M <= 0.0 or dt <= 0.0 or drive <= 0.0:
        return 0.0, 0.0, False
    if rho_per <= 0.0 or kback <= 0.0:
        raise RuntimeError(
            "persistent-site emission requires positive backstress coupling"
        )
    sigma0 = drive - kback * math.sqrt(rho0)
    if sigma0 <= 0.0:
        return 0.0, 0.0, False
    rate0 = max(float(rate_function(sigma0)), 0.0)
    if not math.isfinite(rate0) or rate0 <= 0.0:
        return 0.0, 0.0, False
    rho_block = (drive / kback) ** 2
    upper = max((rho_block - rho0) / rho_per, 0.0)
    if upper <= 0.0:
        return 0.0, 0.0, False

    def residual(value: float) -> float:
        rho = rho0 + rho_per * max(float(value), 0.0)
        sigma_eff = drive - kback * math.sqrt(max(rho, 0.0))
        rate = (
            0.0
            if sigma_eff <= 0.0
            else max(float(rate_function(sigma_eff)), 0.0)
        )
        if not math.isfinite(rate):
            rate = 0.0
        return float(value) - M * rate * dt

    scale = max(abs(upper), 1.0)
    retreat = max(tol * scale, 16.0 * math.ulp(upper))
    hi_inside = max(upper - retreat, 0.0)
    if hi_inside <= 0.0:
        return upper, upper, True
    r_hi_inside = residual(hi_inside)
    if r_hi_inside <= 0.0:
        return upper, upper, True

    lo = 0.0
    hi = min(hi_inside, M * rate0 * dt)
    if residual(hi) < 0.0:
        hi = hi_inside
    hi_scale = max(abs(hi), 1.0)
    if abs(residual(hi)) <= tol * hi_scale:
        return min(max(hi, 0.0), upper), upper, False
    for _ in range(int(max_iterations)):
        mid = 0.5 * (lo + hi)
        value = residual(mid)
        root_scale = max(abs(mid), 1.0)
        if abs(value) <= tol * root_scale:
            return min(max(mid, 0.0), upper), upper, False
        if value > 0.0:
            hi = mid
        else:
            lo = mid
        interval_scale = max(abs(lo), abs(hi), 1.0)
        if (hi - lo) <= tol * interval_scale:
            root = 0.5 * (lo + hi)
            return min(max(root, 0.0), upper), upper, False
    root = 0.5 * (lo + hi)
    return min(max(root, 0.0), upper), upper, False


__all__ = ["solve_backstress_limited_activations"]
