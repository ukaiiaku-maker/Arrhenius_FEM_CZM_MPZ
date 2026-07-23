"""PF v10.2.22 persistent-site signed moving-process-zone support.

This module contains the mechanically derived signed-kernel contract, the
mesh-independent along-front width, moving-frame translation helpers, and the
audited backstress-complementarity emission solver.
"""
from __future__ import annotations

from dataclasses import dataclass
import json
import math
from pathlib import Path
from typing import Any, Callable

import numpy as np

MODEL_ID = "FEM_CZM_persistent_signed_MPZ_PF_v10_2_22_parity_v10_0_5_14"
KERNEL_SCHEMA = "v10.0.5.14_signed_shielding_kernel"


@dataclass(frozen=True)
class SignedShieldingKernelV100514:
    active_kernel_Pa_sqrt_m_per_signed_line: np.ndarray
    wake_kernel_Pa_sqrt_m_per_signed_line: np.ndarray
    activation_to_line_content_by_system: np.ndarray
    metadata: dict[str, Any]
    source_path: str = ""
    active_x_m: np.ndarray | None = None
    wake_x_m: np.ndarray | None = None

    @classmethod
    def from_json(cls, path: str | Path) -> "SignedShieldingKernelV100514":
        source = Path(path).expanduser().resolve()
        payload = json.loads(source.read_text())
        schema = payload.get("schema")
        if schema not in {KERNEL_SCHEMA, "v10.2.5_2d_unit_signed_shielding_kernel"}:
            raise ValueError(f"unsupported signed-kernel schema: {schema!r}")
        active = np.asarray(
            payload["active_kernel_Pa_sqrt_m_per_signed_line"], dtype=float
        )
        wake = np.asarray(
            payload.get("wake_kernel_Pa_sqrt_m_per_signed_line", []), dtype=float
        )
        conversion = np.asarray(
            payload["activation_to_line_content_by_system"], dtype=float
        ).reshape(-1)
        if active.ndim != 2 or not np.all(np.isfinite(active)):
            raise ValueError("active signed shielding kernel must be a finite 2-D array")
        if wake.size == 0:
            wake = np.zeros((active.shape[0], active.shape[1]), dtype=float)
        if wake.ndim != 2 or not np.all(np.isfinite(wake)):
            raise ValueError("wake signed shielding kernel must be a finite 2-D array")
        if conversion.shape != (active.shape[0],) or np.any(conversion <= 0.0):
            raise ValueError(
                "one positive activation-to-line conversion is required per system"
            )
        required = {
            "candidate_independent": True,
            "counts_are_signed_burgers_lines": True,
            "normalization_is_mechanically_derived": True,
        }
        for key, expected in required.items():
            if payload.get(key) is not expected:
                raise ValueError(f"signed kernel requires {key}={expected}")
        active_x = payload.get("active_x_m")
        wake_x = payload.get("wake_x_m")
        active_x_array = (
            None if active_x is None else np.asarray(active_x, dtype=float).reshape(-1)
        )
        wake_x_array = (
            None if wake_x is None else np.asarray(wake_x, dtype=float).reshape(-1)
        )
        metadata = {
            key: value
            for key, value in payload.items()
            if key
            not in {
                "active_kernel_Pa_sqrt_m_per_signed_line",
                "wake_kernel_Pa_sqrt_m_per_signed_line",
                "activation_to_line_content_by_system",
                "active_x_m",
                "wake_x_m",
            }
        }
        return cls(
            active,
            wake,
            conversion,
            metadata,
            str(source),
            active_x_array,
            wake_x_array,
        )

    def validate(self, n_systems: int, n_bins: int) -> None:
        shape = (int(n_systems), int(n_bins))
        if self.active_kernel_Pa_sqrt_m_per_signed_line.shape != shape:
            raise ValueError(
                "active signed shielding kernel shape "
                f"{self.active_kernel_Pa_sqrt_m_per_signed_line.shape} != {shape}"
            )
        if self.wake_kernel_Pa_sqrt_m_per_signed_line.shape[0] != int(n_systems):
            raise ValueError("wake signed shielding kernel system count mismatch")
        if self.activation_to_line_content_by_system.shape != (int(n_systems),):
            raise ValueError("activation-to-line conversion system count mismatch")
        if self.active_x_m is not None:
            active_x = np.asarray(self.active_x_m, dtype=float).reshape(-1)
            if active_x.shape != (int(n_bins),) or not np.all(np.isfinite(active_x)):
                raise ValueError("active_x_m must match the active MPZ bins")
        if self.wake_x_m is not None:
            wake_x = np.asarray(self.wake_x_m, dtype=float).reshape(-1)
            if wake_x.shape != (
                self.wake_kernel_Pa_sqrt_m_per_signed_line.shape[1],
            ):
                raise ValueError("wake_x_m must match the wake-kernel bins")
            if not np.all(np.isfinite(wake_x)):
                raise ValueError("wake_x_m must be finite")


def effective_front_width_m(
    rho_unsigned_m2: float,
    *,
    reference_width_m: float,
    reference_density_m2: float,
    minimum_physical_width_m: float,
    burgers_m: float,
    maximum_width_m: float,
) -> float:
    """PF v10.2.22 along-front width; ahead-of-tip dx is absent."""
    rho = max(float(rho_unsigned_m2), float(reference_density_m2), 1.0)
    width = float(reference_width_m) * math.sqrt(float(reference_density_m2) / rho)
    lower = max(float(minimum_physical_width_m), abs(float(burgers_m)), 1.0e-30)
    upper = max(float(maximum_width_m), lower)
    return min(max(width, lower), upper)


def persistent_site_multiplicity(
    rho_site0_m2: float,
    active_arc_factor: float,
    tip_radius_m: float,
    front_width_m: float,
) -> float:
    area = (
        max(float(active_arc_factor), 0.0)
        * max(float(tip_radius_m), 0.0)
        * max(float(front_width_m), 0.0)
    )
    return max(float(rho_site0_m2), 0.0) * area


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
    """Return activation increment, mechanical block increment, and block flag."""
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

    hi_inside = math.nextafter(upper, 0.0)
    if hi_inside <= 0.0:
        return upper, upper, True
    r_hi_inside = residual(hi_inside)
    if r_hi_inside <= 0.0:
        return upper, upper, True
    lo = 0.0
    hi = min(hi_inside, M * rate0 * dt)
    if residual(hi) < 0.0:
        hi = hi_inside
    scale = max(abs(hi), 1.0)
    if abs(residual(hi)) <= tol * scale:
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


def _translate_toward_tip(
    field: np.ndarray,
    distance_m: float,
    dx: float,
    wake_bins: int,
    wake_dx: float,
) -> tuple[np.ndarray, np.ndarray, float]:
    source = np.asarray(field, dtype=float)
    active = np.zeros_like(source)
    wake = np.zeros((source.shape[0], int(wake_bins)), dtype=float)
    shift = max(float(distance_m), 0.0) / max(float(dx), 1.0e-30)
    for i in range(source.shape[1]):
        dest = float(i) - shift
        j0 = math.floor(dest)
        frac = dest - j0
        for j, weight in ((j0, 1.0 - frac), (j0 + 1, frac)):
            if weight <= 0.0:
                continue
            mass = source[:, i] * weight
            if 0 <= j < source.shape[1]:
                active[:, j] += mass
            elif j < 0:
                y = (-float(j) - 0.5) * float(dx)
                k = max(int(y / max(float(wake_dx), 1.0e-30)), 0)
                if k < int(wake_bins):
                    wake[:, k] += mass
    total = float(np.sum(source))
    discarded = max(total - float(np.sum(active)) - float(np.sum(wake)), 0.0)
    return active, wake, discarded


def _shift_wake_forward(
    field: np.ndarray, distance_m: float, dx: float
) -> tuple[np.ndarray, float]:
    source = np.asarray(field, dtype=float)
    out = np.zeros_like(source)
    shift = max(float(distance_m), 0.0) / max(float(dx), 1.0e-30)
    for i in range(source.shape[1]):
        dest = float(i) + shift
        j0 = math.floor(dest)
        frac = dest - j0
        for j, weight in ((j0, 1.0 - frac), (j0 + 1, frac)):
            if weight > 0.0 and 0 <= j < source.shape[1]:
                out[:, j] += source[:, i] * weight
    return out, max(float(np.sum(source)) - float(np.sum(out)), 0.0)


__all__ = [
    "MODEL_ID",
    "KERNEL_SCHEMA",
    "SignedShieldingKernelV100514",
    "effective_front_width_m",
    "persistent_site_multiplicity",
    "solve_backstress_limited_activations",
    "_translate_toward_tip",
    "_shift_wake_forward",
]
