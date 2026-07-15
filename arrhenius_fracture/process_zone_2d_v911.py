"""Extract dimensionless 2-D crack-tip profiles for the v9.11 MPZ coupling.

The FEM supplies two kinds of information:

1. Plastic redistribution is already contained in the domain-integral K/J drive.
2. The finite-radius stress and scalar forest-density fields supply only spatial
   *shape* and local forest density to the unresolved 1-D process-zone model.

Raw finite-radius FEM stress is never substituted for the calibrated sharp-tip
magnitude. Scalar rho is never converted into signed shielding or backstress;
that requires slip-system-resolved GND/Nye-tensor state that this 2-D model does
not currently carry.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class ProcessZone2DProfile:
    x_m: np.ndarray
    forest_density_m2: np.ndarray
    stress_shape: np.ndarray
    reliable: bool
    coverage_fraction: float
    selected_elements: int
    reason: str

    def diagnostics(self) -> dict[str, float | bool | str]:
        return {
            "mpz_2d_profile_reliable": bool(self.reliable),
            "mpz_2d_profile_coverage_fraction": float(self.coverage_fraction),
            "mpz_2d_profile_selected_elements": int(self.selected_elements),
            "mpz_2d_profile_reason": str(self.reason),
            "mpz_2d_rho_min_m2": float(np.min(self.forest_density_m2)),
            "mpz_2d_rho_median_m2": float(np.median(self.forest_density_m2)),
            "mpz_2d_rho_max_m2": float(np.max(self.forest_density_m2)),
            "mpz_2d_stress_shape_min": float(np.min(self.stress_shape)),
            "mpz_2d_stress_shape_max": float(np.max(self.stress_shape)),
            "bulk_scalar_rho_used_for_signed_shielding": False,
        }


def _unit(v) -> np.ndarray:
    q = np.asarray(v, dtype=float).reshape(2)
    n = float(np.linalg.norm(q))
    return np.array([1.0, 0.0]) if n <= 1.0e-30 else q / n


def _fill_bins(values: np.ndarray, valid: np.ndarray, fallback: float) -> np.ndarray:
    values = np.asarray(values, float).copy()
    idx = np.arange(values.size)
    if np.count_nonzero(valid) == 0:
        values[:] = fallback
    elif np.count_nonzero(valid) == 1:
        values[:] = values[valid][0]
    else:
        values[~valid] = np.interp(idx[~valid], idx[valid], values[valid])
    return values


def sample_process_zone_profile(
    mesh,
    sigma_gp: np.ndarray,
    rho_gp: np.ndarray | None,
    damage_nodal: np.ndarray,
    crack_tip,
    crack_direction,
    *,
    length_m: float,
    n_bins: int,
    forest_floor_m2: float = 5.0e12,
    sector_half_angle_deg: float = 45.0,
    damage_cutoff: float = 0.85,
    min_elements: int = 8,
    poisson: float = 0.28,
) -> ProcessZone2DProfile:
    n_bins = max(int(n_bins), 4)
    length_m = max(float(length_m), 1.0e-12)
    xgrid = (np.arange(n_bins, dtype=float) + 0.5) * length_m / n_bins
    fallback_rho = max(float(forest_floor_m2), 1.0)
    fallback_shape = np.sqrt(max(xgrid[0], 1.0e-30) / np.maximum(xgrid, xgrid[0]))

    cent = np.asarray(mesh.nodes, float)[np.asarray(mesh.elems, int)].mean(axis=1)
    area = np.maximum(np.asarray(mesh.area_e, float), 1.0e-30)
    tip = np.asarray(crack_tip, float).reshape(2)
    t = _unit(crack_direction)
    n = np.array([-t[1], t[0]])
    rel = cent - tip
    x = rel @ t
    y = rel @ n
    angle = np.degrees(np.arctan2(np.abs(y), np.maximum(x, 1.0e-30)))

    d = np.asarray(damage_nodal, float)
    de = d[np.asarray(mesh.elems, int)].mean(axis=1)
    selected = (
        (x > 0.0) & (x <= length_m) &
        (angle <= float(sector_half_angle_deg)) &
        (de < float(damage_cutoff))
    )
    ids = np.flatnonzero(selected)
    if ids.size < int(min_elements):
        return ProcessZone2DProfile(
            xgrid, np.full(n_bins, fallback_rho), fallback_shape,
            False, 0.0, int(ids.size), "insufficient_forward_undamaged_elements",
        )

    sigma = np.asarray(sigma_gp, float)
    if sigma.ndim != 2 or sigma.shape[0] < 3 or sigma.shape[1] != len(cent):
        raise ValueError("sigma_gp must have shape (>=3, ne)")
    sxx, syy, sxy = sigma[0], sigma[1], sigma[2]
    szz = float(poisson) * (sxx + syy)
    mean = (sxx + syy + szz) / 3.0
    dxx, dyy, dzz = sxx - mean, syy - mean, szz - mean
    seq2 = np.sqrt(np.maximum(1.5 * (dxx*dxx + dyy*dyy + dzz*dzz + 2.0*sxy*sxy), 0.0))

    if rho_gp is None:
        rho = np.full(len(cent), fallback_rho)
    else:
        rho = np.asarray(rho_gp, float).reshape(-1)
        if rho.size != len(cent):
            raise ValueError(f"rho_gp length {rho.size} != element count {len(cent)}")
        rho = np.maximum(rho, fallback_rho)

    edges = np.linspace(0.0, length_m, n_bins + 1)
    ibin = np.clip(np.searchsorted(edges, x[ids], side="right") - 1, 0, n_bins - 1)
    rho_bin = np.zeros(n_bins)
    stress_bin = np.zeros(n_bins)
    weight_bin = np.zeros(n_bins)
    for e, b in zip(ids, ibin):
        w = area[e]
        rho_bin[b] += w * rho[e]
        stress_bin[b] += w * max(seq2[e], 0.0)
        weight_bin[b] += w
    valid = weight_bin > 0.0
    rho_bin[valid] /= weight_bin[valid]
    stress_bin[valid] /= weight_bin[valid]
    rho_bin = _fill_bins(rho_bin, valid, fallback_rho)
    stress_bin = _fill_bins(stress_bin, valid, 1.0)

    near = np.flatnonzero(valid)[: max(1, min(5, np.count_nonzero(valid)))]
    ref = float(np.max(stress_bin[near])) if near.size else float(stress_bin[0])
    if not np.isfinite(ref) or ref <= 1.0e-30:
        stress_shape = fallback_shape
        reliable = False
        reason = "nonpositive_fem_stress_shape_reference"
    else:
        stress_shape = np.clip(stress_bin / ref, 0.02, 1.0)
        reliable = True
        reason = "ok"
    coverage = float(np.count_nonzero(valid) / n_bins)
    return ProcessZone2DProfile(
        xgrid,
        np.maximum(rho_bin, fallback_rho),
        stress_shape,
        bool(reliable),
        coverage,
        int(ids.size),
        reason,
    )


__all__ = ["ProcessZone2DProfile", "sample_process_zone_profile"]
