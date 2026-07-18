"""Tensor-resolved FEM drive capture and cohesive-stepper coupling for v10.0.5.

The wrapper observes the same stress state used by the live J-integral evaluation.
It converts the finite-radius FEM tensor only into dimensionless channel shapes;
the sharp-tip K amplitude remains the calibrated absolute stress scale.
"""
from __future__ import annotations

import copy
from dataclasses import asdict, dataclass
import math
from typing import Any, Callable, Mapping

import numpy as np

from .crystal import bcc_slip_traces
from .kinetic_cohesive_stepper import KineticCohesiveStepper


@dataclass(frozen=True)
class TensorResolvedDriveConfig:
    crystal_theta_deg: float = 45.0
    probe_radius_m: float = 10.0e-6
    sector_half_angle_deg: float = 25.0
    damage_cutoff: float = 0.85
    min_elements: int = 3
    schmid_reference: float = 0.5

    def validate(self) -> "TensorResolvedDriveConfig":
        if float(self.probe_radius_m) <= 0.0:
            raise ValueError("probe_radius_m must be positive")
        if not (0.0 < float(self.sector_half_angle_deg) < 90.0):
            raise ValueError("sector_half_angle_deg must lie in (0, 90)")
        if not (0.0 <= float(self.damage_cutoff) <= 1.0):
            raise ValueError("damage_cutoff must lie in [0, 1]")
        if int(self.min_elements) < 1:
            raise ValueError("min_elements must be positive")
        if float(self.schmid_reference) <= 0.0:
            raise ValueError("schmid_reference must be positive")
        return self


_RUNTIME: dict[str, Any] = {}


def reset_tensor_drive_runtime(
    config: TensorResolvedDriveConfig | None = None,
) -> None:
    cfg = (config or TensorResolvedDriveConfig()).validate()
    _RUNTIME.clear()
    _RUNTIME.update({
        "schema": "tensor_resolved_parallel_opening_emission_v10_0_5",
        "config": asdict(cfg),
        "capture_count": 0,
        "finite_capture_count": 0,
        "nonzero_emission_drive_capture_count": 0,
        "max_emission_drive_factor": 0.0,
        "min_emission_drive_factor": math.inf,
        "latest": None,
        "failures": [],
        "resolved_tensor_drives_active": True,
        "drive_factor_normalization_or_clipping_active": False,
        "directional_multiplier_applied_after_hazard": False,
        "fit_derived_shielding_cap_active": False,
    })


def _config() -> TensorResolvedDriveConfig:
    if not _RUNTIME:
        reset_tensor_drive_runtime()
    return TensorResolvedDriveConfig(**_RUNTIME["config"]).validate()


def _unit(vector: np.ndarray) -> np.ndarray:
    value = np.asarray(vector, dtype=float).reshape(2)
    norm = float(np.linalg.norm(value))
    if norm <= 1.0e-30:
        return np.array([1.0, 0.0])
    return value / norm


def _stress_tensor_on_ray(
    mesh: Any,
    sigma_gp: np.ndarray,
    damage: np.ndarray,
    tip: np.ndarray,
    ray_direction: np.ndarray,
    config: TensorResolvedDriveConfig,
) -> tuple[np.ndarray, dict[str, Any]]:
    nodes = np.asarray(mesh.nodes, dtype=float)
    elems = np.asarray(mesh.elems, dtype=int)
    centroids = nodes[elems].mean(axis=1)
    rel = centroids - np.asarray(tip, dtype=float).reshape(1, 2)
    direction = _unit(ray_direction)
    normal = np.array([-direction[1], direction[0]])
    longitudinal = rel @ direction
    transverse = rel @ normal
    radius = float(config.probe_radius_m)
    angle = np.degrees(np.arctan2(np.abs(transverse), np.maximum(longitudinal, 1.0e-30)))
    r = np.sqrt(np.sum(rel * rel, axis=1))
    damage_e = np.asarray(damage, dtype=float)[elems].mean(axis=1)
    valid_damage = damage_e < float(config.damage_cutoff)

    selected = np.zeros(elems.shape[0], dtype=bool)
    expansion_used = None
    for expansion in (1.0, 1.5, 2.0, 3.0):
        selected = (
            valid_damage
            & (longitudinal > 0.0)
            & (r >= 0.25 * radius / expansion)
            & (r <= 1.75 * radius * expansion)
            & (angle <= min(float(config.sector_half_angle_deg) * expansion, 85.0))
        )
        if int(np.count_nonzero(selected)) >= int(config.min_elements):
            expansion_used = expansion
            break

    if int(np.count_nonzero(selected)) < int(config.min_elements):
        candidates = np.flatnonzero(valid_damage & (longitudinal > 0.0))
        if candidates.size < int(config.min_elements):
            raise RuntimeError(
                "tensor-resolved drive probe has too few undamaged forward elements"
            )
        target = np.asarray(tip, dtype=float) + radius * direction
        distance = np.linalg.norm(centroids[candidates] - target[None, :], axis=1)
        keep = candidates[np.argsort(distance)[: int(config.min_elements)]]
        selected = np.zeros(elems.shape[0], dtype=bool)
        selected[keep] = True
        expansion_used = "nearest_forward_fallback"

    idx = np.flatnonzero(selected)
    area = np.asarray(getattr(mesh, "area_e", np.ones(elems.shape[0])), dtype=float)[idx]
    area = np.maximum(area, 1.0e-30)
    weights = area / float(np.sum(area))
    sigma = np.asarray(sigma_gp, dtype=float)
    sxx = float(weights @ sigma[0, idx])
    syy = float(weights @ sigma[1, idx])
    sxy = float(weights @ sigma[2, idx])
    tensor = np.array([[sxx, sxy], [sxy, syy]], dtype=float)
    if not np.all(np.isfinite(tensor)):
        raise RuntimeError("tensor-resolved drive probe produced nonfinite stress")
    return tensor, {
        "n_elements": int(idx.size),
        "probe_radius_m": radius,
        "expansion": expansion_used,
        "direction": direction.tolist(),
    }


def capture_tensor_resolved_drives(
    *,
    mesh: Any,
    sigma_gp: np.ndarray,
    damage: np.ndarray,
    crack_tip: np.ndarray,
    crack_direction: np.ndarray,
    KJ_Pa_sqrt_m: float,
) -> dict[str, Any]:
    config = _config()
    crack_t = _unit(crack_direction)
    crack_n = np.array([-crack_t[1], crack_t[0]])
    opening_tensor, opening_probe = _stress_tensor_on_ray(
        mesh,
        sigma_gp,
        damage,
        crack_tip,
        crack_t,
        config,
    )
    principal = np.linalg.eigvalsh(opening_tensor)
    sigma1 = max(float(principal[-1]), 0.0)
    sigma_nn = float(crack_n @ opening_tensor @ crack_n)
    amplitude = max(sigma1, max(sigma_nn, 0.0), 1.0e-30)
    opening_shape = max(sigma_nn, 0.0) / amplitude

    names: list[str] = []
    signed_tau: list[float] = []
    abs_tau: list[float] = []
    factors: list[float] = []
    probe_meta: list[dict[str, Any]] = []
    for system in bcc_slip_traces(float(config.crystal_theta_deg)):
        direction = _unit(np.asarray(system["t"], dtype=float))
        if float(direction @ crack_t) < 0.0:
            direction = -direction
        normal = _unit(np.asarray(system["n"], dtype=float))
        tensor, meta = _stress_tensor_on_ray(
            mesh,
            sigma_gp,
            damage,
            crack_tip,
            direction,
            config,
        )
        tau = float(direction @ tensor @ normal)
        equivalent_shape = abs(tau) / (
            float(config.schmid_reference) * amplitude
        )
        names.append(str(system["name"]))
        signed_tau.append(tau)
        abs_tau.append(abs(tau))
        factors.append(equivalent_shape)
        probe_meta.append(meta)

    factor_array = np.asarray(factors, dtype=float)
    if factor_array.size == 0 or not np.all(np.isfinite(factor_array)):
        raise RuntimeError("tensor-resolved emission factors are missing or nonfinite")
    if np.any(factor_array < 0.0):
        raise RuntimeError("tensor-resolved emission factors must be nonnegative")

    record = {
        "KJ_Pa_sqrt_m": float(KJ_Pa_sqrt_m),
        "crack_tip_m": np.asarray(crack_tip, dtype=float).tolist(),
        "crack_direction": crack_t.tolist(),
        "opening_probe_stress_tensor_Pa": opening_tensor.tolist(),
        "opening_probe_sigma1_Pa": sigma1,
        "opening_probe_sigma_nn_Pa": sigma_nn,
        "opening_shape_factor": opening_shape,
        "slip_system_names": names,
        "slip_system_tau_signed_Pa": signed_tau,
        "slip_system_tau_abs_Pa": abs_tau,
        "slip_system_drive_factors": factors,
        "opening_probe": opening_probe,
        "slip_system_probes": probe_meta,
        "tensor_resolved_drive_active": True,
        "drive_factor_normalization_or_clipping_active": False,
    }
    _RUNTIME["capture_count"] += 1
    _RUNTIME["finite_capture_count"] += 1
    if np.any(factor_array > 0.0):
        _RUNTIME["nonzero_emission_drive_capture_count"] += 1
    _RUNTIME["max_emission_drive_factor"] = max(
        float(_RUNTIME["max_emission_drive_factor"]),
        float(np.max(factor_array)),
    )
    _RUNTIME["min_emission_drive_factor"] = min(
        float(_RUNTIME["min_emission_drive_factor"]),
        float(np.min(factor_array)),
    )
    _RUNTIME["latest"] = copy.deepcopy(record)
    return record


def make_tensor_resolved_J_wrapper(
    original_compute_J: Callable[..., tuple[float, float, dict]],
) -> Callable[..., tuple[float, float, dict]]:
    def wrapped(
        mesh,
        u,
        sigma_gp,
        psi_e_gp,
        d,
        crack_tip,
        crack_direction,
        mat,
        ell,
        cfg=None,
        crack_segments=None,
        exclude_radius=0.0,
    ):
        result = original_compute_J(
            mesh,
            u,
            sigma_gp,
            psi_e_gp,
            d,
            crack_tip,
            crack_direction,
            mat,
            ell,
            cfg=cfg,
            crack_segments=crack_segments,
            exclude_radius=exclude_radius,
        )
        J, KJ, info = result
        drive = capture_tensor_resolved_drives(
            mesh=mesh,
            sigma_gp=sigma_gp,
            damage=d,
            crack_tip=crack_tip,
            crack_direction=crack_direction,
            KJ_Pa_sqrt_m=KJ,
        )
        enriched = dict(info)
        enriched.update({
            "tensor_resolved_drive_active": True,
            "opening_shape_factor": drive["opening_shape_factor"],
            "slip_system_drive_factors": drive["slip_system_drive_factors"],
            "slip_system_tau_signed_Pa": drive["slip_system_tau_signed_Pa"],
        })
        return J, KJ, enriched

    wrapped.__name__ = getattr(original_compute_J, "__name__", "compute_J_integral")
    wrapped.__qualname__ = getattr(original_compute_J, "__qualname__", wrapped.__name__)
    wrapped.__doc__ = getattr(original_compute_J, "__doc__", None)
    wrapped._v1005_tensor_resolved_J_wrapper = True
    wrapped._v1005_original = original_compute_J
    return wrapped


def latest_tensor_drive() -> dict[str, Any]:
    latest = _RUNTIME.get("latest")
    if latest is None:
        raise RuntimeError(
            "tensor-resolved cohesive step requested before a live FEM/J drive capture"
        )
    return copy.deepcopy(latest)


class TensorResolvedKineticCohesiveStepper(KineticCohesiveStepper):
    """Use the latest live FEM tensor projection instead of unit system weights."""

    tensor_resolved_parallel_coupling = True

    @staticmethod
    def _drives(
        mechanics: Mapping[str, Any],
    ) -> tuple[float, float, np.ndarray | None]:
        K_open, K_cleave, _legacy_weights = KineticCohesiveStepper._drives(mechanics)
        drive = latest_tensor_drive()
        KJ = float(drive["KJ_Pa_sqrt_m"])
        tolerance = max(1.0e-6, 1.0e-10 * max(abs(K_open), abs(KJ), 1.0))
        if abs(K_open - KJ) > tolerance:
            raise RuntimeError(
                "tensor-resolved drive cache is stale relative to the mechanics solve: "
                f"Kopen={K_open:.16g}, cached KJ={KJ:.16g}"
            )
        factors = np.asarray(drive["slip_system_drive_factors"], dtype=float)
        if factors.ndim != 1 or factors.size < 1 or not np.all(np.isfinite(factors)):
            raise RuntimeError("invalid tensor-resolved slip-system drive factors")
        if isinstance(mechanics, dict):
            mechanics.update(copy.deepcopy(drive))
        return K_open, K_cleave, factors

    def audit_payload(self) -> dict[str, Any]:
        payload = super().audit_payload()
        payload.update({
            "schema": "kinetic_cohesive_stepper_v10_0_5",
            "tensor_resolved_parallel_coupling": True,
            "legacy_unit_slip_weights_ignored": True,
            "directional_multiplier_applied_after_hazard": False,
        })
        return payload


def tensor_drive_runtime_payload() -> dict[str, Any]:
    payload = copy.deepcopy(_RUNTIME)
    if math.isinf(float(payload.get("min_emission_drive_factor", math.inf))):
        payload["min_emission_drive_factor"] = None
    payload["implementation_certified"] = bool(
        int(payload.get("capture_count", 0)) > 0
        and int(payload.get("finite_capture_count", 0))
        == int(payload.get("capture_count", 0))
        and not payload.get("failures")
        and payload.get("resolved_tensor_drives_active", False)
        and not payload.get("drive_factor_normalization_or_clipping_active", True)
        and not payload.get("directional_multiplier_applied_after_hazard", True)
        and not payload.get("fit_derived_shielding_cap_active", True)
    )
    return payload


__all__ = [
    "TensorResolvedDriveConfig",
    "reset_tensor_drive_runtime",
    "capture_tensor_resolved_drives",
    "make_tensor_resolved_J_wrapper",
    "latest_tensor_drive",
    "TensorResolvedKineticCohesiveStepper",
    "tensor_drive_runtime_payload",
]
