"""PF v10.2.14 active-only signed shielding atlas for FEM/CZM v10.0.5.14.1.

The PF production artifact is a candidate-independent family of mechanically
measured signed-Burgers kernels indexed by cumulative crack-path extension.  It
is not a single frozen 2 x N matrix.  This loader reproduces the family-level
activation normalization, crack-extension interpolation, and spatial projection
onto the runtime moving-process-zone grid.
"""
from __future__ import annotations

from dataclasses import dataclass
import copy
import json
from pathlib import Path
from typing import Any

import numpy as np

from .persistent_site_signed_support_v100514 import SignedShieldingKernelV100514

FAMILY_SCHEMA = "v10.2.14_active_only_real_signed_2d_shielding_atlas"
POINT_RELEASE = "10.0.5.14.1"


def _finite_vector(value: Any, name: str) -> np.ndarray:
    array = np.asarray(value, dtype=float).reshape(-1)
    if array.size == 0 or not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be a nonempty finite vector")
    return array


def _finite_matrix(value: Any, name: str) -> np.ndarray:
    array = np.asarray(value, dtype=float)
    if array.ndim != 2 or not np.all(np.isfinite(array)):
        raise ValueError(f"{name} must be a finite 2-D array")
    return array


def _strictly_increasing(array: np.ndarray, name: str) -> None:
    if np.any(np.diff(array) <= 0.0):
        raise ValueError(f"{name} must be strictly increasing")


@dataclass(frozen=True)
class SignedKernelFamilyStateV1005141:
    state_id: str
    crack_extension_m: float
    active_I: np.ndarray
    active_II: np.ndarray
    wake_I: np.ndarray
    wake_II: np.ndarray
    r_eff_over_r0: float
    opening_strength_fraction: float


@dataclass(frozen=True)
class SignedShieldingKernelFamilyV1005141:
    active_x_m: np.ndarray
    wake_x_m: np.ndarray
    activation_to_line_content_by_system: np.ndarray
    states: tuple[SignedKernelFamilyStateV1005141, ...]
    interpolation_method: str
    interpolation_neighbors: int
    interpolation_power: float
    envelope_relative_tolerance: float
    extrapolation_allowed: bool
    metadata: dict[str, Any]
    source_path: str

    @classmethod
    def from_json(
        cls, path: str | Path
    ) -> "SignedShieldingKernelFamilyV1005141":
        source = Path(path).expanduser().resolve()
        payload = json.loads(source.read_text())
        if payload.get("schema") != FAMILY_SCHEMA:
            raise ValueError(
                f"unsupported signed-kernel family schema: {payload.get('schema')!r}"
            )
        required_true = (
            "candidate_independent",
            "counts_are_signed_burgers_lines",
            "normalization_is_mechanically_derived",
            "active_kernel_mechanically_measured",
            "kernel_from_signed_interaction_integral",
            "signed_burgers_population_required",
            "production_parameterization_allowed",
        )
        for name in required_true:
            if payload.get(name) is not True:
                raise ValueError(f"PF kernel family requires {name}=true")
        if payload.get("constitutive_K_shield_cap") is not False:
            raise ValueError("PF parity requires no constitutive shielding cap")
        if payload.get("constitutive_K_shield_cap_present") is not False:
            raise ValueError("PF parity requires no shielding-cap field")
        if payload.get("wake_kernel_forced_zero") is not True:
            raise ValueError("v10.0.5.14.1 supports the active-only zero-wake atlas")
        if payload.get("wake_shielding_supported") is not False:
            raise ValueError("active-only PF family must report wake shielding unsupported")
        if payload.get("crack_extension_m_semantics") != (
            "cumulative_crack_path_extension_m"
        ):
            raise ValueError("kernel family must use cumulative crack-path extension")
        semantics = dict(payload.get("v10_2_13_state_semantics", {}))
        if semantics.get("cumulative_crack_path_extension_used_for_interpolation") is not True:
            raise ValueError("crack extension must be the active interpolation coordinate")
        if semantics.get("crack_extension_extrapolation_allowed") is not False:
            raise ValueError("PF kernel-family extrapolation must be disabled")
        if semantics.get("analytical_r_eff_used_for_interpolation") is not False:
            raise ValueError("r_eff must remain a disabled compatibility coordinate")
        if semantics.get("opening_strength_fraction_used_for_interpolation") is not False:
            raise ValueError("opening strength must remain a disabled compatibility coordinate")

        active_x = _finite_vector(payload.get("active_x_m"), "active_x_m")
        wake_x = _finite_vector(payload.get("wake_x_m"), "wake_x_m")
        _strictly_increasing(active_x, "active_x_m")
        _strictly_increasing(wake_x, "wake_x_m")
        conversion = _finite_vector(
            payload.get("activation_to_line_content_by_system"),
            "activation_to_line_content_by_system",
        )
        if np.any(conversion <= 0.0):
            raise ValueError("activation-to-line conversion must be positive")

        raw_states = payload.get("states")
        if not isinstance(raw_states, list) or len(raw_states) < 2:
            raise ValueError("kernel family requires at least two measured states")
        states: list[SignedKernelFamilyStateV1005141] = []
        for index, raw in enumerate(raw_states):
            active_I = _finite_matrix(
                raw.get("active_kernel_I_Pa_sqrt_m_per_signed_line"),
                f"states[{index}].active_kernel_I",
            )
            active_II = _finite_matrix(
                raw.get("active_kernel_II_Pa_sqrt_m_per_signed_line"),
                f"states[{index}].active_kernel_II",
            )
            wake_I = _finite_matrix(
                raw.get("wake_kernel_I_Pa_sqrt_m_per_signed_line"),
                f"states[{index}].wake_kernel_I",
            )
            wake_II = _finite_matrix(
                raw.get("wake_kernel_II_Pa_sqrt_m_per_signed_line"),
                f"states[{index}].wake_kernel_II",
            )
            expected_active = (conversion.size, active_x.size)
            expected_wake = (conversion.size, wake_x.size)
            if active_I.shape != expected_active or active_II.shape != expected_active:
                raise ValueError(
                    f"state {raw.get('state_id')} active kernels do not match {expected_active}"
                )
            if wake_I.shape != expected_wake or wake_II.shape != expected_wake:
                raise ValueError(
                    f"state {raw.get('state_id')} wake kernels do not match {expected_wake}"
                )
            if not np.allclose(wake_I, 0.0) or not np.allclose(wake_II, 0.0):
                raise ValueError("active-only PF family must contain exactly zero wake kernels")
            states.append(
                SignedKernelFamilyStateV1005141(
                    state_id=str(raw.get("state_id", f"state_{index}")),
                    crack_extension_m=float(raw["crack_extension_m"]),
                    active_I=active_I,
                    active_II=active_II,
                    wake_I=wake_I,
                    wake_II=wake_II,
                    r_eff_over_r0=float(raw.get("r_eff_over_r0", 1.0)),
                    opening_strength_fraction=float(
                        raw.get("opening_strength_fraction", 0.0)
                    ),
                )
            )
        states.sort(key=lambda state: state.crack_extension_m)
        extensions = np.asarray([state.crack_extension_m for state in states])
        _strictly_increasing(extensions, "state crack extensions")
        if len({state.r_eff_over_r0 for state in states}) != 1:
            raise ValueError("this parity loader requires the constant r_eff compatibility axis")
        if len({state.opening_strength_fraction for state in states}) != 1:
            raise ValueError("this parity loader requires the constant opening compatibility axis")

        interpolation = dict(payload.get("interpolation", {}))
        method = str(interpolation.get("method", ""))
        if method != "inverse_distance":
            raise ValueError(f"unsupported PF family interpolation method: {method!r}")
        neighbors = int(interpolation.get("neighbors", len(states)))
        power = float(interpolation.get("power", 2.0))
        tolerance = float(interpolation.get("envelope_relative_tolerance", 1.0e-10))
        extrapolation = bool(interpolation.get("extrapolation_allowed", False))
        if neighbors < 1 or power <= 0.0 or tolerance < 0.0:
            raise ValueError("invalid inverse-distance interpolation configuration")
        if extrapolation:
            raise ValueError("PF production atlas forbids crack-extension extrapolation")

        excluded = {
            "active_x_m",
            "wake_x_m",
            "activation_to_line_content_by_system",
            "states",
        }
        metadata = {
            key: copy.deepcopy(value)
            for key, value in payload.items()
            if key not in excluded
        }
        return cls(
            active_x_m=active_x,
            wake_x_m=wake_x,
            activation_to_line_content_by_system=conversion,
            states=tuple(states),
            interpolation_method=method,
            interpolation_neighbors=neighbors,
            interpolation_power=power,
            envelope_relative_tolerance=tolerance,
            extrapolation_allowed=extrapolation,
            metadata=metadata,
            source_path=str(source),
        )

    @property
    def schema(self) -> str:
        return FAMILY_SCHEMA

    @property
    def extension_levels_m(self) -> np.ndarray:
        return np.asarray([state.crack_extension_m for state in self.states])

    def validate(
        self,
        n_systems: int,
        n_bins: int,
        active_x_m: np.ndarray | None = None,
        wake_x_m: np.ndarray | None = None,
    ) -> None:
        if self.activation_to_line_content_by_system.shape != (int(n_systems),):
            raise ValueError("kernel-family system count does not match the MPZ")
        if active_x_m is not None:
            target = _finite_vector(active_x_m, "runtime active_x_m")
            if target.shape != (int(n_bins),):
                raise ValueError("runtime active grid does not match the MPZ bin count")
            _strictly_increasing(target, "runtime active_x_m")
            if target[0] < 0.0:
                raise ValueError("runtime active coordinates must be nonnegative")
            if target[-1] > self.active_x_m[-1] + 1.0e-15:
                raise ValueError(
                    "runtime MPZ extends beyond the mechanically measured active kernel"
                )
        if wake_x_m is not None:
            target_wake = _finite_vector(wake_x_m, "runtime wake_x_m")
            _strictly_increasing(target_wake, "runtime wake_x_m")

    def _weights(self, cumulative_crack_path_extension_m: float) -> np.ndarray:
        extensions = self.extension_levels_m
        target = float(cumulative_crack_path_extension_m)
        span = max(float(extensions[-1] - extensions[0]), 1.0e-30)
        envelope_tolerance = max(
            self.envelope_relative_tolerance * span, 1.0e-15
        )
        if target < extensions[0] - envelope_tolerance:
            raise ValueError("crack extension lies below the PF kernel-family envelope")
        if target > extensions[-1] + envelope_tolerance:
            raise ValueError(
                "crack extension exceeds the PF kernel-family envelope; extrapolation is disabled"
            )
        target = min(max(target, float(extensions[0])), float(extensions[-1]))
        distances = np.abs(extensions - target)
        exact = int(np.argmin(distances))
        if distances[exact] <= envelope_tolerance:
            weights = np.zeros_like(extensions)
            weights[exact] = 1.0
            return weights
        neighbor_count = min(self.interpolation_neighbors, extensions.size)
        selected = np.argsort(distances)[:neighbor_count]
        normalized = distances[selected] / span
        raw = 1.0 / np.power(normalized, self.interpolation_power)
        weights = np.zeros_like(extensions)
        weights[selected] = raw / float(np.sum(raw))
        return weights

    @staticmethod
    def _project_spatial(
        source_x_m: np.ndarray,
        source_values: np.ndarray,
        target_x_m: np.ndarray,
    ) -> np.ndarray:
        target = np.asarray(target_x_m, dtype=float).reshape(-1)
        projected = np.empty((source_values.shape[0], target.size), dtype=float)
        for system, row in enumerate(source_values):
            # PF v10.2.14 uses piecewise-linear physical projection with exact
            # endpoint coverage.  Runtime points inside the endpoint half-cells
            # inherit the measured endpoint coefficient rather than extrapolating.
            projected[system] = np.interp(
                target,
                source_x_m,
                row,
                left=float(row[0]),
                right=float(row[-1]),
            )
        return projected

    def snapshot(
        self,
        cumulative_crack_path_extension_m: float,
        runtime_active_x_m: np.ndarray,
        runtime_wake_x_m: np.ndarray,
    ) -> SignedShieldingKernelV100514:
        self.validate(
            self.activation_to_line_content_by_system.size,
            len(runtime_active_x_m),
            active_x_m=runtime_active_x_m,
            wake_x_m=runtime_wake_x_m,
        )
        weights = self._weights(cumulative_crack_path_extension_m)
        active_I_source = np.zeros(
            (self.activation_to_line_content_by_system.size, self.active_x_m.size)
        )
        active_II_source = np.zeros_like(active_I_source)
        for weight, state in zip(weights, self.states):
            active_I_source += float(weight) * state.active_I
            active_II_source += float(weight) * state.active_II
        active_I = self._project_spatial(
            self.active_x_m,
            active_I_source,
            np.asarray(runtime_active_x_m, dtype=float),
        )
        active_II = self._project_spatial(
            self.active_x_m,
            active_II_source,
            np.asarray(runtime_active_x_m, dtype=float),
        )
        wake = np.zeros(
            (
                self.activation_to_line_content_by_system.size,
                len(runtime_wake_x_m),
            ),
            dtype=float,
        )
        state_weights = {
            state.state_id: float(weight)
            for state, weight in zip(self.states, weights)
            if weight > 0.0
        }
        metadata = {
            **copy.deepcopy(self.metadata),
            "schema": FAMILY_SCHEMA,
            "kernel_artifact_kind": "crack_extension_family",
            "cumulative_crack_path_extension_m": float(
                cumulative_crack_path_extension_m
            ),
            "state_weights": state_weights,
            "state_extension_levels_m": self.extension_levels_m.tolist(),
            "spatial_projection": "piecewise_linear_with_endpoint_hold",
            "runtime_active_x_m": np.asarray(runtime_active_x_m).tolist(),
            "runtime_active_kernel_II_Pa_sqrt_m_per_signed_line": (
                active_II.tolist()
            ),
            "wake_kernel_forced_zero": True,
            "constitutive_K_shield_cap": False,
        }
        return SignedShieldingKernelV100514(
            active_kernel_Pa_sqrt_m_per_signed_line=active_I,
            wake_kernel_Pa_sqrt_m_per_signed_line=wake,
            activation_to_line_content_by_system=(
                self.activation_to_line_content_by_system.copy()
            ),
            metadata=metadata,
            source_path=self.source_path,
            active_x_m=np.asarray(runtime_active_x_m, dtype=float).copy(),
            wake_x_m=np.asarray(runtime_wake_x_m, dtype=float).copy(),
        )

    def audit_payload(self) -> dict[str, Any]:
        return {
            "schema": FAMILY_SCHEMA,
            "point_release": POINT_RELEASE,
            "source_path": self.source_path,
            "artifact_kind": "crack_extension_kernel_family",
            "state_ids": [state.state_id for state in self.states],
            "crack_extension_levels_m": self.extension_levels_m.tolist(),
            "source_active_grid_points": int(self.active_x_m.size),
            "interpolation": {
                "method": self.interpolation_method,
                "neighbors": self.interpolation_neighbors,
                "power": self.interpolation_power,
                "extrapolation_allowed": self.extrapolation_allowed,
            },
            "spatial_projection": "piecewise_linear_with_endpoint_hold",
            "activation_to_line_content_by_system": (
                self.activation_to_line_content_by_system.tolist()
            ),
            "wake_kernel_forced_zero": True,
            "constitutive_K_shield_cap": False,
            "candidate_independent": True,
            "counts_are_signed_burgers_lines": True,
            "normalization_is_mechanically_derived": True,
        }


def load_signed_shielding_artifact_v1005141(
    path: str | Path,
) -> SignedShieldingKernelV100514 | SignedShieldingKernelFamilyV1005141:
    source = Path(path).expanduser().resolve()
    payload = json.loads(source.read_text())
    if payload.get("schema") == FAMILY_SCHEMA:
        return SignedShieldingKernelFamilyV1005141.from_json(source)
    return SignedShieldingKernelV100514.from_json(source)


__all__ = [
    "FAMILY_SCHEMA",
    "POINT_RELEASE",
    "SignedKernelFamilyStateV1005141",
    "SignedShieldingKernelFamilyV1005141",
    "load_signed_shielding_artifact_v1005141",
]
