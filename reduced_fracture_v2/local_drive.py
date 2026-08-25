"""Typed local mechanics drives for the V2 reduced crack-front model.

Global fracture metrics and local kinetic fields deliberately live in separate
types.  In particular, a scalar KJ can never satisfy the PF signed-emission
contract.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Mapping, Sequence
from pathlib import Path

import numpy as np
import pandas as pd


class DriveQualificationError(RuntimeError):
    """Raised when a provider cannot supply a source-qualified drive."""


def _vec2(value: Sequence[float], name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if out.shape != (2,) or not np.all(np.isfinite(out)):
        raise DriveQualificationError(f"{name} must be a finite 2-vector")
    return out


def _sym2(value: Sequence[Sequence[float]], name: str) -> np.ndarray:
    out = np.asarray(value, dtype=float)
    if out.shape != (2, 2) or not np.all(np.isfinite(out)):
        raise DriveQualificationError(f"{name} must be a finite 2x2 tensor")
    if not np.array_equal(out, out.T):
        raise DriveQualificationError(f"{name} must be exactly symmetric")
    return out


@dataclass(frozen=True)
class GlobalStructuralDrive:
    opening_m: float
    reaction_N: float
    native_J_J_m2: float
    native_KJ_Pa_sqrt_m: float
    qualified_G_J_m2: float | None = None
    qualified_KG_Pa_sqrt_m: float | None = None


@dataclass(frozen=True)
class LocalTensorSample:
    crack_local_s_m: float
    crack_local_n_m: float
    stress_Pa: tuple[tuple[float, float], tuple[float, float]]
    strain: tuple[tuple[float, float], tuple[float, float]] | None = None

    def stress_array(self) -> np.ndarray:
        return _sym2(self.stress_Pa, "stress_Pa")


@dataclass(frozen=True)
class EmissionDrive:
    opening_tensor_Pa: tuple[tuple[float, float], tuple[float, float]]
    channel_tensors_Pa: tuple[
        tuple[tuple[float, float], tuple[float, float]], ...
    ]
    trace_directions: tuple[tuple[float, float], ...]
    trace_normals: tuple[tuple[float, float], ...]
    resolved_tau_signed_Pa: tuple[float, ...]
    drive_factors: tuple[float, ...]
    reliable: bool
    mechanics_serial: int

    def validate(self) -> "EmissionDrive":
        _sym2(self.opening_tensor_Pa, "opening_tensor_Pa")
        n = len(self.channel_tensors_Pa)
        if n != 2 or not all(len(x) == n for x in (
            self.trace_directions, self.trace_normals,
            self.resolved_tau_signed_Pa, self.drive_factors,
        )):
            raise DriveQualificationError("PF production requires exactly two channels")
        for i, tensor in enumerate(self.channel_tensors_Pa):
            _sym2(tensor, f"channel_tensors_Pa[{i}]")
        for i, (t_raw, n_raw) in enumerate(zip(self.trace_directions, self.trace_normals)):
            t, normal = _vec2(t_raw, f"trace_direction[{i}]"), _vec2(n_raw, f"trace_normal[{i}]")
            if not np.isclose(np.linalg.norm(t), 1.0, rtol=0.0, atol=1e-12):
                raise DriveQualificationError("trace direction is not unit length")
            if not np.isclose(np.linalg.norm(normal), 1.0, rtol=0.0, atol=1e-12):
                raise DriveQualificationError("trace normal is not unit length")
            if not np.isclose(t @ normal, 0.0, rtol=0.0, atol=1e-12):
                raise DriveQualificationError("trace basis is not orthogonal")
            tau = float(t @ _sym2(self.channel_tensors_Pa[i], "channel") @ normal)
            if not np.isclose(tau, self.resolved_tau_signed_Pa[i], rtol=2e-14, atol=1e-6):
                raise DriveQualificationError("resolved signed shear disagrees with tensor projection")
        if not self.reliable:
            raise DriveQualificationError("PF persistent emission requires a reliable 2-D tensor drive")
        return self


@dataclass(frozen=True)
class NativeDriveBundle:
    global_drive: GlobalStructuralDrive
    crack_tip_xy_m: tuple[float, float]
    crack_tangent: tuple[float, float]
    crack_normal: tuple[float, float]
    samples: tuple[LocalTensorSample, ...]
    emission: EmissionDrive | None
    tip_radius_m: float
    front_width_m: float
    mesh_wake_topology_fingerprint: str
    source_commit: str
    map_sha256: str
    extension_bounds_m: tuple[float, float]
    sampling_bounds_m: tuple[float, float, float, float]
    selected_emission_system: str | None = None
    probe_support: Mapping[str, tuple[int, ...]] | None = None
    probe_weights: Mapping[str, tuple[float, ...]] | None = None
    process_zone_metadata: Mapping[str, Any] | None = None
    field_snapshot_hash: str = ""
    probe_query_hash: str = ""

    def require_pf_signed_emission(self) -> EmissionDrive:
        if self.emission is None:
            raise DriveQualificationError(
                "scalar KJ is insufficient: PF persistent signed emission requires local tensors"
            )
        return self.emission.validate()

    def validate_basis(self) -> None:
        tangent = _vec2(self.crack_tangent, "crack_tangent")
        normal = _vec2(self.crack_normal, "crack_normal")
        if not np.isclose(np.linalg.norm(tangent), 1.0, atol=1e-12, rtol=0.0):
            raise DriveQualificationError("crack tangent is not unit length")
        if not np.isclose(np.linalg.norm(normal), 1.0, atol=1e-12, rtol=0.0):
            raise DriveQualificationError("crack normal is not unit length")
        if not np.isclose(tangent @ normal, 0.0, atol=1e-12, rtol=0.0):
            raise DriveQualificationError("crack-local basis is not orthogonal")

    def require_radius_in_domain(self) -> None:
        s0, s1, n0, n1 = self.sampling_bounds_m
        if not (s0 <= self.tip_radius_m <= s1 and n0 <= 0.0 <= n1):
            raise DriveQualificationError("dynamic tip radius is outside the qualified profile domain")


@dataclass(frozen=True)
class PFProbeRegimeKey:
    """Discrete identity of a PF sharp-wake tensor-probe state.

    Continuous tensor values are deliberately excluded.  A change in any
    member means that interpolation would mix different discrete operators.
    Signed-system activity is categorical and is included because a sign or
    system switch is a hard kinetic boundary, not a small interpolation error.
    """

    wake_topology_fingerprint: str
    front_basis_fingerprint: str
    probe_support_fingerprint: str
    probe_weight_fingerprint: str
    selected_system_signature: str
    resolved_sign_signature: str
    reliable: bool

    @classmethod
    def from_row(cls, row) -> "PFProbeRegimeKey":
        required = (
            "wake_fingerprint", "front_basis_fingerprint",
            "probe_support_fingerprint", "probe_weight_fingerprint",
            "selected_system_signature", "resolved_sign_signature", "reliable",
        )
        missing = [name for name in required if name not in row]
        if missing:
            raise DriveQualificationError(
                "PF regime metadata missing: " + ", ".join(missing)
            )
        reliable = row["reliable"]
        if isinstance(reliable, str):
            if reliable.strip().lower() not in {"true", "false"}:
                raise DriveQualificationError("PF regime reliability is not boolean")
            reliable = reliable.strip().lower() == "true"
        return cls(
            *(str(row[name]) for name in required[:-1]),
            bool(reliable),
        )

    def fingerprint(self) -> str:
        import hashlib
        import json
        payload = json.dumps(self.__dict__, sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(payload.encode()).hexdigest()


def interpolate_symmetric_tensor(
    x: float, x0: float, x1: float,
    y0: Sequence[Sequence[float]], y1: Sequence[Sequence[float]],
) -> np.ndarray:
    """Bounded piecewise-linear interpolation retaining sign and symmetry."""
    if not x0 <= x <= x1 or x1 <= x0:
        raise DriveQualificationError("tensor interpolation is outside qualified bounds")
    a, b = _sym2(y0, "y0"), _sym2(y1, "y1")
    weight = (x - x0) / (x1 - x0)
    out = (1.0 - weight) * a + weight * b
    return 0.5 * (out + out.T)


class PFTensorDriveProvider:
    """Bounded PF tensor-node provider with exact source projections.

    Sparse-map qualification is intentionally external to this evaluator.  A
    manifest marked unqualified must not be admitted by ``from_artifacts``.
    """
    def __init__(
        self,
        nodes: pd.DataFrame,
        manifest: dict,
        *,
        exact_query: Callable[[float, float, PFProbeRegimeKey | None, str], EmissionDrive] | None = None,
    ):
        self.nodes = nodes.sort_values("actual_extension_um").drop_duplicates(
            "actual_extension_um", keep="last"
        ).reset_index(drop=True)
        self.manifest = dict(manifest)
        self.exact_query = exact_query
        self._exact_cache: dict[tuple[float, float, str], EmissionDrive] = {}

    @classmethod
    def from_artifacts(cls, summary_csv: str | Path, manifest_json: str | Path,
                       *, allow_source_nodes_only: bool = False,
                       exact_query=None):
        import json
        manifest = json.loads(Path(manifest_json).read_text())
        status = str(manifest.get("status", ""))
        if status != "PF_TENSOR_DRIVE_MAP_QUALIFIED" and not (
            allow_source_nodes_only and status == "PF_RAW_TENSOR_PROFILE_MAP_QUALIFIED_AT_SOURCE_NODES"
        ):
            raise DriveQualificationError("PF tensor map interpolation is not qualified")
        return cls(pd.read_csv(summary_csv), manifest, exact_query=exact_query)

    @staticmethod
    def _values(row) -> dict[str, float]:
        return {
            column: float(row[column])
            for column in row.index
            if column.endswith("_Pa") or column.startswith("factor")
        }

    @staticmethod
    def _sign_signature(values: dict[str, float]) -> tuple[int, int]:
        return tuple(int(np.sign(values[f"tau{i}_signed_Pa"])) for i in range(2))

    def _exact_fallback(
        self, extension_m: float, opening_m: float, regime_key: PFProbeRegimeKey | None,
        geometry_state_fingerprint: str | None,
    ) -> EmissionDrive:
        if self.exact_query is None:
            raise DriveQualificationError(
                "PF probe/topology transition requires an exact deterministic query"
            )
        if regime_key is None and not geometry_state_fingerprint:
            raise DriveQualificationError(
                "exact deterministic query requires a realized geometry fingerprint or regime key"
            )
        key_fingerprint = (
            str(geometry_state_fingerprint)
            if geometry_state_fingerprint
            else regime_key.fingerprint()
        )
        cache_key = (float(extension_m), float(opening_m), key_fingerprint)
        if cache_key not in self._exact_cache:
            result = self.exact_query(
                extension_m, opening_m, regime_key, key_fingerprint
            )
            self._exact_cache[cache_key] = result.validate()
        return self._exact_cache[cache_key]

    def evaluate(
        self, extension_m: float, opening_m: float, *,
        regime_key: PFProbeRegimeKey | None = None,
        geometry_state_fingerprint: str | None = None,
        exact_node_only: bool = False,
    ) -> EmissionDrive:
        x = float(extension_m) * 1e6
        grid = self.nodes.actual_extension_um.to_numpy(float)
        if x < grid[0] or x > grid[-1]:
            raise DriveQualificationError("extension is outside the PF tensor map")
        exact = np.flatnonzero(np.isclose(grid, x, rtol=0.0, atol=1e-12))
        if exact.size:
            row = self.nodes.iloc[int(exact[0])]
            source_key = PFProbeRegimeKey.from_row(row)
            geometry_mismatch = (
                geometry_state_fingerprint is not None
                and "event_geometry_fingerprint" in row
                and str(geometry_state_fingerprint) != str(row["event_geometry_fingerprint"])
            )
            if (regime_key is not None and regime_key != source_key) or geometry_mismatch:
                return self._exact_fallback(
                    extension_m, opening_m, regime_key, geometry_state_fingerprint
                )
            values = self._values(row)
        else:
            if exact_node_only:
                raise DriveQualificationError("sparse PF tensor interpolation is unqualified")
            right = int(np.searchsorted(grid, x, side="right"))
            left = right - 1
            left_row, right_row = self.nodes.iloc[left], self.nodes.iloc[right]
            left_key = PFProbeRegimeKey.from_row(left_row)
            right_key = PFProbeRegimeKey.from_row(right_row)
            if left_key != right_key or (regime_key is not None and regime_key != left_key):
                return self._exact_fallback(
                    extension_m, opening_m, regime_key, geometry_state_fingerprint
                )
            values = {
                column: float(np.interp(x, grid[[left, right]], self.nodes[column].to_numpy(float)[[left, right]]))
                for column in self.nodes.columns
                if column.endswith("_Pa") or column.startswith("factor")
            }
            endpoint_signs = {self._sign_signature(self._values(left_row)), self._sign_signature(self._values(right_row))}
            if len(endpoint_signs) != 1 or self._sign_signature(values) not in endpoint_signs:
                raise DriveQualificationError("signed emission-system transition cannot be interpolated")
        scale = float(opening_m) / float(self.manifest["reference_opening_m"])
        tensor=lambda prefix: ((values[f"{prefix}_sigma_xx_Pa"]*scale, values[f"{prefix}_sigma_xy_Pa"]*scale),
                               (values[f"{prefix}_sigma_xy_Pa"]*scale, values[f"{prefix}_sigma_yy_Pa"]*scale))
        root2=float(np.sqrt(2.0))
        return EmissionDrive(
            tensor("opening"), (tensor("channel0"), tensor("channel1")),
            ((1/root2,1/root2),(1/root2,-1/root2)),
            ((-1/root2,1/root2),(1/root2,1/root2)),
            (values["tau0_signed_Pa"]*scale,values["tau1_signed_Pa"]*scale),
            (values["factor0"],values["factor1"]), True, -1,
        ).validate()
