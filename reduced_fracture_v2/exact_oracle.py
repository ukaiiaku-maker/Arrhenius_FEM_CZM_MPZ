"""Exact deterministic elastic-field oracle contracts for V2.

The field and probe layers are intentionally separate.  An elastic snapshot is
candidate independent; a process-zone state may change a probe query but can
never change the field-cache key.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable, Mapping, Protocol, Sequence
import hashlib
import json

import numpy as np

from .local_drive import DriveQualificationError, NativeDriveBundle


def _digest(value: Any) -> str:
    def default(item: Any):
        if isinstance(item, np.ndarray):
            array = np.ascontiguousarray(item)
            return {
                "dtype": str(array.dtype), "shape": array.shape,
                "sha256": hashlib.sha256(array.tobytes()).hexdigest(),
            }
        raise TypeError(type(item).__name__)
    raw = json.dumps(value, sort_keys=True, separators=(",", ":"), default=default)
    return hashlib.sha256(raw.encode()).hexdigest()


def _readonly_array(value: Any, shape_tail: tuple[int, ...], name: str) -> np.ndarray:
    array = np.array(value, dtype=float, copy=True)
    if shape_tail and (array.ndim < len(shape_tail) or tuple(array.shape[-len(shape_tail):]) != shape_tail):
        raise DriveQualificationError(f"{name} has invalid shape {array.shape}")
    if not np.all(np.isfinite(array)):
        raise DriveQualificationError(f"{name} must be finite")
    array.setflags(write=False)
    return array


@dataclass(frozen=True)
class ElasticFieldCacheKey:
    backend_source_hash: str
    mechanics_factory_hash: str
    geometry_fingerprint: str
    mesh_topology_fingerprint: str
    crack_basis_fingerprint: str
    reference_opening_m: float

    def fingerprint(self) -> str:
        return _digest(self.__dict__)


@dataclass(frozen=True)
class ProbeQueryKey:
    field_snapshot_hash: str
    tip_radius_m: float
    front_width_m: float
    source_area_m2: float
    source_geometry_fingerprint: str
    source_positions_fingerprint: str
    slip_system_identity: str
    probe_convention: str

    def fingerprint(self) -> str:
        return _digest(self.__dict__)


@dataclass(frozen=True)
class NormalizedElasticFieldSnapshot:
    backend: str
    source_repository: str
    source_commit: str
    factory_wrapper_identity: str
    mechanics_factory_hash: str
    bulk_constraint: str
    crack_representation: str
    failed_region_state: str
    geometry_fingerprint: str
    mesh_fingerprint: str
    topology_wake_fingerprint: str
    crack_tip_xy_m: tuple[float, float]
    crack_tangent: tuple[float, float]
    crack_normal: tuple[float, float]
    nodes_m: np.ndarray
    elements: np.ndarray
    element_centroids_m: np.ndarray
    element_areas_m2: np.ndarray
    stress_per_opening_Pa_per_m: np.ndarray
    strain_per_opening_per_m: np.ndarray | None
    damage_or_interface_state: np.ndarray
    reaction_per_opening_N_per_m2: float
    elastic_energy_per_opening2_J_per_m3: float
    native_J_per_opening2_J_per_m4: float
    native_KJ_per_opening_Pa_sqrt_m_per_m: float
    native_J_domain_metadata: Mapping[str, Any]
    qualified_structural_G_per_opening2_J_per_m4: float | None = None
    qualified_structural_KG_per_opening_Pa_sqrt_m_per_m: float | None = None
    qualified_structural_G_metadata: Mapping[str, Any] = field(default_factory=dict)
    validity_flags: tuple[str, ...] = ()
    reference_opening_m: float = 1.0

    def __post_init__(self) -> None:
        if self.reference_opening_m <= 0 or not np.isfinite(self.reference_opening_m):
            raise DriveQualificationError("reference opening must be positive and finite")
        for name in ("candidate_specific", "process_zone", "mobile_count", "retained_count", "backstress", "shielding"):
            if name in " ".join(self.validity_flags).lower():
                raise DriveQualificationError("candidate/process-zone state is forbidden in elastic snapshots")
        object.__setattr__(self, "nodes_m", _readonly_array(self.nodes_m, (2,), "nodes_m"))
        elements = np.array(self.elements, dtype=np.int64, copy=True)
        if elements.ndim != 2 or elements.shape[1] != 3:
            raise DriveQualificationError("elements must have shape (n,3)")
        elements.setflags(write=False); object.__setattr__(self, "elements", elements)
        object.__setattr__(self, "element_centroids_m", _readonly_array(self.element_centroids_m, (2,), "element_centroids_m"))
        object.__setattr__(self, "element_areas_m2", _readonly_array(self.element_areas_m2, (), "element_areas_m2"))
        stress = _readonly_array(self.stress_per_opening_Pa_per_m, (2, 2), "stress_per_opening_Pa_per_m")
        if not np.array_equal(stress, np.swapaxes(stress, -1, -2)):
            raise DriveQualificationError("normalized stress tensors must be exactly symmetric")
        object.__setattr__(self, "stress_per_opening_Pa_per_m", stress)
        if self.strain_per_opening_per_m is not None:
            strain = _readonly_array(self.strain_per_opening_per_m, (2, 2), "strain_per_opening_per_m")
            object.__setattr__(self, "strain_per_opening_per_m", strain)
        state = np.array(self.damage_or_interface_state, dtype=float, copy=True)
        if not np.all(np.isfinite(state)):
            raise DriveQualificationError("damage/interface state must be finite")
        state.setflags(write=False); object.__setattr__(self, "damage_or_interface_state", state)
        counts = (len(self.elements), len(self.element_centroids_m), len(self.element_areas_m2), len(stress))
        if len(set(counts)) != 1:
            raise DriveQualificationError(f"element field lengths disagree: {counts}")

    @property
    def cache_key(self) -> ElasticFieldCacheKey:
        basis = _digest({"tangent": self.crack_tangent, "normal": self.crack_normal})
        return ElasticFieldCacheKey(
            self.source_commit, self.mechanics_factory_hash, self.geometry_fingerprint,
            _digest({"mesh": self.mesh_fingerprint, "topology": self.topology_wake_fingerprint}),
            basis, self.reference_opening_m,
        )

    @property
    def snapshot_hash(self) -> str:
        return _digest({
            "key": self.cache_key.__dict__, "nodes": self.nodes_m,
            "elements": self.elements, "stress": self.stress_per_opening_Pa_per_m,
            "state": self.damage_or_interface_state,
        })

    def stress_at_opening(self, opening_m: float) -> np.ndarray:
        if not np.isfinite(opening_m):
            raise DriveQualificationError("opening must be finite")
        return self.stress_per_opening_Pa_per_m * float(opening_m)


class BackendElasticFieldOracle(Protocol):
    solve_count: int
    def solve(self, geometry: Mapping[str, Any]) -> NormalizedElasticFieldSnapshot: ...


class BackendSourceProbeOperator(Protocol):
    query_count: int
    def evaluate(self, snapshot: NormalizedElasticFieldSnapshot,
                 process_zone_state: Mapping[str, Any], opening_m: float) -> NativeDriveBundle: ...


class BackendNativeStateClosure(Protocol):
    def evaluate_state(self, geometry: Mapping[str, Any], process_zone_state: Mapping[str, Any],
                       opening_m: float) -> NativeDriveBundle: ...


class ExactOracleCache:
    """Two independent deterministic caches with explicit counters."""
    def __init__(self):
        self.fields: dict[str, NormalizedElasticFieldSnapshot] = {}
        self.probes: dict[str, NativeDriveBundle] = {}
        self.field_hits = self.field_misses = 0
        self.probe_hits = self.probe_misses = 0

    def field(self, key: ElasticFieldCacheKey,
              solve: Callable[[], NormalizedElasticFieldSnapshot]) -> NormalizedElasticFieldSnapshot:
        token = key.fingerprint()
        if token in self.fields:
            self.field_hits += 1; return self.fields[token]
        self.field_misses += 1
        value = solve()
        if value.cache_key != key:
            raise DriveQualificationError("oracle returned a snapshot for a different field-cache key")
        self.fields[token] = value
        return value

    def probe(self, key: ProbeQueryKey, evaluate: Callable[[], NativeDriveBundle]) -> NativeDriveBundle:
        token = key.fingerprint()
        if token in self.probes:
            self.probe_hits += 1; return self.probes[token]
        self.probe_misses += 1
        value = evaluate(); self.probes[token] = value
        return value


class OptionalRegimeAwareSurrogate:
    """Acceleration layer that remains locked until exact Level 3 qualifies."""
    def __init__(self, exact_level3_status: str):
        self.exact_level3_status = exact_level3_status

    def evaluate(self, *args, **kwargs):
        if self.exact_level3_status != "NATIVE_CLOSURE_QUALIFIED":
            raise DriveQualificationError("surrogate use is forbidden before exact Level-3 qualification")
        raise NotImplementedError("surrogate construction is intentionally deferred")
