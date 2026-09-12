"""Versioned state for the aligned one-dimensional V3 voiding reduction."""
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from enum import Enum
import hashlib
import json
import math
from typing import Any, Mapping

from reduced_fracture_v2.state import ReducedState


SCHEMA = "oneD.v3.stateful-voiding/1"


class VoidPhase(str, Enum):
    """The production V5 lifecycle, retained without a damage projection."""

    AVAILABLE_SITE = "AVAILABLE_SITE"
    EMBRYO = "EMBRYO"
    HEALED_SITE = "HEALED_SITE"
    STABLE_SUBGRID_VOID = "STABLE_SUBGRID_VOID"
    RESOLVED_VOID = "RESOLVED_VOID"
    CONNECTED_VOID = "CONNECTED_VOID"
    DOWNSTREAM_FRONT_ACTIVE = "DOWNSTREAM_FRONT_ACTIVE"
    MERGED_OR_CONSUMED = "MERGED_OR_CONSUMED"


class MechanicsRepresentation(str, Enum):
    """Numerical representation, separate from physical lifecycle state."""

    INACTIVE = "INACTIVE"
    SUBGRID = "SUBGRID"
    SURROGATE_ACTIVE = "SURROGATE_ACTIVE"
    CONNECTED_CAVITY = "CONNECTED_CAVITY"
    CHILD_ACTIVE = "CHILD_ACTIVE"
    CONSUMED = "CONSUMED"


@dataclass(frozen=True)
class HazardClock:
    accumulated: float = 0.0
    threshold: float = 1.0

    def __post_init__(self) -> None:
        if not math.isfinite(self.accumulated) or self.accumulated < 0.0:
            raise ValueError("hazard accumulation must be finite and nonnegative")
        if not math.isfinite(self.threshold) or self.threshold <= 0.0:
            raise ValueError("hazard threshold must be finite and positive")
        if self.accumulated > self.threshold:
            raise ValueError("hazard accumulation cannot exceed its threshold")

    @property
    def remaining(self) -> float:
        return max(self.threshold - self.accumulated, 0.0)


@dataclass(frozen=True)
class LengthLedgers:
    """Separate newly fractured, pre-existing free-span, and front ledgers."""

    fractured_material_m: float = 0.0
    free_span_m: float = 0.0
    front_advance_m: float = 0.0
    connected_free_surface_extent_m: float = 0.0

    def __post_init__(self) -> None:
        for name, value in asdict(self).items():
            if not math.isfinite(value) or value < 0.0:
                raise ValueError(f"{name} must be finite and nonnegative")


@dataclass(frozen=True)
class OneDVoidState:
    site_id: str
    x_void_m: float
    phase: VoidPhase = VoidPhase.AVAILABLE_SITE
    mechanics_representation: MechanicsRepresentation = MechanicsRepresentation.INACTIVE
    radius_m: float | None = None
    area_m2: float = 0.0
    hit_count: int = 0
    required_hit_count: int = 2
    candidate_weight: float = 1.0
    birth: HazardClock = field(default_factory=HazardClock)
    stabilization: HazardClock = field(default_factory=HazardClock)
    healing: HazardClock = field(default_factory=HazardClock)
    ligament: HazardClock = field(default_factory=HazardClock)
    downstream: HazardClock = field(default_factory=HazardClock)
    available_inventory_area_m2: float = 0.0
    consumed_inventory_area_m2: float = 0.0
    lineage: tuple[str, ...] = ()
    rng_state: Mapping[str, Any] = field(default_factory=dict)
    length_ledgers: LengthLedgers = field(default_factory=LengthLedgers)
    active_tip_position_m: float | None = None
    source_v5_state_sha256: str | None = None
    schema: str = SCHEMA

    def __post_init__(self) -> None:
        if self.schema != SCHEMA:
            raise ValueError(f"unsupported V3 state schema: {self.schema}")
        if not self.site_id:
            raise ValueError("site_id is required")
        if not math.isfinite(self.x_void_m):
            raise ValueError("void position must be finite")
        if self.hit_count < 0 or self.required_hit_count <= 0:
            raise ValueError("hit counts must be nonnegative with a positive target")
        if self.hit_count > self.required_hit_count:
            raise ValueError("hit count exceeds its target")
        if not math.isfinite(self.candidate_weight) or self.candidate_weight <= 0.0:
            raise ValueError("candidate weight must be finite and positive")
        for name in ("available_inventory_area_m2", "consumed_inventory_area_m2"):
            value = float(getattr(self, name))
            if not math.isfinite(value) or value < 0.0:
                raise ValueError(f"{name} must be finite and nonnegative")

        cavity_phases = {
            VoidPhase.RESOLVED_VOID,
            VoidPhase.CONNECTED_VOID,
            VoidPhase.DOWNSTREAM_FRONT_ACTIVE,
            VoidPhase.MERGED_OR_CONSUMED,
        }
        if self.phase in cavity_phases:
            if self.radius_m is None or not math.isfinite(self.radius_m) or self.radius_m <= 0.0:
                raise ValueError("cavity phases require a finite positive radius")
            expected = math.pi * self.radius_m**2
            if not math.isclose(self.area_m2, expected, rel_tol=1.0e-12, abs_tol=1.0e-24):
                raise ValueError("one-dimensional void area retains the 2-D pi*R^2 convention")
        elif self.phase == VoidPhase.STABLE_SUBGRID_VOID:
            if self.radius_m is None:
                if self.area_m2 != 0.0:
                    raise ValueError("an unseeded stable site cannot own cavity area")
            else:
                if not math.isfinite(self.radius_m) or self.radius_m <= 0.0:
                    raise ValueError("seeded cavity radius must be finite and positive")
                expected = math.pi * self.radius_m**2
                if not math.isclose(self.area_m2, expected, rel_tol=1.0e-12, abs_tol=1.0e-24):
                    raise ValueError("one-dimensional void area retains the 2-D pi*R^2 convention")
        elif self.radius_m is not None or self.area_m2 != 0.0:
            raise ValueError("pre-cavity phases cannot own cavity geometry")

        allowed = {
            VoidPhase.AVAILABLE_SITE: {MechanicsRepresentation.INACTIVE},
            VoidPhase.EMBRYO: {MechanicsRepresentation.INACTIVE},
            VoidPhase.HEALED_SITE: {MechanicsRepresentation.INACTIVE},
            VoidPhase.STABLE_SUBGRID_VOID: {MechanicsRepresentation.SUBGRID},
            VoidPhase.RESOLVED_VOID: {MechanicsRepresentation.SURROGATE_ACTIVE},
            VoidPhase.CONNECTED_VOID: {MechanicsRepresentation.CONNECTED_CAVITY},
            VoidPhase.DOWNSTREAM_FRONT_ACTIVE: {MechanicsRepresentation.CHILD_ACTIVE},
            VoidPhase.MERGED_OR_CONSUMED: {MechanicsRepresentation.CONSUMED},
        }
        if self.mechanics_representation not in allowed[self.phase]:
            raise ValueError("physical phase and mechanics representation are inconsistent")
        if self.phase == VoidPhase.DOWNSTREAM_FRONT_ACTIVE:
            if self.active_tip_position_m is None or not math.isfinite(self.active_tip_position_m):
                raise ValueError("a downstream child requires an active tip position")
        elif self.active_tip_position_m is not None:
            raise ValueError("only an active downstream child can own a V3 tip position")

    @property
    def near_surface_m(self) -> float | None:
        return None if self.radius_m is None else self.x_void_m - self.radius_m

    @property
    def far_surface_m(self) -> float | None:
        return None if self.radius_m is None else self.x_void_m + self.radius_m


@dataclass(frozen=True)
class OneDStateV3:
    """Composition keeps fracture, void, and later fatigue ownership separate."""

    fracture_state: ReducedState
    void_state: OneDVoidState | None = None
    fatigue_state: Mapping[str, Any] | None = None
    schema: str = SCHEMA

    def __post_init__(self) -> None:
        if self.schema != SCHEMA:
            raise ValueError(f"unsupported V3 composite schema: {self.schema}")


def no_void_upgrade(fracture_state: ReducedState) -> OneDStateV3:
    """Wrap V2 without copying or translating any V2-owned state."""

    return OneDStateV3(fracture_state=fracture_state)


def no_void_downgrade(state: OneDStateV3) -> ReducedState:
    if state.void_state is not None or state.fatigue_state is not None:
        raise ValueError("only a no-void monotonic V3 state can be downgraded to V2")
    return state.fracture_state


def _encode(value: Any) -> Any:
    if isinstance(value, Enum):
        return value.value
    if hasattr(value, "__dataclass_fields__"):
        return {name: _encode(getattr(value, name)) for name in value.__dataclass_fields__}
    if isinstance(value, Mapping):
        return {str(key): _encode(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_encode(item) for item in value]
    return value


def serialize(state: OneDStateV3) -> str:
    return json.dumps(_encode(state), sort_keys=True, separators=(",", ":"), allow_nan=False)


def fingerprint(state: OneDStateV3) -> str:
    return hashlib.sha256(serialize(state).encode()).hexdigest()


def deserialize(payload: str) -> OneDStateV3:
    raw = json.loads(payload)
    if raw.get("schema") != SCHEMA:
        raise ValueError("unsupported V3 checkpoint schema")
    fracture = ReducedState(**raw["fracture_state"])
    void_raw = raw.get("void_state")
    void = None
    if void_raw is not None:
        void = OneDVoidState(
            **{
                **void_raw,
                "phase": VoidPhase(void_raw["phase"]),
                "mechanics_representation": MechanicsRepresentation(
                    void_raw["mechanics_representation"]
                ),
                "birth": HazardClock(**void_raw["birth"]),
                "stabilization": HazardClock(**void_raw["stabilization"]),
                "healing": HazardClock(**void_raw["healing"]),
                "ligament": HazardClock(**void_raw["ligament"]),
                "downstream": HazardClock(**void_raw["downstream"]),
                "lineage": tuple(void_raw["lineage"]),
                "length_ledgers": LengthLedgers(**void_raw["length_ledgers"]),
            }
        )
    return OneDStateV3(
        fracture_state=fracture,
        void_state=void,
        fatigue_state=raw.get("fatigue_state"),
        schema=raw["schema"],
    )


__all__ = [
    "HazardClock",
    "LengthLedgers",
    "MechanicsRepresentation",
    "OneDStateV3",
    "OneDVoidState",
    "SCHEMA",
    "VoidPhase",
    "deserialize",
    "fingerprint",
    "no_void_downgrade",
    "no_void_upgrade",
    "serialize",
]
