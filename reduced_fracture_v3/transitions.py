"""Atomic monotonic transitions for the aligned V3 void lifecycle."""
from __future__ import annotations

from dataclasses import replace
import math
from typing import Any, Mapping

from reduced_fracture_v2.state import ReducedState

from .state import (
    HazardClock,
    LengthLedgers,
    MechanicsRepresentation,
    OneDStateV3,
    OneDVoidState,
    VoidPhase,
)


def _require_void(state: OneDStateV3) -> OneDVoidState:
    if state.void_state is None:
        raise ValueError("transition requires a V3 void state")
    return state.void_state


def _advance_clock(clock: HazardClock, increment: float) -> tuple[HazardClock, bool]:
    value = float(increment)
    if not math.isfinite(value) or value < 0.0:
        raise ValueError("hazard increment must be finite and nonnegative")
    accumulated = min(clock.accumulated + value, clock.threshold)
    return replace(clock, accumulated=accumulated), accumulated == clock.threshold


def advance_birth(
    state: OneDStateV3,
    hazard_increment: float,
    *,
    renewed_threshold: float | None = None,
    rng_state_after_renewal: Mapping[str, Any] | None = None,
) -> tuple[OneDStateV3, str | None]:
    void = _require_void(state)
    if void.phase != VoidPhase.AVAILABLE_SITE:
        raise ValueError("birth clock is owned only by an available site")
    clock, crossed = _advance_clock(void.birth, hazard_increment * void.candidate_weight)
    if not crossed:
        return replace(state, void_state=replace(void, birth=clock)), None
    hits = void.hit_count + 1
    if hits < void.required_hit_count:
        if renewed_threshold is None or rng_state_after_renewal is None:
            raise ValueError("an intermediate birth hit requires owned threshold/RNG renewal")
        renewed = HazardClock(0.0, float(renewed_threshold))
        updated = replace(
            void,
            hit_count=hits,
            birth=renewed,
            rng_state=dict(rng_state_after_renewal),
            lineage=void.lineage + ("BIRTH_HIT",),
        )
        return replace(state, void_state=updated), "BIRTH_HIT"
    updated = replace(
        void,
        hit_count=hits,
        birth=clock,
        phase=VoidPhase.EMBRYO,
        lineage=void.lineage + ("EMBRYO",),
    )
    return replace(state, void_state=updated), "EMBRYO"


def advance_embryo(
    state: OneDStateV3,
    *,
    stabilization_increment: float,
    healing_increment: float,
) -> tuple[OneDStateV3, str | None]:
    void = _require_void(state)
    if void.phase != VoidPhase.EMBRYO:
        raise ValueError("embryo clocks require the embryo phase")
    stable_fraction = void.stabilization.remaining / max(stabilization_increment, 1.0e-300)
    heal_fraction = void.healing.remaining / max(healing_increment, 1.0e-300)
    stabilization, stable = _advance_clock(void.stabilization, stabilization_increment)
    healing, healed = _advance_clock(void.healing, healing_increment)
    if stable and healed:
        healed = heal_fraction <= stable_fraction
        stable = not healed
    if healed:
        updated = replace(
            void,
            stabilization=stabilization,
            healing=healing,
            phase=VoidPhase.HEALED_SITE,
            lineage=void.lineage + ("HEALED",),
        )
        return replace(state, void_state=updated), "HEALED"
    if stable:
        updated = replace(
            void,
            stabilization=stabilization,
            healing=healing,
            phase=VoidPhase.STABLE_SUBGRID_VOID,
            mechanics_representation=MechanicsRepresentation.SUBGRID,
            lineage=void.lineage + ("STABILIZED",),
        )
        # Geometry is created by seed_cavity; keep the phase transition atomic
        # to the same call so an invalid radius cannot publish a partial state.
        return replace(state, void_state=updated), "STABILIZED"
    return replace(
        state,
        void_state=replace(void, stabilization=stabilization, healing=healing),
    ), None


def seed_cavity(state: OneDStateV3, radius_m: float) -> OneDStateV3:
    void = _require_void(state)
    if void.phase != VoidPhase.STABLE_SUBGRID_VOID or void.radius_m is not None:
        raise ValueError("only an unseeded stabilized site can create a cavity")
    radius = float(radius_m)
    if not math.isfinite(radius) or radius <= 0.0:
        raise ValueError("seed radius must be finite and positive")
    area = math.pi * radius**2
    if area > void.available_inventory_area_m2:
        raise ValueError("cavity seed exceeds available inventory")
    updated = replace(
        void,
        radius_m=radius,
        area_m2=area,
        available_inventory_area_m2=void.available_inventory_area_m2 - area,
        consumed_inventory_area_m2=void.consumed_inventory_area_m2 + area,
        lineage=void.lineage + ("INITIAL_CAVITY_SEED_INVENTORY_DEBIT",),
    )
    return replace(state, void_state=updated)


def grow_cavity(state: OneDStateV3, delta_radius_m: float) -> OneDStateV3:
    void = _require_void(state)
    if void.phase not in (VoidPhase.STABLE_SUBGRID_VOID, VoidPhase.RESOLVED_VOID):
        raise ValueError("growth requires an unconnected cavity")
    if void.radius_m is None:
        raise ValueError("growth requires seeded cavity geometry")
    radius = void.radius_m + float(delta_radius_m)
    if not math.isfinite(radius) or radius <= 0.0:
        raise ValueError("growth would consume or invalidate the cavity")
    area = math.pi * radius**2
    delta_area = area - void.area_m2
    if delta_area > void.available_inventory_area_m2:
        raise ValueError("growth exceeds available defect inventory")
    if -delta_area > void.consumed_inventory_area_m2:
        raise ValueError("shrinkage would return unowned inventory")
    updated = replace(
        void,
        radius_m=radius,
        area_m2=area,
        available_inventory_area_m2=void.available_inventory_area_m2 - delta_area,
        consumed_inventory_area_m2=void.consumed_inventory_area_m2 + delta_area,
    )
    return replace(state, void_state=updated)


def promote_cavity(state: OneDStateV3) -> OneDStateV3:
    """Activate the resolved map without changing clocks, geometry, or ledgers."""

    void = _require_void(state)
    if void.phase != VoidPhase.STABLE_SUBGRID_VOID or void.radius_m is None:
        raise ValueError("promotion requires a seeded stable subgrid cavity")
    updated = replace(
        void,
        phase=VoidPhase.RESOLVED_VOID,
        mechanics_representation=MechanicsRepresentation.SURROGATE_ACTIVE,
        lineage=void.lineage + ("MECHANICS_REPRESENTATION_PROMOTION",),
    )
    return replace(state, void_state=updated)


def connect_ligament(
    state: OneDStateV3,
    *,
    active_tip_position_m: float,
    candidate_id: str,
) -> OneDStateV3:
    void = _require_void(state)
    if void.phase != VoidPhase.RESOLVED_VOID or void.near_surface_m is None:
        raise ValueError("ligament rupture requires a resolved preconnection cavity")
    tip = float(active_tip_position_m)
    ligament = void.near_surface_m - tip
    if not math.isfinite(tip) or ligament <= 0.0:
        raise ValueError("ligament must have finite positive material length")
    ledgers = replace(
        void.length_ledgers,
        fractured_material_m=void.length_ledgers.fractured_material_m + ligament,
        front_advance_m=void.length_ledgers.front_advance_m + ligament,
        connected_free_surface_extent_m=2.0 * void.radius_m,
    )
    updated = replace(
        void,
        phase=VoidPhase.CONNECTED_VOID,
        mechanics_representation=MechanicsRepresentation.CONNECTED_CAVITY,
        active_tip_position_m=None,
        ligament=replace(void.ligament, accumulated=void.ligament.threshold),
        length_ledgers=ledgers,
        lineage=void.lineage + (f"CRACK_TO_VOID_LIGAMENT:{candidate_id}",),
    )
    return replace(state, void_state=updated)


def nucleate_downstream_child(
    state: OneDStateV3,
    *,
    child_length_m: float,
    renewed_fracture_state: ReducedState,
    candidate_id: str,
) -> OneDStateV3:
    void = _require_void(state)
    if void.phase != VoidPhase.CONNECTED_VOID or void.far_surface_m is None:
        raise ValueError("downstream nucleation requires a connected dormant cavity")
    child = float(child_length_m)
    if not math.isfinite(child) or child <= 0.0:
        raise ValueError("downstream nucleation requires a finite child segment")
    if math.isclose(renewed_fracture_state.tip_radius_m, void.radius_m, rel_tol=0.0, abs_tol=0.0):
        raise ValueError("r_tip must remain distinct from R_void")
    free_span = 2.0 * void.radius_m
    ledgers = replace(
        void.length_ledgers,
        fractured_material_m=void.length_ledgers.fractured_material_m + child,
        free_span_m=void.length_ledgers.free_span_m + free_span,
        front_advance_m=void.length_ledgers.front_advance_m + free_span + child,
    )
    updated = replace(
        void,
        phase=VoidPhase.DOWNSTREAM_FRONT_ACTIVE,
        mechanics_representation=MechanicsRepresentation.CHILD_ACTIVE,
        active_tip_position_m=void.far_surface_m + child,
        downstream=replace(void.downstream, accumulated=void.downstream.threshold),
        length_ledgers=ledgers,
        lineage=void.lineage + (f"DOWNSTREAM_FIRST_PASSAGE:{candidate_id}",),
    )
    return OneDStateV3(
        fracture_state=renewed_fracture_state,
        void_state=updated,
        fatigue_state=state.fatigue_state,
    )


def continue_downstream_child(state: OneDStateV3, advance_m: float) -> OneDStateV3:
    void = _require_void(state)
    if void.phase != VoidPhase.DOWNSTREAM_FRONT_ACTIVE or void.active_tip_position_m is None:
        raise ValueError("child continuation requires an active downstream sharp tip")
    advance = float(advance_m)
    if not math.isfinite(advance) or advance <= 0.0:
        raise ValueError("child advance must be finite and positive")
    ledgers = replace(
        void.length_ledgers,
        fractured_material_m=void.length_ledgers.fractured_material_m + advance,
        front_advance_m=void.length_ledgers.front_advance_m + advance,
    )
    updated_void = replace(
        void,
        active_tip_position_m=void.active_tip_position_m + advance,
        length_ledgers=ledgers,
        lineage=void.lineage + ("CHILD_CONTINUATION",),
    )
    updated_fracture = replace(
        state.fracture_state,
        crack_extension_m=state.fracture_state.crack_extension_m + advance,
    )
    return replace(state, fracture_state=updated_fracture, void_state=updated_void)


__all__ = [
    "advance_birth",
    "advance_embryo",
    "connect_ligament",
    "continue_downstream_child",
    "grow_cavity",
    "nucleate_downstream_child",
    "promote_cavity",
    "seed_cavity",
]
