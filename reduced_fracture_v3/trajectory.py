"""Small deterministic monotonic V3 state/event prototype."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

from reduced_fracture_v2.state import ReducedState

from .state import OneDStateV3, fingerprint
from .transitions import (
    advance_birth,
    advance_embryo,
    connect_ligament,
    continue_downstream_child,
    grow_cavity,
    nucleate_downstream_child,
    promote_cavity,
    seed_cavity,
)


@dataclass(frozen=True)
class MonotonicPrototypePlan:
    renewed_birth_threshold: float
    rng_state_after_birth_renewal: Mapping[str, Any]
    seed_radius_m: float
    promotion_radius_m: float
    connection_radius_m: float
    root_tip_position_m: float
    ligament_candidate_id: str
    downstream_candidate_id: str
    child_length_m: float
    child_continuation_m: float
    renewed_child_fracture_state: ReducedState


def _record(state: OneDStateV3, operation: str) -> dict[str, Any]:
    void = state.void_state
    return {
        "operation": operation,
        "phase": None if void is None else void.phase.value,
        "mechanics_representation": (
            None if void is None else void.mechanics_representation.value
        ),
        "state_fingerprint": fingerprint(state),
        "active_tip_position_m": None if void is None else void.active_tip_position_m,
        "length_ledgers": None if void is None else {
            "fractured_material_m": void.length_ledgers.fractured_material_m,
            "free_span_m": void.length_ledgers.free_span_m,
            "front_advance_m": void.length_ledgers.front_advance_m,
            "connected_free_surface_extent_m": (
                void.length_ledgers.connected_free_surface_extent_m
            ),
        },
    }


def run_monotonic_prototype(
    initial: OneDStateV3, plan: MonotonicPrototypePlan
) -> tuple[OneDStateV3, tuple[dict[str, Any], ...]]:
    """Commit one bounded event chain from available site through continuation."""

    if initial.void_state is None:
        raise ValueError("prototype requires one available V3 void site")
    state = initial
    trace = [_record(state, "INITIAL")]

    increment = state.void_state.birth.remaining / state.void_state.candidate_weight
    state, event = advance_birth(
        state,
        increment,
        renewed_threshold=plan.renewed_birth_threshold,
        rng_state_after_renewal=plan.rng_state_after_birth_renewal,
    )
    if event != "BIRTH_HIT":
        raise RuntimeError("first localized birth crossing did not create the first hit")
    trace.append(_record(state, event))

    increment = state.void_state.birth.remaining / state.void_state.candidate_weight
    state, event = advance_birth(state, increment)
    if event != "EMBRYO":
        raise RuntimeError("second localized birth crossing did not create an embryo")
    trace.append(_record(state, event))

    state, event = advance_embryo(
        state,
        stabilization_increment=state.void_state.stabilization.remaining,
        healing_increment=0.0,
    )
    if event != "STABILIZED":
        raise RuntimeError("localized embryo competition did not stabilize")
    trace.append(_record(state, event))

    state = seed_cavity(state, plan.seed_radius_m)
    trace.append(_record(state, "INITIAL_CAVITY_SEED"))
    state = grow_cavity(state, plan.promotion_radius_m - plan.seed_radius_m)
    trace.append(_record(state, "SUBGRID_GROWTH"))
    state = promote_cavity(state)
    trace.append(_record(state, "SURROGATE_PROMOTION"))
    state = grow_cavity(state, plan.connection_radius_m - plan.promotion_radius_m)
    trace.append(_record(state, "RESOLVED_GROWTH"))
    state = connect_ligament(
        state,
        active_tip_position_m=plan.root_tip_position_m,
        candidate_id=plan.ligament_candidate_id,
    )
    trace.append(_record(state, "LIGAMENT_FIRST_PASSAGE"))
    state = nucleate_downstream_child(
        state,
        child_length_m=plan.child_length_m,
        renewed_fracture_state=plan.renewed_child_fracture_state,
        candidate_id=plan.downstream_candidate_id,
    )
    trace.append(_record(state, "DOWNSTREAM_FIRST_PASSAGE"))
    state = continue_downstream_child(state, plan.child_continuation_m)
    trace.append(_record(state, "CHILD_CONTINUATION"))
    return state, tuple(trace)


__all__ = ["MonotonicPrototypePlan", "run_monotonic_prototype"]
