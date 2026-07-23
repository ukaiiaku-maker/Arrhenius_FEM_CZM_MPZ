"""Autonomous v9.13 persistent-site one-dimensional R-curve driver.

The constitutive replay in :mod:`emergent_gnd_campaign_v913` accepts an
already-computed ``K, dt, da`` history.  That is useful for checking state
transfer, but it cannot predict a toughness curve.  This module closes that
loop with the two shared pieces of mechanics that are not candidate
parameters:

* a displacement-to-K geometry map, indexed by accepted crack event; and
* a common-random-number sequence of cleavage thresholds and event lengths.

The candidate cleavage, emission, Peierls, Taylor, source-density,
correlation, and blunting parameters are never modified by this driver.
"""

from __future__ import annotations

import copy
from dataclasses import asdict, dataclass, field
import math
from typing import Any, Iterable, Sequence

import numpy as np

from .emergent_gnd_state_v913 import EmergentGNDState
from .emergent_gnd_types_v913 import CandidateParameters, CommonPhysics


@dataclass(frozen=True)
class RCurveLoadingMap:
    """Shared reduced mechanics and stochastic event clock.

    ``K_per_U`` is the Mode-I stress-intensity response per applied
    displacement for the geometry that precedes each accepted event.  Path
    advance convects the one-dimensional process-zone state.  Projected
    advance is the horizontal R-curve abscissa used by the two-dimensional
    campaign.
    """

    K_per_U_MPa_sqrt_m_per_m: tuple[float, ...]
    threshold_actions: tuple[float, ...]
    path_advances_m: tuple[float, ...]
    projected_advances_m: tuple[float, ...]
    nominal_dU_m: float
    nominal_dt_s: float
    seed: int
    reference_candidate_id: str
    reference_temperature_K: float
    reference_event_K_MPa_sqrt_m: tuple[float, ...] = ()
    reference_event_U_m: tuple[float, ...] = ()
    provenance: dict[str, Any] = field(default_factory=dict)

    def validate(self) -> None:
        lengths = {
            len(self.K_per_U_MPa_sqrt_m_per_m),
            len(self.threshold_actions),
            len(self.path_advances_m),
            len(self.projected_advances_m),
        }
        if len(lengths) != 1 or not lengths or next(iter(lengths)) < 1:
            raise ValueError(
                "R-curve loading-map arrays must have equal nonzero length"
            )
        for values, name in (
            (self.K_per_U_MPa_sqrt_m_per_m, "K_per_U"),
            (self.threshold_actions, "threshold_actions"),
            (self.path_advances_m, "path_advances_m"),
            (self.projected_advances_m, "projected_advances_m"),
        ):
            data = np.asarray(values, dtype=float)
            if np.any(~np.isfinite(data)) or np.any(data <= 0.0):
                raise ValueError(f"{name} must contain positive finite values")
        if not math.isfinite(self.nominal_dU_m) or self.nominal_dU_m <= 0.0:
            raise ValueError("nominal_dU_m must be positive and finite")
        if not math.isfinite(self.nominal_dt_s) or self.nominal_dt_s <= 0.0:
            raise ValueError("nominal_dt_s must be positive and finite")
        for values, name in (
            (self.reference_event_K_MPa_sqrt_m, "reference_event_K"),
            (self.reference_event_U_m, "reference_event_U"),
        ):
            if values and len(values) != len(self.threshold_actions):
                raise ValueError(f"{name} must be empty or have one value per event")

    @property
    def displacement_rate_m_s(self) -> float:
        return float(self.nominal_dU_m) / float(self.nominal_dt_s)

    @property
    def n_events(self) -> int:
        return len(self.threshold_actions)

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> "RCurveLoadingMap":
        tuple_fields = (
            "K_per_U_MPa_sqrt_m_per_m",
            "threshold_actions",
            "path_advances_m",
            "projected_advances_m",
            "reference_event_K_MPa_sqrt_m",
            "reference_event_U_m",
        )
        data = dict(payload)
        for name in tuple_fields:
            if name in data:
                data[name] = tuple(data[name])
        out = cls(**data)
        out.validate()
        return out


@dataclass
class RCurveEvent:
    event_index: int
    threshold_action: float
    applied_displacement_m: float
    elapsed_time_s: float
    K_MPa_sqrt_m: float
    path_advance_m: float
    projected_advance_m: float
    cumulative_path_extension_m: float
    cumulative_projected_extension_m: float
    tip_radius_pre_advance_m: float
    tip_radius_post_advance_m: float
    front_width_pre_advance_m: float
    backstress_pre_advance_Pa: float
    source_multiplicity_pre_advance: float
    cumulative_source_activations: float
    cumulative_line_content: float
    integration_substeps: int

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass
class RCurveResult:
    candidate_id: str
    temperature_K: float
    status: str
    seed: int
    target_projected_extension_m: float
    events: list[RCurveEvent] = field(default_factory=list)
    max_tip_radius_m: float = 0.0
    min_front_width_m: float = float("inf")
    max_backstress_Pa: float = 0.0
    max_source_multiplicity: float = 0.0
    final_applied_displacement_m: float = 0.0
    final_elapsed_time_s: float = 0.0
    numerical_integration: dict[str, Any] = field(default_factory=dict)

    @property
    def achieved_projected_extension_m(self) -> float:
        if not self.events:
            return 0.0
        return float(self.events[-1].cumulative_projected_extension_m)

    def checkpoint_K(self, extension_m: float) -> float:
        """Return the first accepted event at or beyond an extension.

        The v10.2.22 campaign reports checkpoint toughness this way; it does
        not interpolate across the finite stochastic crack jumps.
        """
        if not self.events:
            return float("nan")
        target = max(float(extension_m), 0.0)
        x = np.asarray(
            [event.cumulative_projected_extension_m for event in self.events],
            dtype=float,
        )
        y = np.asarray([event.K_MPa_sqrt_m for event in self.events], dtype=float)
        index = int(np.searchsorted(x, target, side="left"))
        return float(y[min(max(index, 0), len(y) - 1)])

    def summary_dict(self) -> dict[str, Any]:
        return {
            "candidate_id": self.candidate_id,
            "temperature_K": self.temperature_K,
            "status": self.status,
            "seed": self.seed,
            "n_events": len(self.events),
            "K_first_MPa_sqrt_m": (
                self.events[0].K_MPa_sqrt_m if self.events else float("nan")
            ),
            "K_10um_MPa_sqrt_m": self.checkpoint_K(10.0e-6),
            "K_25um_MPa_sqrt_m": self.checkpoint_K(25.0e-6),
            "K_50um_MPa_sqrt_m": self.checkpoint_K(50.0e-6),
            "achieved_projected_extension_um": (
                self.achieved_projected_extension_m * 1.0e6
            ),
            "max_backstress_GPa": self.max_backstress_Pa * 1.0e-9,
            "min_front_width_um": self.min_front_width_m * 1.0e6,
            "max_tip_radius_um": self.max_tip_radius_m * 1.0e6,
            "max_source_multiplicity": self.max_source_multiplicity,
            "final_applied_displacement_m": self.final_applied_displacement_m,
            "final_elapsed_time_s": self.final_elapsed_time_s,
        }

    def as_dict(self) -> dict[str, Any]:
        data = self.summary_dict()
        data["target_projected_extension_m"] = self.target_projected_extension_m
        data["events"] = [event.as_dict() for event in self.events]
        data["numerical_integration"] = dict(self.numerical_integration)
        return data


def _state_extrema(state: EmergentGNDState) -> tuple[float, float, float, float]:
    geometry = state.source_geometry()
    backstress = state.backstress_state()[2]
    return (
        float(geometry["tip_radius_m"]),
        float(geometry["front_width_m"]),
        # The 2-D archive reports its equivalent reduced emission backstress,
        # which is the two-channel mean used by the constitutive replay.
        float(np.mean(backstress)),
        float(geometry["multiplicity_per_system"]),
    )


def _update_extrema(result: RCurveResult, state: EmergentGNDState) -> None:
    radius, width, backstress, multiplicity = _state_extrema(state)
    result.max_tip_radius_m = max(result.max_tip_radius_m, radius)
    result.min_front_width_m = min(result.min_front_width_m, width)
    result.max_backstress_Pa = max(result.max_backstress_Pa, backstress)
    result.max_source_multiplicity = max(
        result.max_source_multiplicity,
        multiplicity,
    )


def _adaptive_displacement_increment(
    *,
    state: EmergentGNDState,
    geometry_factor: float,
    displacement_m: float,
    temperature_K: float,
    displacement_rate_m_s: float,
    nominal_dU_m: float,
    max_hazard_increment: float,
) -> float:
    """Choose a load increment using the current-state cleavage clock."""
    dU = float(nominal_dU_m)
    K0 = geometry_factor * displacement_m
    K1 = geometry_factor * (displacement_m + dU)
    rate0 = state.cleavage_rate_s(K0, temperature_K)
    rate1 = state.cleavage_rate_s(K1, temperature_K)
    rate_bound = max(rate0, rate1, 0.0)
    if rate_bound > 0.0:
        dU = min(
            dU,
            max_hazard_increment * displacement_rate_m_s / rate_bound,
        )
    return max(float(dU), 1.0e-14 * float(nominal_dU_m))


def _state_advance_is_resolved(
    state: EmergentGNDState,
    K_MPa_sqrt_m: float,
    temperature_K: float,
    duration_s: float,
    *,
    minimum_expected_activations: float,
) -> bool:
    """Return whether a zero state needs explicit plastic integration.

    At low temperature the expected aggregate emission can be many orders of
    magnitude below one activation over the entire R-curve.  Calling the full
    spatial integrator for thousands of cleavage-clock quadrature points then
    changes no resolved state and dominates runtime.  Once any state exists,
    or the expected activation count is resolvable, the full integrator is
    always used.
    """
    if (
        np.any(state.mobile_m2 > 0.0)
        or np.any(state.retained_m2 > 0.0)
        or np.any(state.accumulated_slip_m2 > 0.0)
    ):
        return True
    rates = state.local_rates(K_MPa_sqrt_m, temperature_K)
    per_site = np.asarray(rates["emission_rate_per_site_s"], dtype=float)
    multiplicity = float(state.source_geometry()["multiplicity_per_system"])
    expected = (
        max(float(duration_s), 0.0)
        * multiplicity
        * float(np.sum(np.maximum(per_site, 0.0)))
    )
    return expected >= float(minimum_expected_activations)


def run_autonomous_rcurve(
    candidate: CandidateParameters,
    physics: CommonPhysics,
    loading_map: RCurveLoadingMap,
    temperature_K: float,
    *,
    target_projected_extension_m: float = 50.0e-6,
    max_hazard_increment: float = 0.05,
    maximum_applied_displacement_m: float = 2.0e-4,
    maximum_integration_substeps: int = 200_000,
    minimum_resolved_source_activations: float = 1.0e-4,
    translation_mode: str = "hazard_coupled",
    translation_action_exponent: float = 1.0,
) -> RCurveResult:
    """Predict one R-curve without supplying a 2-D K/dt/da history.

    Loading is displacement controlled.  The current event geometry converts
    applied displacement to K.  The persistent-site state evolves during each
    load increment, cleavage hazard is integrated to the next CRN threshold,
    and the resulting event length translates the process-zone state.
    """
    loading_map.validate()
    physics.validate()
    if target_projected_extension_m <= 0.0:
        raise ValueError("target_projected_extension_m must be positive")
    if not 0.0 < max_hazard_increment <= 0.25:
        raise ValueError("max_hazard_increment must lie in (0, 0.25]")
    if translation_mode not in ("hazard_coupled", "event_commit"):
        raise ValueError("translation_mode must be 'hazard_coupled' or 'event_commit'")
    if (
        not math.isfinite(float(translation_action_exponent))
        or translation_action_exponent <= 0.0
    ):
        raise ValueError("translation_action_exponent must be positive and finite")
    state = EmergentGNDState(candidate, physics)
    result = RCurveResult(
        candidate_id=candidate.candidate_id,
        temperature_K=float(temperature_K),
        status="running",
        seed=int(loading_map.seed),
        target_projected_extension_m=float(target_projected_extension_m),
        numerical_integration={
            **state.integration_metadata(),
            "driver": "v9.13_autonomous_displacement_controlled_first_passage",
            "loading_map_reference_candidate": (loading_map.reference_candidate_id),
            "loading_map_reference_temperature_K": (
                loading_map.reference_temperature_K
            ),
            "max_hazard_increment": float(max_hazard_increment),
            "minimum_resolved_source_activations": float(
                minimum_resolved_source_activations
            ),
            "translation_mode": translation_mode,
            "translation_action_exponent": float(translation_action_exponent),
            "candidate_parameters_modified_by_driver": False,
        },
    )
    _update_extrema(result, state)

    displacement = 0.0
    elapsed = 0.0
    path_extension = 0.0
    projected_extension = 0.0
    displacement_rate = loading_map.displacement_rate_m_s
    total_substeps = 0

    for event_index in range(loading_map.n_events):
        if projected_extension >= target_projected_extension_m:
            break
        threshold = float(loading_map.threshold_actions[event_index])
        geometry_factor = float(loading_map.K_per_U_MPa_sqrt_m_per_m[event_index])
        event_path_advance = float(loading_map.path_advances_m[event_index])
        event_geometry_extension = float(state.extension_m)
        path_advance_committed = 0.0
        accumulated_hazard = 0.0
        event_substeps = 0
        event_translation_coupled = translation_mode == "hazard_coupled"

        while accumulated_hazard < threshold:
            if displacement >= maximum_applied_displacement_m:
                result.status = "right_censored_maximum_displacement"
                result.final_applied_displacement_m = displacement
                result.final_elapsed_time_s = elapsed
                return result
            total_substeps += 1
            event_substeps += 1
            if total_substeps > maximum_integration_substeps:
                result.status = "right_censored_maximum_substeps"
                result.final_applied_displacement_m = displacement
                result.final_elapsed_time_s = elapsed
                return result

            K0 = geometry_factor * displacement
            rate0 = state.cleavage_rate_s(K0, temperature_K)
            dU = _adaptive_displacement_increment(
                state=state,
                geometry_factor=geometry_factor,
                displacement_m=displacement,
                temperature_K=temperature_K,
                displacement_rate_m_s=displacement_rate,
                nominal_dU_m=loading_map.nominal_dU_m,
                max_hazard_increment=max_hazard_increment,
            )
            dt = dU / displacement_rate
            K1 = geometry_factor * (displacement + dU)
            K_mid = 0.5 * (K0 + K1)
            rate1_predictor = state.cleavage_rate_s(K1, temperature_K)
            dH_predictor = 0.5 * (rate0 + rate1_predictor) * dt
            da_step = 0.0
            if event_translation_coupled:
                action_end = min(
                    accumulated_hazard + dH_predictor,
                    threshold,
                )
                desired_path_fraction = max(action_end / threshold, 0.0) ** float(
                    translation_action_exponent
                )
                da_step = min(
                    event_path_advance - path_advance_committed,
                    max(
                        event_path_advance * desired_path_fraction
                        - path_advance_committed,
                        0.0,
                    ),
                )
            resolved_state_advance = _state_advance_is_resolved(
                state,
                K_mid,
                temperature_K,
                dt,
                minimum_expected_activations=(minimum_resolved_source_activations),
            )
            state_before = (
                copy.deepcopy(state)
                if resolved_state_advance or da_step > 0.0
                else None
            )
            if resolved_state_advance and da_step > 0.0:
                state.advance_coupled_segment(
                    duration_s=dt,
                    da_m=da_step,
                    K_start_MPa_sqrt_m=K0,
                    K_end_MPa_sqrt_m=K1,
                    T_K=temperature_K,
                    geometry_extension_override_m=event_geometry_extension,
                )
            elif resolved_state_advance:
                state.advance_time(dt, K_mid, temperature_K)
            elif da_step > 0.0:
                state.time_s += dt
                state.translate_tip(da_step)
            rate1 = state.cleavage_rate_s(K1, temperature_K)
            dH = 0.5 * (rate0 + rate1) * dt

            if accumulated_hazard + dH >= threshold:
                fraction = float(
                    np.clip(
                        (threshold - accumulated_hazard) / max(dH, 1.0e-300),
                        0.0,
                        1.0,
                    )
                )
                if state_before is not None:
                    state = state_before
                dU *= fraction
                dt *= fraction
                K1 = geometry_factor * (displacement + dU)
                da_final = (
                    event_path_advance - path_advance_committed
                    if event_translation_coupled
                    else 0.0
                )
                if resolved_state_advance and da_final > 0.0 and dt > 0.0:
                    state.advance_coupled_segment(
                        duration_s=dt,
                        da_m=da_final,
                        K_start_MPa_sqrt_m=K0,
                        K_end_MPa_sqrt_m=K1,
                        T_K=temperature_K,
                        geometry_extension_override_m=event_geometry_extension,
                    )
                elif resolved_state_advance and dt > 0.0:
                    state.advance_time(
                        dt,
                        0.5 * (K0 + K1),
                        temperature_K,
                    )
                elif da_final > 0.0:
                    state.time_s += dt
                    state.translate_tip(da_final)
                elif dt > 0.0:
                    state.time_s += dt
                path_advance_committed += da_final
                displacement += dU
                elapsed += dt
                accumulated_hazard = threshold
            else:
                path_advance_committed += da_step
                if not resolved_state_advance:
                    if da_step <= 0.0:
                        state.time_s += dt
                displacement += dU
                elapsed += dt
                accumulated_hazard += dH
            _update_extrema(result, state)

        geometry_pre = state.source_geometry()
        backstress_pre = float(np.mean(state.backstress_state()[2]))
        radius_pre = float(geometry_pre["tip_radius_m"])
        path_advance = event_path_advance
        projected_advance = float(loading_map.projected_advances_m[event_index])
        event_K = geometry_factor * displacement

        if not event_translation_coupled:
            state.translate_tip(path_advance)
        elif path_advance_committed < path_advance:
            state.translate_tip(path_advance - path_advance_committed)
        path_extension += path_advance
        projected_extension += projected_advance
        _update_extrema(result, state)

        result.events.append(
            RCurveEvent(
                event_index=event_index,
                threshold_action=threshold,
                applied_displacement_m=displacement,
                elapsed_time_s=elapsed,
                K_MPa_sqrt_m=event_K,
                path_advance_m=path_advance,
                projected_advance_m=projected_advance,
                cumulative_path_extension_m=path_extension,
                cumulative_projected_extension_m=projected_extension,
                tip_radius_pre_advance_m=radius_pre,
                tip_radius_post_advance_m=state.tip_radius_m(),
                front_width_pre_advance_m=float(geometry_pre["front_width_m"]),
                backstress_pre_advance_Pa=backstress_pre,
                source_multiplicity_pre_advance=float(
                    geometry_pre["multiplicity_per_system"]
                ),
                cumulative_source_activations=float(
                    np.sum(state.cumulative_source_activations)
                ),
                cumulative_line_content=float(np.sum(state.cumulative_line_content)),
                integration_substeps=event_substeps,
            )
        )

    result.final_applied_displacement_m = displacement
    result.final_elapsed_time_s = elapsed
    result.status = (
        "complete"
        if projected_extension >= target_projected_extension_m
        else "right_censored_loading_map_exhausted"
    )
    return result


def event_rows(results: Iterable[RCurveResult]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for result in results:
        for event in result.events:
            rows.append(
                {
                    "candidate_id": result.candidate_id,
                    "temperature_K": result.temperature_K,
                    "status": result.status,
                    "seed": result.seed,
                    **event.as_dict(),
                }
            )
    return rows


def checkpoint_names() -> Sequence[tuple[str, float]]:
    return (
        ("K_first_MPa_sqrt_m", 0.0),
        ("K_10um_MPa_sqrt_m", 10.0e-6),
        ("K_25um_MPa_sqrt_m", 25.0e-6),
        ("K_50um_MPa_sqrt_m", 50.0e-6),
    )


__all__ = [
    "RCurveEvent",
    "RCurveLoadingMap",
    "RCurveResult",
    "checkpoint_names",
    "event_rows",
    "run_autonomous_rcurve",
]
