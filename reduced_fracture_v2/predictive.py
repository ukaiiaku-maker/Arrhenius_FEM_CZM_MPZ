"""Provider-native predictive V2 mechanics and fast process-zone closure."""
from __future__ import annotations

from dataclasses import dataclass
import copy
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd

from arrhenius_fracture.emergent_gnd_rcurve_v913 import RCurveLoadingMap
from arrhenius_fracture.emergent_gnd_types_v913 import CandidateParameters, CommonPhysics
from arrhenius_fracture.zero_d_persistent_v913 import (
    ZeroDState,
    _advance_state,
    _state_geometry,
    _translation_retention_factor,
    reduction_geometry,
)


class SourceDriveMapBoundsError(ValueError):
    """A source-drive query lies outside the deterministically probed map."""


@dataclass(frozen=True)
class ProviderMechanicsMap:
    backend: str
    extension_um: np.ndarray
    native_K_per_opening_MPa_sqrt_m_per_m: np.ndarray
    native_J_per_opening2_J_per_m4: np.ndarray
    qualified_K_per_opening_MPa_sqrt_m_per_m: np.ndarray | None = None
    qualified_G_per_opening2_J_per_m4: np.ndarray | None = None

    @classmethod
    def from_csv(cls, backend: str, native: str | Path, qualified: str | Path | None = None):
        frame = pd.read_csv(native)
        x, indices = np.unique(frame["actual_extension_um"].to_numpy(float), return_index=True)
        k = frame["KJ_native_over_U"].to_numpy(float)[indices]
        j = frame["J_native_over_U2"].to_numpy(float)[indices]
        qk = qg = None
        if qualified is not None:
            q = pd.read_csv(qualified)
            qx, qi = np.unique(q["actual_extension_um"].to_numpy(float), return_index=True)
            if not np.array_equal(x, qx):
                raise ValueError("native and qualified mechanics grids differ")
            qk = q["KG_over_U"].to_numpy(float)[qi]
            qg = q["G_over_U2"].to_numpy(float)[qi]
        return cls(str(backend), x, k, j, qk, qg)

    def _at(self, values: np.ndarray, extension_m: float) -> float:
        x = float(extension_m) * 1.0e6
        if x < self.extension_um[0] - 1.0e-12 or x > self.extension_um[-1] + 1.0e-12:
            raise ValueError("mechanics query outside qualified map")
        return float(np.interp(x, self.extension_um, values))

    def native_K_per_opening(self, extension_m: float) -> float:
        return self._at(self.native_K_per_opening_MPa_sqrt_m_per_m, extension_m)

    def metrics(self, extension_m: float, opening_m: float) -> dict[str, float | None]:
        U = float(opening_m)
        out: dict[str, float | None] = {
            "native_KJ_MPa_sqrt_m": self.native_K_per_opening(extension_m) * U,
            "native_J_J_per_m2": self._at(self.native_J_per_opening2_J_per_m4, extension_m) * U * U,
            "qualified_KG_MPa_sqrt_m": None,
            "qualified_G_J_per_m2": None,
        }
        if self.qualified_K_per_opening_MPa_sqrt_m_per_m is not None:
            out["qualified_KG_MPa_sqrt_m"] = self._at(
                self.qualified_K_per_opening_MPa_sqrt_m_per_m, extension_m
            ) * U
            out["qualified_G_J_per_m2"] = self._at(
                self.qualified_G_per_opening2_J_per_m4, extension_m
            ) * U * U
        return out


@dataclass(frozen=True)
class SourceDriveMap:
    backend: str
    extension_um: np.ndarray
    radius_um: np.ndarray
    system_factors: tuple[np.ndarray, ...]

    @classmethod
    def from_csv(cls, path: str | Path):
        frame = pd.read_csv(path)
        extensions = np.sort(frame.extension_um.unique().astype(float))
        radii = np.sort(frame.tip_radius_um.unique().astype(float))
        systems = []
        index = frame.pivot(index="extension_um", columns="tip_radius_um")
        number = 0
        while f"drive_factor_system_{number}" in frame.columns:
            systems.append(
                index[f"drive_factor_system_{number}"].loc[extensions, radii].to_numpy(float)
            )
            number += 1
        if not systems:
            raise ValueError("source drive map contains no systems")
        if not frame.reliable.astype(bool).all():
            raise ValueError("source drive map contains unreliable probes")
        return cls(str(frame.backend.iloc[0]), extensions, radii, tuple(systems))

    def evaluate(self, extension_m: float, radius_m: float) -> tuple[float, ...]:
        x = float(extension_m * 1.0e6)
        r = float(radius_m * 1.0e6)
        if (
            x < self.extension_um[0] - 1.0e-12
            or x > self.extension_um[-1] + 1.0e-12
            or r < self.radius_um[0] - 1.0e-12
            or r > self.radius_um[-1] + 1.0e-12
        ):
            raise SourceDriveMapBoundsError(
                "source-drive query outside qualified map: "
                f"extension_um={x:.12g}, radius_um={r:.12g}"
            )
        values = []
        for grid in self.system_factors:
            along_radius = np.asarray([
                np.interp(r, self.radius_um, row) for row in grid
            ])
            values.append(float(np.interp(x, self.extension_um, along_radius)))
        return tuple(values)


def stochastic_loading_map(
    mechanics: ProviderMechanicsMap,
    *,
    seed: int,
    target_extension_m: float,
    nominal_advance_m: float = 5.0e-6,
    nominal_dU_m: float = 0.2e-6,
    nominal_dt_s: float = 8.4,
) -> RCurveLoadingMap:
    rng = np.random.default_rng(int(seed))
    thresholds: list[float] = []
    advances: list[float] = []
    coefficients: list[float] = []
    extension = 0.0
    mean_clip = 0.5 + math.exp(-0.5) - math.exp(-4.0)
    while extension < target_extension_m:
        xi = max(-math.log(max(float(rng.random()), float(np.finfo(float).tiny))), 1.0e-12)
        factor = min(max(xi, 0.5), 4.0) / mean_clip
        advance = nominal_advance_m * factor
        thresholds.append(xi)
        advances.append(advance)
        coefficients.append(mechanics.native_K_per_opening(extension))
        extension += advance
    return RCurveLoadingMap(
        tuple(coefficients), tuple(thresholds), tuple(advances), tuple(advances),
        nominal_dU_m, nominal_dt_s, int(seed), "candidate_independent", 1000.0,
        provenance={"event_length": "threshold_scaled_mean_preserving", "mechanics_backend": mechanics.backend},
    )


def run_zero_d_predictive(
    candidate: CandidateParameters,
    physics: CommonPhysics,
    mechanics: ProviderMechanicsMap,
    drive: SourceDriveMap,
    loading: RCurveLoadingMap,
    temperature_K: float,
    *,
    target_extension_m: float,
    maximum_opening_m: float = 500.0e-6,
    maximum_intervals: int = 500_000,
    maximum_hazard_increment: float = 0.05,
) -> dict[str, Any]:
    """Fast adaptive first-passage screen with provider-owned mechanics.

    This is the repository's qualified zero-D homolog of the spatial
    persistent-site equation.  It uses exact source barriers and the same
    implicit backstress root, with mechanics and normalized source-drive maps
    supplied independently by the selected backend.
    """
    reduced = reduction_geometry(physics)
    nsys = int(physics.n_systems)
    state = ZeroDState(
        mobile_m2=np.zeros(nsys), retained_m2=np.zeros(nsys),
        local_slip_count_by_system=np.zeros(nsys),
        cumulative_activations=np.zeros(nsys), max_tip_radius_m=float(physics.r0_m),
    )
    opening = time_s = extension = 0.0
    opening_rate = loading.displacement_rate_m_s
    events: list[dict[str, Any]] = []
    intervals = 0
    for event_index, threshold in enumerate(loading.threshold_actions):
        if extension >= target_extension_m - 1.0e-15:
            break
        action = 0.0
        while action < threshold:
            intervals += 1
            if intervals > maximum_intervals or opening >= maximum_opening_m:
                return _predictive_result(candidate, mechanics, temperature_K, target_extension_m,
                                          events, state, opening, time_s, extension, intervals,
                                          "RIGHT_CENSORED_NUMERICAL_BOUND")
            geometry = _state_geometry(candidate, physics, reduced, state)
            try:
                factors = drive.evaluate(extension, float(geometry["tip_radius_m"]))
            except SourceDriveMapBoundsError:
                return _predictive_result(
                    candidate, mechanics, temperature_K, target_extension_m,
                    events, state, opening, time_s, extension, intervals,
                    "RIGHT_CENSORED_DRIVE_MAP_BOUND",
                )
            coefficient = mechanics.native_K_per_opening(extension)
            K0 = coefficient * opening
            rate0 = _advance_state(candidate, physics, reduced, state, duration_s=0.0,
                                   K_MPa_sqrt_m=K0, temperature_K=temperature_K,
                                   drive_factors_override=factors)
            dU = float(loading.nominal_dU_m)
            if rate0 > 0.0:
                dU = min(dU, maximum_hazard_increment * opening_rate / rate0)
            dU = max(dU, loading.nominal_dU_m * 1.0e-12)
            dt = dU / opening_rate
            backup = copy.deepcopy(state)
            K1 = coefficient * (opening + dU)
            rate1 = _advance_state(candidate, physics, reduced, state, duration_s=dt,
                                   K_MPa_sqrt_m=0.5 * (K0 + K1), temperature_K=temperature_K,
                                   drive_factors_override=factors)
            dH = 0.5 * (rate0 + rate1) * dt
            if action + dH >= threshold:
                fraction = float(np.clip((threshold - action) / max(dH, 1.0e-300), 0.0, 1.0))
                state = backup
                dU *= fraction
                dt *= fraction
                K1 = coefficient * (opening + dU)
                _advance_state(candidate, physics, reduced, state, duration_s=dt,
                               K_MPa_sqrt_m=0.5 * (K0 + K1), temperature_K=temperature_K,
                               drive_factors_override=factors)
                opening += dU
                time_s += dt
                action = float(threshold)
            else:
                opening += dU
                time_s += dt
                action += dH
        geometry = _state_geometry(candidate, physics, reduced, state)
        proposed = float(loading.path_advances_m[event_index])
        committed = min(proposed, max(target_extension_m - extension, 0.0))
        metrics = mechanics.metrics(extension, opening)
        prior_opening = None if not events else float(events[-1]["event_opening_m"])
        reload_increment = None if prior_opening is None else opening - prior_opening
        reload = bool(reload_increment is not None and reload_increment >= loading.nominal_dU_m * (1.0 - 1.0e-9))
        avalanche = 0 if not events else int(events[-1]["physical_avalanche_index"]) + int(reload)
        try:
            event_drive_factors = drive.evaluate(
                extension, float(geometry["tip_radius_m"])
            )
        except SourceDriveMapBoundsError:
            return _predictive_result(
                candidate, mechanics, temperature_K, target_extension_m,
                events, state, opening, time_s, extension, intervals,
                "RIGHT_CENSORED_DRIVE_MAP_BOUND",
            )
        events.append({
            "event_index": event_index, "physical_avalanche_index": avalanche,
            "event_time_s": time_s, "event_opening_m": opening,
            "reload_opening_m": reload_increment, "reload_separated": reload,
            "extension_before_m": extension, "extension_after_m": extension + committed,
            "proposed_event_length_m": proposed, "event_length_m": committed,
            **metrics,
            "tip_radius_m": float(geometry["tip_radius_m"]),
            "front_width_m": float(geometry["front_width_m"]),
            "source_multiplicity": float(geometry["multiplicity_per_system"]),
            "backstress_Pa": float(np.mean(geometry["sigma_back_by_system_Pa"])),
            "mobile_density_m2": float(np.sum(state.mobile_m2)),
            "retained_density_m2": float(np.sum(state.retained_m2)),
            "cumulative_source_activations": float(np.sum(state.cumulative_activations)),
            "drive_factors": list(event_drive_factors),
            "right_censored_at_target": extension + committed >= target_extension_m - 1.0e-15,
        })
        keep = _translation_retention_factor(committed, physics, reduced)
        state.mobile_m2 *= keep
        state.retained_m2 *= keep
        state.local_slip_count_by_system *= keep
        state.extension_m += committed
        extension += committed
    status = "TARGET_RIGHT_CENSORED" if extension >= target_extension_m - 1.0e-15 else "LOADING_MAP_EXHAUSTED"
    return _predictive_result(candidate, mechanics, temperature_K, target_extension_m,
                              events, state, opening, time_s, extension, intervals, status)


def _predictive_result(candidate, mechanics, temperature_K, target, events, state,
                       opening, time_s, extension, intervals, status):
    counts: dict[int, int] = {}
    for event in events:
        key = int(event["physical_avalanche_index"])
        counts[key] = counts.get(key, 0) + 1
    largest = max(counts.values(), default=0)
    return {
        "schema": "oneD_v2_provider_native_zeroD_predictive_v1",
        "backend": mechanics.backend, "candidate_id": candidate.candidate_id,
        "temperature_K": float(temperature_K), "target_extension_m": float(target),
        "status": status, "event_count": len(events),
        "physical_avalanche_count": len(counts),
        "precursor_reinitiation_count": sum(bool(x["reload_separated"]) for x in events),
        "largest_avalanche_fraction": largest / max(len(events), 1),
        "first_event_opening_m": None if not events else events[0]["event_opening_m"],
        "first_event_native_KJ_MPa_sqrt_m": None if not events else events[0]["native_KJ_MPa_sqrt_m"],
        "terminal_opening_m": opening, "terminal_time_s": time_s,
        "terminal_extension_m": extension, "interval_count": intervals,
        "max_tip_radius_m": state.max_tip_radius_m,
        "max_backstress_Pa": state.max_backstress_Pa,
        "min_front_width_m": state.min_front_width_m,
        "max_source_multiplicity": state.max_multiplicity,
        "events": events,
        "model_contract": "SOURCE_NORMALIZED_DRIVE+V913_ZEROD_IMPLICIT_BACKSTRESS_V1",
    }


__all__ = [
    "ProviderMechanicsMap", "SourceDriveMap", "SourceDriveMapBoundsError",
    "run_zero_d_predictive",
    "stochastic_loading_map",
]
