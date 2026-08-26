"""Reload-separated resistance descriptors for predictive fracture histories.

Only the first pre-event state and later pre-event states separated by a
qualified reload are resistance candidates.  Event-wise drive values inside a
physical avalanche are retained as model-native trajectory data, never as
R-curve points.
"""
from __future__ import annotations

from collections import defaultdict
import json
from typing import Any, Mapping, Sequence

import numpy as np


def reload_separated_onsets(result: Mapping[str, Any]) -> list[dict[str, Any]]:
    """Return pre-event resistance candidates, one per physical avalanche."""
    events = list(result.get("events", ()))
    selected = [
        event for event in events
        if int(event.get("event_index", -1)) == 0
        or bool(event.get("reload_separated", False))
    ]
    output: list[dict[str, Any]] = []
    for ordinal, event in enumerate(selected):
        onset = {
            "onset_ordinal": ordinal,
            "event_index": int(event["event_index"]),
            "physical_avalanche_index": int(event["physical_avalanche_index"]),
            "event_time_s": float(event["event_time_s"]),
            "event_opening_m": float(event["event_opening_m"]),
            "extension_before_m": float(event["extension_before_m"]),
            "native_KJ_MPa_sqrt_m": float(event["native_KJ_MPa_sqrt_m"]),
            "qualified_G_J_per_m2": (
                None if event.get("qualified_G_J_per_m2") is None
                else float(event["qualified_G_J_per_m2"])
            ),
            "reload_opening_m": (
                None if event.get("reload_opening_m") is None
                else float(event["reload_opening_m"])
            ),
            "tip_radius_m": float(event["tip_radius_m"]),
            "front_width_m": float(event["front_width_m"]),
            "mobile_density_m2": float(event["mobile_density_m2"]),
            "retained_density_m2": float(event["retained_density_m2"]),
            "source_multiplicity": float(event["source_multiplicity"]),
            "backstress_Pa": float(event["backstress_Pa"]),
        }
        # These are source-owned diagnostics added by the Taylor/Peierls audit.
        # Keep the postprocessor backward compatible with archived events that
        # predate them, while preserving the complete decomposition whenever
        # the producer emitted it.
        optional_scalars = (
            "K_app_or_common_reference_MPa_sqrt_m",
            "K_native_MPa_sqrt_m",
            "K_shield_MPa_sqrt_m",
            "K_effective_local_equivalent_MPa_sqrt_m",
            "source_opening_stress_Pa",
            "transport_distance_m",
            "transport_time_s_min",
            "retention_encounter_time_s_min",
            "taylor_completion_time_s_min",
            "chi_ret_max",
            "chi_taylor_completion_max",
            "retained_equilibrium_fraction_max",
        )
        optional_vectors = (
            "resolved_emission_drive_Pa_by_system",
            "peierls_rate_s_by_system",
            "peierls_velocity_m_s_by_system",
            "taylor_completion_rate_s_by_system",
            "encounter_rate_s_by_system",
        )
        optional_labels = ("K_shield_status", "timescale_source")
        for field in optional_scalars:
            if field in event:
                onset[field] = float(event[field])
        for field in optional_vectors:
            if field in event:
                onset[field] = [float(value) for value in event[field]]
        for field in optional_labels:
            if field in event:
                onset[field] = str(event[field])
        output.append(onset)
    return output


def _sequence_monotonicity(values: Sequence[float]) -> str:
    if len(values) < 2:
        return "NOT_APPLICABLE"
    increments = np.diff(np.asarray(values, dtype=float))
    tolerance = 1.0e-12 * max(float(np.max(np.abs(values))), 1.0)
    if bool(np.all(increments >= -tolerance)):
        return "NONDECREASING"
    if bool(np.all(increments <= tolerance)):
        return "NONINCREASING"
    return "MIXED"


def summarize_rcurve_propensity(result: Mapping[str, Any]) -> dict[str, Any]:
    """Summarize reload-separated resistance and avalanche topology."""
    events = list(result.get("events", ()))
    onsets = reload_separated_onsets(result)
    native = [float(item["native_KJ_MPa_sqrt_m"]) for item in onsets]
    qualified = [
        float(item["qualified_G_J_per_m2"])
        for item in onsets if item["qualified_G_J_per_m2"] is not None
    ]
    by_avalanche: dict[int, float] = defaultdict(float)
    for event in events:
        by_avalanche[int(event["physical_avalanche_index"])] += (
            float(event["event_length_m"]) * 1.0e6
        )
    avalanche_sizes = [by_avalanche[key] for key in sorted(by_avalanche)]
    event_sizes = [float(event["event_length_m"]) * 1.0e6 for event in events]
    native_in_avalanche = [float(event["native_KJ_MPa_sqrt_m"]) for event in events]
    reload_increments = [
        float(item["reload_opening_m"]) * 1.0e6
        for item in onsets[1:] if item["reload_opening_m"] is not None
    ]
    initial_native = native[0] if native else np.nan
    initial_qualified = qualified[0] if qualified else np.nan
    return {
        "resistance_candidate_policy": "PRE_EVENT_RELOAD_SEPARATED_ONLY",
        "target_termination_semantics": "RIGHT_CENSORED_NOT_PHYSICAL_ARREST",
        "initial_onset_native_KJ_MPa_sqrt_m": initial_native,
        "maximum_onset_native_KJ_MPa_sqrt_m": max(native, default=np.nan),
        "deltaK_reinit_MPa_sqrt_m": (
            max(native) - initial_native if native else np.nan
        ),
        "relative_deltaK_reinit": (
            (max(native) - initial_native) / max(abs(initial_native), 1.0e-12)
            if native else np.nan
        ),
        "initial_onset_qualified_G_J_per_m2": initial_qualified,
        "maximum_onset_qualified_G_J_per_m2": max(qualified, default=np.nan),
        "deltaG_reinit_J_per_m2": (
            max(qualified) - initial_qualified if qualified else np.nan
        ),
        "relative_deltaG_reinit": (
            (max(qualified) - initial_qualified) / max(abs(initial_qualified), 1.0e-12)
            if qualified else np.nan
        ),
        "N_reinit": max(len(onsets) - 1, 0),
        "physical_avalanche_count": int(result.get("physical_avalanche_count", 0)),
        "largest_avalanche_fraction": float(result.get("largest_avalanche_fraction", 0.0)),
        "onset_sequence_monotonicity": _sequence_monotonicity(native),
        "reload_increment_min_um": min(reload_increments, default=np.nan),
        "reload_increment_max_um": max(reload_increments, default=np.nan),
        "mean_event_size_um": float(np.mean(event_sizes)) if event_sizes else np.nan,
        "median_event_size_um": float(np.median(event_sizes)) if event_sizes else np.nan,
        "maximum_event_size_um": max(event_sizes, default=np.nan),
        "final_avalanche_extension_um": avalanche_sizes[-1] if avalanche_sizes else 0.0,
        "onset_candidates_json": json.dumps(onsets, sort_keys=True, separators=(",", ":")),
        "avalanche_extensions_um_json": json.dumps(avalanche_sizes, separators=(",", ":")),
        "event_sizes_um_json": json.dumps(event_sizes, separators=(",", ":")),
        "in_avalanche_native_drive_json": json.dumps(
            native_in_avalanche, separators=(",", ":")
        ),
        "in_avalanche_drive_interpretation": "MODEL_NATIVE_TRAJECTORY_NOT_RESISTANCE",
    }


__all__ = ["reload_separated_onsets", "summarize_rcurve_propensity"]
