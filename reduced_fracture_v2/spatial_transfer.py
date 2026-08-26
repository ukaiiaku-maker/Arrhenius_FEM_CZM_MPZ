"""Neutral spatial-state observers for the Taylor/Peierls transfer audit.

The observer records the existing v9.13 state and an observational ledger of
line content that leaves the moving ahead-tip window.  Neither the snapshots
nor the ledger are read by the constitutive equations.
"""
from __future__ import annotations

import math
from typing import Any, Mapping

import numpy as np

from arrhenius_fracture.emergent_gnd_state_v913 import EmergentGNDState


CHECKPOINTS_UM = (0.0, 25.0, 50.0, 100.0, 200.0, 300.0)


class ObservedEmergentGNDState(EmergentGNDState):
    """Production spatial state plus a strictly read-only wake-departure ledger."""

    def __init__(self, candidate, physics):
        super().__init__(candidate, physics)
        self.observer_mobile_departed_line_content = 0.0
        self.observer_retained_departed_line_content = 0.0
        self.observer_retained_departed_first_moment = 0.0
        self.observer_retained_departed_second_moment = 0.0

    def translate_tip(self, da_m: float) -> None:
        da = max(float(da_m), 0.0)
        if da <= 0.0:
            return
        old_extension = float(self.extension_m)
        mobile_before = float(np.sum(self.mobile_m2) * self.cell_area_m2)
        retained_before = float(np.sum(self.retained_m2) * self.cell_area_m2)
        super().translate_tip(da)
        mobile_after = float(np.sum(self.mobile_m2) * self.cell_area_m2)
        retained_after = float(np.sum(self.retained_m2) * self.cell_area_m2)
        mobile_departed = max(mobile_before - mobile_after, 0.0)
        retained_departed = max(retained_before - retained_after, 0.0)
        lab_position = old_extension + 0.5 * da
        self.observer_mobile_departed_line_content += mobile_departed
        self.observer_retained_departed_line_content += retained_departed
        self.observer_retained_departed_first_moment += retained_departed * lab_position
        self.observer_retained_departed_second_moment += (
            retained_departed * lab_position * lab_position
        )


def _moments(values: np.ndarray, positions: np.ndarray) -> tuple[float, float]:
    weights = np.maximum(np.asarray(values, dtype=float), 0.0)
    total = float(np.sum(weights))
    if total <= 0.0:
        return 0.0, 0.0
    centroid = float(np.sum(weights * positions) / total)
    width = math.sqrt(max(float(np.sum(weights * (positions - centroid) ** 2) / total), 0.0))
    return centroid, width


def state_snapshot(
    phase: str,
    context: Mapping[str, Any],
    state: EmergentGNDState,
    *,
    temperature_K: float,
) -> dict[str, Any]:
    """Aggregate one existing spatial state without mutating it."""
    mobile_by_bin = np.sum(np.maximum(state.mobile_m2, 0.0), axis=(0, 1))
    retained_by_bin = np.sum(np.maximum(state.retained_m2, 0.0), axis=(0, 1))
    mobile_line = mobile_by_bin * state.cell_area_m2
    retained_line = retained_by_bin * state.cell_area_m2
    mobile_total = float(np.sum(mobile_line))
    retained_total = float(np.sum(retained_line))
    total = mobile_total + retained_total
    mobile_centroid, mobile_width = _moments(mobile_line, state.x)
    retained_centroid, retained_width = _moments(retained_line, state.x)
    near = np.asarray(state.x) <= max(float(state.c.source_zone_length_m), float(state.dx))
    geometry = state.source_geometry()
    backstress = state.backstress_state()[2]
    K = float(context.get("K_MPa_sqrt_m", 0.0))
    rates = state.local_rates(K, float(temperature_K))
    peierls_velocity = np.abs(np.asarray(rates["peierls_velocity_m_s"], dtype=float))
    encounter = np.asarray(rates["encounter_s"], dtype=float)
    taylor = np.asarray(rates["taylor_completion_s"], dtype=float)
    wake_retained = float(getattr(state, "observer_retained_departed_line_content", 0.0))
    wake_first = float(getattr(state, "observer_retained_departed_first_moment", 0.0))
    wake_second = float(getattr(state, "observer_retained_departed_second_moment", 0.0))
    wake_centroid = wake_first / wake_retained if wake_retained > 0.0 else 0.0
    wake_width = (
        math.sqrt(max(wake_second / wake_retained - wake_centroid**2, 0.0))
        if wake_retained > 0.0 else 0.0
    )
    return {
        "phase": str(phase),
        **{key: value for key, value in context.items()},
        "temperature_K": float(temperature_K),
        "state_extension_m": float(state.extension_m),
        "mobile_line_content": mobile_total,
        "retained_line_content": retained_total,
        "retained_fraction": retained_total / max(total, 1.0e-300),
        "mobile_centroid_tip_relative_m": mobile_centroid,
        "mobile_width_m": mobile_width,
        "retained_centroid_tip_relative_m": retained_centroid,
        "retained_width_m": retained_width,
        "mobile_centroid_laboratory_m": float(state.extension_m) + mobile_centroid,
        "retained_centroid_laboratory_m": float(state.extension_m) + retained_centroid,
        "near_tip_mobile_line_content": float(np.sum(mobile_line[near])),
        "near_tip_retained_line_content": float(np.sum(retained_line[near])),
        "wake_mobile_departed_line_content": float(
            getattr(state, "observer_mobile_departed_line_content", 0.0)
        ),
        "wake_retained_departed_line_content": wake_retained,
        "wake_retained_centroid_laboratory_m": wake_centroid,
        "wake_retained_width_m": wake_width,
        "tip_radius_m": float(geometry["tip_radius_m"]),
        "front_width_m": float(geometry["front_width_m"]),
        "backstress_Pa": float(np.mean(backstress)),
        "K_shield_MPa_sqrt_m": float(state.K_shield_MPa_sqrt_m()),
        "source_multiplicity": float(geometry["multiplicity_per_system"]),
        "peierls_transport_velocity_abs_max_m_s": float(np.max(peierls_velocity)),
        "peierls_transport_velocity_abs_mean_m_s": float(np.mean(peierls_velocity)),
        "encounter_retention_rate_max_s": float(np.max(encounter)),
        "encounter_retention_rate_mean_s": float(np.mean(encounter)),
        "taylor_completion_rate_max_s": float(np.max(taylor)),
        "taylor_completion_rate_mean_s": float(np.mean(taylor)),
        "observer_feedback": False,
        "observer_wake_ledger_semantics": (
            "CUMULATIVE_LINE_CONTENT_DEPARTING_AHEAD_TIP_WINDOW_NO_FEEDBACK"
        ),
    }


def physical_onset_event_indices(
    events: list[Any], *, reload_gap_threshold_m: float
) -> tuple[list[int], list[int]]:
    """Apply the PF V2 unresolved-pair latch and return onset/avalanche indices."""
    if not events:
        return [], []
    onsets = [0]
    avalanches = [0]
    latched = False
    avalanche = 0
    for index in range(1, len(events)):
        gap = float(events[index].applied_displacement_m) - float(
            events[index - 1].applied_displacement_m
        )
        reload = (not latched) and gap >= float(reload_gap_threshold_m)
        if reload:
            avalanche += 1
            onsets.append(index)
        else:
            latched = True
        avalanches.append(avalanche)
    return onsets, avalanches


def sampled_history(
    snapshots: list[dict[str, Any]], onset_indices: list[int]
) -> list[dict[str, Any]]:
    """Select initial, finite-event checkpoint crossings, and physical onsets."""
    output: list[dict[str, Any]] = []
    initial = next(row for row in snapshots if row["phase"] == "INITIAL")
    output.append({**initial, "sample_role": "EXTENSION_CHECKPOINT", "requested_extension_um": 0.0})
    posts = [row for row in snapshots if row["phase"] == "POST_EVENT"]
    for checkpoint in CHECKPOINTS_UM[1:]:
        candidates = [
            row for row in posts
            if float(row["projected_extension_m"]) * 1.0e6 >= checkpoint - 1.0e-12
        ]
        if not candidates:
            raise RuntimeError(f"missing spatial checkpoint at {checkpoint:g} um")
        output.append({
            **candidates[0],
            "sample_role": "EXTENSION_CHECKPOINT",
            "requested_extension_um": checkpoint,
        })
    pre = {int(row["event_index"]): row for row in snapshots if row["phase"] == "PRE_EVENT"}
    for ordinal, event_index in enumerate(onset_indices):
        output.append({
            **pre[event_index],
            "sample_role": "PHYSICAL_ONSET",
            "physical_onset_ordinal": ordinal,
            "requested_extension_um": float(pre[event_index]["projected_extension_m"]) * 1.0e6,
        })
    return output


PROFILE_FIELDS = (
    "mobile_m2", "retained_m2", "accumulated_slip_m2",
)


def profile_rows(
    candidate_id: str,
    sample: Mapping[str, Any],
    state_payload: Mapping[str, Any],
) -> list[dict[str, Any]]:
    """Expand a saved active-window payload into tidy system/sign/bin rows."""
    rows: list[dict[str, Any]] = []
    x = np.asarray(state_payload["x_m"], dtype=float)
    mobile = np.asarray(state_payload["mobile_m2"], dtype=float)
    retained = np.asarray(state_payload["retained_m2"], dtype=float)
    slip = np.asarray(state_payload["accumulated_slip_m2"], dtype=float)
    for system in range(mobile.shape[0]):
        for sign_index, sign in enumerate((-1, 1)):
            for bin_index, distance in enumerate(x):
                rows.append({
                    "candidate_id": candidate_id,
                    "sample_role": sample["sample_role"],
                    "requested_extension_um": sample["requested_extension_um"],
                    "actual_projected_extension_um": float(sample["projected_extension_m"]) * 1.0e6,
                    "physical_onset_ordinal": sample.get("physical_onset_ordinal"),
                    "system": system,
                    "sign": sign,
                    "bin_index": bin_index,
                    "tip_relative_x_m": float(distance),
                    "laboratory_x_m": float(sample["state_extension_m"]) + float(distance),
                    "mobile_density_m2": float(mobile[system, sign_index, bin_index]),
                    "retained_density_m2": float(retained[system, sign_index, bin_index]),
                    "accumulated_slip_density_m2": float(slip[system, sign_index, bin_index]),
                    "observer_feedback": False,
                })
    return rows


__all__ = [
    "CHECKPOINTS_UM", "ObservedEmergentGNDState", "PROFILE_FIELDS",
    "physical_onset_event_indices", "profile_rows", "sampled_history",
    "state_snapshot",
]
