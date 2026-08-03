"""Adaptive regularity search for v10.0.5.18.3.9 atomic corridors.

The first full archived-seed qualification established that the corridor
architecture committed seven exact physical renewals, including an 18.38 um
event, before reaching a cavity for which all candidates were generated with
only 1, 2, or 5 degree Triangle minimum-angle requests.  The best final mesh
passed the production q floor but reached a local support-area ratio of 0.0657,
just below the unchanged 0.08 gate.

This module completes the bounded search over standard shape-regular Triangle
angle targets.  It does not relax either acceptance floor.  Earlier low-angle
candidates remain first, so already-feasible events retain their prior mesh.
Only a rejected cavity proceeds to 10, 15, 20, 25, and 28 degree requests.
"""
from __future__ import annotations

import copy
import math
from typing import Any

import numpy as np

from . import atomic_path_corridor_czm_v10051839 as _base
from .atomic_path_corridor_local_scale_v10051839 import (
    CertifiedAtomicPathCorridorCZMBackendV10051839,
)


MODEL_ID = "atomic_path_corridor_adaptive_regular_quality_search_v10_0_5_18_3_9"
ANGLE_SCHEDULE_DEG = (1.0, 2.0, 5.0, 10.0, 15.0, 20.0, 25.0, 28.0)
WIDTH_SCHEDULE = (0.75, 1.0, 1.5, 2.25, 3.25)
RING_SCHEDULE = (0, 1, 2)


class AdaptiveQualityAtomicPathCorridorCZMBackendV10051839(
    CertifiedAtomicPathCorridorCZMBackendV10051839
):
    """Certified corridor backend with a complete bounded regularity search."""

    name = "adaptive_czm_v10051839_atomic_path_corridor_adaptive_quality"

    def _build_corridor(
        self,
        mesh: Any,
        damage: np.ndarray,
        displacement: np.ndarray,
        p0: np.ndarray,
        p1: np.ndarray,
        direction: np.ndarray,
        front_id: int,
    ):
        tip_ids = self._tip_geometric_node_ids(mesh, p0, int(front_id))
        incident = (
            self._incident_elements(mesh.elems, set(map(int, tip_ids)))
            if tip_ids else np.zeros(0, int)
        )
        lengths = _base._unique_edge_lengths(mesh, incident)
        h_local = (
            float(np.median(lengths))
            if len(lengths) else float(max(mesh.hbar_tip, 1.0e-12))
        )
        length = float(np.linalg.norm(p1 - p0))
        base_n = max(
            1,
            int(math.ceil(length / max(0.75 * h_local, 1.0e-12))),
        )
        max_path = _base._int_env(
            "ARRHENIUS_CORRIDOR_MAX_PATH_SUBSEGMENTS", 16
        )
        n_values: list[int] = []
        for value in (
            1,
            base_n,
            base_n + 1,
            2 * base_n,
            3 * base_n,
            max_path,
        ):
            value = min(max(int(value), 1), max_path)
            if value not in n_values:
                n_values.append(value)

        attempts: list[dict[str, Any]] = []
        for nseg in n_values:
            for width in WIDTH_SCHEDULE:
                for rings in RING_SCHEDULE:
                    for angle in ANGLE_SCHEDULE_DEG:
                        state, record = _base._corridor_attempt(
                            self,
                            mesh,
                            damage,
                            displacement,
                            p0,
                            p1,
                            direction,
                            front_id,
                            width,
                            rings,
                            angle,
                            nseg,
                        )
                        record.update(
                            attempt_index=int(len(attempts)),
                            accepted=state is not None,
                            adaptive_quality_search_active=True,
                            quality_angle_schedule_deg=list(ANGLE_SCHEDULE_DEG),
                        )
                        attempts.append(copy.deepcopy(record))
                        if state is not None:
                            record["accepted_quality_angle_deg"] = float(angle)
                            return state, record, attempts

        return None, {
            "reason": "no_feasible_atomic_path_corridor_after_adaptive_quality_search",
            "attempt_count": int(len(attempts)),
            "quality_angle_schedule_deg": list(ANGLE_SCHEDULE_DEG),
            "path_subsegment_schedule": n_values,
            "width_schedule": list(WIDTH_SCHEDULE),
            "ring_schedule": list(RING_SCHEDULE),
        }, attempts


__all__ = [
    "AdaptiveQualityAtomicPathCorridorCZMBackendV10051839",
    "ANGLE_SCHEDULE_DEG",
    "MODEL_ID",
    "RING_SCHEDULE",
    "WIDTH_SCHEDULE",
]
