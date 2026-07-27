"""Numerical resilience helpers for FEM/CZM v10.0.5.18.3.

These helpers do not alter fracture or emission kinetics.  They make two
numerical operations retryable:

* a physical crack event may be represented by several collinear cohesive
  topology subsegments whose total length equals the original stochastic event;
* a failed finite-radius tensor probe is repeated with a modestly enlarged
  physical annulus and relaxed damage exclusion before the bounded last-reliable
  drive fallback is used.
"""
from __future__ import annotations

import copy
from typing import Callable

import numpy as np

from .crack_backend import AdaptiveCZMBackend, CrackAdvanceResult

MODEL_ID = "FEM_CZM_numerical_resilience_v10_0_5_18_3"
CZM_RETRY_SCHEMA = "v10.0.5.18.3_atomic_collinear_czm_partition_retry"
PROBE_RETRY_SCHEMA = "v10.0.5.18.3_expanded_physical_annulus_tensor_probe_retry"


class RetryingAdaptiveCZMBackendV1005183(AdaptiveCZMBackend):
    """Represent one physical event by atomic collinear topology subsegments."""

    name = "adaptive_czm_v1005183_retry"
    partition_counts = (2, 4, 8)

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.partition_retry_attempts = 0
        self.partition_retry_successes = 0
        self.partition_retry_failures = 0
        self.partition_retry_log: list[dict] = []

    def advance(self, **kwargs) -> CrackAdvanceResult:
        # Recursive exact-ray marching performed by AdaptiveCZMBackend must not
        # recursively start a new partition policy.
        if int(kwargs.get("_subdepth", 0) or 0) > 0 or bool(
            kwargs.get("_partition_retry_active", False)
        ):
            return super().advance(**kwargs)

        first = super().advance(**kwargs)
        if first.inserted:
            return first

        mesh0 = kwargs["mesh"]
        boundary0 = kwargs["boundary"]
        damage0 = kwargs["damage"]
        displacement0 = kwargs["displacement"]
        p0 = np.asarray(kwargs["p0"], dtype=float)
        p1 = np.asarray(kwargs["p1"], dtype=float)
        direction = np.asarray(kwargs["direction"], dtype=float)
        front_id = int(kwargs.get("front_id", 0))
        norm = float(np.linalg.norm(direction))
        requested = float(np.linalg.norm(p1 - p0))
        if norm <= 0.0 or requested <= 0.0:
            return first
        direction = direction / norm

        base_kwargs = {
            key: value
            for key, value in kwargs.items()
            if key not in {
                "mesh",
                "boundary",
                "damage",
                "displacement",
                "p0",
                "p1",
                "direction",
                "front_id",
                "_subdepth",
                "_partition_retry_active",
            }
        }

        for partitions in self.partition_counts:
            self.partition_retry_attempts += 1
            backend_snapshot = self._transaction_snapshot()
            mesh = mesh0
            boundary = boundary0
            damage = damage0
            displacement = displacement0
            point = p0.copy()
            parent_map = np.arange(mesh0.ne, dtype=int)
            moved_total = 0.0
            max_angle = 0.0
            failure_reason = "none"
            log_start = len(self.advance_log)
            succeeded = True

            for index in range(partitions):
                target = p0 + requested * float(index + 1) / float(partitions) * direction
                result = super().advance(
                    mesh=mesh,
                    boundary=boundary,
                    damage=damage,
                    displacement=displacement,
                    p0=point,
                    p1=target,
                    direction=direction,
                    front_id=front_id,
                    _partition_retry_active=True,
                    **base_kwargs,
                )
                if not result.inserted or result.moved <= 0.0:
                    failure_reason = str(result.reason)
                    succeeded = False
                    break
                if result.elem_parent_map is not None:
                    parent_map = parent_map[
                        np.asarray(result.elem_parent_map, dtype=int)
                    ]
                mesh = result.mesh
                boundary = result.boundary
                damage = result.damage
                displacement = result.displacement
                moved_total += float(result.moved)
                max_angle = max(max_angle, abs(float(result.angle_error_deg)))
                if front_id in self.tip_nodes:
                    point = np.asarray(self.tip_nodes[front_id][2], dtype=float).copy()
                elif self.advance_log:
                    row = self.advance_log[-1]
                    point = np.array([row["x1"], row["y1"]], dtype=float)
                else:
                    point = target.copy()

            length_ok = abs(moved_total - requested) <= max(
                1.0e-12, 1.0e-6 * requested
            )
            if succeeded and length_ok:
                self.partition_retry_successes += 1
                for row in self.advance_log[log_start:]:
                    row["physical_event_partition_retry"] = True
                    row["physical_event_partition_count"] = int(partitions)
                    row["physical_event_requested_length_m"] = requested
                    row["physical_event_retry_schema"] = CZM_RETRY_SCHEMA
                self.partition_retry_log.append(
                    {
                        "front_id": front_id,
                        "partitions": int(partitions),
                        "requested_length_m": requested,
                        "moved_length_m": moved_total,
                        "initial_failure_reason": str(first.reason),
                        "success": True,
                    }
                )
                return CrackAdvanceResult(
                    mesh,
                    boundary,
                    damage,
                    displacement,
                    moved_total,
                    True,
                    angle_error_deg=max_angle,
                    selected_edge_length=moved_total,
                    reason=f"partition_retry_{partitions}",
                    elem_parent_map=parent_map,
                )

            self._transaction_rollback(backend_snapshot)
            self.partition_retry_failures += 1
            self.partition_retry_log.append(
                {
                    "front_id": front_id,
                    "partitions": int(partitions),
                    "requested_length_m": requested,
                    "moved_length_m": moved_total,
                    "initial_failure_reason": str(first.reason),
                    "partition_failure_reason": failure_reason,
                    "success": False,
                }
            )

        return CrackAdvanceResult(
            mesh0,
            boundary0,
            damage0,
            displacement0,
            0.0,
            False,
            angle_error_deg=float(first.angle_error_deg),
            selected_edge_length=0.0,
            reason=(
                f"{first.reason};partition_retry_exhausted:"
                + ",".join(str(value) for value in self.partition_counts)
            ),
            elem_parent_map=None,
        )

    def write_diagnostics(self, out_dir: str) -> None:
        super().write_diagnostics(out_dir)
        import json
        from pathlib import Path

        out = Path(out_dir)
        payload = {
            "schema": CZM_RETRY_SCHEMA,
            "attempts": int(self.partition_retry_attempts),
            "successes": int(self.partition_retry_successes),
            "failures": int(self.partition_retry_failures),
            "records": copy.deepcopy(self.partition_retry_log),
        }
        (out / "adaptive_czm_partition_retry_v10_0_5_18_3.json").write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n"
        )


def make_robust_process_zone_traction_probe(
    original: Callable,
) -> Callable:
    """Retry the tensor probe without changing its stress normalization."""

    def robust(
        mesh,
        sigma_gp,
        d,
        tip,
        crack_direction,
        radius_m,
        annulus_half_width=0.45,
        sector_half_angle_deg=40.0,
        damage_cutoff=0.85,
        min_elements=4,
    ):
        attempts = (
            (1.0, float(damage_cutoff), int(min_elements)),
            (1.25, max(float(damage_cutoff), 0.95), max(int(min_elements), 4)),
            (1.75, max(float(damage_cutoff), 0.99), max(int(min_elements) - 1, 3)),
            (2.5, max(float(damage_cutoff), 0.999), max(int(min_elements) - 2, 2)),
        )
        best = None
        for level, (radius_factor, cutoff, required) in enumerate(attempts):
            result = original(
                mesh,
                sigma_gp,
                d,
                tip,
                crack_direction,
                float(radius_m) * radius_factor,
                annulus_half_width=annulus_half_width,
                sector_half_angle_deg=sector_half_angle_deg,
                damage_cutoff=cutoff,
                min_elements=required,
            )
            result = dict(result)
            result["probe_retry_schema"] = PROBE_RETRY_SCHEMA
            result["probe_retry_level"] = int(level)
            result["probe_retry_radius_factor"] = float(radius_factor)
            result["probe_retry_damage_cutoff"] = float(cutoff)
            result["probe_retry_min_elements"] = int(required)
            if bool(result.get("reliable", False)):
                result["probe_retry_active"] = bool(level > 0)
                return result
            if best is None or int(result.get("n_elements", 0)) > int(
                best.get("n_elements", 0)
            ):
                best = result
        best = {} if best is None else best
        best["probe_retry_active"] = True
        best["probe_retry_exhausted"] = True
        return best

    return robust


__all__ = [
    "CZM_RETRY_SCHEMA",
    "MODEL_ID",
    "PROBE_RETRY_SCHEMA",
    "RetryingAdaptiveCZMBackendV1005183",
    "make_robust_process_zone_traction_probe",
]
