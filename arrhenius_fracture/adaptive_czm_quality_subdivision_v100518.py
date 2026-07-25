"""Quality-safe atomic subdivision for arbitrary-length adaptive-CZM events.

The v9.18.5.6 post-transaction gate correctly rejects poor local triangles, but
its original implementation applies the child/parent area-ratio floor again to
a complete multi-subsegment event. A composed parent map then compares a
later-generation child with the event-origin element, even though each
immediate topology operation already passed the same floor. It also counts one
veto once per recursive stack frame.

This additive repair is installed only by the standalone four-class FEM path.
It keeps every production quality floor unchanged. It evaluates composed
multi-subsegment transactions using immediate-parent ratios from per-step logs,
counts a rejected physical event once, and retries a rejected arbitrary-length
event as an atomic sequence of collinear quality-checked subsegments. No
constitutive, hazard, loading, or material law is changed.
"""
from __future__ import annotations

from contextlib import contextmanager
import math
import os
from typing import Any, Iterator

import numpy as np

from . import crack_backend as _cb
from . import mode_i_first_passage_v9_18_5 as _v9185
from . import mode_i_first_passage_v9_18_5_4 as _v91854
from . import mode_i_first_passage_v9_18_5_6 as _v91856

MODEL_ID = "adaptive_CZM_atomic_quality_subdivision_v10_0_5_18"


def _float_env(name: str, default: float) -> float:
    try:
        value = float(os.environ.get(name, default))
    except (TypeError, ValueError):
        value = float(default)
    return value if math.isfinite(value) else float(default)


def _subdivision_counts() -> tuple[int, ...]:
    raw = os.environ.get("ARRHENIUS_QUALITY_SUBDIVISION_COUNTS", "2,3,4,6,8")
    values: list[int] = []
    for token in raw.replace(";", ",").split(","):
        token = token.strip()
        if not token:
            continue
        try:
            value = int(token)
        except ValueError:
            continue
        if value >= 2 and value not in values:
            values.append(value)
    return tuple(values or (2, 3, 4, 6, 8))


def _audit_lengths() -> dict[str, int]:
    return {
        "accepted": len(_v91856._AUDIT.get("accepted_events", [])),
        "warnings": len(_v91856._AUDIT.get("resolution_warnings", [])),
        "vetoes": len(_v91856._AUDIT.get("quality_vetoes", [])),
        "runtime_vetoes": len(_v9185._RUNTIME.get("quality_vetoes", [])),
    }


def _truncate_audits(lengths: dict[str, int]) -> None:
    del _v91856._AUDIT.setdefault("accepted_events", [])[lengths["accepted"] :]
    del _v91856._AUDIT.setdefault("resolution_warnings", [])[lengths["warnings"] :]
    del _v91856._AUDIT.setdefault("quality_vetoes", [])[lengths["vetoes"] :]
    del _v9185._RUNTIME.setdefault("quality_vetoes", [])[lengths["runtime_vetoes"] :]
    _v91856._AUDIT["consecutive_veto_abort"] = None


def _counter_snapshot(self: Any) -> tuple[int, Any]:
    return (
        int(getattr(self, "_v91856_consecutive_geometry_vetoes", 0)),
        getattr(self, "_v91856_last_veto_reason", None),
    )


def _restore_counter(self: Any, snapshot: tuple[int, Any]) -> None:
    self._v91856_consecutive_geometry_vetoes = int(snapshot[0])
    self._v91856_last_veto_reason = snapshot[1]


def _local_area_ratios_from_logs(self: Any, start: int) -> list[float]:
    logs = list(getattr(self, "advance_log", []) or [])[int(start) :]
    ratios: list[float] = []
    for row in logs:
        if not isinstance(row, dict):
            continue
        for key in (
            "v100518_min_child_area_ratio",
            "v91856_min_child_area_ratio",
            "v9185_min_child_area_ratio",
            "min_area_ratio",
        ):
            value = row.get(key)
            if value is None:
                continue
            try:
                number = float(value)
            except (TypeError, ValueError):
                continue
            if math.isfinite(number) and number > 0.0:
                ratios.append(number)
                break
    return ratios


def _assess(
    self: Any,
    old_mesh: Any,
    result: Any,
    kwargs: dict[str, Any],
    log_start: int,
):
    new_mesh = result.mesh
    affected = _v9185._affected_elements(old_mesh, new_mesh)
    quality = self._triangle_quality(new_mesh.nodes, new_mesh.elems[affected])
    qmin = float(np.min(quality)) if quality.size else 1.0

    if not hasattr(self, "_v91856_production_q_floor"):
        self._v91856_production_q_floor = float(self.min_triangle_quality)
        self._v91856_production_area_floor = float(self.min_area_ratio)
    qfloor = _float_env(
        "ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY",
        float(self._v91856_production_q_floor),
    )
    afloor = _float_env(
        "ARRHENIUS_MIN_ACCEPTED_CHILD_AREA_RATIO",
        float(self._v91856_production_area_floor),
    )

    new_logs = list(getattr(self, "advance_log", []) or [])[int(log_start) :]
    local_ratios = _local_area_ratios_from_logs(self, log_start)
    multigeneration = len(new_logs) > 1
    if multigeneration and local_ratios:
        amin = float(min(local_ratios))
        ratio_reference = "minimum_immediate_topology_subtransaction"
    else:
        parent_map = getattr(result, "elem_parent_map", None)
        ratios: list[float] = []
        if parent_map is not None:
            pm = np.asarray(parent_map, dtype=int)
            for e in affected:
                if e < len(pm) and 0 <= pm[e] < int(old_mesh.ne):
                    ratios.append(
                        float(new_mesh.area_e[e])
                        / max(float(old_mesh.area_e[pm[e]]), 1.0e-300)
                    )
        else:
            for e in affected:
                if e < int(old_mesh.ne):
                    ratios.append(
                        float(new_mesh.area_e[e])
                        / max(float(old_mesh.area_e[e]), 1.0e-300)
                    )
        amin = float(min(ratios)) if ratios else 1.0
        ratio_reference = "direct_parent_map_for_single_topology_transaction"

    p0 = np.asarray(kwargs.get("p0"), float)
    p1 = np.asarray(kwargs.get("p1"), float)
    da = max(float(np.linalg.norm(p1 - p0)), 1.0e-300)
    endpoint = _v91854._active_endpoint(self, result, kwargs)
    local = _v91854._active_tip_one_ring_resolution(new_mesh, endpoint)
    tip_ratio = float(local["active_tip_h_mean_m"]) / da
    requested_tip_ratio = _float_env("ARRHENIUS_MAX_TIP_H_OVER_DA", 0.75)

    incidence = np.bincount(
        np.asarray(new_mesh.elems, int).ravel(), minlength=int(new_mesh.nn)
    )
    orphan = np.where(incidence <= 0)[0]
    bad_endpoint: list[int] = []
    for elem in self.cohesive_network.elements:
        for nid in elem.nodes4:
            if (
                int(nid) < 0
                or int(nid) >= int(new_mesh.nn)
                or incidence[int(nid)] <= 0
            ):
                bad_endpoint.append(int(nid))

    issues: list[str] = []
    if not np.all(np.isfinite(new_mesh.area_e)) or np.any(new_mesh.area_e <= 0.0):
        issues.append("nonpositive_or_nonfinite_area")
    if qmin < qfloor:
        issues.append(f"triangle_quality={qmin:.6e}<{qfloor:.6e}")
    if amin < afloor:
        issues.append(f"child_area_ratio={amin:.6e}<{afloor:.6e}")
    if orphan.size:
        issues.append(f"orphan_bulk_nodes={orphan[:10].tolist()}")
    if bad_endpoint:
        issues.append(f"unsupported_cohesive_endpoints={bad_endpoint[:10]}")

    warning = bool(math.isfinite(tip_ratio) and tip_ratio > requested_tip_ratio)
    row = {
        "schema": MODEL_ID,
        "front_id": int(kwargs.get("front_id", -1)),
        "min_triangle_quality": qmin,
        "triangle_quality_floor": qfloor,
        "min_child_area_ratio": amin,
        "child_area_ratio_floor": afloor,
        "child_area_ratio_reference": ratio_reference,
        "multi_generation_atomic_event": multigeneration,
        "topology_subtransaction_count": int(len(new_logs)),
        "active_tip_h_over_da": tip_ratio,
        "requested_tip_h_over_da": requested_tip_ratio,
        "tip_h_over_da_enforced_as_veto": False,
        "resolution_warning": warning,
        "requested_da_m": da,
        "affected_element_count": int(len(affected)),
        "accepted": not bool(issues),
        "issues": issues,
        "resolution_metric": "active_endpoint_one_ring_unique_nonzero_edge_mean",
        **local,
    }
    return issues, row, warning


def _accept(self: Any, result: Any, row: dict[str, Any], warning: bool) -> Any:
    result.mesh.hbar_tip = float(row["active_tip_h_mean_m"])
    if getattr(self, "advance_log", None):
        self.advance_log[-1].update(
            {
                "v100518_quality_gate_passed": True,
                "v100518_min_triangle_quality": row["min_triangle_quality"],
                "v100518_min_child_area_ratio": row["min_child_area_ratio"],
                "v100518_child_area_ratio_reference": row[
                    "child_area_ratio_reference"
                ],
                "v100518_active_tip_h_over_da": row["active_tip_h_over_da"],
                "v100518_resolution_warning": warning,
                "v100518_active_tip_h_mean_m": row["active_tip_h_mean_m"],
            }
        )
    _v91856._AUDIT.setdefault("accepted_events", []).append(row)
    if warning:
        _v91856._AUDIT.setdefault("resolution_warnings", []).append(row.copy())
    self._v91856_consecutive_geometry_vetoes = 0
    self._v91856_last_veto_reason = None
    return result


def _veto_result(
    self: Any,
    old_mesh: Any,
    kwargs: dict[str, Any],
    result: Any,
    row: dict[str, Any],
):
    _v91856._AUDIT.setdefault("quality_vetoes", []).append(row)
    _v9185._RUNTIME.setdefault("quality_vetoes", []).append(row.copy())
    return _cb.CrackAdvanceResult(
        mesh=old_mesh,
        boundary=kwargs["boundary"],
        damage=kwargs["damage"],
        displacement=kwargs["displacement"],
        moved=0.0,
        inserted=False,
        angle_error_deg=float(getattr(result, "angle_error_deg", 0.0)),
        reason="v100518_quality_veto:" + ";".join(row.get("issues", [])),
        elem_parent_map=None,
    )


def _compose_parent_map(
    previous: np.ndarray,
    result: Any,
    previous_ne: int,
) -> np.ndarray:
    parent_map = getattr(result, "elem_parent_map", None)
    if parent_map is None:
        if int(result.mesh.ne) != int(previous_ne):
            raise RuntimeError(
                "quality subdivision cannot compose a missing parent map across a "
                "topology-changing segment"
            )
        return previous.copy()
    pm = np.asarray(parent_map, dtype=int)
    if np.any(pm < 0) or np.any(pm >= len(previous)):
        raise RuntimeError("quality subdivision received an invalid element parent map")
    return previous[pm]


def _retryable(reason: str) -> bool:
    text = str(reason).lower()
    return any(
        token in text
        for token in (
            "quality_veto",
            "triangle_quality",
            "child_area_ratio",
            "no_quality_safe_exact_ray_move",
            "local_hrefine",
            "exit_at_target",
            "no_ray_exit",
        )
    )


def _try_subdivisions(
    self: Any,
    original: Any,
    args: tuple[Any, ...],
    kwargs: dict[str, Any],
    outer_snapshot: Any,
    audit_snapshot: dict[str, int],
    counter_snapshot: tuple[int, Any],
):
    if args:
        return None
    p0 = np.asarray(kwargs.get("p0"), float).reshape(2)
    p1 = np.asarray(kwargs.get("p1"), float).reshape(2)
    direction = np.asarray(kwargs.get("direction"), float).reshape(2)
    length = float(np.linalg.norm(p1 - p0))
    norm = float(np.linalg.norm(direction))
    if not math.isfinite(length) or length <= 0.0 or norm <= 1.0e-30:
        return None
    direction = direction / norm

    old_mesh = kwargs["mesh"]
    for count in _subdivision_counts():
        self._transaction_rollback(outer_snapshot)
        _truncate_audits(audit_snapshot)
        _restore_counter(self, counter_snapshot)

        mesh = kwargs["mesh"]
        boundary = kwargs["boundary"]
        damage = kwargs["damage"]
        displacement = kwargs["displacement"]
        current_p0 = p0.copy()
        parent_to_origin = np.arange(int(old_mesh.ne), dtype=int)
        moved_total = 0.0
        success = True
        last_result = None

        for index in range(1, count + 1):
            target = p0 + direction * (length * float(index) / float(count))
            local = dict(kwargs)
            local.update(
                {
                    "mesh": mesh,
                    "boundary": boundary,
                    "damage": damage,
                    "displacement": displacement,
                    "p0": current_p0,
                    "p1": target,
                    "direction": direction,
                    "_subdepth": 1,
                }
            )
            log_start = len(getattr(self, "advance_log", []) or [])
            previous_ne = int(mesh.ne)
            result = original(self, **local)
            if not bool(getattr(result, "inserted", False)):
                success = False
                break
            issues, row, warning = _assess(self, mesh, result, local, log_start)
            row.update(
                {
                    "quality_subdivision_active": True,
                    "quality_subdivision_count": int(count),
                    "quality_subdivision_index": int(index),
                    "physical_event_total_length_m": length,
                }
            )
            if issues:
                success = False
                break
            _accept(self, result, row, warning)
            parent_to_origin = _compose_parent_map(
                parent_to_origin,
                result,
                previous_ne,
            )
            moved_total += float(result.moved)
            mesh = result.mesh
            boundary = result.boundary
            damage = result.damage
            displacement = result.displacement
            current_p0 = target
            last_result = result

        if success and last_result is not None and math.isclose(
            moved_total,
            length,
            rel_tol=1.0e-8,
            abs_tol=5.0e-10,
        ):
            if getattr(self, "advance_log", None):
                self.advance_log[-1].update(
                    {
                        "v100518_atomic_quality_subdivision_recovery": True,
                        "v100518_subdivision_count": int(count),
                        "v100518_physical_event_total_length_m": length,
                    }
                )
            self._v91856_consecutive_geometry_vetoes = 0
            self._v91856_last_veto_reason = None
            return _cb.CrackAdvanceResult(
                mesh=mesh,
                boundary=boundary,
                damage=damage,
                displacement=displacement,
                moved=moved_total,
                inserted=True,
                angle_error_deg=0.0,
                selected_edge_length=moved_total,
                reason=f"ok_atomic_quality_subdivision_{count}",
                elem_parent_map=parent_to_origin,
            )

    self._transaction_rollback(outer_snapshot)
    _truncate_audits(audit_snapshot)
    _restore_counter(self, counter_snapshot)
    return None


def _quality_subdividing_advance_v100518(self: Any, *args, **kwargs):
    original = _quality_subdividing_advance_v100518._original
    depth = int(kwargs.get("_subdepth", 0) or 0)
    old_mesh = kwargs["mesh"]
    outer_snapshot = self._transaction_snapshot()
    audit_snapshot = _audit_lengths()
    counter = _counter_snapshot(self)
    log_start = len(getattr(self, "advance_log", []) or [])

    result = original(self, *args, **kwargs)
    initial_row: dict[str, Any] | None = None
    if bool(getattr(result, "inserted", False)):
        issues, row, warning = _assess(self, old_mesh, result, kwargs, log_start)
        if not issues:
            return _accept(self, result, row, warning)
        initial_row = row
        self._transaction_rollback(outer_snapshot)
        _truncate_audits(audit_snapshot)
        _restore_counter(self, counter)
        result = _veto_result(self, old_mesh, kwargs, result, row)

    if depth > 0:
        return result

    reason = str(getattr(result, "reason", ""))
    if initial_row is not None or _retryable(reason):
        self._transaction_rollback(outer_snapshot)
        _truncate_audits(audit_snapshot)
        _restore_counter(self, counter)
        recovered = _try_subdivisions(
            self,
            original,
            args,
            kwargs,
            outer_snapshot,
            audit_snapshot,
            counter,
        )
        if recovered is not None:
            if initial_row is not None:
                recovery = dict(initial_row)
                recovery.update(
                    {
                        "accepted": False,
                        "recovered_by_atomic_quality_subdivision": True,
                        "constitutive_event_consumed_once": True,
                    }
                )
                _v91856._AUDIT.setdefault("quality_vetoes", []).append(recovery)
                _v9185._RUNTIME.setdefault("quality_vetoes", []).append(
                    recovery.copy()
                )
            return recovered

    if initial_row is not None:
        self._transaction_rollback(outer_snapshot)
        _truncate_audits(audit_snapshot)
        _restore_counter(self, counter)
        result = _veto_result(self, old_mesh, kwargs, result, initial_row)
    return _v91856._record_or_raise(self, kwargs, result)


@contextmanager
def installed_quality_subdivision_v100518() -> Iterator[None]:
    original = _v91856._strict_quality_advance_v91856
    _v91856._strict_quality_advance_v91856 = _quality_subdividing_advance_v100518
    try:
        yield
    finally:
        _v91856._strict_quality_advance_v91856 = original


__all__ = [
    "MODEL_ID",
    "installed_quality_subdivision_v100518",
]
