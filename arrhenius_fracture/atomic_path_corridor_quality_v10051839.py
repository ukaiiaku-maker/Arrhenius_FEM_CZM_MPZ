"""Strict production quality wrapper for v10.0.5.18.3.9 atomic remeshing.

The triangle-shape, finite-area, orphan-node, cohesive-support, and tip-resolution
checks are inherited unchanged.  For a certified whole-cavity remesh, the 0.08
area threshold is evaluated using the remesher's local path-support cell rather
than pretending the final mesh was produced by one immediate split of a single
old parent.  Old-parent mappings remain authoritative for elastic state
inheritance and their direct ratios remain in the certificate as diagnostics.
"""
from __future__ import annotations

import math
import os
from typing import Any

import numpy as np

from . import crack_backend as _cb
from . import mode_i_first_passage_v9_18_5 as _v9185
from . import mode_i_first_passage_v9_18_5_4 as _v91854
from . import mode_i_first_passage_v9_18_5_6 as _v91856


MODEL_ID = "certified_atomic_path_corridor_quality_wrapper_v10_0_5_18_3_9"


def _float_env(name: str, default: float) -> float:
    try:
        value = float(os.environ.get(name, default))
    except (TypeError, ValueError):
        value = float(default)
    return value if math.isfinite(value) else float(default)


def strict_quality_advance_v10051839(self: Any, *args, **kwargs):
    original = strict_quality_advance_v10051839._original
    old_mesh = kwargs["mesh"]
    snap = self._transaction_snapshot()
    result = original(self, *args, **kwargs)
    if not bool(getattr(result, "inserted", False)):
        return _v91856._record_or_raise(self, kwargs, result)

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

    certificate = dict(
        getattr(result, "v10051839_quality_certificate", {}) or {}
    )
    certified = bool(
        certificate.get("certificate_schema")
        == "v10.0.5.18.3.9_local_support_area_quality_certificate"
    )

    parent_map = getattr(result, "elem_parent_map", None)
    old_parent_ratios: list[float] = []
    if parent_map is not None:
        pm = np.asarray(parent_map, dtype=int)
        for e in affected:
            if e < len(pm) and 0 <= pm[e] < int(old_mesh.ne):
                old_parent_ratios.append(
                    float(new_mesh.area_e[e])
                    / max(float(old_mesh.area_e[pm[e]]), 1.0e-300)
                )
    else:
        for e in affected:
            if e < int(old_mesh.ne):
                old_parent_ratios.append(
                    float(new_mesh.area_e[e])
                    / max(float(old_mesh.area_e[e]), 1.0e-300)
                )
    min_old_parent_ratio = (
        float(min(old_parent_ratios)) if old_parent_ratios else 1.0
    )

    if certified:
        amin = float(certificate.get("min_local_support_area_ratio", -np.inf))
        area_metric = str(certificate.get("area_gate_metric"))
    else:
        amin = min_old_parent_ratio
        area_metric = "new_element_area_over_mapped_old_parent_area"

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
        "front_id": int(kwargs.get("front_id", -1)),
        "min_triangle_quality": qmin,
        "triangle_quality_floor": qfloor,
        "min_child_area_ratio": amin,
        "child_area_ratio_floor": afloor,
        "area_gate_metric": area_metric,
        "atomic_path_corridor_certificate_active": certified,
        "min_old_parent_area_ratio_diagnostic": min_old_parent_ratio,
        "active_tip_h_over_da": tip_ratio,
        "requested_tip_h_over_da": requested_tip_ratio,
        "tip_h_over_da_enforced_as_veto": False,
        "resolution_warning": warning,
        "requested_da_m": da,
        "affected_element_count": int(len(affected)),
        "accepted": not bool(issues),
        "issues": issues,
        "resolution_metric": "active_endpoint_one_ring_unique_nonzero_edge_mean",
        "quality_certificate": certificate if certified else None,
        **local,
    }

    if not issues:
        new_mesh.hbar_tip = float(local["active_tip_h_mean_m"])
        if getattr(self, "advance_log", None):
            for log_row in self.advance_log:
                if log_row.get("physical_event_index") == len(
                    getattr(self, "advance_log", [])
                ):
                    continue
            self.advance_log[-1].update(
                {
                    "v91856_quality_gate_passed": True,
                    "v91856_min_triangle_quality": qmin,
                    "v91856_min_child_area_ratio": amin,
                    "v91856_area_gate_metric": area_metric,
                    "v10051839_quality_certificate_active": certified,
                    "v91856_active_tip_h_over_da": tip_ratio,
                    "v91856_resolution_warning": warning,
                    "v91856_active_tip_h_mean_m": local["active_tip_h_mean_m"],
                    "v91856_legacy_stored_hbar_tip_m": local[
                        "legacy_stored_hbar_tip_m"
                    ],
                }
            )
        _v91856._AUDIT["accepted_events"].append(row)
        if warning:
            _v91856._AUDIT["resolution_warnings"].append(row.copy())
        return _v91856._record_or_raise(self, kwargs, result)

    _v91856._AUDIT["quality_vetoes"].append(row)
    _v9185._RUNTIME["quality_vetoes"].append(row)
    self._transaction_rollback(snap)
    veto = _cb.CrackAdvanceResult(
        mesh=old_mesh,
        boundary=kwargs["boundary"],
        damage=kwargs["damage"],
        displacement=kwargs["displacement"],
        moved=0.0,
        inserted=False,
        angle_error_deg=float(getattr(result, "angle_error_deg", 0.0)),
        reason="v10051839_quality_veto:" + ";".join(issues),
        elem_parent_map=None,
    )
    return _v91856._record_or_raise(self, kwargs, veto)


__all__ = ["MODEL_ID", "strict_quality_advance_v10051839"]
