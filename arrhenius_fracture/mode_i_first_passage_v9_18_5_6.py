"""v9.18.5.6: explicit production-quality wrapper chain.

v9.18.5.5 attempted to make h_tip/da audit-only by wrapping the v9.18.5.4
quality function.  The v9.18.5 runtime subsequently assigned ``_original`` on
whichever function occupied its module-global quality slot.  That reassignment
caused the nested v9.18.5.5 wrapper to call the raw adaptive-CZM backend rather
than the inherited strict-quality gate.  The run therefore completed with zero
accepted quality-gate records.

This revision installs one self-contained function directly into the v9.18.5
quality slot before entering the v9.18.5.3 corridor/startup chain.  v9.18.5 then
binds that function directly to the raw backend advance implementation.  The
function itself performs every production validity check and treats h_tip/da as
an audit warning only.

No barrier, hazard, cohesive-opening, MPZ, wake, shielding, loading, or material
law is changed.
"""
from __future__ import annotations

import json
import math
import os
from pathlib import Path
import sys
from typing import Any

import numpy as np

from . import crack_backend as _cb
from . import mode_i_first_passage_v9_18_5 as _v9185
from . import mode_i_first_passage_v9_18_5_3 as _v91853
from . import mode_i_first_passage_v9_18_5_4 as _v91854


_AUDIT: dict[str, Any] = {
    "accepted_events": [],
    "resolution_warnings": [],
    "quality_vetoes": [],
    "consecutive_veto_abort": None,
}


def _float_env(name: str, default: float) -> float:
    try:
        value = float(os.environ.get(name, default))
    except (TypeError, ValueError):
        value = float(default)
    return value if math.isfinite(value) else float(default)


def _record_or_raise(self: Any, kwargs: dict[str, Any], result: Any) -> Any:
    """Abort any uninterrupted sequence of rejected geometry transactions."""
    if bool(getattr(result, "inserted", False)):
        self._v91856_consecutive_geometry_vetoes = 0
        self._v91856_last_veto_reason = None
        return result

    reason = str(getattr(result, "reason", "unknown"))
    count = int(getattr(self, "_v91856_consecutive_geometry_vetoes", 0)) + 1
    self._v91856_consecutive_geometry_vetoes = count
    self._v91856_last_veto_reason = reason
    limit = max(int(os.environ.get("ARRHENIUS_MAX_IDENTICAL_GEOMETRY_VETOES", "12")), 1)
    if count >= limit:
        payload = {
            "front_id": int(kwargs.get("front_id", -1)),
            "count": count,
            "limit": limit,
            "reason": reason,
            "p0_m": np.asarray(kwargs.get("p0", []), float).tolist(),
            "p1_m": np.asarray(kwargs.get("p1", []), float).tolist(),
        }
        _AUDIT["consecutive_veto_abort"] = payload
        raise RuntimeError(
            "v9.18.5.6 consecutive geometry veto limit reached; physical renewal "
            f"remains unconsumed: front={payload['front_id']} "
            f"count={count}/{limit} reason={reason}"
        )
    return result


def _strict_quality_advance_v91856(self: Any, *args, **kwargs):
    """Raw backend advance plus the complete production geometry gate.

    ``v9.18.5.main`` assigns ``_original`` on this exact function before it is
    installed as ``AdaptiveCZMBackend.advance``.  There is therefore no nested
    wrapper whose ``_original`` attribute can be rebound to the wrong callable.
    """
    original = _strict_quality_advance_v91856._original
    old_mesh = kwargs["mesh"]
    snap = self._transaction_snapshot()
    result = original(self, *args, **kwargs)
    if not bool(getattr(result, "inserted", False)):
        return _record_or_raise(self, kwargs, result)

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

    parent_map = getattr(result, "elem_parent_map", None)
    ratios: list[float] = []
    if parent_map is not None:
        pm = np.asarray(parent_map, dtype=int)
        for e in affected:
            if e < len(pm) and 0 <= pm[e] < int(old_mesh.ne):
                ratios.append(
                    float(new_mesh.area_e[e]) /
                    max(float(old_mesh.area_e[pm[e]]), 1.0e-300)
                )
    else:
        for e in affected:
            if e < int(old_mesh.ne):
                ratios.append(
                    float(new_mesh.area_e[e]) /
                    max(float(old_mesh.area_e[e]), 1.0e-300)
                )
    amin = float(min(ratios)) if ratios else 1.0

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
            if int(nid) < 0 or int(nid) >= int(new_mesh.nn) or incidence[int(nid)] <= 0:
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

    if not issues:
        new_mesh.hbar_tip = float(local["active_tip_h_mean_m"])
        if getattr(self, "advance_log", None):
            self.advance_log[-1].update({
                "v91856_quality_gate_passed": True,
                "v91856_min_triangle_quality": qmin,
                "v91856_min_child_area_ratio": amin,
                "v91856_active_tip_h_over_da": tip_ratio,
                "v91856_resolution_warning": warning,
                "v91856_active_tip_h_mean_m": local["active_tip_h_mean_m"],
                "v91856_legacy_stored_hbar_tip_m": local["legacy_stored_hbar_tip_m"],
            })
        _AUDIT["accepted_events"].append(row)
        if warning:
            _AUDIT["resolution_warnings"].append(row.copy())
        return _record_or_raise(self, kwargs, result)

    _AUDIT["quality_vetoes"].append(row)
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
        reason="v91856_quality_veto:" + ";".join(issues),
        elem_parent_map=None,
    )
    return _record_or_raise(self, kwargs, veto)


def _option_value(argv: list[str], name: str) -> str | None:
    for i, token in enumerate(argv):
        if token == name and i + 1 < len(argv):
            return argv[i + 1]
        if token.startswith(name + "="):
            return token.split("=", 1)[1]
    return None


def _write_audit(argv: list[str], error: BaseException | None) -> None:
    out_value = _option_value(argv, "--out")
    if out_value is None:
        return
    out = Path(out_value)
    out.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema": "explicit_quality_wrapper_chain_v91856_v1",
        "single_self_contained_quality_wrapper_installed": True,
        "nested_original_rebinding_bypass_removed": True,
        "tip_h_over_da_role": "audit_warning_only",
        "tip_h_over_da_requested_threshold": _float_env(
            "ARRHENIUS_MAX_TIP_H_OVER_DA", 0.75
        ),
        "triangle_quality_veto_retained": True,
        "child_area_ratio_veto_retained": True,
        "finite_positive_area_veto_retained": True,
        "orphan_node_veto_retained": True,
        "cohesive_endpoint_support_veto_retained": True,
        "max_consecutive_geometry_vetoes": int(
            os.environ.get("ARRHENIUS_MAX_IDENTICAL_GEOMETRY_VETOES", "12")
        ),
        "run_completed_without_exception": error is None,
        "runtime_error_type": None if error is None else type(error).__name__,
        "runtime_error": None if error is None else str(error),
        "constitutive_physics_changed": False,
        **_AUDIT,
    }
    (out / "explicit_quality_wrapper_chain_v91856.json").write_text(
        json.dumps(payload, indent=2, default=str)
    )


def main(argv=None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    _AUDIT["accepted_events"] = []
    _AUDIT["resolution_warnings"] = []
    _AUDIT["quality_vetoes"] = []
    _AUDIT["consecutive_veto_abort"] = None

    original = _v9185._strict_quality_advance
    _v9185._strict_quality_advance = _strict_quality_advance_v91856
    error: BaseException | None = None
    try:
        # Enter below v9.18.5.4/v9.18.5.5 so no nested quality wrapper can have
        # its ``_original`` pointer rebound by v9.18.5.main.
        return _v91853.main(user_args)
    except BaseException as exc:
        error = exc
        raise
    finally:
        _write_audit(user_args, error)
        _v9185._strict_quality_advance = original


if __name__ == "__main__":
    main()
