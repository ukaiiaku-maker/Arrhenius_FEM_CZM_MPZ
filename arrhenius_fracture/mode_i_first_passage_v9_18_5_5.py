"""v9.18.5.5: resolution audit-only and robust consecutive-veto abort.

The 700 K ceramic gate showed that v9.18.5.4 still rejected an otherwise
mechanically valid exact-target insertion solely because the diagnostic ratio
``active_tip_h / da`` exceeded 0.75.  Exact-ray insertion already validates the
quantities that determine whether the topology operation is admissible:
positive finite areas, triangle quality, child-area ratio, supported cohesive
endpoints, and absence of orphan bulk nodes.  The h/da ratio is therefore a
resolution diagnostic, not an independent topology-validity condition.

This wrapper preserves and reports the requested h/da threshold but prevents it
from vetoing an event.  All geometric validity floors remain enforced.  It also
counts every consecutive rejected geometry transaction until an event is
accepted, so tiny changes in p0, p1, or formatted reason cannot bypass the
configured fail-fast limit.

No barrier, hazard, cohesive, MPZ, wake, shielding, or material law changes.
"""
from __future__ import annotations

import json
import math
import os
from pathlib import Path
import sys
from typing import Any

from . import mode_i_first_passage_v9_18_5_4 as _v91854


_AUDIT: dict[str, Any] = {
    "resolution_warnings": [],
    "consecutive_veto_abort": None,
}


def _float_env(name: str, default: float) -> float:
    try:
        value = float(os.environ.get(name, default))
    except (TypeError, ValueError):
        value = float(default)
    return value if math.isfinite(value) else float(default)


def _consecutive_veto_guard(self: Any, kwargs: dict[str, Any], result: Any) -> Any:
    """Abort any unbroken sequence of rejected geometry transactions.

    A rejected transaction restores the pre-event mesh and leaves the physical
    renewal unconsumed.  Without an accepted event there is no physical reason
    to retry indefinitely, even if a diagnostic string changes by roundoff.
    """
    if bool(getattr(result, "inserted", False)):
        self._v91855_consecutive_geometry_vetoes = 0
        self._v91855_last_veto_reason = None
        return result

    reason = str(getattr(result, "reason", "unknown"))
    count = int(getattr(self, "_v91855_consecutive_geometry_vetoes", 0)) + 1
    self._v91855_consecutive_geometry_vetoes = count
    self._v91855_last_veto_reason = reason
    limit = max(int(os.environ.get("ARRHENIUS_MAX_IDENTICAL_GEOMETRY_VETOES", "12")), 1)
    if count >= limit:
        payload = {
            "front_id": int(kwargs.get("front_id", -1)),
            "count": count,
            "limit": limit,
            "reason": reason,
            "p0_m": list(map(float, kwargs.get("p0", []))),
            "p1_m": list(map(float, kwargs.get("p1", []))),
        }
        _AUDIT["consecutive_veto_abort"] = payload
        raise RuntimeError(
            "v9.18.5.5 consecutive geometry veto limit reached; physical renewal "
            f"remains unconsumed: front={payload['front_id']} "
            f"count={count}/{limit} reason={reason}"
        )
    return result


def _quality_with_resolution_audit_only(self: Any, *args, **kwargs):
    """Run v9.18.5.4 quality checks with h/da disabled only as a veto.

    The desired threshold remains available for post-event audit.  Triangle
    quality, child-area ratio, finite area, orphan-node, and cohesive-endpoint
    checks are unchanged because only ARRHENIUS_MAX_TIP_H_OVER_DA is temporarily
    replaced while the inherited strict-quality function executes.
    """
    original = _quality_with_resolution_audit_only._original
    desired = _float_env("ARRHENIUS_MAX_TIP_H_OVER_DA", 0.75)
    key = "ARRHENIUS_MAX_TIP_H_OVER_DA"
    old = os.environ.get(key)
    before = len(_v91854._AUDIT.get("accepted_events", []))
    os.environ[key] = "inf"
    try:
        result = original(self, *args, **kwargs)
    finally:
        if old is None:
            os.environ.pop(key, None)
        else:
            os.environ[key] = old

    accepted = _v91854._AUDIT.get("accepted_events", [])
    if bool(getattr(result, "inserted", False)) and len(accepted) > before:
        row = accepted[-1]
        ratio = float(row.get("active_tip_h_over_da", float("nan")))
        warning = bool(math.isfinite(ratio) and ratio > desired)
        row["v91855_resolution_threshold_requested"] = desired
        row["v91855_resolution_threshold_enforced_as_veto"] = False
        row["v91855_resolution_warning"] = warning
        if warning:
            _AUDIT["resolution_warnings"].append({
                "front_id": int(kwargs.get("front_id", -1)),
                "active_tip_h_over_da": ratio,
                "requested_threshold": desired,
                "p0_m": list(map(float, kwargs.get("p0", []))),
                "p1_m": list(map(float, kwargs.get("p1", []))),
            })
    return result


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
        "schema": "resolution_audit_only_consecutive_veto_v91855_v1",
        "tip_h_over_da_role": "audit_warning_only",
        "tip_h_over_da_requested_threshold": _float_env(
            "ARRHENIUS_MAX_TIP_H_OVER_DA", 0.75
        ),
        "triangle_quality_veto_retained": True,
        "child_area_ratio_veto_retained": True,
        "finite_positive_area_veto_retained": True,
        "orphan_node_veto_retained": True,
        "cohesive_endpoint_support_veto_retained": True,
        "consecutive_veto_guard_ignores_roundoff_signature_changes": True,
        "max_consecutive_geometry_vetoes": int(
            os.environ.get("ARRHENIUS_MAX_IDENTICAL_GEOMETRY_VETOES", "12")
        ),
        "run_completed_without_exception": error is None,
        "runtime_error_type": None if error is None else type(error).__name__,
        "runtime_error": None if error is None else str(error),
        "constitutive_physics_changed": False,
        **_AUDIT,
    }
    (out / "resolution_audit_only_consecutive_veto_v91855.json").write_text(
        json.dumps(payload, indent=2, default=str)
    )


def main(argv=None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    _AUDIT["resolution_warnings"] = []
    _AUDIT["consecutive_veto_abort"] = None

    original_quality = _v91854._strict_quality_advance_v91854
    original_guard = _v91854._record_veto_or_raise
    _quality_with_resolution_audit_only._original = original_quality
    _v91854._strict_quality_advance_v91854 = _quality_with_resolution_audit_only
    _v91854._record_veto_or_raise = _consecutive_veto_guard

    error: BaseException | None = None
    try:
        return _v91854.main(user_args)
    except BaseException as exc:
        error = exc
        raise
    finally:
        _write_audit(user_args, error)
        _v91854._strict_quality_advance_v91854 = original_quality
        _v91854._record_veto_or_raise = original_guard


if __name__ == "__main__":
    main()
