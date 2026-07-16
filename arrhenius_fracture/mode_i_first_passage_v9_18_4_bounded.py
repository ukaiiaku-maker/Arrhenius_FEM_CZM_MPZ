"""Bounded runtime wrapper for v9.18.4 geometry candidates.

This module keeps the v9.18.4 constitutive and geometry helpers but ensures that
candidate retry count is finite and that backend state is rolled back if the
underlying adaptive advance raises unexpectedly.
"""
from __future__ import annotations

import json
import math
import os
from pathlib import Path
import sys

import numpy as np

from . import crack_backend as _cb
from . import fem as _fem
from . import mode_i_first_passage_v9_18_1 as _v9181
from . import mode_i_first_passage_v9_18_4 as _v9184


def _bounded_insert(self, mesh, displacement, p0, target, front_id):
    result = _v9184._mechanically_regularized_insert_target(
        self, mesh, displacement, p0, target, front_id
    )
    meta = result[4] if len(result) > 4 and isinstance(result[4], dict) else {}
    if bool(meta.get("edge_front_regularization_used", False)):
        self._v9184_regularization_candidate_count = int(
            meta.get("regularization_candidate_count", 1) or 1
        )
    return result


def _bounded_advance(self, *args, **kwargs):
    original = _bounded_advance._original
    snap = self._transaction_snapshot()
    try:
        result = original(self, *args, **kwargs)
    except Exception:
        self._transaction_rollback(snap)
        raise

    if bool(getattr(result, "inserted", False)):
        issues = _v9184._mechanical_topology_issues(self, result)
        if not issues:
            self._v9184_regularization_attempt = 0
            self._v9184_regularization_candidate_count = 0
            self._v9184_last_veto_signature = None
            self._v9184_identical_veto_count = 0
            return result

        self._transaction_rollback(snap)
        attempt = int(getattr(self, "_v9184_regularization_attempt", 0)) + 1
        self._v9184_regularization_attempt = attempt
        candidate_count = int(getattr(self, "_v9184_regularization_candidate_count", 0))
        reason = "mechanical_topology_veto:" + ";".join(issues)
        if candidate_count > 0 and attempt >= candidate_count:
            raise RuntimeError(
                "v9.18.4 exhausted all same-length front-regularization "
                f"candidates without a mechanically valid topology: "
                f"attempts={attempt}/{candidate_count} reason={reason}"
            )
        return _cb.CrackAdvanceResult(
            mesh=kwargs["mesh"],
            boundary=kwargs["boundary"],
            damage=kwargs["damage"],
            displacement=kwargs["displacement"],
            moved=0.0,
            inserted=False,
            angle_error_deg=float(getattr(result, "angle_error_deg", 0.0)),
            reason=reason,
            elem_parent_map=None,
        )

    reason = str(getattr(result, "reason", "unknown"))
    front_id = int(kwargs.get("front_id", -1))
    p0 = np.asarray(kwargs.get("p0", [math.nan, math.nan]), float)
    p1 = np.asarray(kwargs.get("p1", [math.nan, math.nan]), float)
    scale = max(float(getattr(kwargs.get("mesh", None), "hbar_tip", 1.0e-12)), 1.0e-12)
    signature = (
        front_id,
        tuple(np.round(p0 / scale, 8)),
        tuple(np.round(p1 / scale, 8)),
        reason,
        int(getattr(self, "_v9184_regularization_attempt", 0)),
    )
    if signature == getattr(self, "_v9184_last_veto_signature", None):
        count = int(getattr(self, "_v9184_identical_veto_count", 0)) + 1
    else:
        count = 1
    self._v9184_last_veto_signature = signature
    self._v9184_identical_veto_count = count
    limit = max(int(os.environ.get("ARRHENIUS_MAX_IDENTICAL_GEOMETRY_VETOES", "12")), 1)
    if count >= limit:
        raise RuntimeError(
            "v9.18.4 repeated identical geometry veto; physical event remains "
            f"unconsumed: front={front_id} count={count}/{limit} reason={reason}"
        )
    return result


def main(argv=None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    original_insert = _cb.AdaptiveCZMBackend._insert_target_in_incident_triangle
    original_advance = _cb.AdaptiveCZMBackend.advance
    original_solve = _fem.solve_dirichlet
    _bounded_advance._original = original_advance
    _v9184._finite_solve_dirichlet._original = original_solve
    _cb.AdaptiveCZMBackend._insert_target_in_incident_triangle = _bounded_insert
    _cb.AdaptiveCZMBackend.advance = _bounded_advance
    _fem.solve_dirichlet = _v9184._finite_solve_dirichlet
    try:
        results = _v9181.main(user_args)
    finally:
        _cb.AdaptiveCZMBackend._insert_target_in_incident_triangle = original_insert
        _cb.AdaptiveCZMBackend.advance = original_advance
        _fem.solve_dirichlet = original_solve

    out_value = _v9184._option_value(user_args, "--out")
    if out_value is not None:
        out = Path(out_value)
        out.mkdir(parents=True, exist_ok=True)
        payload = {
            "schema": "mechanically_valid_front_regularization_v9184_v2",
            "same_length_interior_regularization_enabled": True,
            "regularization_angles_deg": _v9184._angle_schedule_deg(),
            "mechanical_topology_validation_enabled": True,
            "candidate_retry_count_bounded": True,
            "unexpected_backend_exception_rolls_back": True,
            "nonfinite_mechanics_fail_fast_enabled": True,
            "constitutive_physics_changed": False,
        }
        (out / "mechanical_topology_validation_v9184.json").write_text(
            json.dumps(payload, indent=2, default=str)
        )
    return results


if __name__ == "__main__":
    main()
