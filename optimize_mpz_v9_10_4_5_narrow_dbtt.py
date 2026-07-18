#!/usr/bin/env python3
"""Crash-safe first-passage search wrapper for v9.10.4.5.

The v9.10.4 search objective intentionally returns a compact result when a
candidate is incomplete or rejected before the transition calculation. The
result-export loop, however, assumes that every detailed evaluation contains
``parameters``, ``temperature_detail``, and ``event_detail``. This wrapper
stabilizes that schema before delegating to the existing search driver.

No objective, constitutive rate, timestep, transition gate, or optimization
algorithm is changed.
"""
from __future__ import annotations

from typing import Any

import numpy as np

import optimize_mpz_v9_10_4_narrow_dbtt as _base


_ORIGINAL_EVALUATE = _base.NarrowDBTTObjective.evaluate


def stabilize_detailed_result(
    result: dict[str, Any],
    x: np.ndarray,
    *,
    details: bool,
) -> dict[str, Any]:
    """Return a stable export schema for every detailed objective result."""
    stable = dict(result)
    if not details:
        return stable

    parameters = stable.get("parameters")
    if not isinstance(parameters, dict):
        try:
            parameters = _base.decode(np.asarray(x, dtype=float))
        except Exception:
            parameters = {}
    stable["parameters"] = parameters

    temperature_detail = stable.get("temperature_detail")
    stable["temperature_detail"] = (
        list(temperature_detail) if isinstance(temperature_detail, list) else []
    )

    event_detail = stable.get("event_detail")
    stable["event_detail"] = list(event_detail) if isinstance(event_detail, list) else []

    if float(stable.get("completion_loss", 0.0)) > 0.0:
        stable.setdefault("evaluation_status", "INCOMPLETE_CANDIDATE")
    elif bool(stable.get("invalid_parameter_vector", False)):
        stable.setdefault("evaluation_status", "INVALID_PARAMETER_VECTOR")
    elif "transition_shelf_ratio" not in stable:
        stable.setdefault("evaluation_status", "EARLY_REJECTED_CANDIDATE")
    else:
        stable.setdefault("evaluation_status", "COMPLETE_CANDIDATE")
    return stable


def _stable_evaluate(self, x: np.ndarray, *, details: bool = False) -> dict[str, Any]:
    result = _ORIGINAL_EVALUATE(self, x, details=details)
    return stabilize_detailed_result(result, np.asarray(x, dtype=float), details=details)


_base.NarrowDBTTObjective.evaluate = _stable_evaluate


if __name__ == "__main__":
    _base.main()
