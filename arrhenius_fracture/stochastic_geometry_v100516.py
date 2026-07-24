"""PF-equivalent stochastic crack-event geometry for FEM/CZM v10.0.5.16.

The constitutive front engine enqueues one realized event length after a stochastic
cleavage renewal completes.  The existing adaptive-CZM backend remains the only
geometry transaction; this adapter changes only the requested endpoint distance
from the nominal checkpoint to the threshold-correlated PF reward.
"""
from __future__ import annotations

from collections import deque
from contextlib import contextmanager
import copy
import json
import math
from pathlib import Path
from typing import Any, Iterator

import numpy as np

from . import crack_backend as _crack_backend

GEOMETRY_MODEL = "PF_v10_1_7_3_threshold_correlated_single_checked_CZM_commit_v10_0_5_16"

_PENDING_EVENTS: deque[dict[str, Any]] = deque()
_REALIZED_EVENTS: list[dict[str, Any]] = []


def clear_stochastic_geometry_state() -> None:
    _PENDING_EVENTS.clear()
    _REALIZED_EVENTS.clear()


def enqueue_stochastic_geometry_event(event: dict[str, Any]) -> None:
    payload = copy.deepcopy(dict(event))
    length = float(payload.get("event_advance_m", 0.0))
    if not np.isfinite(length) or length <= 0.0:
        raise ValueError("stochastic geometry event requires a positive finite length")
    _PENDING_EVENTS.append(payload)


def pending_stochastic_geometry_events() -> int:
    return len(_PENDING_EVENTS)


def realized_stochastic_geometry_events() -> list[dict[str, Any]]:
    return copy.deepcopy(_REALIZED_EVENTS)


class _StochasticLengthBackend:
    """Delegate to the validated backend after replacing only the event endpoint."""

    def __init__(self, delegate: Any):
        self._delegate = delegate

    def __getattr__(self, name: str) -> Any:
        return getattr(self._delegate, name)

    @property
    def name(self) -> str:
        return str(self._delegate.name)

    def advance(self, **kwargs):
        if not _PENDING_EVENTS:
            raise RuntimeError(
                "v10.0.5.16 cohesive geometry advance has no pending stochastic "
                "cleavage event; constitutive and geometry transactions are desynchronized"
            )
        event = _PENDING_EVENTS.popleft()
        p0 = np.asarray(kwargs["p0"], dtype=float).reshape(2)
        direction = np.asarray(kwargs["direction"], dtype=float).reshape(2)
        norm = float(np.linalg.norm(direction))
        if not np.isfinite(norm) or norm <= 1.0e-30:
            raise RuntimeError("stochastic geometry event has an invalid crack direction")
        direction = direction / norm
        requested_p1 = np.asarray(kwargs["p1"], dtype=float).reshape(2)
        nominal_length = float(np.linalg.norm(requested_p1 - p0))
        event_length = float(event["event_advance_m"])
        realized_request = p0 + event_length * direction

        local = dict(kwargs)
        local["p1"] = realized_request
        local["direction"] = direction
        result = self._delegate.advance(**local)
        moved = float(getattr(result, "moved", 0.0))
        inserted = bool(getattr(result, "inserted", False))
        synchronized = (not inserted) or math.isclose(
            moved,
            event_length,
            rel_tol=1.0e-6,
            abs_tol=5.0e-10,
        )

        record = {
            "schema": GEOMETRY_MODEL,
            **copy.deepcopy(event),
            "front_id": int(kwargs.get("front_id", 0)),
            "requested_fixed_length_m": nominal_length,
            "requested_stochastic_length_m": event_length,
            "p0_m": p0.tolist(),
            "nominal_p1_m": requested_p1.tolist(),
            "stochastic_p1_requested_m": realized_request.tolist(),
            "inserted": inserted,
            "moved_m": moved,
            "event_length_mismatch_m": moved - event_length,
            "microstructure_geometry_length_synchronized": synchronized,
            "reason": str(getattr(result, "reason", "")),
            "angle_error_deg": float(getattr(result, "angle_error_deg", 0.0)),
        }
        logs = getattr(self._delegate, "advance_log", None)
        if logs:
            last = logs[-1]
            record["actual_p1_m"] = [float(last["x1"]), float(last["y1"])]
        _REALIZED_EVENTS.append(record)
        if inserted and not synchronized:
            raise RuntimeError(
                "adaptive-CZM stochastic event length does not match the already "
                "translated MPZ: "
                f"requested={event_length:.9e} m moved={moved:.9e} m"
            )
        return result


@contextmanager
def installed_stochastic_geometry_v100516() -> Iterator[None]:
    """Install the endpoint adapter while preserving backend identity and quality gates."""

    original = _crack_backend.build_crack_backend
    clear_stochastic_geometry_state()

    def build(local_args, geometry):
        return _StochasticLengthBackend(original(local_args, geometry))

    _crack_backend.build_crack_backend = build
    try:
        yield
    finally:
        _crack_backend.build_crack_backend = original


def write_stochastic_geometry_events(out: str | Path) -> Path:
    target = Path(out).expanduser().resolve() / "stochastic_geometry_events_v10_0_5_16.json"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(_REALIZED_EVENTS, indent=2, sort_keys=True) + "\n")
    return target


__all__ = [
    "GEOMETRY_MODEL",
    "clear_stochastic_geometry_state",
    "enqueue_stochastic_geometry_event",
    "installed_stochastic_geometry_v100516",
    "pending_stochastic_geometry_events",
    "realized_stochastic_geometry_events",
    "write_stochastic_geometry_events",
]
