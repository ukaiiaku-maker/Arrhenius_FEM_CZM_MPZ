"""v9.14 direct Mode-I entry using event-driven adaptive-CZM refinement.

The direct Mode-I runner is routed through the mature directional
``_advance_polyline`` path while the admissible cleavage inventory is constrained
to one exactly forward plane.  A same-load correction is scheduled once per
accepted *physical Arrhenius renewal*, not once per internal exact-ray CZM
subsegment.  Geometry-veto rollback cancels the corresponding pending correction.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys
from typing import Any

import numpy as np

from . import crack_backend as _crack_backend
from . import crystal as _crystal
from . import mode_i_first_passage_v9_13 as _base
from . import sharp_front as _sharp_front
from .mpz_front_engine import MovingProcessZoneFrontEngine


class _PostEventEquilibriumController:
    """Bookkeep one zero-time/zero-load correction per physical renewal."""

    def __init__(self) -> None:
        self.pending_event_id: int | None = None
        self.dt_seen = False
        self.next_event_id = 1
        self.events_scheduled = 0
        self.corrections_consumed = 0
        self.events_cancelled = 0
        self.duplicate_schedule_calls = 0
        self.scheduled_event_ids: list[int] = []
        self.corrected_event_ids: list[int] = []
        self.cancelled_event_ids: list[int] = []

    @property
    def pending(self) -> bool:
        return self.pending_event_id is not None

    def schedule_physical_event(self) -> int:
        # A single-front accepted solve can contain several adaptive-CZM
        # subsegments, but only one Arrhenius renewal.  Do not multiply the
        # correction count when a lower-level routine is called repeatedly.
        if self.pending_event_id is not None:
            self.duplicate_schedule_calls += 1
            return int(self.pending_event_id)
        event_id = int(self.next_event_id)
        self.next_event_id += 1
        self.pending_event_id = event_id
        self.dt_seen = False
        self.events_scheduled += 1
        self.scheduled_event_ids.append(event_id)
        return event_id

    # Backward-compatible test/helper name.
    def schedule(self) -> int:
        return self.schedule_physical_event()

    def cancel_pending_event(self) -> int | None:
        event_id = self.pending_event_id
        if event_id is not None:
            self.events_cancelled += 1
            self.cancelled_event_ids.append(int(event_id))
        self.pending_event_id = None
        self.dt_seen = False
        return event_id

    def value(self, base: float, role: str, factor: Any) -> float:
        factor_f = float(factor)
        if self.pending_event_id is None:
            return float(base) * factor_f
        if role == "dt":
            self.dt_seen = True
            return 0.0
        if role == "dU":
            # The solver evaluates dt before dU in each trial. Consume the
            # correction only after both have been replaced by zero.
            if self.dt_seen:
                event_id = int(self.pending_event_id)
                self.pending_event_id = None
                self.dt_seen = False
                self.corrections_consumed += 1
                self.corrected_event_ids.append(event_id)
            return 0.0
        return float(base) * factor_f


class _ControlledScalar(float):
    def __new__(cls, value: float, controller: _PostEventEquilibriumController, role: str):
        obj = float.__new__(cls, float(value))
        obj.controller = controller
        obj.role = role
        return obj

    def __mul__(self, other):
        return self.controller.value(float(self), self.role, other)

    def __rmul__(self, other):
        return self.controller.value(float(self), self.role, other)


class _LoadingProxy:
    """Delegate loading settings while wrapping assigned dt and dU scalars."""

    def __init__(self, base: Any, controller: _PostEventEquilibriumController):
        object.__setattr__(self, "_base", base)
        object.__setattr__(self, "_controller", controller)

    def __getattr__(self, name: str):
        return getattr(self._base, name)

    def __setattr__(self, name: str, value: Any) -> None:
        if name == "dt":
            setattr(self._base, name, _ControlledScalar(value, self._controller, "dt"))
        elif name == "dU_top":
            setattr(self._base, name, _ControlledScalar(value, self._controller, "dU"))
        else:
            setattr(self._base, name, value)


def _forward_mode_i_plane(theta_deg: float, *args: Any, **kwargs: Any):
    del theta_deg, args, kwargs
    return [{
        "name": "v914_forward_mode_I",
        "family": "mode_I_forward",
        "angle_deg": 0.0,
        "t": np.array([1.0, 0.0], dtype=float),
        "n": np.array([0.0, 1.0], dtype=float),
        "gamma_rel": 1.0,
    }]


def _inject_once(argv: list[str], name: str, value: str | None = None) -> None:
    if any(token == name or token.startswith(name + "=") for token in argv):
        return
    argv.append(name)
    if value is not None:
        argv.append(value)


def _option_value(argv: list[str], name: str) -> str | None:
    for i, token in enumerate(argv):
        if token == name and i + 1 < len(argv):
            return argv[i + 1]
        if token.startswith(name + "="):
            return token.split("=", 1)[1]
    return None


def main(argv=None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    _inject_once(user_args, "--crystal-aniso")
    # Branching is disabled by the inherited parser default.
    _inject_once(user_args, "--crack-backend", "adaptive_czm")
    _inject_once(user_args, "--adaptive-events")
    _inject_once(user_args, "--min-global-forward", "0.999999")

    controller = _PostEventEquilibriumController()
    original_planes = _crystal.bcc_cleavage_traces
    original_config = _sharp_front.make_emergent_config
    original_backend_advance = _crack_backend.AdaptiveCZMBackend.advance
    original_engine_step = MovingProcessZoneFrontEngine.step
    original_geometry_veto = MovingProcessZoneFrontEngine.restore_geometry_veto

    def configured_model():
        cfg = original_config()
        cfg.loading = _LoadingProxy(cfg.loading, controller)
        return cfg

    def step_and_schedule(self, *args, **kwargs):
        out = original_engine_step(self, *args, **kwargs)
        if int(out.get("n_fire", 0) or 0) > 0:
            out["v914_physical_event_id"] = controller.schedule_physical_event()
        return out

    def veto_and_cancel(self, n_restore: int):
        result = original_geometry_veto(self, n_restore)
        controller.cancel_pending_event()
        return result

    def advance_and_tag(self, *args, **kwargs):
        before = len(getattr(self, "advance_log", []))
        result = original_backend_advance(self, *args, **kwargs)
        event_id = controller.pending_event_id
        new_rows = getattr(self, "advance_log", [])[before:]
        for subsegment_index, row in enumerate(new_rows):
            row["physical_event_id"] = int(event_id) if event_id is not None else -1
            row["physical_subsegment_index"] = int(subsegment_index)
            row["physical_subsegment_count"] = int(len(new_rows))
        return result

    _crystal.bcc_cleavage_traces = _forward_mode_i_plane
    _sharp_front.make_emergent_config = configured_model
    MovingProcessZoneFrontEngine.step = step_and_schedule
    MovingProcessZoneFrontEngine.restore_geometry_veto = veto_and_cancel
    _crack_backend.AdaptiveCZMBackend.advance = advance_and_tag
    try:
        results = _base.main(user_args)
    finally:
        _crystal.bcc_cleavage_traces = original_planes
        _sharp_front.make_emergent_config = original_config
        MovingProcessZoneFrontEngine.step = original_engine_step
        MovingProcessZoneFrontEngine.restore_geometry_veto = original_geometry_veto
        _crack_backend.AdaptiveCZMBackend.advance = original_backend_advance

    out_value = _option_value(user_args, "--out")
    if out_value is not None:
        out = Path(out_value)
        payload = {
            "schema": "post_event_same_load_equilibrium_v914",
            "events_scheduled": controller.events_scheduled,
            "corrections_consumed": controller.corrections_consumed,
            "events_cancelled": controller.events_cancelled,
            "duplicate_schedule_calls": controller.duplicate_schedule_calls,
            "scheduled_event_ids": controller.scheduled_event_ids,
            "corrected_event_ids": controller.corrected_event_ids,
            "cancelled_event_ids": controller.cancelled_event_ids,
            "pending_event_id": controller.pending_event_id,
            "pending_at_exit": controller.pending,
            "correction_definition": (
                "one accepted FEM/J solve with dt=0 and dU=0 after each physical "
                "Arrhenius renewal; adaptive-CZM subsegments share one event id"
            ),
            "all_scheduled_corrections_consumed": (
                controller.events_scheduled
                == controller.corrections_consumed + controller.events_cancelled
                and not controller.pending
            ),
        }
        out.mkdir(parents=True, exist_ok=True)
        (out / "post_event_equilibrium_audit_v914.json").write_text(
            json.dumps(payload, indent=2)
        )
    return results


if __name__ == "__main__":
    main()
