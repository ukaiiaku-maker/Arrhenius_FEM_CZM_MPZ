"""v9.14 direct Mode-I entry using event-driven adaptive-CZM refinement.

The direct Mode-I runner is routed through the mature directional
``_advance_polyline`` path while the admissible cleavage inventory is constrained
to one exactly forward plane.  After every successful cohesive insertion a
single correction solve is scheduled with both ``dt=0`` and ``dU=0``.  The
refined mesh is therefore re-equilibrated at the same physical event time and
remote displacement before loading and Arrhenius clocks resume.
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


class _PostEventEquilibriumController:
    def __init__(self) -> None:
        self.pending = False
        self.dt_seen = False
        self.events_scheduled = 0
        self.corrections_consumed = 0

    def schedule(self) -> None:
        self.pending = True
        self.dt_seen = False
        self.events_scheduled += 1

    def value(self, base: float, role: str, factor: Any) -> float:
        factor_f = float(factor)
        if not self.pending:
            return float(base) * factor_f
        if role == "dt":
            self.dt_seen = True
            return 0.0
        if role == "dU":
            # The solver evaluates dt before dU in each trial.  Consume the
            # correction only after both have been replaced by zero.
            if self.dt_seen:
                self.pending = False
                self.dt_seen = False
                self.corrections_consumed += 1
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
    _inject_once(user_args, "--no-crystal-branch")
    _inject_once(user_args, "--crack-backend", "adaptive_czm")
    _inject_once(user_args, "--adaptive-events")
    _inject_once(user_args, "--min-global-forward", "0.999999")

    controller = _PostEventEquilibriumController()
    original_planes = _crystal.bcc_cleavage_traces
    original_config = _sharp_front.make_emergent_config
    original_advance = _crack_backend.AdaptiveCZMBackend.advance

    def configured_model():
        cfg = original_config()
        cfg.loading = _LoadingProxy(cfg.loading, controller)
        return cfg

    def advance_and_schedule(self, *args, **kwargs):
        result = original_advance(self, *args, **kwargs)
        if bool(getattr(result, "inserted", False)) and float(getattr(result, "moved", 0.0)) > 0.0:
            controller.schedule()
        return result

    _crystal.bcc_cleavage_traces = _forward_mode_i_plane
    _sharp_front.make_emergent_config = configured_model
    _crack_backend.AdaptiveCZMBackend.advance = advance_and_schedule
    try:
        results = _base.main(user_args)
    finally:
        _crystal.bcc_cleavage_traces = original_planes
        _sharp_front.make_emergent_config = original_config
        _crack_backend.AdaptiveCZMBackend.advance = original_advance

    out_value = _option_value(user_args, "--out")
    if out_value is not None:
        out = Path(out_value)
        payload = {
            "schema": "post_event_same_load_equilibrium_v914",
            "events_scheduled": controller.events_scheduled,
            "corrections_consumed": controller.corrections_consumed,
            "pending_at_exit": controller.pending,
            "correction_definition": "one accepted FEM/J solve with dt=0 and dU=0 after each successful adaptive-CZM insertion",
            "all_scheduled_corrections_consumed": (
                controller.events_scheduled == controller.corrections_consumed and not controller.pending
            ),
        }
        out.mkdir(parents=True, exist_ok=True)
        (out / "post_event_equilibrium_audit_v914.json").write_text(json.dumps(payload, indent=2))
    return results


if __name__ == "__main__":
    main()
