"""v9.15 direct Mode-I entry with coupled cohesive-opening/MPZ relaxation.

A completed Arrhenius cleavage renewal nucleates one adaptive-CZM segment.  The
segment is then opened over a finite physical event window at fixed remote
opening.  During that window the FEM equilibrium, J-derived driving force, and
moving-process-zone emission/transport state are advanced together.  Cleavage
renewal is held only while the already-nucleated event is opening.

This is an experimental event-resolution protocol.  The Arrhenius clock remains
the fracture nucleation criterion; no independent critical traction or Gc rule
is introduced.
"""
from __future__ import annotations

import json
import math
import os
from pathlib import Path
import sys
from typing import Any

import numpy as np

from . import crack_backend as _crack_backend
from . import crystal as _crystal
from . import mode_i_first_passage_v9_13 as _base
from . import sharp_front as _sharp_front
from .mpz_front_engine_v911 import MovingProcessZone2DFrontEngine


class CoupledEventRelaxationController:
    """Resolve one nucleated crack quantum over finite time at fixed load."""

    def __init__(self) -> None:
        self.tau_event_s = max(
            float(os.environ.get("ARRHENIUS_EVENT_RELAXATION_TIME_S", "1e-6")),
            1.0e-15,
        )
        self.n_substeps = max(
            int(os.environ.get("ARRHENIUS_EVENT_RELAXATION_SUBSTEPS", "12")), 2
        )
        self.initial_damage = float(np.clip(
            float(os.environ.get("ARRHENIUS_EVENT_INITIAL_DAMAGE", "0")), 0.0, 0.999999
        ))
        self.next_event_id = 1
        self.active_event_id: int | None = None
        self.active_elements: list[Any] = []
        self.active_log_rows: list[dict[str, Any]] = []
        self.substep_index = 0
        self.prepared = False
        self.events: list[dict[str, Any]] = []
        self._active_record: dict[str, Any] | None = None
        self.geometry_retry_attempts = 0
        self.geometry_retry_successes = 0

    @property
    def active(self) -> bool:
        return self.active_event_id is not None

    @property
    def ready_for_relaxation(self) -> bool:
        return self.active and bool(self.active_elements)

    @property
    def dt_substep_s(self) -> float:
        return self.tau_event_s / float(self.n_substeps)

    @staticmethod
    def _smooth_progress(s: float) -> float:
        s = float(np.clip(s, 0.0, 1.0))
        return s * s * (3.0 - 2.0 * s)

    def schedule_event(self, out: dict[str, Any]) -> int:
        if self.active_event_id is not None:
            return int(self.active_event_id)
        event_id = int(self.next_event_id)
        self.next_event_id += 1
        self.active_event_id = event_id
        self.active_elements = []
        self.active_log_rows = []
        self.substep_index = 0
        self.prepared = False
        rec = {
            "event_id": event_id,
            "status": "nucleated_geometry_pending",
            "tau_event_s": float(self.tau_event_s),
            "n_relaxation_substeps_requested": int(self.n_substeps),
            "initial_damage": float(self.initial_damage),
            "KJ_nucleation_Pa_sqrt_m": _finite_or_none(
                out.get("anisotropic_KJ_Pa_sqrt_m")
            ),
            "Kshield_nucleation_Pa_sqrt_m": _finite_or_none(
                out.get("mpz_K_shield_pre_renewal_Pa_sqrt_m")
            ),
            "emitted_total_nucleation": _finite_or_none(
                out.get("mpz_emitted_total_pre_renewal")
            ),
            "substeps": [],
            "geometry_retry_level": 0,
        }
        self.events.append(rec)
        self._active_record = rec
        return event_id

    def register_geometry(
        self,
        elements: list[Any],
        log_rows: list[dict[str, Any]],
        retry_level: int = 0,
    ) -> None:
        if self.active_event_id is None:
            return
        self.active_elements = list(elements)
        self.active_log_rows = list(log_rows)
        for elem in self.active_elements:
            elem.damage = float(self.initial_damage)
            elem.metadata["v915_physical_event_id"] = int(self.active_event_id)
            elem.metadata["v915_coupled_relaxation"] = True
        for i, row in enumerate(self.active_log_rows):
            row["physical_event_id"] = int(self.active_event_id)
            row["physical_subsegment_index"] = int(i)
            row["physical_subsegment_count"] = int(len(self.active_log_rows))
            row["v915_initial_damage"] = float(self.initial_damage)
            row["v915_geometry_retry_level"] = int(retry_level)
        if self._active_record is not None:
            self._active_record.update({
                "status": "relaxation_pending",
                "n_cohesive_subsegments": int(len(self.active_elements)),
                "cohesive_length_m": float(sum(
                    max(float(getattr(e, "length", 0.0)), 0.0)
                    for e in self.active_elements
                )),
                "geometry_retry_level": int(retry_level),
            })

    def cancel_event(self, reason: str = "geometry_veto") -> None:
        if self._active_record is not None:
            self._active_record.update({
                "status": "cancelled",
                "cancel_reason": str(reason),
                "relaxation_completed": False,
            })
        self.active_event_id = None
        self.active_elements = []
        self.active_log_rows = []
        self.substep_index = 0
        self.prepared = False
        self._active_record = None

    def prepare_substep(self) -> float:
        if not self.ready_for_relaxation:
            return 0.0
        if self.prepared:
            return self.dt_substep_s
        self.substep_index += 1
        s = min(float(self.substep_index) / float(self.n_substeps), 1.0)
        q = self.initial_damage + (1.0 - self.initial_damage) * self._smooth_progress(s)
        for elem in self.active_elements:
            elem.damage = float(np.clip(q, 0.0, 1.0))
            elem.clock = float(self.substep_index * self.dt_substep_s)
        self.prepared = True
        return self.dt_substep_s

    def finish_substep(self, out: dict[str, Any]) -> None:
        if not self.ready_for_relaxation or not self.prepared:
            return
        q = max((float(getattr(e, "damage", 0.0)) for e in self.active_elements), default=0.0)
        row = {
            "substep": int(self.substep_index),
            "time_s": float(self.substep_index * self.dt_substep_s),
            "dt_s": float(self.dt_substep_s),
            "damage_progress": float(q),
            "remote_displacement_increment_m": 0.0,
            "KJ_Pa_sqrt_m": _finite_or_none(out.get("anisotropic_KJ_Pa_sqrt_m")),
            "Kshield_Pa_sqrt_m": _finite_or_none(out.get("mpz_K_shield_pre_renewal_Pa_sqrt_m")),
            "sigma_emit_tip_Pa": _finite_or_none(out.get("sigma_emit_tip")),
            "dN_emit": _finite_or_none(out.get("dN_emit")),
            "dN_emit_raw": _finite_or_none(out.get("dN_emit_raw")),
            "mobile_count": _finite_or_none(out.get("mpz_mobile_count")),
            "retained_count": _finite_or_none(out.get("mpz_retained_count")),
            "emitted_total": _finite_or_none(out.get("mpz_emitted_total")),
        }
        if self._active_record is not None:
            self._active_record["substeps"].append(row)
        self.prepared = False
        if self.substep_index >= self.n_substeps:
            for elem in self.active_elements:
                elem.damage = 1.0
                elem.clock = float(self.tau_event_s)
            if self._active_record is not None:
                sub = self._active_record.get("substeps", [])
                emitted = [x.get("dN_emit") for x in sub if x.get("dN_emit") is not None]
                Kvals = [x.get("KJ_Pa_sqrt_m") for x in sub if x.get("KJ_Pa_sqrt_m") is not None]
                self._active_record.update({
                    "status": "complete",
                    "relaxation_completed": True,
                    "n_relaxation_substeps_completed": int(self.substep_index),
                    "final_damage": 1.0,
                    "relaxation_time_s": float(self.tau_event_s),
                    "dN_emit_relaxation": float(sum(emitted)) if emitted else 0.0,
                    "KJ_relaxation_start_Pa_sqrt_m": float(Kvals[0]) if Kvals else None,
                    "KJ_relaxation_end_Pa_sqrt_m": float(Kvals[-1]) if Kvals else None,
                })
            self.active_event_id = None
            self.active_elements = []
            self.active_log_rows = []
            self.substep_index = 0
            self._active_record = None

    def controlled_value(self, base: float, role: str, factor: Any) -> float:
        factor_f = float(factor)
        if not self.ready_for_relaxation:
            return float(base) * factor_f
        if role == "dt":
            return float(self.prepare_substep())
        if role == "dU":
            return 0.0
        return float(base) * factor_f

    def payload(self) -> dict[str, Any]:
        return {
            "schema": "coupled_cohesive_mpz_event_relaxation_v915_v1",
            "event_relaxation_time_s": float(self.tau_event_s),
            "event_relaxation_substeps": int(self.n_substeps),
            "event_initial_damage": float(self.initial_damage),
            "active_event_id_at_exit": self.active_event_id,
            "geometry_retry_attempts": int(self.geometry_retry_attempts),
            "geometry_retry_successes": int(self.geometry_retry_successes),
            "events": self.events,
        }


def _finite_or_none(value: Any) -> float | None:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if math.isfinite(out) else None


class _ControlledScalar(float):
    def __new__(cls, value: float, controller: CoupledEventRelaxationController, role: str):
        obj = float.__new__(cls, float(value))
        obj.controller = controller
        obj.role = role
        return obj

    def __mul__(self, other):
        return self.controller.controlled_value(float(self), self.role, other)

    def __rmul__(self, other):
        return self.controller.controlled_value(float(self), self.role, other)


class _LoadingProxy:
    def __init__(self, base: Any, controller: CoupledEventRelaxationController):
        object.__setattr__(self, "_base", base)
        object.__setattr__(self, "_controller", controller)

    def __getattr__(self, name: str):
        value = getattr(self._base, name)
        if name == "dt":
            return _ControlledScalar(value, self._controller, "dt")
        if name == "dU_top":
            return _ControlledScalar(value, self._controller, "dU")
        return value

    def __setattr__(self, name: str, value: Any) -> None:
        if name in {"dt", "dU_top"}:
            setattr(self._base, name, float(value))
        else:
            setattr(self._base, name, value)


def _forward_mode_i_plane(theta_deg: float, *args: Any, **kwargs: Any):
    del theta_deg, args, kwargs
    return [{
        "name": "v915_forward_mode_I",
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
    _inject_once(user_args, "--crack-backend", "adaptive_czm")
    _inject_once(user_args, "--adaptive-events")
    _inject_once(user_args, "--min-global-forward", "0.999999")
    # Give exact-ray marching more room before declaring a geometry veto.  The
    # physical 5 um renewal is still atomic; retries must realize its full length.
    _inject_once(user_args, "--czm-max-hrefine-subsegments", "4096")

    controller = CoupledEventRelaxationController()
    original_planes = _crystal.bcc_cleavage_traces
    original_config = _sharp_front.make_emergent_config
    original_backend_advance = _crack_backend.AdaptiveCZMBackend.advance
    original_predict = MovingProcessZone2DFrontEngine.predict_clock_increment_drives
    original_step_drives = MovingProcessZone2DFrontEngine.step_drives
    original_geometry_veto = MovingProcessZone2DFrontEngine.restore_geometry_veto

    def configured_model():
        cfg = original_config()
        cfg.loading = _LoadingProxy(cfg.loading, controller)
        return cfg

    def predict_with_relaxation(self, K_cleave, K_emit, T, dt):
        if controller.ready_for_relaxation:
            return 0.0
        return original_predict(self, K_cleave, K_emit, T, dt)

    def step_with_relaxation(self, K_cleave, K_emit, T, dt, metadata=None):
        if controller.ready_for_relaxation:
            saved_B = float(self.B)
            saved_nu0 = float(self.f.nu0_c)
            saved_lam = getattr(self, "_lambda_c_prev", None)
            saved_Kc = getattr(self, "_K_cleave_prev", None)
            self.B = 0.0
            self.f.nu0_c = 0.0
            try:
                out = original_step_drives(
                    self, K_cleave, K_emit, T, dt, metadata=metadata
                )
            finally:
                self.f.nu0_c = saved_nu0
            self.B = saved_B
            self._lambda_c_prev = saved_lam
            self._K_cleave_prev = saved_Kc
            out.update({
                "fired": False,
                "n_fire": 0,
                "n_fire_available": 0,
                "v_crack": 0.0,
                "B": float(self.B),
                "lambda_c": 0.0,
                "lambda_c_raw": 0.0,
                "v915_event_relaxation_active": 1.0,
                "v915_event_relaxation_id": float(controller.active_event_id or -1),
                "v915_event_relaxation_substep": float(controller.substep_index),
            })
            controller.finish_substep(out)
            return out

        out = original_step_drives(self, K_cleave, K_emit, T, dt, metadata=metadata)
        if int(out.get("n_fire", 0) or 0) > 0:
            out["v915_physical_event_id"] = controller.schedule_event(out)
        out["v915_event_relaxation_active"] = 0.0
        return out

    def veto_and_cancel(self, n_restore: int):
        result = original_geometry_veto(self, n_restore)
        controller.cancel_event("geometry_veto")
        return result

    def advance_and_register(self, *args, **kwargs):
        before_elems = len(getattr(self, "cohesive_network").elements)
        before_log = len(getattr(self, "advance_log", []))
        result = original_backend_advance(self, *args, **kwargs)
        retry_level = 0
        if not bool(getattr(result, "inserted", False)) and controller.active:
            # Retry the same complete physical increment with progressively less
            # restrictive local quality thresholds.  No partial physical advance is
            # accepted and the Arrhenius renewal is not redrawn.
            original_values = (
                self.min_area_ratio,
                self.min_triangle_quality,
                self.max_node_move_factor,
                self.max_hrefine_subsegments,
            )
            for retry_level in range(1, 4):
                controller.geometry_retry_attempts += 1
                self.min_area_ratio = max(original_values[0] * (0.5 ** retry_level), 0.01)
                self.min_triangle_quality = max(original_values[1] * (0.5 ** retry_level), 0.005)
                self.max_node_move_factor = original_values[2] * (1.0 + 0.5 * retry_level)
                self.max_hrefine_subsegments = max(
                    original_values[3], 4096 * (2 ** (retry_level - 1))
                )
                result = original_backend_advance(self, *args, **kwargs)
                if bool(getattr(result, "inserted", False)):
                    controller.geometry_retry_successes += 1
                    break
            (
                self.min_area_ratio,
                self.min_triangle_quality,
                self.max_node_move_factor,
                self.max_hrefine_subsegments,
            ) = original_values

        if bool(getattr(result, "inserted", False)) and controller.active:
            new_elements = self.cohesive_network.elements[before_elems:]
            new_rows = self.advance_log[before_log:]
            controller.register_geometry(new_elements, new_rows, retry_level=retry_level)
        return result

    _crystal.bcc_cleavage_traces = _forward_mode_i_plane
    _sharp_front.make_emergent_config = configured_model
    MovingProcessZone2DFrontEngine.predict_clock_increment_drives = predict_with_relaxation
    MovingProcessZone2DFrontEngine.step_drives = step_with_relaxation
    MovingProcessZone2DFrontEngine.restore_geometry_veto = veto_and_cancel
    _crack_backend.AdaptiveCZMBackend.advance = advance_and_register
    try:
        results = _base.main(user_args)
    finally:
        _crystal.bcc_cleavage_traces = original_planes
        _sharp_front.make_emergent_config = original_config
        MovingProcessZone2DFrontEngine.predict_clock_increment_drives = original_predict
        MovingProcessZone2DFrontEngine.step_drives = original_step_drives
        MovingProcessZone2DFrontEngine.restore_geometry_veto = original_geometry_veto
        _crack_backend.AdaptiveCZMBackend.advance = original_backend_advance

    out_value = _option_value(user_args, "--out")
    if out_value is not None:
        out = Path(out_value)
        out.mkdir(parents=True, exist_ok=True)
        (out / "coupled_event_relaxation_v915.json").write_text(
            json.dumps(controller.payload(), indent=2, default=str)
        )
    return results


if __name__ == "__main__":
    main()
