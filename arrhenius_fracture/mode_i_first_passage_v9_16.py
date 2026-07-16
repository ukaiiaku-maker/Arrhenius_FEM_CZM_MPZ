"""v9.16 kinetic trial-cohesive opening with deferred MPZ renewal.

A completed Arrhenius cleavage clock creates a *trial* adaptive-CZM segment.
The segment then opens at fixed remote displacement with a progress rate tied to
its instantaneous cleavage hazard.  Emission/transport evolve during the same
physical time.  The moving-process-zone renewal is committed only after the
cohesive event reaches unit progress.  A sufficiently suppressed event is
paused at partial damage and may resume after later loading.

The Arrhenius clock remains the sole nucleation criterion.  No critical
traction, critical opening, Griffith threshold, or empirical Gc criterion is
introduced.
"""
from __future__ import annotations

import json
import math
import os
from pathlib import Path
import sys
from typing import Any

import numpy as np

from . import cohesive as _cohesive
from . import crack_backend as _crack_backend
from . import crystal as _crystal
from . import mode_i_first_passage_v9_13 as _base
from . import sharp_front as _sharp_front
from .mode_i_first_passage_v9_15 import (
    _LoadingProxy,
    _finite_or_none,
    _forward_mode_i_plane,
    _inject_once,
    _option_value,
)
from .mpz_front_engine_v911 import MovingProcessZone2DFrontEngine


class KineticTrialEventController:
    """Resolve one nucleated fracture quantum as a kinetic trial interface."""

    def __init__(self) -> None:
        self.tau_event_s = max(
            float(os.environ.get("ARRHENIUS_EVENT_RELAXATION_TIME_S", "1e-6")),
            1.0e-15,
        )
        self.substeps_per_tau = max(
            int(os.environ.get("ARRHENIUS_EVENT_RELAXATION_SUBSTEPS", "24")), 2
        )
        self.initial_damage = float(np.clip(
            float(os.environ.get("ARRHENIUS_EVENT_INITIAL_DAMAGE", "0")),
            0.0,
            0.999999,
        ))
        self.rate_exponent = max(
            float(os.environ.get("ARRHENIUS_EVENT_RATE_EXPONENT", "1")), 0.0
        )
        self.min_rate_ratio = max(
            float(os.environ.get("ARRHENIUS_EVENT_MIN_RATE_RATIO", "1e-3")), 0.0
        )
        self.resume_rate_ratio = max(
            float(os.environ.get("ARRHENIUS_EVENT_RESUME_RATE_RATIO", "2e-3")),
            self.min_rate_ratio,
        )
        self.arrest_substeps = max(
            int(os.environ.get("ARRHENIUS_EVENT_ARREST_SUBSTEPS", "4")), 1
        )
        self.max_time_multiplier = max(
            float(os.environ.get("ARRHENIUS_EVENT_MAX_TIME_MULTIPLIER", "20")), 1.0
        )
        self.completion_tolerance = max(
            float(os.environ.get("ARRHENIUS_EVENT_COMPLETION_TOL", "1e-10")), 0.0
        )

        self.next_event_id = 1
        self.active_event_id: int | None = None
        self.active_elements: list[Any] = []
        self.active_log_rows: list[dict[str, Any]] = []
        self.active_engine: Any | None = None
        self.pending_nfire = 0
        self.pending_distance_m = 0.0
        self.pending_wake_preview: dict[str, float] = {}
        self.progress = self.initial_damage
        self.substep_index = 0
        self.physical_time_s = 0.0
        self.low_rate_count = 0
        self.prepared = False
        self.paused = False
        self.events: list[dict[str, Any]] = []
        self._active_record: dict[str, Any] | None = None
        self._last_W_emit: float | None = None
        self.geometry_retry_attempts = 0
        self.geometry_retry_successes = 0
        self.total_committed_distance_m = 0.0
        self.total_arrests = 0
        self.total_resumes = 0

    @property
    def active(self) -> bool:
        return self.active_event_id is not None

    @property
    def ready_for_relaxation(self) -> bool:
        return self.active and bool(self.active_elements) and not self.paused

    @property
    def dt_substep_s(self) -> float:
        return self.tau_event_s / float(self.substeps_per_tau)

    def defer_engine_renewal(
        self,
        engine: Any,
        nfire: int,
        distance_m: float,
        wake_preview: dict[str, float] | None = None,
    ) -> None:
        """Register a consumed cleavage renewal whose MPZ translation is deferred."""
        if self.active:
            raise RuntimeError("cannot defer a second renewal while a trial event is active")
        self.active_engine = engine
        self.pending_nfire = max(int(nfire), 0)
        self.pending_distance_m = max(float(distance_m), 0.0)
        self.pending_wake_preview = dict(wake_preview or {})

    def schedule_event(self, out: dict[str, Any]) -> int:
        if self.active_event_id is not None:
            return int(self.active_event_id)
        event_id = int(self.next_event_id)
        self.next_event_id += 1
        self.active_event_id = event_id
        self.active_elements = []
        self.active_log_rows = []
        self.progress = self.initial_damage
        self.substep_index = 0
        self.physical_time_s = 0.0
        self.low_rate_count = 0
        self.prepared = False
        self.paused = False
        lam_nuc = max(float(out.get("lambda_c", 0.0) or 0.0), 0.0)
        lam_raw_nuc = max(float(out.get("lambda_c_raw", 0.0) or 0.0), 0.0)
        self._last_W_emit = _finite_or_none(out.get("W_emit"))
        rec = {
            "event_id": event_id,
            "status": "nucleated_trial_geometry_pending",
            "tau_event_s": float(self.tau_event_s),
            "substeps_per_tau": int(self.substeps_per_tau),
            "max_time_multiplier": float(self.max_time_multiplier),
            "initial_damage": float(self.initial_damage),
            "progress_law": "dq_dt=(lambda_c/lambda_c_nucleation)^p/tau_event",
            "progress_rate_exponent": float(self.rate_exponent),
            "min_rate_ratio": float(self.min_rate_ratio),
            "resume_rate_ratio": float(self.resume_rate_ratio),
            "arrest_substeps": int(self.arrest_substeps),
            "lambda_c_nucleation_s-1": float(lam_nuc),
            "lambda_c_raw_nucleation_s-1": float(lam_raw_nuc),
            "KJ_nucleation_Pa_sqrt_m": _finite_or_none(
                out.get("anisotropic_KJ_Pa_sqrt_m")
            ),
            "Kshield_nucleation_Pa_sqrt_m": _finite_or_none(
                out.get("mpz_K_shield_pre_renewal_Pa_sqrt_m")
            ),
            "emitted_total_nucleation": _finite_or_none(
                out.get("mpz_emitted_total_pre_renewal")
            ),
            "trial_geometry_inserted": False,
            "mpz_renewal_deferred": bool(self.active_engine is not None),
            "mpz_renewal_committed": False,
            "pending_nfire": int(self.pending_nfire),
            "pending_distance_m": float(self.pending_distance_m),
            "substeps": [],
            "reload_probes": [],
            "geometry_retry_level": 0,
            "cohesive_work_J_per_m": 0.0,
            "tip_emission_work_J_per_m": 0.0,
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
            elem.metadata["v916_physical_event_id"] = int(self.active_event_id)
            elem.metadata["v916_trial_interface"] = True
            elem.metadata["v916_committed"] = False
        for i, row in enumerate(self.active_log_rows):
            row["physical_event_id"] = int(self.active_event_id)
            row["physical_subsegment_index"] = int(i)
            row["physical_subsegment_count"] = int(len(self.active_log_rows))
            row["v916_initial_damage"] = float(self.initial_damage)
            row["v916_geometry_retry_level"] = int(retry_level)
            row["v916_trial_interface"] = True
        if self._active_record is not None:
            self._active_record.update({
                "status": "kinetic_relaxation_pending",
                "trial_geometry_inserted": True,
                "n_cohesive_subsegments": int(len(self.active_elements)),
                "cohesive_length_m": float(sum(
                    max(float(getattr(e, "length", 0.0)), 0.0)
                    for e in self.active_elements
                )),
                "geometry_retry_level": int(retry_level),
            })

    def _clear_active(self) -> None:
        self.active_event_id = None
        self.active_elements = []
        self.active_log_rows = []
        self.active_engine = None
        self.pending_nfire = 0
        self.pending_distance_m = 0.0
        self.pending_wake_preview = {}
        self.progress = self.initial_damage
        self.substep_index = 0
        self.physical_time_s = 0.0
        self.low_rate_count = 0
        self.prepared = False
        self.paused = False
        self._active_record = None
        self._last_W_emit = None

    def cancel_event(self, reason: str = "geometry_veto") -> None:
        if self._active_record is not None:
            self._active_record.update({
                "status": "cancelled",
                "cancel_reason": str(reason),
                "relaxation_completed": False,
                "mpz_renewal_committed": False,
            })
        self._clear_active()

    def prepare_substep(self) -> float:
        if not self.ready_for_relaxation:
            return 0.0
        if self.prepared:
            return self.dt_substep_s
        self.substep_index += 1
        for elem in self.active_elements:
            elem.damage = float(np.clip(self.progress, 0.0, 1.0))
            elem.clock = float(self.physical_time_s)
        self.prepared = True
        return self.dt_substep_s

    def current_rate_ratio(self, lambda_c: float) -> float:
        if self._active_record is None:
            return 0.0
        ref = max(float(self._active_record.get("lambda_c_nucleation_s-1", 0.0) or 0.0), 1.0e-300)
        return max(float(lambda_c), 0.0) / ref

    def _cohesive_energy_state(self) -> tuple[float, float, float]:
        recoverable = 0.0
        opening_max = 0.0
        traction_max = 0.0
        for elem in self.active_elements:
            md = getattr(elem, "metadata", {})
            recoverable += max(float(md.get("v916_recoverable_energy_J_per_m", 0.0) or 0.0), 0.0)
            opening_max = max(opening_max, abs(float(md.get("v916_opening_jump_m", 0.0) or 0.0)))
            traction_max = max(traction_max, abs(float(md.get("v916_normal_traction_Pa", 0.0) or 0.0)))
        return recoverable, opening_max, traction_max

    def _cohesive_damage_work(self, dq: float) -> float:
        if dq <= 0.0:
            return 0.0
        work = 0.0
        for elem in self.active_elements:
            md = getattr(elem, "metadata", {})
            intact_energy = max(float(md.get("v916_intact_reference_energy_J_per_m", 0.0) or 0.0), 0.0)
            work += dq * intact_energy
        return float(work)

    def _commit_deferred_renewal(self) -> dict[str, float]:
        wake = {
            "wake_mobile": 0.0,
            "wake_retained": 0.0,
            "wake_slip": 0.0,
            "source_sites_refreshed": 0.0,
        }
        eng = self.active_engine
        if eng is not None and self.pending_nfire > 0 and self.pending_distance_m > 0.0:
            wake = eng.mpz_state.advance(self.pending_distance_m)
            eng.a_adv += self.pending_distance_m
            eng.n_adv += self.pending_nfire
            eng._sync_compat()
        self.total_committed_distance_m += self.pending_distance_m
        for elem in self.active_elements:
            elem.damage = 1.0
            elem.clock = float(self.physical_time_s)
            elem.metadata["v916_trial_interface"] = False
            elem.metadata["v916_committed"] = True
        for row in self.active_log_rows:
            row["v916_trial_interface"] = False
            row["v916_committed"] = True
        return {k: float(v) for k, v in wake.items()}

    def finish_substep(
        self,
        out: dict[str, Any],
        *,
        lambda_c_current: float,
        lambda_c_raw_current: float,
        K_cleave: float,
        T: float,
    ) -> None:
        if not self.ready_for_relaxation or not self.prepared:
            return
        dt = self.dt_substep_s
        q0 = float(self.progress)
        ratio = self.current_rate_ratio(lambda_c_current)
        rate_factor = ratio ** self.rate_exponent if self.rate_exponent > 0.0 else 1.0
        dq = max(dt / self.tau_event_s * rate_factor, 0.0)
        q1 = min(q0 + dq, 1.0)
        dq_used = q1 - q0
        self.progress = q1
        self.physical_time_s += dt

        if ratio < self.min_rate_ratio:
            self.low_rate_count += 1
        else:
            self.low_rate_count = 0

        cohesive_recoverable, opening_max, traction_max = self._cohesive_energy_state()
        dWcoh = self._cohesive_damage_work(dq_used)
        W_emit = _finite_or_none(out.get("W_emit"))
        dWemit = 0.0
        if W_emit is not None and self._last_W_emit is not None:
            dWemit = max(W_emit - self._last_W_emit, 0.0)
        if W_emit is not None:
            self._last_W_emit = W_emit

        row = {
            "substep": int(self.substep_index),
            "time_s": float(self.physical_time_s),
            "dt_s": float(dt),
            "damage_progress_before": q0,
            "damage_progress": q1,
            "damage_increment": dq_used,
            "remote_displacement_increment_m": 0.0,
            "KJ_Pa_sqrt_m": _finite_or_none(out.get("anisotropic_KJ_Pa_sqrt_m")),
            "K_cleave_drive_Pa_sqrt_m": float(K_cleave),
            "Kshield_Pa_sqrt_m": _finite_or_none(out.get("mpz_K_shield_pre_renewal_Pa_sqrt_m")),
            "sigma_emit_tip_Pa": _finite_or_none(out.get("sigma_emit_tip")),
            "lambda_c_current_s-1": float(max(lambda_c_current, 0.0)),
            "lambda_c_raw_current_s-1": float(max(lambda_c_raw_current, 0.0)),
            "lambda_c_ratio_to_nucleation": float(ratio),
            "temperature_K": float(T),
            "dN_emit": _finite_or_none(out.get("dN_emit")),
            "dN_emit_raw": _finite_or_none(out.get("dN_emit_raw")),
            "mobile_count": _finite_or_none(out.get("mpz_mobile_count")),
            "retained_count": _finite_or_none(out.get("mpz_retained_count")),
            "emitted_total": _finite_or_none(out.get("mpz_emitted_total")),
            "cohesive_recoverable_energy_J_per_m": cohesive_recoverable,
            "cohesive_damage_work_increment_J_per_m": dWcoh,
            "tip_emission_work_increment_J_per_m": dWemit,
            "opening_jump_max_m": opening_max,
            "normal_traction_max_Pa": traction_max,
        }
        if self._active_record is not None:
            self._active_record["substeps"].append(row)
            self._active_record["cohesive_work_J_per_m"] = float(
                self._active_record.get("cohesive_work_J_per_m", 0.0) + dWcoh
            )
            self._active_record["tip_emission_work_J_per_m"] = float(
                self._active_record.get("tip_emission_work_J_per_m", 0.0) + dWemit
            )
        self.prepared = False

        completed = self.progress >= 1.0 - self.completion_tolerance
        max_time = self.physical_time_s >= self.tau_event_s * self.max_time_multiplier
        arrested = self.low_rate_count >= self.arrest_substeps or max_time
        if completed:
            self.progress = 1.0
            wake = self._commit_deferred_renewal()
            if self._active_record is not None:
                sub = self._active_record.get("substeps", [])
                emitted = [x.get("dN_emit") for x in sub if x.get("dN_emit") is not None]
                Kvals = [x.get("KJ_Pa_sqrt_m") for x in sub if x.get("KJ_Pa_sqrt_m") is not None]
                ratios = [x.get("lambda_c_ratio_to_nucleation") for x in sub]
                self._active_record.update({
                    "status": "complete_committed",
                    "relaxation_completed": True,
                    "n_relaxation_substeps_completed": int(self.substep_index),
                    "final_damage": 1.0,
                    "relaxation_time_s": float(self.physical_time_s),
                    "dN_emit_relaxation": float(sum(emitted)) if emitted else 0.0,
                    "KJ_relaxation_start_Pa_sqrt_m": float(Kvals[0]) if Kvals else None,
                    "KJ_relaxation_end_Pa_sqrt_m": float(Kvals[-1]) if Kvals else None,
                    "min_lambda_c_ratio": float(min(ratios)) if ratios else None,
                    "max_lambda_c_ratio": float(max(ratios)) if ratios else None,
                    "mpz_renewal_committed": True,
                    "committed_distance_m": float(self.pending_distance_m),
                    "committed_nfire": int(self.pending_nfire),
                    "wake_on_commit": wake,
                })
            self._clear_active()
        elif arrested:
            self.paused = True
            self.total_arrests += 1
            if self._active_record is not None:
                self._active_record.update({
                    "status": "arrested_pending_reload",
                    "relaxation_completed": False,
                    "arrested": True,
                    "arrest_reason": "max_event_time" if max_time else "cleavage_rate_suppressed",
                    "arrest_damage": float(self.progress),
                    "arrest_time_s": float(self.physical_time_s),
                    "mpz_renewal_committed": False,
                })

    def note_reload_probe(
        self,
        *,
        lambda_c_current: float,
        K_cleave: float,
        KJ: Any,
        T: float,
    ) -> bool:
        if not self.active or not self.paused:
            return False
        ratio = self.current_rate_ratio(lambda_c_current)
        probe = {
            "K_cleave_drive_Pa_sqrt_m": float(K_cleave),
            "KJ_Pa_sqrt_m": _finite_or_none(KJ),
            "temperature_K": float(T),
            "lambda_c_current_s-1": float(max(lambda_c_current, 0.0)),
            "lambda_c_ratio_to_nucleation": float(ratio),
            "damage_progress": float(self.progress),
        }
        if self._active_record is not None:
            self._active_record.setdefault("reload_probes", []).append(probe)
        if ratio >= self.resume_rate_ratio:
            self.paused = False
            self.low_rate_count = 0
            self.total_resumes += 1
            if self._active_record is not None:
                self._active_record.update({
                    "status": "kinetic_relaxation_resumed",
                    "resume_count": int(self._active_record.get("resume_count", 0) + 1),
                })
            return True
        return False

    def controlled_value(self, base: float, role: str, factor: Any) -> float:
        factor_f = float(factor)
        if role == "dt" and self.ready_for_relaxation:
            return float(self.prepare_substep())
        if role == "dU" and self.ready_for_relaxation:
            return 0.0
        return float(base) * factor_f

    def payload(self) -> dict[str, Any]:
        return {
            "schema": "kinetic_trial_cohesive_mpz_event_v916_v1",
            "event_relaxation_time_s": float(self.tau_event_s),
            "event_substeps_per_tau": int(self.substeps_per_tau),
            "event_rate_exponent": float(self.rate_exponent),
            "event_min_rate_ratio": float(self.min_rate_ratio),
            "event_resume_rate_ratio": float(self.resume_rate_ratio),
            "event_arrest_substeps": int(self.arrest_substeps),
            "event_max_time_multiplier": float(self.max_time_multiplier),
            "active_event_id_at_exit": self.active_event_id,
            "active_event_paused_at_exit": bool(self.paused),
            "active_event_progress_at_exit": float(self.progress) if self.active else None,
            "geometry_retry_attempts": int(self.geometry_retry_attempts),
            "geometry_retry_successes": int(self.geometry_retry_successes),
            "total_committed_distance_m": float(self.total_committed_distance_m),
            "total_arrests": int(self.total_arrests),
            "total_resumes": int(self.total_resumes),
            "events": self.events,
        }


def _install_cohesive_state_audit():
    """Wrap cohesive assembly so each trial element exposes jump/energy state."""
    original = _cohesive.cohesive_contribution

    def wrapped(network, u, ndof):
        K, R = original(network, u, ndof)
        if network is None:
            return K, R
        kn = float(getattr(network, "penalty_normal_Pa_per_m", 0.0))
        kt = float(getattr(network, "penalty_tangent_Pa_per_m", 0.0))
        for elem in getattr(network, "elements", []):
            try:
                p0, p1 = elem.plus_nodes
                m0, m1 = elem.minus_nodes
                up = 0.5 * (u[2 * p0:2 * p0 + 2] + u[2 * p1:2 * p1 + 2])
                um = 0.5 * (u[2 * m0:2 * m0 + 2] + u[2 * m1:2 * m1 + 2])
                jump = np.asarray(up - um, dtype=float)
                dn = float(np.asarray(elem.normal) @ jump)
                dtan = float(np.asarray(elem.tangent) @ jump)
                intact = max(1.0 - float(elem.damage), 0.0)
                dn_tension = max(dn, 0.0)
                L = max(float(elem.length), 0.0)
                intact_ref = 0.5 * (kn * dn_tension * dn_tension + kt * dtan * dtan) * L
                recoverable = intact * intact_ref
                elem.metadata.update({
                    "v916_opening_jump_m": dn,
                    "v916_sliding_jump_m": dtan,
                    "v916_normal_traction_Pa": kn * intact * dn_tension,
                    "v916_tangential_traction_Pa": kt * intact * dtan,
                    "v916_intact_reference_energy_J_per_m": intact_ref,
                    "v916_recoverable_energy_J_per_m": recoverable,
                })
            except Exception:
                continue
        return K, R

    _cohesive.cohesive_contribution = wrapped
    return original


def main(argv=None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    _inject_once(user_args, "--crystal-aniso")
    _inject_once(user_args, "--crack-backend", "adaptive_czm")
    _inject_once(user_args, "--adaptive-events")
    _inject_once(user_args, "--min-global-forward", "0.999999")
    _inject_once(user_args, "--czm-max-hrefine-subsegments", "4096")

    controller = KineticTrialEventController()
    original_planes = _crystal.bcc_cleavage_traces
    original_config = _sharp_front.make_emergent_config
    original_backend_advance = _crack_backend.AdaptiveCZMBackend.advance
    original_predict = MovingProcessZone2DFrontEngine.predict_clock_increment_drives
    original_step_drives = MovingProcessZone2DFrontEngine.step_drives
    original_renew = MovingProcessZone2DFrontEngine._renew
    original_geometry_veto = MovingProcessZone2DFrontEngine.restore_geometry_veto
    original_cohesive = _install_cohesive_state_audit()

    def configured_model():
        cfg = original_config()
        cfg.loading = _LoadingProxy(cfg.loading, controller)
        return cfg

    def predict_with_trial(self, K_cleave, K_emit, T, dt):
        if controller.active:
            return 0.0
        return original_predict(self, K_cleave, K_emit, T, dt)

    def renew_deferred(self, dt):
        pre_a = float(self.a_adv)
        pre_n = int(self.n_adv)
        out = original_renew(self, dt)
        nfire = int(out.get("n_fire", 0) or 0)
        if nfire > 0:
            distance = float(self.f.da) * nfire
            wake_preview = {
                "wake_mobile": float(out.get("mpz_wake_mobile_block", 0.0) or 0.0),
                "wake_retained": float(out.get("mpz_wake_retained_block", 0.0) or 0.0),
                "wake_slip": float(out.get("mpz_wake_slip_block", 0.0) or 0.0),
                "source_sites_refreshed": float(out.get("mpz_source_sites_refreshed_on_advance", 0.0) or 0.0),
            }
            if getattr(self, "_last_pre_renewal_state", None) is not None:
                self.mpz_state = self._last_pre_renewal_state.copy()
            self.a_adv = pre_a
            self.n_adv = pre_n
            self._sync_compat()
            controller.defer_engine_renewal(self, nfire, distance, wake_preview)
            out.update({
                "v_crack": 0.0,
                "N_em_retained": float(self.N_em),
                "N_em_shed_to_wake": 0.0,
                "mpz_wake_mobile_block": 0.0,
                "mpz_wake_retained_block": 0.0,
                "mpz_wake_slip_block": 0.0,
                "mpz_source_sites_refreshed_on_advance": 0.0,
                "v916_mpz_renewal_deferred": 1.0,
                "v916_pending_advance_m": distance,
            })
        else:
            out["v916_mpz_renewal_deferred"] = 0.0
        return out

    def emission_only_step(self, K_cleave, K_emit, T, dt, metadata=None):
        saved_B = float(self.B)
        saved_nu0 = float(self.f.nu0_c)
        saved_lam = getattr(self, "_lambda_c_prev", None)
        saved_Kc = getattr(self, "_K_cleave_prev", None)
        self.B = 0.0
        self.f.nu0_c = 0.0
        try:
            out = original_step_drives(self, K_cleave, K_emit, T, dt, metadata=metadata)
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
        })
        return out

    def step_with_trial(self, K_cleave, K_emit, T, dt, metadata=None):
        if controller.ready_for_relaxation:
            out = emission_only_step(self, K_cleave, K_emit, T, dt, metadata=metadata)
            sig_c = self.sigma_tip(K_cleave)
            lam_c, lam_raw, _ = self.lambda_cleave(sig_c, T)
            out.update({
                "v916_event_relaxation_active": 1.0,
                "v916_event_id": float(controller.active_event_id or -1),
                "v916_event_substep": float(controller.substep_index),
                "v916_event_progress": float(controller.progress),
                "v916_lambda_c_probe": float(lam_c),
            })
            controller.finish_substep(
                out,
                lambda_c_current=lam_c,
                lambda_c_raw_current=lam_raw,
                K_cleave=float(K_cleave),
                T=float(T),
            )
            return out

        if controller.active and controller.paused:
            out = emission_only_step(self, K_cleave, K_emit, T, dt, metadata=metadata)
            sig_c = self.sigma_tip(K_cleave)
            lam_c, lam_raw, _ = self.lambda_cleave(sig_c, T)
            controller.note_reload_probe(
                lambda_c_current=lam_c,
                K_cleave=float(K_cleave),
                KJ=out.get("anisotropic_KJ_Pa_sqrt_m"),
                T=float(T),
            )
            out.update({
                "v916_event_relaxation_active": 0.0,
                "v916_event_paused": 1.0,
                "v916_event_id": float(controller.active_event_id or -1),
                "v916_event_progress": float(controller.progress),
                "v916_lambda_c_probe": float(lam_c),
                "v916_lambda_c_raw_probe": float(lam_raw),
            })
            return out

        out = original_step_drives(self, K_cleave, K_emit, T, dt, metadata=metadata)
        if int(out.get("n_fire", 0) or 0) > 0:
            out["v916_physical_event_id"] = controller.schedule_event(out)
        out["v916_event_relaxation_active"] = 0.0
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
    MovingProcessZone2DFrontEngine.predict_clock_increment_drives = predict_with_trial
    MovingProcessZone2DFrontEngine.step_drives = step_with_trial
    MovingProcessZone2DFrontEngine._renew = renew_deferred
    MovingProcessZone2DFrontEngine.restore_geometry_veto = veto_and_cancel
    _crack_backend.AdaptiveCZMBackend.advance = advance_and_register
    try:
        results = _base.main(user_args)
    finally:
        _crystal.bcc_cleavage_traces = original_planes
        _sharp_front.make_emergent_config = original_config
        MovingProcessZone2DFrontEngine.predict_clock_increment_drives = original_predict
        MovingProcessZone2DFrontEngine.step_drives = original_step_drives
        MovingProcessZone2DFrontEngine._renew = original_renew
        MovingProcessZone2DFrontEngine.restore_geometry_veto = original_geometry_veto
        _crack_backend.AdaptiveCZMBackend.advance = original_backend_advance
        _cohesive.cohesive_contribution = original_cohesive

    out_value = _option_value(user_args, "--out")
    if out_value is not None:
        out = Path(out_value)
        out.mkdir(parents=True, exist_ok=True)
        (out / "kinetic_trial_event_relaxation_v916.json").write_text(
            json.dumps(controller.payload(), indent=2, default=str)
        )
    return results


if __name__ == "__main__":
    main()
