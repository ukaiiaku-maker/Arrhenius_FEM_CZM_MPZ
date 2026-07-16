"""v9.17 absolute-hazard trial opening with opening-proportional source refresh.

One calibrated cleavage renewal creates one physical crack quantum.  The trial
cohesive segment then consumes one additional unit of the same *absolute*
cleavage hazard,

    dq/dt = lambda_c(K, S, T),  q in [0, 1].

No prescribed event time, normalized rate, rate-ratio arrest threshold, or
independent cohesive fracture criterion is introduced.  Fresh source sites are
created incrementally with opened surface area, while moving-frame translation
is still deferred until q reaches one.
"""
from __future__ import annotations

import json
import math
import os
from pathlib import Path
import sys
from typing import Any

import numpy as np

from . import mode_i_first_passage_v9_16 as _v916


def _finite_or_none(value: Any) -> float | None:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if math.isfinite(out) else None


class HazardClockTrialEventController(_v916.KineticTrialEventController):
    """Advance cohesive progress by the calibrated absolute cleavage hazard."""

    def __init__(self) -> None:
        super().__init__()
        self.target_hazard_increment = float(np.clip(
            float(os.environ.get("ARRHENIUS_EVENT_TARGET_DQ", "0.05")),
            1.0e-6,
            0.5,
        ))
        self.min_event_dt_s = max(
            float(os.environ.get("ARRHENIUS_EVENT_MIN_DT_S", "1e-12")),
            1.0e-18,
        )
        raw_max = float(os.environ.get("ARRHENIUS_EVENT_MAX_FIXED_HOLD_S", "inf"))
        self.max_fixed_hold_s = raw_max if math.isfinite(raw_max) and raw_max > 0 else math.inf
        self.external_dt_s: float | None = None
        self._last_lambda_c = 0.0
        self._prepared_dt_s = 0.0
        self._prepared_progress = self.progress
        self._prepared_hazard_increment = 0.0
        self._prepared_source_refresh = 0.0
        self._prepared_hit_external_cap = False
        self._source_refresh_entitlement_applied: np.ndarray | None = None
        self._source_refresh_total = 0.0
        self._source_available_at_nucleation: np.ndarray | None = None
        self.loading_resume_requests = 0

    @property
    def dt_substep_s(self) -> float:
        return float(self._prepared_dt_s)

    def defer_engine_renewal(
        self,
        engine: Any,
        nfire: int,
        distance_m: float,
        wake_preview: dict[str, float] | None = None,
    ) -> None:
        if int(nfire) != 1:
            raise RuntimeError(
                "v9.17 requires exactly one renewal per trial-opening event; "
                f"received nfire={nfire}. Run with --max-advances-per-step 1."
            )
        super().defer_engine_renewal(engine, nfire, distance_m, wake_preview)

    def schedule_event(self, out: dict[str, Any]) -> int:
        event_id = super().schedule_event(out)
        self._last_lambda_c = max(float(out.get("lambda_c", 0.0) or 0.0), 0.0)
        self._prepared_dt_s = 0.0
        self._prepared_progress = self.progress
        self._prepared_hazard_increment = 0.0
        self._prepared_source_refresh = 0.0
        self._prepared_hit_external_cap = False
        self._source_refresh_total = 0.0
        state = getattr(self.active_engine, "mpz_state", None)
        if state is not None:
            self._source_available_at_nucleation = np.asarray(
                state.available_sites, dtype=float
            ).copy()
            self._source_refresh_entitlement_applied = np.zeros_like(
                self._source_available_at_nucleation
            )
        else:
            self._source_available_at_nucleation = None
            self._source_refresh_entitlement_applied = None

        if self._active_record is not None:
            for key in (
                "tau_event_s", "substeps_per_tau", "max_time_multiplier",
                "progress_rate_exponent", "min_rate_ratio", "resume_rate_ratio",
                "arrest_substeps",
            ):
                self._active_record.pop(key, None)
            lam = max(self._last_lambda_c, 1.0e-300)
            self._active_record.update({
                "status": "absolute_hazard_opening_pending",
                "progress_law": "dq_dt=lambda_c_absolute",
                "opening_hazard_target": 1.0,
                "target_hazard_increment": float(self.target_hazard_increment),
                "lambda_c_nucleation_s-1": float(self._last_lambda_c),
                "predicted_opening_time_at_nucleation_s": float(1.0 / lam),
                "source_refresh_law": (
                    "incremental_missing_inventory_times_"
                    "min(q*pending_distance/source_refresh_length,1)"
                ),
                "source_sites_available_at_nucleation": (
                    self._source_available_at_nucleation.tolist()
                    if self._source_available_at_nucleation is not None else None
                ),
                "source_sites_refreshed_during_opening": 0.0,
                "loading_resume_requests": 0,
            })
        return event_id

    def _incremental_source_refresh(self, q_target: float) -> float:
        state = getattr(self.active_engine, "mpz_state", None)
        baseline = self._source_available_at_nucleation
        applied = self._source_refresh_entitlement_applied
        if state is None or baseline is None or applied is None:
            return 0.0
        capacity = np.asarray(state.site_capacity, dtype=float)
        missing_at_nucleation = np.maximum(capacity - baseline, 0.0)
        Lrefresh = max(float(state.cfg.source_refresh_length_m), float(state.dx))
        opened_distance = max(float(q_target), 0.0) * max(self.pending_distance_m, 0.0)
        fraction = min(opened_distance / Lrefresh, 1.0)
        entitlement = missing_at_nucleation * fraction
        increment = np.maximum(entitlement - applied, 0.0)
        if np.any(increment > 0.0):
            state.available_sites = np.minimum(
                np.asarray(state.available_sites, dtype=float) + increment,
                capacity,
            )
            self._source_refresh_entitlement_applied = entitlement
            added = float(np.sum(increment))
            self._source_refresh_total += added
            if self.active_engine is not None:
                self.active_engine._sync_compat()
            return added
        self._source_refresh_entitlement_applied = np.maximum(applied, entitlement)
        return 0.0

    def _hold_cap_s(self) -> float:
        caps = [self.max_fixed_hold_s]
        if self.external_dt_s is not None and self.external_dt_s > 0.0:
            caps.append(float(self.external_dt_s))
        finite = [x for x in caps if math.isfinite(x) and x > 0.0]
        return min(finite) if finite else math.inf

    def prepare_substep(self) -> float:
        if not self.ready_for_relaxation:
            return 0.0
        if self.prepared:
            return float(self._prepared_dt_s)

        self.substep_index += 1
        remaining = max(1.0 - float(self.progress), 0.0)
        dq_target = min(self.target_hazard_increment, remaining)
        lam = max(float(self._last_lambda_c), 0.0)
        dt_to_target = dq_target / lam if lam > 0.0 else math.inf
        hold_cap = self._hold_cap_s()
        dt = min(dt_to_target, hold_cap)
        if not math.isfinite(dt) or dt <= 0.0:
            dt = hold_cap if math.isfinite(hold_cap) and hold_cap > 0.0 else self.min_event_dt_s
        dt = max(float(dt), self.min_event_dt_s)
        dq = min(lam * dt, dq_target, remaining)
        q_trial = min(float(self.progress) + dq, 1.0)

        self._prepared_dt_s = dt
        self._prepared_progress = q_trial
        self._prepared_hazard_increment = max(q_trial - float(self.progress), 0.0)
        self._prepared_hit_external_cap = bool(
            math.isfinite(hold_cap)
            and (not math.isfinite(dt_to_target) or dt_to_target > hold_cap * (1.0 + 1.0e-12))
        )
        self._prepared_source_refresh = self._incremental_source_refresh(q_trial)

        for elem in self.active_elements:
            elem.damage = float(q_trial)
            elem.clock = float(self.physical_time_s + dt)
        self.prepared = True
        return float(dt)

    def _commit_deferred_renewal(self) -> dict[str, float]:
        wake = {
            "wake_mobile": 0.0,
            "wake_retained": 0.0,
            "wake_slip": 0.0,
            "source_sites_refreshed": 0.0,
            "source_sites_refreshed_during_opening": float(self._source_refresh_total),
        }
        eng = self.active_engine
        if eng is not None and self.pending_nfire == 1 and self.pending_distance_m > 0.0:
            available_before = np.asarray(eng.mpz_state.available_sites, dtype=float).copy()
            translated = eng.mpz_state.advance(self.pending_distance_m)
            eng.mpz_state.available_sites = np.minimum(
                available_before, np.asarray(eng.mpz_state.site_capacity, dtype=float)
            )
            translated["source_sites_refreshed"] = 0.0
            wake.update({k: float(v) for k, v in translated.items()})
            wake["source_sites_refreshed_during_opening"] = float(self._source_refresh_total)
            eng.a_adv += self.pending_distance_m
            eng.n_adv += 1
            eng._sync_compat()
        self.total_committed_distance_m += self.pending_distance_m
        for elem in self.active_elements:
            elem.damage = 1.0
            elem.clock = float(self.physical_time_s)
            elem.metadata["v916_trial_interface"] = False
            elem.metadata["v916_committed"] = True
            elem.metadata["v917_absolute_hazard_opening"] = True
        for row in self.active_log_rows:
            row["v916_trial_interface"] = False
            row["v916_committed"] = True
            row["v917_absolute_hazard_opening"] = True
        return wake

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
        dt = float(self._prepared_dt_s)
        q0 = float(self.progress)
        q1 = float(self._prepared_progress)
        dq = max(q1 - q0, 0.0)
        self.progress = q1
        self.physical_time_s += dt

        cohesive_recoverable, opening_max, traction_max = self._cohesive_energy_state()
        dWcoh = self._cohesive_damage_work(dq)
        W_emit = _finite_or_none(out.get("W_emit"))
        dWemit = 0.0
        if W_emit is not None and self._last_W_emit is not None:
            dWemit = max(W_emit - self._last_W_emit, 0.0)
        if W_emit is not None:
            self._last_W_emit = W_emit

        lam_now = max(float(lambda_c_current), 0.0)
        remaining = max(1.0 - q1, 0.0)
        predicted_remaining = remaining / lam_now if lam_now > 0.0 else math.inf
        row = {
            "substep": int(self.substep_index),
            "time_s": float(self.physical_time_s),
            "dt_s": dt,
            "damage_progress_before": q0,
            "damage_progress": q1,
            "damage_increment": dq,
            "cleavage_hazard_increment": dq,
            "cleavage_hazard_integral": q1,
            "remote_displacement_increment_m": 0.0,
            "KJ_Pa_sqrt_m": _finite_or_none(out.get("anisotropic_KJ_Pa_sqrt_m")),
            "K_cleave_drive_Pa_sqrt_m": float(K_cleave),
            "Kshield_Pa_sqrt_m": _finite_or_none(out.get("mpz_K_shield_pre_renewal_Pa_sqrt_m")),
            "sigma_emit_tip_Pa": _finite_or_none(out.get("sigma_emit_tip")),
            "lambda_e_current_s-1": _finite_or_none(out.get("lambda_e")),
            "lambda_c_previous_s-1": float(self._last_lambda_c),
            "lambda_c_current_s-1": lam_now,
            "lambda_c_raw_current_s-1": float(max(lambda_c_raw_current, 0.0)),
            "lambda_c_ratio_to_nucleation": float(self.current_rate_ratio(lam_now)),
            "predicted_remaining_opening_time_s": float(predicted_remaining),
            "temperature_K": float(T),
            "dN_emit": _finite_or_none(out.get("dN_emit")),
            "dN_emit_raw": _finite_or_none(out.get("dN_emit_raw")),
            "mobile_count": _finite_or_none(out.get("mpz_mobile_count")),
            "retained_count": _finite_or_none(out.get("mpz_retained_count")),
            "emitted_total": _finite_or_none(out.get("mpz_emitted_total")),
            "available_site_fraction": _finite_or_none(out.get("mpz_available_site_fraction")),
            "source_sites_refreshed_before_substep": float(self._prepared_source_refresh),
            "source_sites_refreshed_cumulative": float(self._source_refresh_total),
            "cohesive_recoverable_energy_J_per_m": cohesive_recoverable,
            "cohesive_damage_work_increment_J_per_m": dWcoh,
            "tip_emission_work_increment_J_per_m": dWemit,
            "opening_jump_max_m": opening_max,
            "normal_traction_max_Pa": traction_max,
            "fixed_hold_capped_by_external_dt": bool(self._prepared_hit_external_cap),
        }
        if self._active_record is not None:
            self._active_record["substeps"].append(row)
            self._active_record["cohesive_work_J_per_m"] = float(
                self._active_record.get("cohesive_work_J_per_m", 0.0) + dWcoh
            )
            self._active_record["tip_emission_work_J_per_m"] = float(
                self._active_record.get("tip_emission_work_J_per_m", 0.0) + dWemit
            )
            self._active_record["source_sites_refreshed_during_opening"] = float(
                self._source_refresh_total
            )

        self._last_lambda_c = lam_now
        self.prepared = False
        completed = self.progress >= 1.0 - self.completion_tolerance
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
                    "opening_hazard_integral": 1.0,
                    "relaxation_time_s": float(self.physical_time_s),
                    "dN_emit_relaxation": float(sum(emitted)) if emitted else 0.0,
                    "KJ_relaxation_start_Pa_sqrt_m": float(Kvals[0]) if Kvals else None,
                    "KJ_relaxation_end_Pa_sqrt_m": float(Kvals[-1]) if Kvals else None,
                    "min_lambda_c_ratio": float(min(ratios)) if ratios else None,
                    "max_lambda_c_ratio": float(max(ratios)) if ratios else None,
                    "mpz_renewal_committed": True,
                    "committed_distance_m": float(self.pending_distance_m),
                    "committed_nfire": 1,
                    "wake_on_commit": wake,
                    "source_sites_refreshed_during_opening": float(self._source_refresh_total),
                })
            self._clear_active()
        elif self._prepared_hit_external_cap:
            self.paused = True
            self.loading_resume_requests += 1
            if self._active_record is not None:
                self._active_record.update({
                    "status": "loading_resumed_hazard_limited",
                    "loading_resume_requests": int(self.loading_resume_requests),
                    "hazard_progress_at_loading_resume": float(self.progress),
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
        lam = max(float(lambda_c_current), 0.0)
        self._last_lambda_c = lam
        remaining_increment = min(
            self.target_hazard_increment,
            max(1.0 - float(self.progress), 0.0),
        )
        dt_to_increment = remaining_increment / lam if lam > 0.0 else math.inf
        hold_cap = self._hold_cap_s()
        resume = bool(dt_to_increment <= hold_cap * (1.0 + 1.0e-12))
        probe = {
            "K_cleave_drive_Pa_sqrt_m": float(K_cleave),
            "KJ_Pa_sqrt_m": _finite_or_none(KJ),
            "temperature_K": float(T),
            "lambda_c_current_s-1": lam,
            "lambda_c_ratio_to_nucleation": float(self.current_rate_ratio(lam)),
            "damage_progress": float(self.progress),
            "time_to_next_hazard_increment_s": float(dt_to_increment),
            "external_loading_interval_s": float(hold_cap),
            "resumed": resume,
        }
        if self._active_record is not None:
            self._active_record.setdefault("reload_probes", []).append(probe)
        if resume:
            self.paused = False
            self.total_resumes += 1
            if self._active_record is not None:
                self._active_record.update({
                    "status": "absolute_hazard_opening_resumed",
                    "resume_count": int(self._active_record.get("resume_count", 0) + 1),
                })
        return resume

    def controlled_value(self, base: float, role: str, factor: Any) -> float:
        factor_f = float(factor)
        normal = float(base) * factor_f
        if role == "dt":
            if normal > 0.0:
                self.external_dt_s = normal
            if self.ready_for_relaxation:
                return float(self.prepare_substep())
            return normal
        if role == "dU" and self.ready_for_relaxation:
            return 0.0
        return normal

    def payload(self) -> dict[str, Any]:
        data = super().payload()
        for key in (
            "event_relaxation_time_s", "event_substeps_per_tau",
            "event_rate_exponent", "event_min_rate_ratio",
            "event_resume_rate_ratio", "event_arrest_substeps",
            "event_max_time_multiplier",
        ):
            data.pop(key, None)
        data.update({
            "schema": "absolute_hazard_trial_cohesive_mpz_event_v917_v1",
            "opening_hazard_target": 1.0,
            "progress_law": "dq_dt=lambda_c_absolute",
            "target_hazard_increment": float(self.target_hazard_increment),
            "min_event_dt_s": float(self.min_event_dt_s),
            "max_fixed_hold_s": (
                float(self.max_fixed_hold_s) if math.isfinite(self.max_fixed_hold_s) else None
            ),
            "latest_external_loading_interval_s": self.external_dt_s,
            "loading_resume_requests": int(self.loading_resume_requests),
            "one_fire_per_event_enforced": True,
            "source_refresh_during_opening": True,
        })
        return data


def main(argv=None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    _v916._inject_once(user_args, "--max-advances-per-step", "1")
    original = _v916.KineticTrialEventController
    _v916.KineticTrialEventController = HazardClockTrialEventController
    try:
        results = _v916.main(user_args)
    finally:
        _v916.KineticTrialEventController = original

    out_value = _v916._option_value(user_args, "--out")
    if out_value is not None:
        out = Path(out_value)
        source = out / "kinetic_trial_event_relaxation_v916.json"
        if source.exists():
            payload = json.loads(source.read_text())
            payload["compatibility_source_filename"] = source.name
            (out / "absolute_hazard_event_relaxation_v917.json").write_text(
                json.dumps(payload, indent=2, default=str)
            )
    return results


if __name__ == "__main__":
    main()
