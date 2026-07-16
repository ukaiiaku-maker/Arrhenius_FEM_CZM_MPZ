from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np


def _read_json(path: Path) -> Any:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text())
    except Exception:
        return {}


def _finite(value: Any) -> float | None:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if math.isfinite(out) else None


def audit_case_v916(
    case_dir: str | Path,
    T_K: float,
    analysis_target_extension_um: float,
) -> dict[str, Any]:
    root = Path(case_dir)
    relax = _read_json(root / "kinetic_trial_event_relaxation_v916.json")
    events = relax.get("events", []) if isinstance(relax, dict) else []
    events = events if isinstance(events, list) else []

    committed: list[dict[str, Any]] = []
    cumulative_um = 0.0
    target = float(analysis_target_extension_um)
    tol = max(1.0e-8, 1.0e-8 * abs(target))
    for event in events:
        if not isinstance(event, dict):
            continue
        distance_um = max(float(event.get("committed_distance_m", 0.0) or 0.0), 0.0) * 1.0e6
        is_committed = bool(event.get("mpz_renewal_committed", False))
        if is_committed and distance_um > 0.0 and cumulative_um + distance_um <= target + tol:
            committed.append(event)
            cumulative_um += distance_um

    completed = []
    positive_time = []
    fixed_load = []
    deferred_then_committed = []
    kinetic_rows_present = []
    energy_rows_present = []
    progress_monotone = []
    progress_rate_responsive = []
    coupling_nonzero = []
    event_ids = []
    event_summaries = []

    for event in committed:
        event_ids.append(int(event.get("event_id", -1) or -1))
        sub = event.get("substeps", [])
        sub = sub if isinstance(sub, list) else []
        completed.append(
            str(event.get("status", "")) == "complete_committed"
            and bool(event.get("relaxation_completed", False))
            and math.isclose(float(event.get("final_damage", np.nan)), 1.0,
                             rel_tol=0.0, abs_tol=1.0e-9)
        )
        positive_time.append(float(event.get("relaxation_time_s", 0.0) or 0.0) > 0.0)
        fixed_load.append(bool(sub) and all(
            abs(float(row.get("remote_displacement_increment_m", np.nan))) <= 1.0e-30
            and float(row.get("dt_s", 0.0) or 0.0) > 0.0
            for row in sub
        ))
        deferred_then_committed.append(
            bool(event.get("mpz_renewal_deferred", False))
            and bool(event.get("mpz_renewal_committed", False))
            and int(event.get("committed_nfire", 0) or 0) > 0
        )
        kinetic_rows_present.append(bool(sub) and all(
            _finite(row.get("lambda_c_current_s-1")) is not None
            and _finite(row.get("lambda_c_ratio_to_nucleation")) is not None
            and _finite(row.get("damage_increment")) is not None
            for row in sub
        ))
        energy_rows_present.append(bool(sub) and all(
            "cohesive_damage_work_increment_J_per_m" in row
            and "cohesive_recoverable_energy_J_per_m" in row
            and "tip_emission_work_increment_J_per_m" in row
            for row in sub
        ))

        q = [float(row.get("damage_progress", np.nan)) for row in sub]
        dq = [float(row.get("damage_increment", np.nan)) for row in sub]
        ratios = [float(row.get("lambda_c_ratio_to_nucleation", np.nan)) for row in sub]
        qfinite = [x for x in q if math.isfinite(x)]
        progress_monotone.append(
            bool(qfinite)
            and all(qfinite[i + 1] + 1.0e-12 >= qfinite[i] for i in range(len(qfinite) - 1))
            and qfinite[-1] >= 1.0 - 1.0e-9
        )
        finite_dq = [x for x in dq if math.isfinite(x)]
        finite_ratios = [x for x in ratios if math.isfinite(x)]
        responsive = False
        if finite_ratios:
            responsive = (max(finite_ratios) - min(finite_ratios) > 1.0e-8)
        if finite_dq:
            responsive = responsive or (max(finite_dq) - min(finite_dq) > 1.0e-12)
        progress_rate_responsive.append(bool(responsive))

        emitted = max(float(event.get("dN_emit_relaxation", 0.0) or 0.0), 0.0)
        emit_work = max(float(event.get("tip_emission_work_J_per_m", 0.0) or 0.0), 0.0)
        Kshield = [
            _finite(row.get("Kshield_Pa_sqrt_m"))
            for row in sub
        ]
        Kshield = [x for x in Kshield if x is not None]
        Kspan = max(Kshield) - min(Kshield) if Kshield else 0.0
        coupling_nonzero.append(
            emitted > 1.0e-8 or emit_work > 1.0e-12 or Kspan > 1.0e3
        )
        event_summaries.append({
            "event_id": int(event.get("event_id", -1) or -1),
            "committed_distance_um": float(event.get("committed_distance_m", 0.0) or 0.0) * 1.0e6,
            "relaxation_time_s": float(event.get("relaxation_time_s", 0.0) or 0.0),
            "n_substeps": len(sub),
            "dN_emit_relaxation": emitted,
            "tip_emission_work_J_per_m": emit_work,
            "cohesive_work_J_per_m": float(event.get("cohesive_work_J_per_m", 0.0) or 0.0),
            "Kshield_span_Pa_sqrt_m": float(Kspan),
            "min_lambda_c_ratio": _finite(event.get("min_lambda_c_ratio")),
            "max_lambda_c_ratio": _finite(event.get("max_lambda_c_ratio")),
            "n_reload_probes": len(event.get("reload_probes", []) or []),
        })

    software_gate = bool(
        committed
        and all(completed)
        and all(positive_time)
        and all(fixed_load)
        and all(deferred_then_committed)
        and all(kinetic_rows_present)
        and all(energy_rows_present)
        and all(progress_monotone)
    )
    kinetic_response_observed = bool(committed and any(progress_rate_responsive))
    physical_coupling_relevant = bool(committed and any(coupling_nonzero))

    active_event_id = relax.get("active_event_id_at_exit") if isinstance(relax, dict) else None
    active_progress = relax.get("active_event_progress_at_exit") if isinstance(relax, dict) else None
    active_paused = bool(relax.get("active_event_paused_at_exit", False)) if isinstance(relax, dict) else False
    target_committed = cumulative_um + tol >= target

    payload = {
        "schema": "kinetic_trial_cohesive_mpz_event_audit_v916_v1",
        "case_dir": str(root),
        "T_K": float(T_K),
        "analysis_target_extension_um": target,
        "committed_extension_um": float(cumulative_um),
        "target_committed": bool(target_committed),
        "n_committed_events_in_analysis": int(len(committed)),
        "committed_event_ids": event_ids,
        "events_complete_and_damage_one": completed,
        "events_positive_physical_time": positive_time,
        "events_fixed_remote_load": fixed_load,
        "events_deferred_then_committed": deferred_then_committed,
        "events_kinetic_rows_present": kinetic_rows_present,
        "events_energy_rows_present": energy_rows_present,
        "events_progress_monotone": progress_monotone,
        "events_progress_rate_responsive": progress_rate_responsive,
        "events_material_coupling_nonzero": coupling_nonzero,
        "all_trial_commit_requirements_passed": software_gate,
        "kinetic_response_observed": kinetic_response_observed,
        "physical_coupling_relevance_observed": physical_coupling_relevant,
        "active_trial_event_at_exit": active_event_id,
        "active_trial_progress_at_exit": active_progress,
        "active_trial_paused_at_exit": active_paused,
        "total_arrests": int(relax.get("total_arrests", 0) or 0) if isinstance(relax, dict) else 0,
        "total_resumes": int(relax.get("total_resumes", 0) or 0) if isinstance(relax, dict) else 0,
        "event_summaries": event_summaries,
        "interpretation": (
            "kinetic_trial_commit_and_material_coupling_observed"
            if software_gate and physical_coupling_relevant
            else "kinetic_trial_commit_observed_but_material_coupling_negligible"
            if software_gate
            else "kinetic_trial_commit_missing_or_incomplete"
        ),
    }
    (root / "kinetic_trial_event_audit_v916.json").write_text(
        json.dumps(payload, indent=2, default=str)
    )
    return payload


__all__ = ["audit_case_v916"]
