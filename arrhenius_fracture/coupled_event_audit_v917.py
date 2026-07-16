from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any


def _read(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    try:
        value = json.loads(path.read_text())
    except Exception:
        return {}
    return value if isinstance(value, dict) else {}


def _finite(value: Any) -> float | None:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if math.isfinite(out) else None


def audit_case_v917(
    case_dir: str | Path,
    T_K: float,
    analysis_target_extension_um: float,
) -> dict[str, Any]:
    root = Path(case_dir)
    data = _read(root / "absolute_hazard_event_relaxation_v917.json")
    events = data.get("events", []) if isinstance(data.get("events", []), list) else []
    target = float(analysis_target_extension_um)
    tol = max(1.0e-8, 1.0e-8 * abs(target))

    committed: list[dict[str, Any]] = []
    extension_um = 0.0
    for event in events:
        if not isinstance(event, dict):
            continue
        d_um = max(float(event.get("committed_distance_m", 0.0) or 0.0), 0.0) * 1.0e6
        if (
            event.get("status") == "complete_committed"
            and bool(event.get("mpz_renewal_committed", False))
            and d_um > 0.0
            and extension_um + d_um <= target + tol
        ):
            committed.append(event)
            extension_um += d_um

    per_event: list[dict[str, Any]] = []
    software_checks: list[bool] = []
    hazard_checks: list[bool] = []
    refresh_checks: list[bool] = []
    coupling_checks: list[bool] = []
    retention_checks: list[bool] = []

    for event in committed:
        rows = event.get("substeps", [])
        rows = rows if isinstance(rows, list) else []
        q = [_finite(row.get("damage_progress")) for row in rows]
        q = [x for x in q if x is not None]
        dq = [_finite(row.get("cleavage_hazard_increment")) for row in rows]
        dq = [x for x in dq if x is not None]
        fixed_load = bool(rows) and all(
            abs(float(row.get("remote_displacement_increment_m", math.nan))) <= 1.0e-30
            and float(row.get("dt_s", 0.0) or 0.0) > 0.0
            for row in rows
        )
        monotone = bool(q) and all(q[i + 1] + 1.0e-12 >= q[i] for i in range(len(q) - 1))
        hazard_integral = float(sum(dq)) if dq else 0.0
        hazard_ok = (
            event.get("progress_law") == "dq_dt=lambda_c_absolute"
            and math.isclose(hazard_integral, 1.0, rel_tol=0.0, abs_tol=2.0e-6)
            and math.isclose(float(event.get("opening_hazard_integral", 0.0) or 0.0), 1.0,
                             rel_tol=0.0, abs_tol=2.0e-6)
        )
        absolute_rates = bool(rows) and all(
            _finite(row.get("lambda_c_current_s-1")) is not None
            and _finite(row.get("lambda_c_raw_current_s-1")) is not None
            and _finite(row.get("lambda_e_current_s-1")) is not None
            for row in rows
        )
        one_fire = int(event.get("committed_nfire", 0) or 0) == 1
        source_refresh = max(
            float(event.get("source_sites_refreshed_during_opening", 0.0) or 0.0), 0.0
        )
        refresh_ok = "source_sites_refreshed_during_opening" in event
        dN_emit = max(float(event.get("dN_emit_relaxation", 0.0) or 0.0), 0.0)
        Kvals = [_finite(row.get("Kshield_Pa_sqrt_m")) for row in rows]
        Kvals = [x for x in Kvals if x is not None]
        Kspan = max(Kvals) - min(Kvals) if Kvals else 0.0
        retained = [_finite(row.get("retained_count")) for row in rows]
        retained = [x for x in retained if x is not None]
        retained_end = retained[-1] if retained else 0.0
        coupling = dN_emit > 1.0e-6 or Kspan > 1.0e3
        retention = retained_end > 1.0e-6

        software = bool(
            event.get("relaxation_completed", False)
            and math.isclose(float(event.get("final_damage", math.nan)), 1.0,
                             rel_tol=0.0, abs_tol=1.0e-9)
            and bool(event.get("mpz_renewal_deferred", False))
            and bool(event.get("mpz_renewal_committed", False))
            and fixed_load
            and monotone
            and absolute_rates
            and one_fire
        )
        software_checks.append(software)
        hazard_checks.append(hazard_ok)
        refresh_checks.append(refresh_ok)
        coupling_checks.append(coupling)
        retention_checks.append(retention)
        per_event.append({
            "event_id": int(event.get("event_id", -1) or -1),
            "committed_distance_um": float(event.get("committed_distance_m", 0.0) or 0.0) * 1.0e6,
            "relaxation_time_s": float(event.get("relaxation_time_s", 0.0) or 0.0),
            "n_substeps": len(rows),
            "opening_hazard_integral": hazard_integral,
            "lambda_c_nucleation_s-1": _finite(event.get("lambda_c_nucleation_s-1")),
            "predicted_opening_time_at_nucleation_s": _finite(
                event.get("predicted_opening_time_at_nucleation_s")
            ),
            "source_sites_refreshed_during_opening": source_refresh,
            "dN_emit_relaxation": dN_emit,
            "retained_count_end": retained_end,
            "Kshield_span_Pa_sqrt_m": Kspan,
            "loading_resume_requests": int(event.get("loading_resume_requests", 0) or 0),
            "n_reload_probes": len(event.get("reload_probes", []) or []),
            "hazard_clock_consistent": hazard_ok,
            "material_coupling_nonzero": coupling,
            "retained_state_nonzero_at_commit": retention,
        })

    target_committed = extension_um + tol >= target
    software_gate = bool(committed and all(software_checks) and all(hazard_checks) and target_committed)
    refresh_observed = bool(committed and any(
        x["source_sites_refreshed_during_opening"] > 0.0 for x in per_event
    ))
    coupling_observed = bool(committed and any(coupling_checks))
    retention_observed = bool(committed and any(retention_checks))

    payload = {
        "schema": "absolute_hazard_trial_opening_audit_v917_v1",
        "case_dir": str(root),
        "T_K": float(T_K),
        "analysis_target_extension_um": target,
        "committed_extension_um": float(extension_um),
        "target_committed": target_committed,
        "n_committed_events_in_analysis": len(committed),
        "one_fire_per_event_required": True,
        "absolute_hazard_clock_required": True,
        "source_refresh_during_opening_required_when_inventory_depleted": True,
        "software_commit_gate_passed": software_gate,
        "hazard_clock_consistency_passed": bool(committed and all(hazard_checks)),
        "source_refresh_during_opening_observed": refresh_observed,
        "physical_coupling_relevance_observed": coupling_observed,
        "retained_state_at_commit_observed": retention_observed,
        "retention_failure_still_possible": not retention_observed,
        "active_trial_event_at_exit": data.get("active_event_id_at_exit"),
        "active_trial_progress_at_exit": data.get("active_event_progress_at_exit"),
        "event_summaries": per_event,
        "interpretation": (
            "hazard_clock_source_refresh_and_retained_coupling_observed"
            if software_gate and refresh_observed and coupling_observed and retention_observed
            else "hazard_clock_and_source_refresh_work_but_retention_is_weak"
            if software_gate and refresh_observed and coupling_observed
            else "hazard_clock_software_passes_but_material_coupling_remains_negligible"
            if software_gate
            else "hazard_clock_or_commit_implementation_failed"
        ),
    }
    (root / "absolute_hazard_event_audit_v917.json").write_text(
        json.dumps(payload, indent=2, default=str)
    )
    return payload


__all__ = ["audit_case_v917"]
