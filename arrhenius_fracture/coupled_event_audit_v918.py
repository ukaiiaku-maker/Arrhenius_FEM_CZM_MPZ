from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

from .coupled_event_audit_v917 import audit_case_v917


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


def audit_case_v918(
    case_dir: str | Path,
    T_K: float,
    analysis_target_extension_um: float,
) -> dict[str, Any]:
    root = Path(case_dir)
    compatibility = audit_case_v917(root, T_K, analysis_target_extension_um)
    data = _read(root / "persistent_wake_event_relaxation_v918.json")
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

    summaries: list[dict[str, Any]] = []
    conservation_checks: list[bool] = []
    wake_commit_checks: list[bool] = []
    wake_shield_checks: list[bool] = []
    carryover_checks: list[bool] = []

    for index, event in enumerate(committed):
        wake = event.get("wake_on_commit", {})
        wake = wake if isinstance(wake, dict) else {}
        active_pre = max(float(wake.get("active_retained_count_precommit", 0.0) or 0.0), 0.0)
        wake_pre = max(float(wake.get("wake_retained_count_precommit", 0.0) or 0.0), 0.0)
        active_post = max(float(wake.get("active_retained_count_postcommit", 0.0) or 0.0), 0.0)
        wake_post = max(float(wake.get("wake_retained_count_postcommit", 0.0) or 0.0), 0.0)
        discarded = max(float(wake.get("wake_retained_discarded", 0.0) or 0.0), 0.0)
        crossing = max(float(wake.get("wake_retained", 0.0) or 0.0), 0.0)
        wake_K_post = max(float(wake.get("wake_K_shield_Pa_sqrt_m_postcommit", 0.0) or 0.0), 0.0)
        total_K_post = max(float(wake.get("total_K_shield_Pa_sqrt_m_postcommit", 0.0) or 0.0), 0.0)

        before_total = active_pre + wake_pre
        after_total = active_post + wake_post + discarded
        scale = max(before_total, after_total, 1.0)
        conserved = abs(before_total - after_total) <= 1.0e-9 * scale
        wake_committed = bool(
            event.get("persistent_wake_state_committed", False)
            and wake_post + discarded + 1.0e-12 >= crossing
        )
        wake_shield = wake_K_post > 1.0e3

        next_wake = None
        carryover = False
        if index + 1 < len(committed):
            next_wake = _finite(committed[index + 1].get("wake_retained_count_at_nucleation"))
            carryover = next_wake is not None and next_wake > 1.0e-8

        conservation_checks.append(conserved)
        wake_commit_checks.append(wake_committed)
        wake_shield_checks.append(wake_shield)
        if index + 1 < len(committed):
            carryover_checks.append(carryover)

        summaries.append({
            "event_id": int(event.get("event_id", -1) or -1),
            "committed_distance_um": float(event.get("committed_distance_m", 0.0) or 0.0) * 1.0e6,
            "active_retained_precommit": active_pre,
            "wake_retained_precommit": wake_pre,
            "active_retained_postcommit": active_post,
            "wake_retained_postcommit": wake_post,
            "wake_retained_crossing": crossing,
            "wake_retained_discarded": discarded,
            "wake_K_shield_postcommit_Pa_sqrt_m": wake_K_post,
            "total_K_shield_postcommit_Pa_sqrt_m": total_K_post,
            "wake_retained_at_next_event_nucleation": next_wake,
            "retained_state_conserved_on_commit": conserved,
            "persistent_wake_committed": wake_committed,
            "wake_shielding_nonzero": wake_shield,
            "wake_carryover_to_next_event": carryover if next_wake is not None else None,
        })

    target_committed = extension_um + tol >= target
    no_uncommitted_trial = data.get("active_event_id_at_exit") in (None, "", -1)
    wake_commit_gate = bool(committed and all(wake_commit_checks) and all(conservation_checks))
    wake_shield_observed = bool(committed and any(wake_shield_checks))
    carryover_observed = bool(carryover_checks and any(carryover_checks))

    payload = {
        "schema": "persistent_plastic_wake_event_audit_v918_v1",
        "case_dir": str(root),
        "T_K": float(T_K),
        "analysis_target_extension_um": target,
        "committed_extension_um": float(extension_um),
        "target_committed": target_committed,
        "n_committed_events_in_analysis": len(committed),
        "v917_software_commit_gate_passed": bool(
            compatibility.get("software_commit_gate_passed", False)
        ),
        "persistent_wake_commit_gate_passed": wake_commit_gate,
        "retained_state_conservation_passed": bool(
            committed and all(conservation_checks)
        ),
        "wake_shielding_observed": wake_shield_observed,
        "wake_carryover_to_later_event_observed": carryover_observed,
        "no_uncommitted_trial_at_exit": no_uncommitted_trial,
        "adaptive_dt_used_as_hold_cap": bool(
            data.get("adaptive_dt_used_as_hold_cap", True)
        ),
        "nominal_loading_dt_s": _finite(data.get("nominal_loading_dt_s")),
        "committed_target_reached_by_controller": bool(
            data.get("committed_target_reached", False)
        ),
        "post_target_renewals_suppressed": int(
            data.get("post_target_renewals_suppressed", 0) or 0
        ),
        "event_summaries": summaries,
        "interpretation": (
            "persistent_wake_conserved_and_shielding_carries_between_events"
            if wake_commit_gate and wake_shield_observed and carryover_observed
            else "persistent_wake_conserved_but_not_mechanically_consequential"
            if wake_commit_gate
            else "persistent_wake_commit_or_conservation_failed"
        ),
    }
    (root / "persistent_wake_event_audit_v918.json").write_text(
        json.dumps(payload, indent=2, default=str)
    )
    return payload


__all__ = ["audit_case_v918"]
