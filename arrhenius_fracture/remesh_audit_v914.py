"""Audit the five v9.14 event-driven remeshing requirements.

The adaptive CZM may realize one physical 5-um Arrhenius renewal through several
collinear cohesive subsegments.  This audit groups those subsegments by their
``physical_event_id`` and evaluates only events inside the requested analysis
window.  A terminal solver-guard event outside that window may remain pending at
exit; it exists solely so the preceding endpoint event receives a same-load
correction solve.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def _read_csv(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path)
    except (pd.errors.EmptyDataError, OSError, ValueError):
        return pd.DataFrame()


def _read_json(path: Path) -> Any:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text())
    except Exception:
        return {}


def _ordered_unique_positive(values) -> list[int]:
    out: list[int] = []
    seen = set()
    for value in values:
        try:
            event_id = int(value)
        except (TypeError, ValueError):
            continue
        if event_id > 0 and event_id not in seen:
            seen.add(event_id)
            out.append(event_id)
    return out


def audit_case(
    case_dir: str | Path,
    T_K: float,
    same_load_rtol: float = 1.0e-10,
    analysis_target_extension_um: float | None = None,
) -> dict[str, Any]:
    root = Path(case_dir)
    tag = f"{int(round(T_K)):04d}K"
    steps = _read_csv(root / f"steps_{tag}.csv")
    diagnostics_dir = root / f"czm_{tag}"
    log_path = diagnostics_dir / "czm_advance_log.json"
    cohesive_path = diagnostics_dir / "cohesive_elements.csv"
    # Backward-compatible fallbacks for early v9.14 experiments.
    if not log_path.exists():
        log_path = root / "czm_advance_log.json"
    if not cohesive_path.exists():
        cohesive_path = root / "cohesive_elements.csv"
    log = _read_json(log_path)
    equilibrium = _read_json(root / "post_event_equilibrium_audit_v914.json")
    subsegments = log if isinstance(log, list) else []

    fire = pd.to_numeric(
        steps.get("n_fire", pd.Series(dtype=float)), errors="coerce"
    ).fillna(0.0)
    all_fire_idx = np.flatnonzero(fire.to_numpy(float) > 0.0)
    extension_um = pd.to_numeric(
        steps.get("crack_extension_m", pd.Series(np.nan, index=steps.index)),
        errors="coerce",
    ).to_numpy(float) * 1.0e6
    target = (
        float(analysis_target_extension_um)
        if analysis_target_extension_um is not None
        else float("inf")
    )
    target_tol = max(1.0e-8, 1.0e-8 * abs(target)) if math.isfinite(target) else 1.0e-8
    analysis_fire_idx = np.asarray([
        int(i) for i in all_fire_idx
        if not math.isfinite(target)
        or (np.isfinite(extension_um[i]) and extension_um[i] <= target + target_tol)
    ], dtype=int)
    guard_fire_idx = np.asarray([
        int(i) for i in all_fire_idx if int(i) not in set(analysis_fire_idx.tolist())
    ], dtype=int)

    one_event = bool(
        len(analysis_fire_idx) == 0
        or np.nanmax(fire.iloc[analysis_fire_idx].to_numpy(float)) <= 1.0
    )
    adaptive = False
    if not steps.empty and "adaptive_frac" in steps:
        vals = pd.to_numeric(steps["adaptive_frac"], errors="coerce").to_numpy(float)
        adaptive = bool(np.any(np.isfinite(vals) & (vals < 1.0)))

    same_load_checks: list[float] = []
    zero_dt_checks: list[bool] = []
    if not steps.empty and "Uapp_m" in steps:
        U = pd.to_numeric(steps["Uapp_m"], errors="coerce").to_numpy(float)
        dt = pd.to_numeric(
            steps.get("dt_cur_s", pd.Series(np.nan, index=steps.index)),
            errors="coerce",
        ).to_numpy(float)
        for i in analysis_fire_idx:
            if i + 1 < len(U) and np.isfinite(U[i]) and np.isfinite(U[i + 1]):
                rel = abs(U[i + 1] - U[i]) / max(abs(U[i]), 1.0e-30)
                same_load_checks.append(float(rel))
                zero_dt_checks.append(
                    bool(np.isfinite(dt[i + 1]) and abs(dt[i + 1]) <= 1.0e-30)
                )

    log_event_ids = _ordered_unique_positive(
        [row.get("physical_event_id", -1) for row in subsegments]
    )
    n_analysis_fire = int(len(analysis_fire_idx))
    analysis_event_ids = log_event_ids[:n_analysis_fire]
    guard_event_ids = log_event_ids[n_analysis_fire:]
    analysis_event_set = set(analysis_event_ids)
    analysis_subsegments = [
        row for row in subsegments
        if int(row.get("physical_event_id", -1) or -1) in analysis_event_set
    ]

    grouped_lengths: dict[int, float] = {}
    for row in analysis_subsegments:
        event_id = int(row.get("physical_event_id", -1) or -1)
        grouped_lengths[event_id] = grouped_lengths.get(event_id, 0.0) + max(
            float(row.get("length_m", 0.0) or 0.0), 0.0
        )
    crack_log_ok = bool(
        n_analysis_fire > 0
        and len(analysis_event_ids) == n_analysis_fire
        and all(grouped_lengths.get(event_id, 0.0) > 0.0 for event_id in analysis_event_ids)
    )

    corrected_ids = set(_ordered_unique_positive(equilibrium.get("corrected_event_ids", [])))
    cancelled_ids = set(_ordered_unique_positive(equilibrium.get("cancelled_event_ids", [])))
    pending_id = equilibrium.get("pending_event_id")
    try:
        pending_id = int(pending_id) if pending_id is not None else None
    except (TypeError, ValueError):
        pending_id = None
    analysis_corrections_ok = bool(
        analysis_event_ids
        and all(event_id in corrected_ids for event_id in analysis_event_ids)
        and not any(event_id in cancelled_ids for event_id in analysis_event_ids)
        and (pending_id is None or pending_id not in analysis_event_set)
    )
    rowwise_same_load_ok = bool(
        len(same_load_checks) == n_analysis_fire
        and max(same_load_checks, default=np.inf) <= float(same_load_rtol)
        and len(zero_dt_checks) == n_analysis_fire
        and all(zero_dt_checks)
    )
    same_load_ok = analysis_corrections_ok and rowwise_same_load_ok

    fields = root / f"field_snapshot_manifest_{int(round(T_K))}K.json"
    state_outputs_ok = fields.exists() and cohesive_path.exists()
    first_four = bool(one_event and adaptive and crack_log_ok and state_outputs_ok)

    payload = {
        "schema": "event_driven_remesh_audit_v914_v2",
        "case_dir": str(root),
        "T_K": float(T_K),
        "analysis_target_extension_um": (
            float(target) if math.isfinite(target) else None
        ),
        "czm_diagnostics_dir": str(diagnostics_dir),
        "czm_advance_log_path": str(log_path),
        "cohesive_elements_path": str(cohesive_path),
        "n_czm_advances": int(len(analysis_event_ids)),
        "n_physical_events_in_analysis": n_analysis_fire,
        "n_czm_subsegments_in_analysis": int(len(analysis_subsegments)),
        "n_czm_subsegments_total": int(len(subsegments)),
        "physical_event_ids_in_analysis": analysis_event_ids,
        "physical_event_subsegment_lengths_m": {
            str(k): float(v) for k, v in grouped_lengths.items()
        },
        "n_fire_rows": n_analysis_fire,
        "n_guard_fire_rows": int(len(guard_fire_idx)),
        "guard_event_ids": guard_event_ids,
        "one_topology_event_per_accepted_solve": one_event,
        "adaptive_event_time_localization_observed": adaptive,
        "cohesive_path_log_valid": crack_log_ok,
        "state_and_field_outputs_present": state_outputs_ok,
        "same_load_post_event_relative_changes": same_load_checks,
        "zero_time_post_event_rows": zero_dt_checks,
        "explicit_correction_audit": equilibrium,
        "analysis_event_ids_corrected": sorted(corrected_ids & analysis_event_set),
        "analysis_events_all_corrected": analysis_corrections_ok,
        "terminal_pending_guard_event_allowed": bool(
            pending_id is not None and pending_id in set(guard_event_ids)
        ),
        "same_load_post_event_reequilibration_observed": same_load_ok,
        "requirements_1_to_4_passed": first_four,
        "requirement_5_passed": same_load_ok,
        "all_five_requirements_passed": bool(first_four and same_load_ok),
        "interpretation": (
            "all_five_event_driven_remesh_requirements_observed"
            if first_four and same_load_ok
            else "post_event_same_load_equilibrium_missing_or_not_demonstrated"
            if first_four
            else "event_localization_or_remesh_transfer_incomplete"
        ),
    }
    (root / "event_driven_remesh_audit_v914.json").write_text(
        json.dumps(payload, indent=2)
    )
    return payload


__all__ = ["audit_case"]
