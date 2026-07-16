from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


def _read_json(path: Path) -> Any:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text())
    except Exception:
        return {}


def _read_csv(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path)
    except (pd.errors.EmptyDataError, OSError, ValueError):
        return pd.DataFrame()


def audit_case_v915(
    case_dir: str | Path,
    T_K: float,
    analysis_target_extension_um: float,
) -> dict[str, Any]:
    root = Path(case_dir)
    tag = f"{int(round(float(T_K))):04d}K"
    steps = _read_csv(root / f"steps_{tag}.csv")
    relax = _read_json(root / "coupled_event_relaxation_v915.json")
    diagnostics_dir = root / f"czm_{tag}"
    log = _read_json(diagnostics_dir / "czm_advance_log.json")
    subsegments = log if isinstance(log, list) else []
    events = relax.get("events", []) if isinstance(relax, dict) else []
    events = events if isinstance(events, list) else []

    fire = pd.to_numeric(
        steps.get("n_fire", pd.Series(dtype=float)), errors="coerce"
    ).fillna(0.0)
    extension_um = pd.to_numeric(
        steps.get("crack_extension_m", pd.Series(np.nan, index=steps.index)),
        errors="coerce",
    ).to_numpy(float) * 1.0e6
    tol = max(1.0e-8, 1.0e-8 * abs(float(analysis_target_extension_um)))
    fire_idx = [
        int(i) for i in np.flatnonzero(fire.to_numpy(float) > 0.0)
        if np.isfinite(extension_um[i])
        and extension_um[i] <= float(analysis_target_extension_um) + tol
    ]
    n_analysis = int(len(fire_idx))

    event_map = {}
    for row in events:
        try:
            event_id = int(row.get("event_id"))
        except (TypeError, ValueError):
            continue
        event_map[event_id] = row
    analysis_ids = list(range(1, n_analysis + 1))
    analysis_events = [event_map.get(i, {}) for i in analysis_ids]

    completed = [
        bool(row.get("relaxation_completed", False))
        and str(row.get("status", "")) == "complete"
        for row in analysis_events
    ]
    positive_time = [
        float(row.get("relaxation_time_s", 0.0) or 0.0) > 0.0
        for row in analysis_events
    ]
    finite_substeps = [
        int(row.get("n_relaxation_substeps_completed", 0) or 0)
        == int(relax.get("event_relaxation_substeps", 0) or 0)
        for row in analysis_events
    ]
    final_damage = [
        math.isclose(float(row.get("final_damage", np.nan)), 1.0, rel_tol=0.0, abs_tol=1e-12)
        for row in analysis_events
    ]
    fixed_load = []
    emission_integrals = []
    K_paths = []
    for row in analysis_events:
        sub = row.get("substeps", []) if isinstance(row, dict) else []
        sub = sub if isinstance(sub, list) else []
        fixed_load.append(bool(sub) and all(
            abs(float(x.get("remote_displacement_increment_m", np.nan))) <= 1e-30
            and float(x.get("dt_s", 0.0) or 0.0) > 0.0
            for x in sub
        ))
        emission_integrals.append(float(row.get("dN_emit_relaxation", 0.0) or 0.0))
        K_paths.append([
            x.get("KJ_Pa_sqrt_m") for x in sub
            if x.get("KJ_Pa_sqrt_m") is not None
        ])

    log_ids = []
    for row in subsegments:
        try:
            event_id = int(row.get("physical_event_id", -1))
        except (TypeError, ValueError):
            continue
        if event_id > 0 and event_id not in log_ids:
            log_ids.append(event_id)
    geometry_ok = bool(
        n_analysis > 0
        and all(i in log_ids for i in analysis_ids)
        and all(
            any(int(seg.get("physical_event_id", -1) or -1) == i
                and float(seg.get("length_m", 0.0) or 0.0) > 0.0
                for seg in subsegments)
            for i in analysis_ids
        )
    )
    fields_ok = bool(
        (root / f"field_snapshot_manifest_{int(round(float(T_K)))}K.json").exists()
        and (diagnostics_dir / "cohesive_elements.csv").exists()
    )
    all_coupled = bool(
        n_analysis > 0
        and len(analysis_events) == n_analysis
        and all(completed)
        and all(positive_time)
        and all(finite_substeps)
        and all(final_damage)
        and all(fixed_load)
        and geometry_ok
        and fields_ok
    )

    payload = {
        "schema": "coupled_cohesive_mpz_event_audit_v915_v1",
        "case_dir": str(root),
        "T_K": float(T_K),
        "analysis_target_extension_um": float(analysis_target_extension_um),
        "n_physical_events_in_analysis": n_analysis,
        "analysis_event_ids": analysis_ids,
        "analysis_events_relaxation_completed": completed,
        "analysis_events_positive_physical_time": positive_time,
        "analysis_events_requested_substeps_completed": finite_substeps,
        "analysis_events_final_damage_one": final_damage,
        "analysis_events_fixed_remote_load": fixed_load,
        "relaxation_emission_counts": emission_integrals,
        "relaxation_KJ_paths_Pa_sqrt_m": K_paths,
        "cohesive_path_log_valid": geometry_ok,
        "state_and_field_outputs_present": fields_ok,
        "geometry_retry_attempts": int(relax.get("geometry_retry_attempts", 0) or 0),
        "geometry_retry_successes": int(relax.get("geometry_retry_successes", 0) or 0),
        "active_guard_event_at_exit": relax.get("active_event_id_at_exit"),
        "all_coupled_event_requirements_passed": all_coupled,
        "interpretation": (
            "coupled_cohesive_mpz_relaxation_observed"
            if all_coupled
            else "coupled_event_relaxation_missing_or_incomplete"
        ),
    }
    (root / "coupled_event_audit_v915.json").write_text(
        json.dumps(payload, indent=2, default=str)
    )
    return payload


__all__ = ["audit_case_v915"]
