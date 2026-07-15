"""Audit the five v9.14 event-driven remeshing requirements."""
from __future__ import annotations

import json
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


def audit_case(case_dir: str | Path, T_K: float, same_load_rtol: float = 1.0e-10) -> dict[str, Any]:
    root = Path(case_dir)
    steps = _read_csv(root / f"steps_{int(round(T_K)):04d}K.csv")
    log = _read_json(root / "czm_advance_log.json")
    equilibrium = _read_json(root / "post_event_equilibrium_audit_v914.json")
    rows = log if isinstance(log, list) else []
    fire = pd.to_numeric(steps.get("n_fire", pd.Series(dtype=float)), errors="coerce").fillna(0.0)
    fire_idx = np.flatnonzero(fire.to_numpy(float) > 0.0)

    one_event = bool(len(fire_idx) == 0 or np.nanmax(fire.iloc[fire_idx].to_numpy(float)) <= 1.0)
    adaptive = False
    if not steps.empty and "adaptive_frac" in steps:
        vals = pd.to_numeric(steps["adaptive_frac"], errors="coerce").to_numpy(float)
        adaptive = bool(np.any(np.isfinite(vals) & (vals < 1.0)))

    same_load_checks = []
    zero_dt_checks = []
    if not steps.empty and "Uapp_m" in steps:
        U = pd.to_numeric(steps["Uapp_m"], errors="coerce").to_numpy(float)
        dt = pd.to_numeric(steps.get("dt_cur_s", pd.Series(np.nan, index=steps.index)), errors="coerce").to_numpy(float)
        for i in fire_idx:
            if i + 1 < len(U) and np.isfinite(U[i]) and np.isfinite(U[i + 1]):
                rel = abs(U[i + 1] - U[i]) / max(abs(U[i]), 1.0e-30)
                same_load_checks.append(float(rel))
                zero_dt_checks.append(bool(np.isfinite(dt[i + 1]) and abs(dt[i + 1]) <= 1.0e-30))

    explicit_corrections_ok = bool(
        equilibrium.get("all_scheduled_corrections_consumed", False)
        and int(equilibrium.get("events_scheduled", -1)) == len(rows)
        and int(equilibrium.get("corrections_consumed", -2)) == len(rows)
    )
    rowwise_same_load_ok = bool(
        same_load_checks
        and max(same_load_checks) <= float(same_load_rtol)
        and all(zero_dt_checks)
    )
    same_load_ok = explicit_corrections_ok and rowwise_same_load_ok

    crack_log_ok = bool(rows and all(float(r.get("length_m", 0.0)) > 0.0 for r in rows))
    fields = root / f"field_snapshot_manifest_{int(round(T_K))}K.json"
    state_outputs_ok = fields.exists() and (root / "cohesive_elements.csv").exists()
    first_four = bool(one_event and adaptive and crack_log_ok and state_outputs_ok)

    payload = {
        "schema": "event_driven_remesh_audit_v914",
        "case_dir": str(root), "T_K": float(T_K),
        "n_czm_advances": len(rows), "n_fire_rows": int(len(fire_idx)),
        "one_topology_event_per_accepted_solve": one_event,
        "adaptive_event_time_localization_observed": adaptive,
        "cohesive_path_log_valid": crack_log_ok,
        "state_and_field_outputs_present": state_outputs_ok,
        "same_load_post_event_relative_changes": same_load_checks,
        "zero_time_post_event_rows": zero_dt_checks,
        "explicit_correction_audit": equilibrium,
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
    (root / "event_driven_remesh_audit_v914.json").write_text(json.dumps(payload, indent=2))
    return payload


__all__ = ["audit_case"]
