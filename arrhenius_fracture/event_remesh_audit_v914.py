"""Numerical and material-transfer audit for v9.14 event remeshing."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

from .material_rcurve_audit_v913 import audit_campaign as audit_material_v913

CLASSES = ("ceramic", "weakT", "DBTT")


def _json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    try:
        value = json.loads(path.read_text())
    except Exception:
        return {}
    return value if isinstance(value, dict) else {}


def _csv(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path)
    except (pd.errors.EmptyDataError, OSError, ValueError):
        return pd.DataFrame()


def _finite(value, default=np.nan) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return float(default)
    return out if np.isfinite(out) else float(default)


def audit_case_v914(case_dir: str | Path, material_class: str, T_K: float) -> dict[str, Any]:
    root = Path(case_dir)
    tag = f"{int(round(T_K)):04d}K"
    remesh_path = root / f"czm_{tag}" / "event_remesh_audit_v914.json"
    equilibrium_path = root / "event_equilibrium_audit_v914.json"
    remesh = _json(remesh_path)
    equilibrium = _json(equilibrium_path)
    steps = _csv(root / f"steps_{tag}.csv")
    cascade = _csv(root / "R_curve_cascade_metrics.csv")
    raw_events = 0
    if not cascade.empty:
        raw_events = int(cascade.iloc[0].get("n_raw_topology_events", 0) or 0)
    max_n_fire = 0.0
    if "n_fire" in steps.columns and not steps.empty:
        vals = pd.to_numeric(steps["n_fire"], errors="coerce").fillna(0.0)
        max_n_fire = float(vals.max())

    n_remesh = int(remesh.get("n_events", 0) or 0)
    n_equilibrium = int(equilibrium.get("n_post_event_equilibria", 0) or 0)
    area_err = _finite(remesh.get("max_parent_relative_area_conservation_error"))
    total_area_err = _finite(remesh.get("max_relative_total_area_error"))
    boundary_drift = _finite(
        equilibrium.get("max_relative_boundary_displacement_drift")
    )
    rho_err = _finite(equilibrium.get("max_relative_rho_area_integral_error"))
    ep_err = _finite(equilibrium.get("max_relative_ep_area_integral_error"))

    event_count_matches = bool(raw_events > 0 and n_remesh == raw_events)
    equilibrium_count_matches = bool(n_remesh > 0 and n_equilibrium == n_remesh)
    one_event_per_solve = bool(max_n_fire <= 1.0 + 1.0e-12)
    parent_maps = bool(remesh.get("all_parent_maps_valid", False))
    cohesive_state = bool(
        remesh.get("all_preexisting_cohesive_states_unchanged", False)
    )
    one_physical_event = bool(
        remesh.get("all_events_one_physical_cohesive_event", False)
    )
    patch_resolution = bool(remesh.get("all_patch_targets_satisfied", False))
    equilibrium_complete = bool(
        remesh.get("all_post_event_equilibria_completed", False)
        and equilibrium.get("all_same_time", False)
        and equilibrium.get("all_zero_hazard_increment", False)
        and equilibrium.get("all_J_recomputed", False)
    )
    no_failed_attempts = int(remesh.get("n_failed_event_attempts", 0) or 0) == 0
    conservation = bool(
        np.isfinite(area_err) and area_err <= 1.0e-10
        and np.isfinite(total_area_err) and total_area_err <= 1.0e-12
        and np.isfinite(rho_err) and rho_err <= 1.0e-10
        and np.isfinite(ep_err) and ep_err <= 1.0e-10
    )
    same_load = bool(np.isfinite(boundary_drift) and boundary_drift <= 1.0e-12)
    numerical_gate = all((
        event_count_matches,
        equilibrium_count_matches,
        one_event_per_solve,
        parent_maps,
        cohesive_state,
        one_physical_event,
        patch_resolution,
        equilibrium_complete,
        no_failed_attempts,
        conservation,
        same_load,
    ))
    return {
        "material_class": str(material_class),
        "case_dir": str(root),
        "remesh_audit_file": str(remesh_path),
        "equilibrium_audit_file": str(equilibrium_path),
        "n_raw_physical_events": raw_events,
        "n_remeshed_events": n_remesh,
        "n_same_load_equilibria": n_equilibrium,
        "max_n_fire_per_accepted_solve": max_n_fire,
        "event_count_matches": event_count_matches,
        "equilibrium_count_matches": equilibrium_count_matches,
        "one_event_per_equilibrium_solve": one_event_per_solve,
        "all_parent_maps_valid": parent_maps,
        "all_preexisting_cohesive_states_unchanged": cohesive_state,
        "all_events_one_physical_cohesive_event": one_physical_event,
        "all_patch_targets_satisfied": patch_resolution,
        "same_time_same_load_equilibrium_complete": equilibrium_complete,
        "no_failed_event_attempts": no_failed_attempts,
        "max_parent_relative_area_error": area_err,
        "max_relative_total_area_error": total_area_err,
        "max_relative_rho_integral_error": rho_err,
        "max_relative_ep_integral_error": ep_err,
        "max_relative_boundary_displacement_drift": boundary_drift,
        "conservative_transfer_gate_passed": conservation,
        "same_load_gate_passed": same_load,
        "numerical_event_remesh_gate_passed": numerical_gate,
    }


def audit_campaign(
    campaign_root: str | Path,
    seed: int,
    T_K: float,
    classes: Iterable[str] = CLASSES,
    bulk_mode: str = "tip_only",
) -> dict[str, Any]:
    root = Path(campaign_root)
    class_list = [str(x) for x in classes]
    material = audit_material_v913(
        root, seed, T_K, classes=class_list, bulk_mode=bulk_mode
    )
    cases = []
    for cls in class_list:
        case = root / f"seed_{int(seed)}" / bulk_mode / cls / f"T{int(round(T_K))}_th45"
        cases.append(audit_case_v914(case, cls, T_K))

    failed_numerical = [
        row["material_class"] for row in cases
        if not row["numerical_event_remesh_gate_passed"]
    ]
    numerical_gate = not failed_numerical and len(cases) == len(class_list) and len(cases) > 0
    material_gate_v913 = bool(material.get("material_transfer_gate_passed", False))
    material_gate_v914 = bool(numerical_gate and material_gate_v913)
    if not numerical_gate:
        interpretation = "event_remesh_or_same_load_equilibrium_failed"
    elif not material_gate_v913:
        interpretation = "numerically_valid_remeshed_run_but_material_differentiation_not_yet_supported"
    else:
        interpretation = "numerically_valid_remeshed_material_differentiation_supported"

    payload = {
        "schema": "event_remesh_campaign_audit_v914",
        "campaign_root": str(root),
        "seed": int(seed),
        "T_K": float(T_K),
        "bulk_mode": bulk_mode,
        "classes": class_list,
        "cases": cases,
        "failed_numerical_remesh_cases": failed_numerical,
        "numerical_event_remesh_gate_passed": numerical_gate,
        "v913_material_transfer_gate_passed": material_gate_v913,
        "material_transfer_gate_passed_v914": material_gate_v914,
        "interpretation": interpretation,
        "material_audit_v913": material,
        "remesh_scope": (
            "event-centered refinement-only recentering; no global coarsening; "
            "piecewise-constant state transfer is exact by parent subdivision"
        ),
    }
    root.mkdir(parents=True, exist_ok=True)
    (root / "event_remesh_campaign_audit_v914.json").write_text(
        json.dumps(payload, indent=2, default=str)
    )
    pd.DataFrame(cases).to_csv(
        root / "event_remesh_case_audit_v914.csv", index=False
    )
    return payload


__all__ = ["audit_case_v914", "audit_campaign"]
