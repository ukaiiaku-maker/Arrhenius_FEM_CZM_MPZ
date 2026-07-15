"""Hardened material-transfer audit for the full 2-D v9.13+ FEM/CZM gate.

The audit deliberately separates execution/completion, field-output, protocol,
and material-differentiation evidence.  A missing comparison, right-censored run,
or failed subprocess can no longer pass vacuously.  Later stacked execution
branches may supply versioned run-config and case-summary files; the newest
recognized metadata is selected first.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

CLASSES = ("ceramic", "weakT", "DBTT")


def _read_csv(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path)
    except (pd.errors.EmptyDataError, OSError, ValueError):
        return pd.DataFrame()


def _read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    try:
        value = json.loads(path.read_text())
    except Exception:
        return {}
    return value if isinstance(value, dict) else {}


def _finite(frame: pd.DataFrame, name: str, which: str, scale: float = 1.0) -> float:
    if name not in frame.columns:
        return float("nan")
    values = pd.to_numeric(frame[name], errors="coerce").to_numpy(float) * float(scale)
    values = values[np.isfinite(values)]
    if not values.size:
        return float("nan")
    return float(values[-1] if which == "last" else np.max(values))


def _observed_kinit(fp: dict[str, Any]) -> float:
    if str(fp.get("control_state", "")).lower() != "first_passage":
        return float("nan")
    for key in ("Kc_first_existing_MPa_sqrt_m", "KJ_reference_first_MPa_sqrt_m"):
        try:
            value = float(fp.get(key))
        except (TypeError, ValueError):
            continue
        if np.isfinite(value):
            return value
    return float("nan")


def _state_dict(front: dict[str, Any]) -> dict[str, Any]:
    state = front.get("state", {})
    if isinstance(state, dict) and isinstance(state.get("state"), dict):
        state = state["state"]
    return state if isinstance(state, dict) else {}


def _max_emitted_total(case_dir: Path, T_K: float) -> float:
    payload = _read_json(case_dir / f"mpz_state_snapshots_{int(round(T_K)):04d}K.json")
    values: list[float] = []
    for snap in payload.get("snapshots", []):
        for front in snap.get("fronts", []):
            try:
                values.append(float(_state_dict(front).get("emitted_total")))
            except (TypeError, ValueError):
                pass
    for front in payload.get("final_fronts", []):
        try:
            values.append(float(_state_dict(front).get("emitted_total")))
        except (TypeError, ValueError):
            pass
    finite = np.asarray([v for v in values if np.isfinite(v)], dtype=float)
    return float(np.max(finite)) if finite.size else float("nan")


def _crack_path(case_dir: Path, T_K: float) -> np.ndarray:
    frame = _read_csv(case_dir / f"crack_path_{int(round(T_K))}K.csv")
    if frame.empty or not {"x_m", "y_m"}.issubset(frame.columns):
        return np.empty((0, 2), dtype=float)
    return frame[["x_m", "y_m"]].to_numpy(float)


def paths_identical(case_a: Path, case_b: Path, T_K: float, atol_m: float = 1.0e-12) -> bool:
    a = _crack_path(case_a, T_K)
    b = _crack_path(case_b, T_K)
    return bool(a.shape == b.shape and a.size > 0 and np.allclose(a, b, rtol=0.0, atol=atol_m))


def _normalized_shape(case_dir: Path, representation: str) -> tuple[np.ndarray, np.ndarray]:
    if representation == "clustered":
        frame = _read_csv(case_dir / "R_curve_load_events_clustered.csv")
        xname, kname = "crack_extension_um", "KJ_MPa_sqrt_m"
    else:
        frame = _read_csv(case_dir / "R_curve_topology_events_raw.csv")
        xname, kname = "crack_extension_after_um", "KJ_MPa_sqrt_m"
    if frame.empty or not {xname, kname}.issubset(frame.columns):
        return np.asarray([], float), np.asarray([], float)
    x = pd.to_numeric(frame[xname], errors="coerce").to_numpy(float)
    k = pd.to_numeric(frame[kname], errors="coerce").to_numpy(float)
    good = np.isfinite(x) & np.isfinite(k)
    x, k = x[good], k[good]
    if x.size < 2 or abs(float(k[0])) < 1.0e-30:
        return np.asarray([], float), np.asarray([], float)
    x = x - x[0]
    return x / max(float(x[-1]), 1.0e-30), k / float(k[0])


def normalized_shape_metrics(case_a: Path, case_b: Path, representation: str = "clustered", n: int = 101) -> dict[str, float]:
    xa, ya = _normalized_shape(case_a, representation)
    xb, yb = _normalized_shape(case_b, representation)
    if min(xa.size, xb.size) < 2:
        return {"correlation": float("nan"), "relative_rmse": float("nan"), "max_relative_difference": float("nan")}
    grid = np.linspace(0.0, 1.0, max(int(n), 3))
    aa = np.interp(grid, xa, ya)
    bb = np.interp(grid, xb, yb)
    corr = float("nan") if min(np.std(aa), np.std(bb)) <= 1.0e-14 else float(np.corrcoef(aa, bb)[0, 1])
    return {
        "correlation": corr,
        "relative_rmse": float(np.sqrt(np.mean((aa - bb) ** 2))),
        "max_relative_difference": float(np.max(np.abs(aa - bb))),
    }


@dataclass
class CaseAudit:
    material_class: str
    case_dir: str
    T_K: float
    subprocess_returncode: int | None
    solver_success: bool
    completion_marker_present: bool
    control_state: str
    initiation_observed: bool
    target_extension_um: float
    final_extension_um: float
    target_completion_fraction: float
    target_reached: bool
    event_statistics: str
    stochastic_emission: bool
    deterministic_mean_protocol: bool
    K_init_MPa_sqrt_m: float
    n_raw_topology_events: int
    n_independent_load_events: int
    n_unstable_same_load_cascades: int
    cascade_event_fraction: float
    largest_same_load_jump_um: float
    max_K_shield_MPa_sqrt_m: float
    max_K_shield_over_K_init: float
    max_retained_count: float
    max_mobile_count: float
    max_local_slip_count: float
    max_emitted_total: float
    full_field_image_present: bool
    tip_zoom_image_present: bool
    field_manifest_present: bool
    response_classification: str


def audit_case(case_dir: str | Path, material_class: str, T_K: float, completion_tolerance: float = 0.98) -> CaseAudit:
    root = Path(case_dir)
    steps = _read_csv(root / f"steps_{int(round(T_K)):04d}K.csv")
    fp = _read_json(root / "anisotropic_calibrated_tip_first_passage_summary.json")
    casc = _read_csv(root / "R_curve_cascade_metrics.csv")
    c0 = casc.iloc[0].to_dict() if not casc.empty else {}
    run_audit = _read_json(root / "rcurve_run_audit.json")
    config = (
        _read_json(root / "v9_14_run_config.json")
        or _read_json(root / "v9_13_run_config.json")
        or _read_json(root / "v9_12_run_config.json")
    )
    summary = (
        _read_json(root / "v9_14_case_summary.json")
        or _read_json(root / "v9_13_case_summary.json")
        or _read_json(root / "v9_12_case_summary.json")
    )

    kinit = _observed_kinit(fp)
    kshield = _finite(steps, "mpz_K_shield_Pa_sqrt_m", "max", 1.0e-6)
    ratio = kshield / kinit if np.isfinite(kshield) and np.isfinite(kinit) and kinit > 0 else float("nan")
    final_extension = _finite(steps, "crack_extension_m", "last", 1.0e6)
    try:
        target = float(run_audit.get("target_extension_um"))
    except (TypeError, ValueError):
        target = float("nan")
    completion_fraction = final_extension / target if np.isfinite(final_extension) and np.isfinite(target) and target > 0 else float("nan")
    marker = (root / ".long_growth_complete").exists()
    target_reached = bool(marker and np.isfinite(completion_fraction) and completion_fraction >= float(completion_tolerance))

    rc_value = summary.get("subprocess_returncode")
    try:
        rc = int(rc_value)
    except (TypeError, ValueError):
        rc = None
    event_statistics = str(config.get("event_statistics", "unknown")).lower()
    stochastic_emission = bool(config.get("stochastic_emission", False))
    deterministic_mean = event_statistics == "deterministic" and not stochastic_emission
    control_state = str(fp.get("control_state", "unknown"))
    initiation = control_state.lower() == "first_passage" and np.isfinite(kinit)
    n_load = int(c0.get("n_independent_load_events", 0) or 0)
    cascade_fraction = float(c0.get("fraction_topology_events_in_cascades", np.nan))
    if not target_reached:
        classification = "incomplete_or_right_censored"
    elif n_load == 0:
        classification = "no_crack_growth"
    elif n_load <= 2 or (np.isfinite(cascade_fraction) and cascade_fraction >= 0.5):
        classification = "unstable_fixed_displacement_propagation"
    else:
        classification = "candidate_stable_resistance_sequence"

    return CaseAudit(
        material_class=str(material_class), case_dir=str(root), T_K=float(T_K),
        subprocess_returncode=rc, solver_success=bool(rc == 0),
        completion_marker_present=marker, control_state=control_state,
        initiation_observed=initiation, target_extension_um=target,
        final_extension_um=final_extension, target_completion_fraction=completion_fraction,
        target_reached=target_reached, event_statistics=event_statistics,
        stochastic_emission=stochastic_emission, deterministic_mean_protocol=deterministic_mean,
        K_init_MPa_sqrt_m=kinit,
        n_raw_topology_events=int(c0.get("n_raw_topology_events", 0) or 0),
        n_independent_load_events=n_load,
        n_unstable_same_load_cascades=int(c0.get("n_unstable_same_load_cascades", 0) or 0),
        cascade_event_fraction=cascade_fraction,
        largest_same_load_jump_um=float(c0.get("largest_same_load_jump_um", np.nan)),
        max_K_shield_MPa_sqrt_m=kshield, max_K_shield_over_K_init=ratio,
        max_retained_count=_finite(steps, "mpz_retained_count", "max"),
        max_mobile_count=_finite(steps, "mpz_mobile_count", "max"),
        max_local_slip_count=_finite(steps, "mpz_local_slip_count", "max"),
        max_emitted_total=_max_emitted_total(root, T_K),
        full_field_image_present=(root / f"field_snapshots_{int(round(T_K))}K.png").exists(),
        tip_zoom_image_present=(root / f"field_snapshots_tip_zoom_{int(round(T_K))}K.png").exists(),
        field_manifest_present=(root / f"field_snapshot_manifest_{int(round(T_K))}K.json").exists(),
        response_classification=classification,
    )


def audit_campaign(campaign_root: str | Path, seed: int, T_K: float, classes: Iterable[str] = CLASSES, bulk_mode: str = "tip_only") -> dict[str, Any]:
    root = Path(campaign_root)
    class_list = [str(x) for x in classes]
    case_dirs = {c: root / f"seed_{int(seed)}" / bulk_mode / c / f"T{int(round(T_K))}_th45" for c in class_list}
    cases = [audit_case(case_dirs[c], c, T_K) for c in class_list]
    pairs: list[dict[str, Any]] = []
    for i, a in enumerate(class_list):
        for b in class_list[i + 1:]:
            clustered = normalized_shape_metrics(case_dirs[a], case_dirs[b], "clustered")
            raw = normalized_shape_metrics(case_dirs[a], case_dirs[b], "raw")
            same_path = paths_identical(case_dirs[a], case_dirs[b], T_K)
            strongly_similar = bool(
                same_path
                and np.isfinite(clustered["correlation"])
                and clustered["correlation"] >= 0.98
                and clustered["relative_rmse"] <= 0.03
            )
            nearly_identical = bool(
                same_path
                and np.isfinite(clustered["correlation"])
                and clustered["correlation"] >= 0.995
                and clustered["relative_rmse"] <= 0.01
            )
            pairs.append({
                "class_a": a, "class_b": b, "crack_path_identical": same_path,
                "clustered_shape_correlation": clustered["correlation"],
                "clustered_relative_rmse": clustered["relative_rmse"],
                "clustered_max_relative_difference": clustered["max_relative_difference"],
                "raw_shape_correlation": raw["correlation"],
                "raw_relative_rmse": raw["relative_rmse"],
                "strongly_similar_geometry_path": strongly_similar,
                "nearly_identical_geometry_scaled_response": nearly_identical,
            })

    expected_pairs = len(class_list) * (len(class_list) - 1) // 2
    failed = [c.material_class for c in cases if not c.solver_success]
    incomplete = [c.material_class for c in cases if not c.target_reached]
    no_init = [c.material_class for c in cases if not c.initiation_observed]
    missing_fields = [c.material_class for c in cases if not (c.full_field_image_present and c.tip_zoom_image_present and c.field_manifest_present)]
    non_mean = [c.material_class for c in cases if not c.deterministic_mean_protocol]
    similar_pairs = [f"{p['class_a']}:{p['class_b']}" for p in pairs if p["strongly_similar_geometry_path"]]
    identical_pairs = [f"{p['class_a']}:{p['class_b']}" for p in pairs if p["nearly_identical_geometry_scaled_response"]]

    execution_gate = not failed
    completion_gate = not incomplete and not no_init
    field_gate = not missing_fields
    coverage_gate = len(class_list) >= 3 and len(pairs) == expected_pairs and expected_pairs >= 3
    protocol_gate = not non_mean
    differentiation_gate = coverage_gate and not similar_pairs
    transfer_gate = all((execution_gate, completion_gate, field_gate, coverage_gate, protocol_gate, differentiation_gate))

    if failed or incomplete or no_init:
        interpretation = "execution_or_completion_failed"
    elif not coverage_gate:
        interpretation = "insufficient_material_pair_coverage"
    elif non_mean:
        interpretation = "stochastic_ensemble_not_deterministic_transfer_gate"
    elif similar_pairs:
        interpretation = "geometry_scaled_or_strongly_similar_material_response"
    else:
        interpretation = "deterministic_material_differentiation_supported_for_this_gate"

    payload = {
        "schema": "material_rcurve_audit_v913",
        "campaign_root": str(root), "seed": int(seed), "T_K": float(T_K), "bulk_mode": bulk_mode,
        "classes_requested": class_list, "n_material_classes": len(class_list),
        "n_pairwise_comparisons": len(pairs), "expected_pairwise_comparisons": expected_pairs,
        "cases": [asdict(c) for c in cases], "pairwise_shape_audit": pairs,
        "failed_solver_cases": failed, "incomplete_cases": incomplete,
        "missing_initiation_cases": no_init, "missing_field_outputs": missing_fields,
        "non_deterministic_mean_cases": non_mean,
        "strongly_similar_pairs": similar_pairs, "nearly_identical_pairs": identical_pairs,
        "execution_gate_passed": execution_gate, "completion_gate_passed": completion_gate,
        "field_output_gate_passed": field_gate, "pairwise_coverage_gate_passed": coverage_gate,
        "deterministic_mean_protocol_gate_passed": protocol_gate,
        "material_differentiation_gate_passed": differentiation_gate,
        "material_transfer_gate_passed": transfer_gate, "interpretation": interpretation,
        "KJ_conversion_note": "J is converted to KJ with the solver isotropic Eprime convention even when anisotropic stiffness is active",
    }
    root.mkdir(parents=True, exist_ok=True)
    (root / "material_rcurve_audit_v913.json").write_text(json.dumps(payload, indent=2, default=str))
    pd.DataFrame(payload["cases"]).to_csv(root / "material_rcurve_case_audit_v913.csv", index=False)
    pd.DataFrame(pairs).to_csv(root / "material_rcurve_pairwise_audit_v913.csv", index=False)
    return payload


__all__ = ["CaseAudit", "audit_case", "audit_campaign", "normalized_shape_metrics", "paths_identical"]
