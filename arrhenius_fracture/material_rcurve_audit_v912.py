"""Material-versus-geometry audit for full 2-D v9.12 FEM/CZM campaigns.

The adaptive cohesive backend may serialize a fixed-displacement instability into
many edge insertions.  A smooth sequence of those topology updates is not, by
itself, a material resistance curve.  This module combines the load-event
postprocessing with the actual MPZ state and full-field output inventory so that
campaigns are classified before their curves are interpreted.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
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
        return json.loads(path.read_text())
    except Exception:
        return {}


def _finite_max(frame: pd.DataFrame, name: str, scale: float = 1.0) -> float:
    if name not in frame.columns:
        return float("nan")
    q = pd.to_numeric(frame[name], errors="coerce").to_numpy(float) * float(scale)
    q = q[np.isfinite(q)]
    return float(np.max(q)) if q.size else float("nan")


def _finite_last(frame: pd.DataFrame, name: str, scale: float = 1.0) -> float:
    if name not in frame.columns:
        return float("nan")
    q = pd.to_numeric(frame[name], errors="coerce").to_numpy(float) * float(scale)
    q = q[np.isfinite(q)]
    return float(q[-1]) if q.size else float("nan")


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


def _crack_path(case_dir: Path, T_K: float) -> np.ndarray:
    path = case_dir / f"crack_path_{int(round(T_K))}K.csv"
    frame = _read_csv(path)
    if frame.empty or not {"x_m", "y_m"}.issubset(frame.columns):
        return np.empty((0, 2), dtype=float)
    return frame[["x_m", "y_m"]].to_numpy(float)


def _normalized_shape(case_dir: Path) -> tuple[np.ndarray, np.ndarray]:
    raw = _read_csv(case_dir / "R_curve_topology_events_raw.csv")
    if raw.empty or not {"crack_extension_after_um", "KJ_MPa_sqrt_m"}.issubset(raw.columns):
        return np.asarray([], float), np.asarray([], float)
    x = pd.to_numeric(raw["crack_extension_after_um"], errors="coerce").to_numpy(float)
    k = pd.to_numeric(raw["KJ_MPa_sqrt_m"], errors="coerce").to_numpy(float)
    good = np.isfinite(x) & np.isfinite(k)
    x, k = x[good], k[good]
    if x.size < 2 or not np.isfinite(k[0]) or abs(k[0]) < 1.0e-30:
        return np.asarray([], float), np.asarray([], float)
    x = x - x[0]
    span = max(float(x[-1]), 1.0e-30)
    return x / span, k / k[0]


def normalized_shape_correlation(case_a: Path, case_b: Path, n: int = 101) -> float:
    xa, ya = _normalized_shape(case_a)
    xb, yb = _normalized_shape(case_b)
    if min(xa.size, xb.size) < 2:
        return float("nan")
    grid = np.linspace(0.0, 1.0, max(int(n), 3))
    aa = np.interp(grid, xa, ya)
    bb = np.interp(grid, xb, yb)
    if np.std(aa) <= 1.0e-14 or np.std(bb) <= 1.0e-14:
        return float("nan")
    return float(np.corrcoef(aa, bb)[0, 1])


def paths_identical(case_a: Path, case_b: Path, T_K: float, atol_m: float = 1.0e-12) -> bool:
    a = _crack_path(case_a, T_K)
    b = _crack_path(case_b, T_K)
    return bool(a.shape == b.shape and a.size > 0 and np.allclose(a, b, rtol=0.0, atol=atol_m))


@dataclass
class CaseAudit:
    material_class: str
    case_dir: str
    T_K: float
    control_state: str
    K_init_MPa_sqrt_m: float
    final_extension_um: float
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
    max_emitted_ledger: float
    full_field_image: str | None
    full_field_image_present: bool
    field_rows: str
    response_classification: str


def audit_case(case_dir: str | Path, material_class: str, T_K: float) -> CaseAudit:
    root = Path(case_dir)
    steps = _read_csv(root / f"steps_{int(round(T_K)):04d}K.csv")
    fp = _read_json(root / "anisotropic_calibrated_tip_first_passage_summary.json")
    casc = _read_csv(root / "R_curve_cascade_metrics.csv")
    c0 = casc.iloc[0].to_dict() if not casc.empty else {}
    kinit = _observed_kinit(fp)
    kshield = _finite_max(steps, "mpz_K_shield_Pa_sqrt_m", 1.0e-6)
    ratio = kshield / kinit if np.isfinite(kshield) and np.isfinite(kinit) and kinit > 0 else float("nan")
    image = root / f"field_snapshots_{int(round(T_K))}K.png"
    n_load = int(c0.get("n_independent_load_events", 0) or 0)
    cascade_fraction = float(c0.get("fraction_topology_events_in_cascades", np.nan))
    if n_load == 0:
        classification = "no_crack_growth"
    elif n_load <= 2 or (np.isfinite(cascade_fraction) and cascade_fraction >= 0.5):
        classification = "unstable_fixed_displacement_propagation"
    else:
        classification = "candidate_stable_resistance_sequence"
    return CaseAudit(
        material_class=str(material_class),
        case_dir=str(root),
        T_K=float(T_K),
        control_state=str(fp.get("control_state", "unknown")),
        K_init_MPa_sqrt_m=kinit,
        final_extension_um=_finite_last(steps, "crack_extension_m", 1.0e6),
        n_raw_topology_events=int(c0.get("n_raw_topology_events", 0) or 0),
        n_independent_load_events=n_load,
        n_unstable_same_load_cascades=int(c0.get("n_unstable_same_load_cascades", 0) or 0),
        cascade_event_fraction=cascade_fraction,
        largest_same_load_jump_um=float(c0.get("largest_same_load_jump_um", np.nan)),
        max_K_shield_MPa_sqrt_m=kshield,
        max_K_shield_over_K_init=ratio,
        max_retained_count=_finite_max(steps, "mpz_retained_count"),
        max_mobile_count=_finite_max(steps, "mpz_mobile_count"),
        max_local_slip_count=_finite_max(steps, "mpz_local_slip_count"),
        max_emitted_ledger=_finite_max(steps, "N_em"),
        full_field_image=str(image) if image.exists() else None,
        full_field_image_present=image.exists(),
        field_rows="damage;log10_rho;sigma1_FEM;equivalent_plastic_strain;crack_path_overlay",
        response_classification=classification,
    )


def audit_campaign(
    campaign_root: str | Path,
    seed: int,
    T_K: float,
    classes: Iterable[str] = CLASSES,
    bulk_mode: str = "tip_only",
) -> dict[str, Any]:
    root = Path(campaign_root)
    class_list = [str(x) for x in classes]
    case_dirs = {
        cls: root / f"seed_{int(seed)}" / bulk_mode / cls / f"T{int(round(T_K))}_th45"
        for cls in class_list
    }
    cases = [audit_case(case_dirs[cls], cls, T_K) for cls in class_list]
    pairs: list[dict[str, Any]] = []
    for i, a in enumerate(class_list):
        for b in class_list[i + 1:]:
            corr = normalized_shape_correlation(case_dirs[a], case_dirs[b])
            same_path = paths_identical(case_dirs[a], case_dirs[b], T_K)
            pairs.append({
                "class_a": a,
                "class_b": b,
                "normalized_raw_shape_correlation": corr,
                "crack_path_identical": same_path,
                "geometry_dominated_similarity": bool(
                    same_path and np.isfinite(corr) and corr >= 0.995
                ),
            })
    geometry_pairs = [p for p in pairs if p["geometry_dominated_similarity"]]
    missing_fields = [c.material_class for c in cases if not c.full_field_image_present]
    payload = {
        "schema": "material_rcurve_audit_v912",
        "campaign_root": str(root),
        "seed": int(seed),
        "T_K": float(T_K),
        "bulk_mode": bulk_mode,
        "cases": [asdict(c) for c in cases],
        "pairwise_shape_audit": pairs,
        "n_geometry_dominated_pairs": len(geometry_pairs),
        "geometry_dominated_pairs": [f"{p['class_a']}:{p['class_b']}" for p in geometry_pairs],
        "missing_full_field_images": missing_fields,
        "material_rcurve_gate_passed": bool(not geometry_pairs and not missing_fields),
        "interpretation": (
            "geometry_or_continuation_dominated_do_not_publish_as_material_R_curves"
            if geometry_pairs
            else "no_pairwise_geometry_dominance_detected"
        ),
    }
    root.mkdir(parents=True, exist_ok=True)
    (root / "material_rcurve_audit_v912.json").write_text(json.dumps(payload, indent=2, default=str))
    pd.DataFrame(payload["cases"]).to_csv(root / "material_rcurve_case_audit_v912.csv", index=False)
    pd.DataFrame(pairs).to_csv(root / "material_rcurve_pairwise_audit_v912.csv", index=False)
    return payload


__all__ = [
    "CaseAudit",
    "audit_case",
    "audit_campaign",
    "normalized_shape_correlation",
    "paths_identical",
]
