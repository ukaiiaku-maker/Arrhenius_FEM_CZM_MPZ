"""Fail-closed v9.12.1 material-versus-geometry publication audit.

The original v9.12 gate only checked for geometry-dominated pairs and missing
images.  Consequently a failed, right-censored, incomplete, or single-class run
could pass vacuously.  This release requires affirmative solver and material
response evidence before pairwise material interpretation is allowed.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

from .material_rcurve_audit_v912 import (
    CLASSES,
    _finite_last,
    _finite_max,
    _observed_kinit,
    _read_csv,
    _read_json,
    normalized_shape_correlation,
    paths_identical,
)

SCHEMA = "material_rcurve_audit_v9121_fail_closed"


@dataclass
class CaseAudit:
    material_class: str
    case_dir: str
    T_K: float
    subprocess_returncode: int | None
    solver_status: str
    control_state: str
    K_init_MPa_sqrt_m: float
    requested_target_extension_um: float
    final_extension_um: float
    solver_invocation_passed: bool
    first_passage_observed: bool
    target_extension_reached: bool
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
    response_classification: str
    stable_resistance_sequence: bool
    case_gate_passed: bool
    case_gate_failures: list[str]


def _coerce_int(value) -> int | None:
    try:
        parsed = int(float(value))
    except (TypeError, ValueError):
        return None
    return parsed


def _coerce_float(value, default=float("nan")) -> float:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return float(default)
    return parsed


def _case_contract(root: Path) -> dict[str, Any]:
    contract = _read_json(root / "v9_12_1_case_contract.json")
    if contract:
        return contract
    return _read_json(root / "v9_12_case_summary.json")


def audit_case(case_dir: str | Path, material_class: str, T_K: float) -> CaseAudit:
    root = Path(case_dir)
    steps = _read_csv(root / f"steps_{int(round(T_K)):04d}K.csv")
    fp = _read_json(root / "anisotropic_calibrated_tip_first_passage_summary.json")
    casc = _read_csv(root / "R_curve_cascade_metrics.csv")
    c0 = casc.iloc[0].to_dict() if not casc.empty else {}
    contract = _case_contract(root)

    returncode = _coerce_int(contract.get("subprocess_returncode"))
    status = str(contract.get("solver_status", contract.get("status", "unknown")))
    control_state = str(fp.get("control_state", contract.get("control_state", "unknown")))
    target_um = _coerce_float(
        contract.get("requested_target_extension_um", contract.get("target_extension_um"))
    )
    final_um = _finite_last(steps, "crack_extension_m", 1.0e6)
    if not np.isfinite(final_um):
        final_um = _coerce_float(contract.get("final_extension_um"))
    kinit = _observed_kinit(fp)
    kshield = _finite_max(steps, "mpz_K_shield_Pa_sqrt_m", 1.0e-6)
    ratio = (
        kshield / kinit
        if np.isfinite(kshield) and np.isfinite(kinit) and kinit > 0.0
        else float("nan")
    )
    image = root / f"field_snapshots_{int(round(T_K))}K.png"
    n_load = int(c0.get("n_independent_load_events", 0) or 0)
    cascade_fraction = float(c0.get("fraction_topology_events_in_cascades", np.nan))
    if n_load == 0:
        classification = "no_crack_growth"
    elif n_load <= 2 or (
        np.isfinite(cascade_fraction) and cascade_fraction >= 0.5
    ):
        classification = "unstable_fixed_displacement_propagation"
    else:
        classification = "candidate_stable_resistance_sequence"

    invocation_passed = returncode == 0
    first_passage = control_state.lower() == "first_passage" and np.isfinite(kinit)
    extension_tolerance = max(1.0e-6, 1.0e-3 * target_um) if np.isfinite(target_um) else 0.0
    target_reached = bool(
        np.isfinite(target_um)
        and target_um > 0.0
        and np.isfinite(final_um)
        and final_um >= target_um - extension_tolerance
    )
    stable = classification == "candidate_stable_resistance_sequence"
    failures: list[str] = []
    if not invocation_passed:
        failures.append("solver_subprocess_failed_or_returncode_missing")
    if not first_passage:
        failures.append("first_passage_not_observed")
    if not target_reached:
        failures.append("requested_target_extension_not_reached")
    if not stable:
        failures.append("no_stable_resistance_sequence")
    if not image.exists():
        failures.append("full_field_snapshot_missing")
    if status.lower() in {
        "right_censored",
        "right_censored_max_blocks",
        "failed",
        "incomplete",
    }:
        failures.append(f"solver_status_{status.lower()}")

    return CaseAudit(
        material_class=str(material_class),
        case_dir=str(root),
        T_K=float(T_K),
        subprocess_returncode=returncode,
        solver_status=status,
        control_state=control_state,
        K_init_MPa_sqrt_m=kinit,
        requested_target_extension_um=target_um,
        final_extension_um=final_um,
        solver_invocation_passed=invocation_passed,
        first_passage_observed=first_passage,
        target_extension_reached=target_reached,
        n_raw_topology_events=int(c0.get("n_raw_topology_events", 0) or 0),
        n_independent_load_events=n_load,
        n_unstable_same_load_cascades=int(
            c0.get("n_unstable_same_load_cascades", 0) or 0
        ),
        cascade_event_fraction=cascade_fraction,
        largest_same_load_jump_um=float(c0.get("largest_same_load_jump_um", np.nan)),
        max_K_shield_MPa_sqrt_m=kshield,
        max_K_shield_over_K_init=ratio,
        max_retained_count=_finite_max(steps, "mpz_retained_count"),
        max_mobile_count=_finite_max(steps, "mpz_mobile_count"),
        max_local_slip_count=_finite_max(steps, "mpz_local_slip_count"),
        max_emitted_ledger=_finite_max(
            steps,
            "mpz_emitted_total" if "mpz_emitted_total" in steps.columns else "N_em",
        ),
        full_field_image=str(image) if image.exists() else None,
        full_field_image_present=image.exists(),
        response_classification=classification,
        stable_resistance_sequence=stable,
        case_gate_passed=not failures,
        case_gate_failures=failures,
    )


def audit_campaign(
    campaign_root: str | Path,
    seed: int,
    T_K: float,
    classes: Iterable[str] = CLASSES,
    bulk_mode: str = "tip_only",
    theta_deg: float = 45.0,
) -> dict[str, Any]:
    root = Path(campaign_root)
    class_list = [str(value) for value in classes]
    theta_token = f"{float(theta_deg):g}"
    case_dirs = {
        cls: root
        / f"seed_{int(seed)}"
        / bulk_mode
        / cls
        / f"T{int(round(T_K))}_th{theta_token}"
        for cls in class_list
    }
    cases = [audit_case(case_dirs[cls], cls, T_K) for cls in class_list]
    pairs: list[dict[str, Any]] = []
    for i, a in enumerate(class_list):
        for b in class_list[i + 1 :]:
            corr = normalized_shape_correlation(case_dirs[a], case_dirs[b])
            same_path = paths_identical(case_dirs[a], case_dirs[b], T_K)
            pairs.append(
                {
                    "class_a": a,
                    "class_b": b,
                    "normalized_raw_shape_correlation": corr,
                    "crack_path_identical": same_path,
                    "geometry_dominated_similarity": bool(
                        same_path and np.isfinite(corr) and corr >= 0.995
                    ),
                }
            )
    geometry_pairs = [row for row in pairs if row["geometry_dominated_similarity"]]
    failed_cases = [case for case in cases if not case.case_gate_passed]
    expected_pairs = len(class_list) * (len(class_list) - 1) // 2
    pairwise_evidence_present = expected_pairs > 0 and len(pairs) == expected_pairs

    failures: list[str] = []
    failures.extend(
        f"case:{case.material_class}:{reason}"
        for case in failed_cases
        for reason in case.case_gate_failures
    )
    if not pairwise_evidence_present:
        failures.append("pairwise_material_comparison_is_vacuous_or_incomplete")
    if geometry_pairs:
        failures.append("geometry_dominated_pair_detected")

    passed = not failures
    if passed:
        interpretation = "complete_stable_material_sequences_with_no_geometry_dominance"
    elif failed_cases:
        interpretation = "incomplete_or_nonpublishable_solver_cases"
    elif not pairwise_evidence_present:
        interpretation = "insufficient_material_classes_for_pairwise_transfer_gate"
    else:
        interpretation = "geometry_or_continuation_dominated_do_not_publish_as_material_R_curves"

    payload = {
        "schema": SCHEMA,
        "campaign_root": str(root),
        "seed": int(seed),
        "T_K": float(T_K),
        "theta_deg": float(theta_deg),
        "bulk_mode": bulk_mode,
        "cases": [asdict(case) for case in cases],
        "n_cases": len(cases),
        "n_case_gate_failures": len(failed_cases),
        "failed_case_classes": [case.material_class for case in failed_cases],
        "pairwise_shape_audit": pairs,
        "n_pairwise_comparisons": len(pairs),
        "expected_pairwise_comparisons": expected_pairs,
        "pairwise_evidence_present": pairwise_evidence_present,
        "n_geometry_dominated_pairs": len(geometry_pairs),
        "geometry_dominated_pairs": [
            f"{row['class_a']}:{row['class_b']}" for row in geometry_pairs
        ],
        "material_rcurve_gate_passed": passed,
        "gate_failures": failures,
        "interpretation": interpretation,
    }
    root.mkdir(parents=True, exist_ok=True)
    (root / "material_rcurve_audit_v9121.json").write_text(
        json.dumps(payload, indent=2, default=str)
    )
    pd.DataFrame(payload["cases"]).to_csv(
        root / "material_rcurve_case_audit_v9121.csv", index=False
    )
    pd.DataFrame(pairs).to_csv(
        root / "material_rcurve_pairwise_audit_v9121.csv", index=False
    )
    return payload


__all__ = [
    "SCHEMA",
    "CaseAudit",
    "audit_case",
    "audit_campaign",
    "normalized_shape_correlation",
    "paths_identical",
]
