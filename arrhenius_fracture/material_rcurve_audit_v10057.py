"""Strict publication gate for v9.12 full-field material R-curve campaigns.

This v10.0.5.7 layer preserves the v9.12 mechanics and constitutive response but
prevents failed, right-censored, non-first-passage, unstable-cascade,
missing-field, or vacuous single-material campaigns from passing the
material-transfer publication gate.
"""
from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

from . import material_rcurve_audit_v912 as _base

POINT_RELEASE = "10.0.5.7"
CLASSES = _base.CLASSES
STABLE_RESPONSE_CLASS = "candidate_stable_resistance_sequence"


def _as_bool(value: Any, default: bool = False) -> bool:
    if value is None:
        return bool(default)
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        if not np.isfinite(float(value)):
            return bool(default)
        return bool(int(value))
    token = str(value).strip().lower()
    if token in {"true", "1", "yes", "y", "complete", "completed"}:
        return True
    if token in {"false", "0", "no", "n", "", "none", "nan"}:
        return False
    return bool(default)


def _as_int(value: Any, default: int = -1) -> int:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return int(default)
    if not np.isfinite(number):
        return int(default)
    return int(number)


def _finite(value: Any) -> bool:
    try:
        return bool(np.isfinite(float(value)))
    except (TypeError, ValueError):
        return False


def audit_case_strict(case_dir: str | Path, material_class: str, T_K: float) -> dict[str, Any]:
    """Return physical, lifecycle, and response-stability audit state."""
    root = Path(case_dir)
    physical = asdict(_base.audit_case(root, material_class, T_K))
    summary_path = root / "v9_12_case_summary.json"
    summary = _base._read_json(summary_path)
    config = _base._read_json(root / "v9_12_run_config.json")

    subprocess_returncode = _as_int(
        summary.get("subprocess_returncode", summary.get("returncode")), -1
    )
    solver_status = str(summary.get("status", "unknown")).strip().lower()
    target_completed = _as_bool(
        summary.get("target_completed"),
        default=(root / ".long_growth_complete").exists(),
    )
    control_state = str(
        summary.get("control_state", physical.get("control_state", "unknown"))
    ).strip().lower()
    target_extension_um = summary.get(
        "target_extension_um", config.get("target_extension_um", np.nan)
    )
    final_extension_um = summary.get(
        "final_extension_um", physical.get("final_extension_um", np.nan)
    )
    first_passage_observed = bool(
        control_state == "first_passage"
        and _finite(summary.get("K_init_MPa_sqrt_m", physical.get("K_init_MPa_sqrt_m")))
    )
    successful_status = solver_status in {"complete", "skipped_complete"}
    case_summary_present = bool(summary_path.exists() and summary)
    response_classification = str(
        physical.get("response_classification", "unknown")
    ).strip().lower()
    response_gate_passed = response_classification == STABLE_RESPONSE_CLASS

    lifecycle_reasons: list[str] = []
    if not case_summary_present:
        lifecycle_reasons.append("missing_v9_12_case_summary")
    if subprocess_returncode != 0:
        lifecycle_reasons.append("subprocess_returncode_nonzero")
    if not successful_status:
        lifecycle_reasons.append("solver_status_not_complete")
    if not target_completed:
        lifecycle_reasons.append("target_extension_not_completed")
    if not first_passage_observed:
        lifecycle_reasons.append("first_passage_not_observed")
    if not bool(physical.get("full_field_image_present", False)):
        lifecycle_reasons.append("missing_full_field_image")

    response_reasons: list[str] = []
    if not response_gate_passed:
        response_reasons.append(
            f"response_not_{STABLE_RESPONSE_CLASS}:{response_classification}"
        )

    out = dict(physical)
    out.update(
        {
            "point_release": POINT_RELEASE,
            "case_summary_present": case_summary_present,
            "subprocess_returncode": subprocess_returncode,
            "solver_status": solver_status,
            "target_extension_um": target_extension_um,
            "final_extension_um": final_extension_um,
            "target_completed": target_completed,
            "control_state": control_state,
            "first_passage_observed": first_passage_observed,
            "solver_gate_passed": not lifecycle_reasons,
            "solver_gate_failure_reasons": ";".join(lifecycle_reasons),
            "response_gate_passed": response_gate_passed,
            "response_gate_failure_reasons": ";".join(response_reasons),
            "case_publication_gate_passed": bool(
                not lifecycle_reasons and response_gate_passed
            ),
        }
    )
    return out


def _interpretation(
    *,
    failed_cases: list[str],
    incomplete_cases: list[str],
    non_first_passage_cases: list[str],
    unstable_response_cases: list[str],
    missing_fields: list[str],
    pairwise_sufficient: bool,
    geometry_pairs: list[dict[str, Any]],
) -> str:
    if failed_cases:
        return "solver_failure_do_not_publish_as_material_R_curves"
    if incomplete_cases:
        return "right_censored_or_incomplete_do_not_publish_as_material_R_curves"
    if non_first_passage_cases:
        return "first_passage_not_observed_do_not_publish_as_material_R_curves"
    if unstable_response_cases:
        return "unstable_fixed_displacement_propagation_do_not_publish_as_material_R_curves"
    if missing_fields:
        return "missing_required_full_field_outputs"
    if not pairwise_sufficient:
        return "insufficient_material_comparisons_for_transfer_gate"
    if geometry_pairs:
        return "geometry_or_continuation_dominated_do_not_publish_as_material_R_curves"
    return "material_transfer_gate_passed"


def audit_campaign(
    campaign_root: str | Path,
    seed: int,
    T_K: float,
    classes: Iterable[str] = CLASSES,
    bulk_mode: str = "tip_only",
) -> dict[str, Any]:
    """Run a fail-closed material-transfer audit.

    Passing requires every case to have a successful subprocess, completed target
    extension, observed first passage, required full-field image, and a stable
    independent-load resistance sequence rather than a serialized same-load
    cascade. At least one pairwise material comparison must exist, and no pair may
    be geometry dominated.
    """
    root = Path(campaign_root)
    class_list = [str(x) for x in classes]
    case_dirs = {
        cls: root / f"seed_{int(seed)}" / bulk_mode / cls / f"T{int(round(T_K))}_th45"
        for cls in class_list
    }
    cases = [audit_case_strict(case_dirs[cls], cls, T_K) for cls in class_list]

    pairs: list[dict[str, Any]] = []
    for i, a in enumerate(class_list):
        for b in class_list[i + 1 :]:
            corr = _base.normalized_shape_correlation(case_dirs[a], case_dirs[b])
            same_path = _base.paths_identical(case_dirs[a], case_dirs[b], T_K)
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

    geometry_pairs = [p for p in pairs if p["geometry_dominated_similarity"]]
    missing_fields = [
        str(c["material_class"]) for c in cases if not c["full_field_image_present"]
    ]
    failed_cases = [
        str(c["material_class"])
        for c in cases
        if int(c["subprocess_returncode"]) != 0
        or str(c["solver_status"]) == "failed"
    ]
    incomplete_cases = [
        str(c["material_class"])
        for c in cases
        if not bool(c["target_completed"])
        or str(c["solver_status"]) not in {"complete", "skipped_complete"}
    ]
    non_first_passage_cases = [
        str(c["material_class"]) for c in cases if not c["first_passage_observed"]
    ]
    unstable_response_cases = [
        str(c["material_class"]) for c in cases if not c["response_gate_passed"]
    ]
    missing_case_summaries = [
        str(c["material_class"]) for c in cases if not c["case_summary_present"]
    ]
    n_pairs = len(pairs)
    pairwise_sufficient = n_pairs > 0
    all_solver_gates_pass = bool(cases) and all(
        bool(c["solver_gate_passed"]) for c in cases
    )
    all_publication_gates_pass = bool(cases) and all(
        bool(c["case_publication_gate_passed"]) for c in cases
    )
    gate_passed = bool(
        all_publication_gates_pass
        and pairwise_sufficient
        and not geometry_pairs
        and not missing_fields
    )
    interpretation = _interpretation(
        failed_cases=failed_cases,
        incomplete_cases=incomplete_cases,
        non_first_passage_cases=non_first_passage_cases,
        unstable_response_cases=unstable_response_cases,
        missing_fields=missing_fields,
        pairwise_sufficient=pairwise_sufficient,
        geometry_pairs=geometry_pairs,
    )

    payload = {
        "schema": "material_rcurve_audit_v10_0_5_7",
        "point_release": POINT_RELEASE,
        "campaign_root": str(root),
        "seed": int(seed),
        "T_K": float(T_K),
        "bulk_mode": bulk_mode,
        "cases": cases,
        "pairwise_shape_audit": pairs,
        "n_pairwise_comparisons": n_pairs,
        "pairwise_comparison_sufficient": pairwise_sufficient,
        "n_geometry_dominated_pairs": len(geometry_pairs),
        "geometry_dominated_pairs": [
            f"{p['class_a']}:{p['class_b']}" for p in geometry_pairs
        ],
        "failed_solver_cases": sorted(set(failed_cases)),
        "incomplete_or_censored_cases": sorted(set(incomplete_cases)),
        "non_first_passage_cases": sorted(set(non_first_passage_cases)),
        "unstable_response_cases": sorted(set(unstable_response_cases)),
        "missing_case_summaries": sorted(set(missing_case_summaries)),
        "missing_full_field_images": sorted(set(missing_fields)),
        "all_case_solver_gates_passed": all_solver_gates_pass,
        "all_case_publication_gates_passed": all_publication_gates_pass,
        "material_rcurve_gate_passed": gate_passed,
        "interpretation": interpretation,
        "constitutive_physics_changed": False,
    }
    root.mkdir(parents=True, exist_ok=True)
    (root / "material_rcurve_audit_v10_0_5_7.json").write_text(
        json.dumps(payload, indent=2, default=str)
    )
    pd.DataFrame(cases).to_csv(
        root / "material_rcurve_case_audit_v10_0_5_7.csv", index=False
    )
    pd.DataFrame(pairs).to_csv(
        root / "material_rcurve_pairwise_audit_v10_0_5_7.csv", index=False
    )
    return payload


__all__ = [
    "POINT_RELEASE",
    "STABLE_RESPONSE_CLASS",
    "audit_case_strict",
    "audit_campaign",
]
