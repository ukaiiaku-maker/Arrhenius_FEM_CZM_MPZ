from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from arrhenius_fracture.material_rcurve_audit_v10057 import audit_campaign
from run_mpz_v10_0_5_7_tip_only_material_rcurve import (
    _authoritative_temperature_summary,
    _persist_repaired_case_summary,
)


def _write_case(
    root: Path,
    cls: str,
    *,
    status: str = "complete",
    returncode: int = 0,
    target_completed: bool = True,
    control_state: str = "first_passage",
    with_field: bool = True,
    y_shift: float = 0.0,
    n_independent_load_events: int = 4,
    cascade_fraction: float = 0.0,
):
    case = root / "seed_1" / "tip_only" / cls / "T700_th45"
    case.mkdir(parents=True, exist_ok=True)
    ext = np.asarray([5.0, 10.0, 15.0, 20.0])
    k0 = {"ceramic": 10.0, "weakT": 20.0, "DBTT": 30.0}[cls]
    shape = {
        "ceramic": np.asarray([1.0, 1.10, 1.22, 1.35]),
        "weakT": np.asarray([1.0, 1.45, 1.18, 1.85]),
        "DBTT": np.asarray([1.0, 1.02, 1.50, 1.12]),
    }[cls]
    k = k0 * shape
    pd.DataFrame(
        {
            "step": [1, 2, 3, 4],
            "crack_extension_m": ext * 1.0e-6,
            "mpz_K_shield_Pa_sqrt_m": [0.0, 0.01e6, 0.02e6, 0.03e6],
            "mpz_retained_count": [0.0, 1.0, 2.0, 3.0],
            "mpz_mobile_count": [0.0, 2.0, 1.0, 0.0],
            "mpz_local_slip_count": [0.0, 1.0, 2.0, 3.0],
            "N_em": [0.0, 1.0, 2.0, 3.0],
        }
    ).to_csv(case / "steps_0700K.csv", index=False)
    pd.DataFrame(
        {
            "raw_event_id": [1, 2, 3, 4],
            "KJ_MPa_sqrt_m": k,
            "crack_extension_after_um": ext,
        }
    ).to_csv(case / "R_curve_topology_events_raw.csv", index=False)
    pd.DataFrame(
        [
            {
                "n_raw_topology_events": 4,
                "n_independent_load_events": n_independent_load_events,
                "n_unstable_same_load_cascades": 1 if cascade_fraction > 0 else 0,
                "largest_same_load_jump_um": 10.0 if cascade_fraction > 0 else 0.0,
                "fraction_topology_events_in_cascades": cascade_fraction,
            }
        ]
    ).to_csv(case / "R_curve_cascade_metrics.csv", index=False)
    (case / "anisotropic_calibrated_tip_first_passage_summary.json").write_text(
        json.dumps(
            {
                "control_state": control_state,
                "Kc_first_existing_MPa_sqrt_m": k0,
            }
        )
    )
    pd.DataFrame(
        {
            "x_m": 0.5e-3 + ext * 1.0e-6,
            "y_m": np.full(4, y_shift),
        }
    ).to_csv(case / "crack_path_700K.csv", index=False)
    if with_field:
        (case / "field_snapshots_700K.png").write_bytes(b"png")
    summary = {
        "class": cls,
        "status": status,
        "subprocess_returncode": returncode,
        "target_completed": target_completed,
        "target_extension_um": 20.0,
        "final_extension_um": 20.0 if target_completed else 5.0,
        "control_state": control_state,
        "K_init_MPa_sqrt_m": k0 if control_state == "first_passage" else None,
    }
    (case / "v9_12_case_summary.json").write_text(json.dumps(summary))
    return case


def test_root_level_temperature_summary_is_selected_and_copied(tmp_path):
    run_root = tmp_path / "seed_1" / "tip_only"
    run_root.mkdir(parents=True)
    pd.DataFrame(
        [
            {
                "class": "DBTT",
                "T_K": 700,
                "status": "complete",
                "target_completed": True,
                "control_state": "first_passage",
            }
        ]
    ).to_csv(run_root / "rcurve_temperature_summary.csv", index=False)
    solver_row, source = _authoritative_temperature_summary(run_root, "DBTT", 700)
    assert source == run_root / "rcurve_temperature_summary.csv"
    assert solver_row["status"] == "complete"

    case = run_root / "DBTT" / "T700_th45"
    case.mkdir(parents=True)
    repaired = _persist_repaired_case_summary(
        legacy_row={
            "class": "DBTT",
            "subprocess_returncode": 0,
            "case_dir": str(case),
        },
        solver_row=solver_row,
        source=source,
        case_dir=case,
    )
    assert repaired["target_completed"] is True
    assert (case / "rcurve_temperature_summary_v9_11.csv").exists()


def test_complete_multimaterial_campaign_passes_strict_gate(tmp_path):
    _write_case(tmp_path, "ceramic", y_shift=0.0)
    _write_case(tmp_path, "weakT", y_shift=1.0e-7)
    _write_case(tmp_path, "DBTT", y_shift=-1.0e-7)
    out = audit_campaign(tmp_path, 1, 700.0)
    assert out["n_pairwise_comparisons"] == 3
    assert out["all_case_solver_gates_passed"]
    assert out["all_case_publication_gates_passed"]
    assert out["material_rcurve_gate_passed"]
    assert out["interpretation"] == "material_transfer_gate_passed"


def test_right_censored_case_fails_publication_gate(tmp_path):
    _write_case(tmp_path, "ceramic", y_shift=0.0)
    _write_case(
        tmp_path,
        "weakT",
        status="right_censored",
        returncode=1,
        target_completed=False,
        y_shift=1.0e-7,
    )
    _write_case(tmp_path, "DBTT", y_shift=-1.0e-7)
    out = audit_campaign(tmp_path, 1, 700.0)
    assert "weakT" in out["failed_solver_cases"]
    assert "weakT" in out["incomplete_or_censored_cases"]
    assert not out["material_rcurve_gate_passed"]
    assert out["interpretation"].startswith("solver_failure")


def test_non_first_passage_case_fails_publication_gate(tmp_path):
    _write_case(tmp_path, "ceramic", y_shift=0.0)
    _write_case(
        tmp_path,
        "weakT",
        control_state="right_censored",
        y_shift=1.0e-7,
    )
    _write_case(tmp_path, "DBTT", y_shift=-1.0e-7)
    out = audit_campaign(tmp_path, 1, 700.0)
    assert "weakT" in out["non_first_passage_cases"]
    assert not out["material_rcurve_gate_passed"]


def test_unstable_fixed_displacement_case_fails_publication_gate(tmp_path):
    _write_case(
        tmp_path,
        "ceramic",
        y_shift=0.0,
        n_independent_load_events=3,
        cascade_fraction=0.75,
    )
    _write_case(tmp_path, "weakT", y_shift=1.0e-7)
    _write_case(tmp_path, "DBTT", y_shift=-1.0e-7)
    out = audit_campaign(tmp_path, 1, 700.0)
    assert out["all_case_solver_gates_passed"]
    assert not out["all_case_publication_gates_passed"]
    assert "ceramic" in out["unstable_response_cases"]
    assert not out["material_rcurve_gate_passed"]
    assert out["interpretation"].startswith("unstable_fixed_displacement")


def test_single_material_campaign_cannot_pass_transfer_gate(tmp_path):
    _write_case(tmp_path, "ceramic")
    out = audit_campaign(tmp_path, 1, 700.0, classes=["ceramic"])
    assert out["n_pairwise_comparisons"] == 0
    assert not out["pairwise_comparison_sufficient"]
    assert not out["material_rcurve_gate_passed"]
    assert out["interpretation"].startswith("insufficient_material_comparisons")
