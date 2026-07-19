from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from arrhenius_fracture.material_rcurve_audit_v9121 import audit_campaign


def _write_case(
    root: Path,
    cls: str,
    k0: float,
    shape,
    *,
    returncode=0,
    control_state="first_passage",
    target_um=20.0,
    final_um=20.0,
    n_load=4,
    cascade_fraction=0.0,
    status="complete",
    with_field=True,
):
    case = root / "seed_1" / "tip_only" / cls / "T700_th45"
    case.mkdir(parents=True, exist_ok=True)
    ext = np.linspace(final_um / 4.0, final_um, 4)
    k = k0 * np.asarray(shape, float)
    pd.DataFrame(
        {
            "step": [1, 2, 3, 4],
            "KJ_Pa_sqrtm": k * 1.0e6,
            "crack_extension_m": ext * 1.0e-6,
            "mpz_K_shield_Pa_sqrt_m": [0.0, 0.01e6, 0.02e6, 0.03e6],
            "mpz_retained_count": [0.0, 1.0, 2.0, 3.0],
            "mpz_mobile_count": [0.0, 2.0, 1.0, 0.0],
            "mpz_local_slip_count": [0.0, 1.0, 2.0, 3.0],
            "mpz_emitted_total": [0.0, 1.0, 2.0, 3.0],
        }
    ).to_csv(case / "steps_0700K.csv", index=False)
    pd.DataFrame(
        {
            "raw_event_id": [1, 2, 3, 4],
            "step": [1, 2, 3, 4],
            "Uapp_m": [1e-4, 1.1e-4, 1.2e-4, 1.3e-4],
            "KJ_MPa_sqrt_m": k,
            "crack_extension_after_um": ext,
            "da_block_um": [final_um / 4.0] * 4,
            "n_fire": [1] * 4,
            "crack_extension_before_um": ext - final_um / 4.0,
        }
    ).to_csv(case / "R_curve_topology_events_raw.csv", index=False)
    pd.DataFrame(
        [
            {
                "n_raw_topology_events": 4,
                "n_independent_load_events": n_load,
                "n_unstable_same_load_cascades": 0,
                "largest_same_load_jump_um": 0.0,
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
    (case / "v9_12_1_case_contract.json").write_text(
        json.dumps(
            {
                "subprocess_returncode": returncode,
                "solver_status": status,
                "requested_target_extension_um": target_um,
                "final_extension_um": final_um,
            }
        )
    )
    pd.DataFrame(
        {"x_m": 0.5e-3 + ext * 1e-6, "y_m": np.zeros(4)}
    ).to_csv(case / "crack_path_700K.csv", index=False)
    if with_field:
        (case / "field_snapshots_700K.png").write_bytes(b"png")
    return case


def _complete_three_class_campaign(root: Path):
    _write_case(root, "ceramic", 10.0, [1.0, 1.1, 1.3, 1.7])
    _write_case(root, "weakT", 20.0, [1.0, 1.4, 1.2, 1.8])
    _write_case(root, "DBTT", 30.0, [1.0, 1.05, 1.4, 1.2])


def test_complete_distinct_three_class_campaign_can_pass(tmp_path):
    _complete_three_class_campaign(tmp_path)
    out = audit_campaign(tmp_path, 1, 700.0)
    assert out["n_pairwise_comparisons"] == 3
    assert out["pairwise_evidence_present"]
    assert out["n_case_gate_failures"] == 0
    assert out["material_rcurve_gate_passed"]


def test_nonzero_returncode_and_right_censor_fail_gate(tmp_path):
    _complete_three_class_campaign(tmp_path)
    _write_case(
        tmp_path,
        "DBTT",
        30.0,
        [1.0, 1.05, 1.4, 1.2],
        returncode=1,
        status="right_censored",
        final_um=5.0,
    )
    out = audit_campaign(tmp_path, 1, 700.0)
    assert not out["material_rcurve_gate_passed"]
    assert "DBTT" in out["failed_case_classes"]
    joined = " ".join(out["gate_failures"])
    assert "solver_subprocess_failed" in joined
    assert "target_extension_not_reached" in joined


def test_missing_first_passage_or_unstable_sequence_fails(tmp_path):
    _complete_three_class_campaign(tmp_path)
    _write_case(
        tmp_path,
        "weakT",
        20.0,
        [1.0, 1.4, 1.2, 1.8],
        control_state="no_first_passage",
        n_load=2,
    )
    out = audit_campaign(tmp_path, 1, 700.0)
    assert not out["material_rcurve_gate_passed"]
    joined = " ".join(out["gate_failures"])
    assert "first_passage_not_observed" in joined
    assert "no_stable_resistance_sequence" in joined


def test_single_class_gate_is_explicitly_nonvacuous(tmp_path):
    _write_case(tmp_path, "ceramic", 10.0, [1.0, 1.1, 1.3, 1.7])
    out = audit_campaign(tmp_path, 1, 700.0, classes=["ceramic"])
    assert out["n_pairwise_comparisons"] == 0
    assert not out["pairwise_evidence_present"]
    assert not out["material_rcurve_gate_passed"]
    assert "pairwise_material_comparison_is_vacuous_or_incomplete" in out["gate_failures"]
