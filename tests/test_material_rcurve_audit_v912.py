from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from arrhenius_fracture.material_rcurve_audit_v912 import (
    audit_campaign,
    normalized_shape_correlation,
)


def _write_case(root: Path, cls: str, k0: float, shape, with_field: bool = True):
    case = root / "seed_1" / "tip_only" / cls / "T700_th45"
    case.mkdir(parents=True, exist_ok=True)
    ext = np.asarray([5.0, 10.0, 15.0, 20.0])
    k = k0 * np.asarray(shape, float)
    pd.DataFrame({
        "step": [1, 2, 3, 4],
        "Uapp_m": [1e-4, 1.1e-4, 1.2e-4, 1.3e-4],
        "KJ_Pa_sqrtm": k * 1e6,
        "a_tip_m": 0.5e-3 + ext * 1e-6,
        "crack_extension_m": ext * 1e-6,
        "n_fire": [1, 1, 1, 1],
        "mpz_K_shield_Pa_sqrt_m": [0.0, 0.01e6, 0.02e6, 0.03e6],
        "mpz_retained_count": [0.0, 1.0, 2.0, 3.0],
        "mpz_mobile_count": [0.0, 2.0, 1.0, 0.0],
        "mpz_local_slip_count": [0.0, 1.0, 2.0, 3.0],
        "N_em": [0.0, 1.0, 2.0, 3.0],
    }).to_csv(case / "steps_0700K.csv", index=False)
    pd.DataFrame({
        "raw_event_id": [1, 2, 3, 4],
        "step": [1, 2, 3, 4],
        "Uapp_m": [1e-4, 1.1e-4, 1.2e-4, 1.3e-4],
        "KJ_MPa_sqrt_m": k,
        "crack_extension_after_um": ext,
        "da_block_um": [5.0] * 4,
        "n_fire": [1] * 4,
        "crack_extension_before_um": ext - 5.0,
    }).to_csv(case / "R_curve_topology_events_raw.csv", index=False)
    pd.DataFrame([{
        "n_raw_topology_events": 4,
        "n_independent_load_events": 4,
        "n_unstable_same_load_cascades": 0,
        "largest_same_load_jump_um": 0.0,
        "fraction_topology_events_in_cascades": 0.0,
        "rcurve_interpretation": "independent_reload_events",
    }]).to_csv(case / "R_curve_cascade_metrics.csv", index=False)
    (case / "anisotropic_calibrated_tip_first_passage_summary.json").write_text(
        json.dumps({
            "control_state": "first_passage",
            "Kc_first_existing_MPa_sqrt_m": k0,
        })
    )
    pd.DataFrame({"x_m": 0.5e-3 + ext * 1e-6, "y_m": np.zeros(4)}).to_csv(
        case / "crack_path_700K.csv", index=False
    )
    if with_field:
        (case / "field_snapshots_700K.png").write_bytes(b"png")
    return case


def test_identical_normalized_shapes_are_flagged_as_geometry_dominated(tmp_path):
    a = _write_case(tmp_path, "ceramic", 10.0, [1.0, 1.1, 1.2, 1.3])
    b = _write_case(tmp_path, "weakT", 20.0, [1.0, 1.1, 1.2, 1.3])
    _write_case(tmp_path, "DBTT", 30.0, [1.0, 1.1, 1.2, 1.3])
    assert normalized_shape_correlation(a, b) > 0.999999
    out = audit_campaign(tmp_path, 1, 700.0)
    assert out["n_geometry_dominated_pairs"] == 3
    assert not out["material_rcurve_gate_passed"]
    assert out["interpretation"].startswith("geometry_or_continuation")


def test_missing_full_field_image_fails_gate(tmp_path):
    _write_case(tmp_path, "ceramic", 10.0, [1.0, 1.1, 1.3, 1.7])
    _write_case(tmp_path, "weakT", 20.0, [1.0, 1.4, 1.2, 1.8])
    _write_case(tmp_path, "DBTT", 30.0, [1.0, 1.05, 1.1, 1.15], with_field=False)
    out = audit_campaign(tmp_path, 1, 700.0)
    assert "DBTT" in out["missing_full_field_images"]
    assert not out["material_rcurve_gate_passed"]
