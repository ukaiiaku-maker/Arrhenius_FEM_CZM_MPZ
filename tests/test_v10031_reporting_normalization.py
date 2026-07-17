from __future__ import annotations

import csv
import json
from pathlib import Path

from normalize_v10_0_3_1_reporting import normalize


def test_reporting_normalization_replaces_stale_legacy_fields(tmp_path: Path):
    (tmp_path / "mode_i_v10_0_3_results.json").write_text(json.dumps([
        {
            "T_K": 700.0,
            "model": "FEM_CZM_Mode_I_kinetic_campaign_czm_v10_0_3_delayed_live_binding_progressive_event_lifecycle",
            "point_release": "10.0.3",
            "front_state_model": "kinetic_campaign_czm",
            "front_state_model_detail": "pf_v10_1_7_1_campaign_calibrated_continuous_tip_reset_safe_v1003",
            "B_final": 0.0,
            "crack_extension_final_m": 5.0e-6,
            "max_N_em": 3.3,
            "source_budget_total": 4.88,
            "source_population_bound": 5.74,
            "progressive_runtime_audit": "kinetic_campaign_czm_progressive_2d_v10_0_3.json",
        }
    ]))
    (tmp_path / "kinetic_campaign_czm_v10_0_3_audit.json").write_text(json.dumps({
        "full_progressive_trial_loop_active": True,
        "live_binding_capture_verified": True,
        "runtime": {
            "source_budget_total": 4.88,
        },
        "result_checks": [{"T_K": 700.0, "B_final": 0.0}],
    }))
    (tmp_path / "anisotropic_calibrated_tip_first_passage_summary.json").write_text(json.dumps({
        "T_K": 700.0,
        "model": "FEM_CZM_mixed_mode_MPZ_v9_11_independent_PT_non_double_counted_shielding",
        "B_final": None,
        "front_state_model_detail": "moving_pz_v911_independent_shapes_2d_profile",
    }))

    audit = normalize(tmp_path)

    normalized = json.loads(
        (tmp_path / "anisotropic_calibrated_tip_first_passage_summary.json").read_text()
    )
    assert normalized["B_final"] == 0.0
    assert normalized["front_state_model"] == "kinetic_campaign_czm"
    assert normalized["front_state_model_detail"].endswith("reset_safe_v1003")
    assert normalized["legacy_wrapper_model"].startswith("FEM_CZM_mixed_mode_MPZ_v9_11")
    assert normalized["source_population_bound"] == 5.74
    assert normalized["reporting_normalization_physics_changed"] is False
    assert audit["physics_recomputed"] is False
    assert audit["after"]["source_population_bound"] == 5.74

    with (tmp_path / "anisotropic_calibrated_tip_first_passage_summary.csv").open() as fp:
        rows = list(csv.DictReader(fp))
    assert len(rows) == 1
    assert rows[0]["front_state_model"] == "kinetic_campaign_czm"
    assert rows[0]["B_final"] == "0.0"
    assert rows[0]["source_population_bound"] == "5.74"

    normalized_results = json.loads(
        (tmp_path / "mode_i_v10_0_3_1_results.json").read_text()
    )
    assert normalized_results[0]["reporting_point_release"] == "10.0.3.1"
    assert normalized_results[0]["reporting_normalization_physics_changed"] is False
