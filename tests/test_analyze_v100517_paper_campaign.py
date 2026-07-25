from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pandas as pd


def _module():
    path = Path(__file__).resolve().parents[1] / "scripts" / "analyze_v100517_paper_parameter_campaign.py"
    spec = importlib.util.spec_from_file_location("analyze_v100517", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _case(root: Path, option: str, temperature: int, seed: int, shift: float) -> None:
    case = root / option / f"T{temperature:04d}K"
    case.mkdir(parents=True)
    pd.DataFrame(
        {
            "step": [10, 20],
            "n_fire": [1, 1],
            "KJ_Pa_sqrtm": [(20.0 + shift) * 1.0e6, (21.0 + shift) * 1.0e6],
            "crack_extension_m": [1.0e-6, 2.1e-6],
        }
    ).to_csv(case / f"steps_{temperature:04d}K.csv", index=False)
    pd.DataFrame(
        {
            "step": [10, 20],
            "front_id": [0, 0],
            "n_fire": [1, 1],
            "J_signed_trial": [2000.0 + 10.0 * shift, 2400.0 + 10.0 * shift],
            "J_effective_trial": [2000.0 + 10.0 * shift, 2400.0 + 10.0 * shift],
            "x_m": [0.501e-3, 0.5021e-3],
            "y_m": [0.5e-6, -0.1e-6],
        }
    ).to_csv(case / f"fronts_{temperature:04d}K.csv", index=False)
    pd.DataFrame(
        {
            "step": [10, 20],
            "angle1_deg": [30.0, -60.0],
        }
    ).to_csv(case / f"branch_diagnostics_{temperature:04d}K.csv", index=False)
    pd.DataFrame(
        {
            "x_m": [0.500e-3, 0.501e-3, 0.5021e-3],
            "y_m": [0.0, 0.5e-6, -0.1e-6],
        }
    ).to_csv(case / f"crack_path_{temperature}K.csv", index=False)
    events = [
        {
            "inserted": True,
            "moved_m": 1.0e-6,
            "threshold_action": 0.7 + shift,
            "hazard_seed": seed,
            "hazard_event_index": 0,
        },
        {
            "inserted": True,
            "moved_m": 1.1e-6,
            "threshold_action": 1.4 + shift,
            "hazard_seed": seed,
            "hazard_event_index": 1,
        },
    ]
    (case / "stochastic_geometry_events_v10_0_5_16.json").write_text(json.dumps(events))
    (case / "persistent_site_production_manifest_v10_0_5_17.json").write_text(
        json.dumps(
            {
                "status": "complete",
                "run_completed_without_exception": True,
                "physics_contract": {"cleavage_hazard_seed": seed},
            }
        )
    )
    (case / "audited_PF_paper_parameter_transfer_v10_0_5_17.json").write_text(
        json.dumps(
            {
                "parameter_option": option,
                "candidate_id": "candidate_test",
                "material_class": "peak",
                "paper_role": "paper primary",
                "parameter_entry": "v10.2.25",
                "selected_row_sha256": f"row-{temperature}",
            }
        )
    )
    (case / f"mpz_state_snapshots_{temperature:04d}K.json").write_text(
        json.dumps({"snapshots": [{"step": 1}, {"step": 2}]})
    )


def test_analyzer_writes_K_J_and_path_outputs(tmp_path):
    module = _module()
    campaign = tmp_path / "campaign"
    option = "v913_paper_peak01_test_persistent_sites"
    _case(campaign, option, 300, 1001, 0.0)
    _case(campaign, option, 400, 1002, 0.2)
    pd.DataFrame(
        [
            ["v10.2.25", option, 300, 1001],
            ["v10.2.25", option, 400, 1002],
        ]
    ).to_csv(campaign / "case_matrix.tsv", sep="\t", index=False, header=False)

    out = tmp_path / "analysis"
    payload = module.analyze(campaign, out, target_um=2.0)
    assert payload["n_parameterizations"] == 1
    assert payload["n_cases"] == 2
    assert payload["analysis_complete"] is True
    assert (out / "K_J_temperature_path_seed_completion_summary.csv").is_file()
    assert (out / option / "K_R_curves_all_temperatures.png").is_file()
    assert (out / option / "J_R_curves_all_temperatures.pdf").is_file()
    assert (out / option / "K_vs_temperature_initial_mean_final.png").is_file()
    assert (out / option / "J_vs_temperature_initial_mean_final.pdf").is_file()
    assert (out / option / "crack_paths_all_temperatures.png").is_file()
    assert (out / option / "committed_crack_direction_vs_extension.pdf").is_file()
    summary = pd.read_csv(out / "K_J_temperature_path_seed_completion_summary.csv")
    assert summary["J_initial_kJ_m2"].iloc[0] == 2.0
    assert summary["path_n_direction_changes_gt_1deg"].min() >= 1


def test_partial_analysis_skips_unfinished_case(tmp_path):
    module = _module()
    campaign = tmp_path / "campaign"
    option = "v913_paper_peak01_test_persistent_sites"
    _case(campaign, option, 300, 1001, 0.0)
    (campaign / option / "T0400K").mkdir(parents=True)
    pd.DataFrame(
        [
            ["v10.2.25", option, 300, 1001],
            ["v10.2.25", option, 400, 1002],
        ]
    ).to_csv(campaign / "case_matrix.tsv", sep="\t", index=False, header=False)
    out = tmp_path / "analysis"
    payload = module.analyze(campaign, out, target_um=2.0, allow_partial=True)
    assert payload["analysis_complete"] is False
    assert payload["n_analyzed_cases"] == 1
    assert payload["n_skipped_cases"] == 1
    assert (out / option / "K_R_curves_all_temperatures.png").is_file()
