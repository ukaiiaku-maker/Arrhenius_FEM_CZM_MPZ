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
            "n_fire": [1, 1],
            "KJ_Pa_sqrtm": [(20.0 + shift) * 1.0e6, (21.0 + shift) * 1.0e6],
            "crack_extension_m": [1.0e-6, 2.1e-6],
        }
    ).to_csv(case / f"steps_{temperature:04d}K.csv", index=False)
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


def test_analyzer_writes_parameter_and_cross_parameter_outputs(tmp_path):
    module = _module()
    campaign = tmp_path / "campaign"
    option = "v913_paper_peak01_test_persistent_sites"
    _case(campaign, option, 300, 1001, 0.0)
    _case(campaign, option, 400, 1002, 0.2)
    out = tmp_path / "analysis"
    payload = module.analyze(campaign, out, target_um=2.0)
    assert payload["n_parameterizations"] == 1
    assert payload["n_cases"] == 2
    assert payload["all_seeds_unique"] is True
    assert payload["all_event_sequences_unique"] is True
    assert (out / "K_temperature_seed_and_completion_summary.csv").is_file()
    assert (out / option / "R_curves_all_temperatures.png").is_file()
    assert (out / option / "K_vs_temperature_initial_mean_max_end.pdf").is_file()
    assert (out / "all_parameterizations_K_initial_vs_temperature.png").is_file()
