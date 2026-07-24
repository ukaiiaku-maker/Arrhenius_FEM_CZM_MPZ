from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pandas as pd
import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "analyze_v100516_stochastic_campaign.py"
SPEC = importlib.util.spec_from_file_location("analyze_v100516", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def write_case(root: Path, temperature: int, seed: int, lengths_um: list[float], k_values: list[float]):
    case = root / f"T{temperature:04d}K"
    case.mkdir(parents=True)
    extension = 0.0
    rows = []
    events = []
    for index, (length_um, kval) in enumerate(zip(lengths_um, k_values)):
        extension += length_um
        rows.append(
            {
                "step": index + 1,
                "n_fire": 1,
                "KJ_Pa_sqrtm": kval * 1.0e6,
                "crack_extension_m": extension * 1.0e-6,
            }
        )
        events.append(
            {
                "inserted": True,
                "moved_m": length_um * 1.0e-6,
                "threshold_action": 0.5 + index + temperature / 10000.0,
                "hazard_seed": seed,
                "hazard_event_index": index,
            }
        )
    pd.DataFrame(rows).to_csv(case / f"steps_{temperature:04d}K.csv", index=False)
    (case / "stochastic_geometry_events_v10_0_5_16.json").write_text(
        json.dumps(events)
    )
    (case / "persistent_site_production_manifest_v10_0_5_16.json").write_text(
        json.dumps(
            {
                "status": "complete",
                "run_completed_without_exception": True,
                "physics_contract": {"cleavage_hazard_seed": seed},
            }
        )
    )
    (case / f"mpz_state_snapshots_{temperature:04d}K.json").write_text(
        json.dumps({"snapshots": [{}, {}, {}]})
    )


def test_campaign_analysis_writes_rcurves_temperature_metrics_and_seed_audit(tmp_path):
    root = tmp_path / "campaign"
    write_case(root, 300, 101, [3.0, 7.0], [10.0, 14.0])
    write_case(root, 400, 202, [4.0, 6.0], [20.0, 25.0])
    out = root / "final_analysis"
    audit = MODULE.analyze_campaign(root, out, target_um=10.0)

    assert audit["unique_temperature_seeds"] is True
    assert audit["unique_event_sequence_fingerprints"] is True
    assert audit["all_cases_reached_target"] is True
    summary = pd.read_csv(out / "K_temperature_and_seed_summary.csv")
    row300 = summary.loc[summary["temperature_K"] == 300].iloc[0]
    assert row300["K_initial_MPa_sqrt_m"] == pytest.approx(10.0)
    assert row300["K_mean_MPa_sqrt_m"] == pytest.approx(12.0)
    assert row300["K_max_MPa_sqrt_m"] == pytest.approx(14.0)
    assert row300["K_end_MPa_sqrt_m"] == pytest.approx(14.0)
    assert row300["state_snapshot_records"] == 3

    for filename in (
        "R_curves_all_temperatures.png",
        "R_curves_all_temperatures.pdf",
        "K_vs_temperature_initial_mean_max_end.png",
        "K_vs_temperature_initial_mean_max_end.pdf",
        "event_advance_vs_crack_extension.png",
        "analysis_audit.json",
        "temperature_seed_map.csv",
    ):
        assert (out / filename).is_file()
    assert (out / "individual_R_curves" / "R_curve_0300K.png").is_file()


def test_campaign_analysis_rejects_duplicate_temperature_seeds(tmp_path):
    root = tmp_path / "campaign"
    write_case(root, 300, 101, [3.0, 7.0], [10.0, 14.0])
    write_case(root, 400, 101, [4.0, 6.0], [20.0, 25.0])
    with pytest.raises(RuntimeError, match="duplicate temperature seeds"):
        MODULE.analyze_campaign(root, root / "final_analysis", target_um=10.0)
