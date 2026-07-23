from __future__ import annotations

from dataclasses import asdict
import csv
import json
from pathlib import Path

import pandas as pd
import pytest

from arrhenius_fracture.emergent_gnd_campaign_v913 import (
    candidate_from_registry_row,
)
from arrhenius_fracture.emergent_gnd_contract_v913 import (
    ACTIVE_CANDIDATE_PARAMETER_FIELDS,
    PERSISTENT_INACTIVE_REGISTRY_FIELDS,
    candidate_feature_record,
    candidate_parameter_fingerprint,
)
from scripts.augment_mpz_v9_12_directional_peak_targets import (
    add_directional_peak_classifications,
    add_trajectory_metrics,
)
from scripts.run_v913_autonomous_dbtt_search import (
    _establish_run_contract,
    _progress_line,
    _progress_payload,
    _select_rows,
    _validate_resumed_payload,
    _write_json_atomic,
)


def _top5_row() -> dict[str, str]:
    path = Path("candidates/v9_13_persistent_sites_top5_registry.csv")
    with path.open(newline="") as stream:
        return next(csv.DictReader(stream))


def test_existing_peak_objective_applies_to_autonomous_K_trajectory():
    table = pd.DataFrame(
        [
            {
                "candidate_id": "peak",
                "y__K25_1d_T700K": 20.0,
                "y__K25_1d_T800K": 30.0,
                "y__K25_1d_T900K": 50.0,
                "y__K25_1d_T1000K": 35.0,
                "y__K25_1d_T1100K": 25.0,
                "y__K25_1d_T1200K": 20.0,
            }
        ]
    )
    add_trajectory_metrics(
        table,
        prefix="y__K25_1d_",
        out_prefix="y__",
        low_max_K=700.0,
        high_min_K=1000.0,
        peak_min_K=800.0,
        peak_max_K=1000.0,
    )
    add_directional_peak_classifications(
        table,
        out_prefix="y__",
        peak_min_K=800.0,
        direction_threshold=5.0,
        peak_threshold=1.0,
    )

    row = table.iloc[0]
    assert row["y__peak_temperature_K"] == 900.0
    assert row["y__peak_rise"] == 30.0
    assert row["y__peak_drop"] == 30.0
    assert row["y__peak_prominence"] == 30.0
    assert bool(row["y__peak_like_1d"]) is True
    assert bool(row["y__direction_correct_1d"]) is False


def test_nested_sobol_prefix_is_balanced_across_selected_peak_parents():
    rows = []
    for family, parent in (
        ("peak", "parent_a"),
        ("peak", "parent_b"),
        ("plateau", "parent_c"),
    ):
        for index in range(4):
            rows.append(
                {
                    "candidate_id": f"{parent}_{index:04d}",
                    "campaign_parent_id": parent,
                    "campaign_parent_family": family,
                }
            )

    selected = _select_rows(
        rows,
        families=("peak",),
        candidate_ids=(),
        per_parent=2,
        parent_offset=1,
    )
    assert [row["candidate_id"] for row in selected] == [
        "parent_a_0001",
        "parent_a_0002",
        "parent_b_0001",
        "parent_b_0002",
    ]


def test_explicit_candidate_selection_preserves_requested_order():
    rows = [
        {"candidate_id": "a", "campaign_parent_family": "peak"},
        {"candidate_id": "b", "campaign_parent_family": "plateau"},
    ]
    selected = _select_rows(
        rows,
        families=("peak",),
        candidate_ids=("b", "a"),
        per_parent=1,
        parent_offset=0,
    )
    assert [row["candidate_id"] for row in selected] == ["b", "a"]


def test_surrogate_features_are_all_and_only_active_candidate_parameters():
    features = candidate_feature_record(_top5_row())
    assert set(features) == {
        f"x_raw__{field}" for field in ACTIVE_CANDIDATE_PARAMETER_FIELDS
    }
    assert not {
        f"x_raw__{field}" for field in PERSISTENT_INACTIVE_REGISTRY_FIELDS
    } & set(features)
    assert "x_raw__Tref_K" in features
    assert "x_raw__cleave_G00_eV" in features
    assert "x_raw__c_blunt" in features


def test_missing_Tref_uses_the_archived_481p33_default_in_fingerprint():
    explicit = _top5_row()
    defaulted = dict(explicit)
    defaulted.pop("Tref_K")
    assert candidate_parameter_fingerprint([explicit]) == (
        candidate_parameter_fingerprint([defaulted])
    )


def test_legacy_refresh_and_recovery_coordinates_are_explicitly_inactive():
    first = _top5_row()
    second = dict(first)
    second.update(
        {
            "source_refresh_length_um": "299.0",
            "recovery_nu0_s": "9.9e99",
            "recovery_H0_eV": "99.0",
            "recovery_activation_entropy_kB": "-99.0",
        }
    )
    parsed_first = candidate_from_registry_row(first)
    parsed_second = candidate_from_registry_row(second)
    assert asdict(parsed_first) == asdict(parsed_second)
    assert parsed_first.source_refresh_length_m == 0.0
    assert parsed_first.recovery_nu0_s == 0.0
    assert parsed_first.recovery_H0_eV == 0.0
    assert parsed_first.recovery_activation_entropy_kB == 0.0


def test_run_contract_refuses_mixed_output_directory(tmp_path: Path):
    case_root = tmp_path / "cases"
    case_root.mkdir()
    path = tmp_path / "run_contract.json"
    first = {"sha256": "first", "contract": {"setting": 1}}
    assert (
        _establish_run_contract(path, first, case_root=case_root)
        == "first"
    )
    assert json.loads(path.read_text()) == first
    with pytest.raises(RuntimeError, match="different autonomous-search contract"):
        _establish_run_contract(
            path,
            {"sha256": "second", "contract": {"setting": 2}},
            case_root=case_root,
        )


def test_resume_payload_requires_exact_contract_case_and_seed():
    payload = {
        "run_contract_sha256": "contract",
        "candidate_id": "candidate",
        "temperature_K": 900.0,
        "seed": 3621,
    }
    _validate_resumed_payload(
        payload,
        candidate_id="candidate",
        temperature_K=900.0,
        contract_sha256="contract",
        loading_map_seed=3621,
    )
    broken = dict(payload, run_contract_sha256="stale")
    with pytest.raises(RuntimeError, match="run_contract_sha256"):
        _validate_resumed_payload(
            broken,
            candidate_id="candidate",
            temperature_K=900.0,
            contract_sha256="contract",
            loading_map_seed=3621,
        )


def test_progress_payload_reports_resume_rate_eta_and_last_case():
    payload = _progress_payload(
        state="running",
        phase="cases",
        started_at_utc="2026-07-23T00:00:00+00:00",
        elapsed_s=1800.0,
        completed_cases=14,
        resumed_cases=4,
        total_cases=100,
        jobs=4,
        contract_sha256="contract",
        last_case={
            "candidate_id": "candidate",
            "temperature_K": 900.0,
        },
    )
    assert payload["newly_completed_cases"] == 10
    assert payload["remaining_cases"] == 86
    assert payload["active_workers_upper_bound"] == 4
    assert payload["progress_fraction"] == pytest.approx(0.14)
    assert payload["new_cases_per_hour"] == pytest.approx(20.0)
    assert payload["eta_s"] == pytest.approx(15480.0)
    assert payload["last_case"]["candidate_id"] == "candidate"
    line = _progress_line(payload)
    assert "completed=14/100" in line
    assert "eta_s=15480.0" in line


def test_progress_json_is_written_atomically(tmp_path: Path):
    path = tmp_path / "autonomous_dbtt_progress.json"
    payload = {"state": "running", "completed_cases": 1}
    _write_json_atomic(path, payload)
    assert json.loads(path.read_text()) == payload
    assert not path.with_suffix(".json.tmp").exists()
