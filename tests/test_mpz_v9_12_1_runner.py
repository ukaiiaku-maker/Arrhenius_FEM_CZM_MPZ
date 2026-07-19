from argparse import Namespace
import json

import pandas as pd

import run_mpz_v9_12_1_tip_only_material_rcurve as runner


def test_shared_temperature_summary_is_captured_before_next_class(monkeypatch, tmp_path):
    args = Namespace(T_K=700.0, crystal_theta_deg=45.0, target_extension_um=20.0)

    def fake_legacy(args, base_seed, class_name, root):
        run_root = root / f"seed_{base_seed}" / "tip_only"
        case = run_root / class_name / "T700_th45"
        case.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(
            [
                {
                    "status": "complete",
                    "final_extension_um": 20.0,
                    "K_init_MPa_sqrt_m": 12.5,
                }
            ]
        ).to_csv(run_root / "rcurve_temperature_summary.csv", index=False)
        (case / "anisotropic_calibrated_tip_first_passage_summary.json").write_text(
            json.dumps({"control_state": "first_passage"})
        )
        return {
            "class": class_name,
            "subprocess_returncode": 0,
            "solver_output_reused": False,
            "status": None,
            "final_extension_um": None,
            "K_init_MPa_sqrt_m": None,
        }

    monkeypatch.setattr(runner, "_legacy_run_case", fake_legacy)
    row = runner.run_case(args, 1, "ceramic", tmp_path)
    case = tmp_path / "seed_1" / "tip_only" / "ceramic" / "T700_th45"
    assert row["status"] == "complete"
    assert row["final_extension_um"] == 20.0
    assert row["K_init_MPa_sqrt_m"] == 12.5
    assert (case / "rcurve_temperature_summary.csv").is_file()
    contract = json.loads((case / "v9_12_1_case_contract.json").read_text())
    assert contract["summary_fresh_for_this_invocation"]
    assert contract["shared_summary_captured_case_locally"]
    assert contract["target_extension_reached"]
    assert contract["subprocess_returncode"] == 0


def test_reused_case_does_not_trust_uncaptured_shared_summary(monkeypatch, tmp_path):
    args = Namespace(T_K=700.0, crystal_theta_deg=45.0, target_extension_um=20.0)

    def fake_reused(args, base_seed, class_name, root):
        run_root = root / f"seed_{base_seed}" / "tip_only"
        case = run_root / class_name / "T700_th45"
        case.mkdir(parents=True, exist_ok=True)
        # This could belong to another material class and must not be imported.
        pd.DataFrame([{"status": "complete", "final_extension_um": 999.0}]).to_csv(
            run_root / "rcurve_temperature_summary.csv", index=False
        )
        return {
            "class": class_name,
            "subprocess_returncode": 0,
            "solver_output_reused": True,
            "status": None,
            "final_extension_um": None,
        }

    monkeypatch.setattr(runner, "_legacy_run_case", fake_reused)
    row = runner.run_case(args, 1, "weakT", tmp_path)
    assert row["status"] is None
    assert row["final_extension_um"] is None
    case = tmp_path / "seed_1" / "tip_only" / "weakT" / "T700_th45"
    contract = json.loads((case / "v9_12_1_case_contract.json").read_text())
    assert contract["summary_source_path"] is None
    assert not contract["summary_fresh_for_this_invocation"]
    assert not contract["target_extension_reached"]


def test_failed_fresh_run_cannot_copy_preexisting_shared_summary(monkeypatch, tmp_path):
    args = Namespace(T_K=700.0, crystal_theta_deg=45.0, target_extension_um=20.0)
    run_root = tmp_path / "seed_1" / "tip_only"
    run_root.mkdir(parents=True)
    shared = run_root / "rcurve_temperature_summary.csv"
    pd.DataFrame(
        [{"status": "complete", "final_extension_um": 999.0}]
    ).to_csv(shared, index=False)

    def fake_failed(args, base_seed, class_name, root):
        case = run_root / class_name / "T700_th45"
        case.mkdir(parents=True, exist_ok=True)
        return {
            "class": class_name,
            "subprocess_returncode": 1,
            "solver_output_reused": False,
            "status": "failed",
            "final_extension_um": None,
        }

    monkeypatch.setattr(runner, "_legacy_run_case", fake_failed)
    row = runner.run_case(args, 1, "DBTT", tmp_path)
    assert row["status"] == "failed"
    assert row["final_extension_um"] is None
    case = run_root / "DBTT" / "T700_th45"
    contract = json.loads((case / "v9_12_1_case_contract.json").read_text())
    assert contract["subprocess_returncode"] == 1
    assert contract["summary_source_path"] is None
    assert not contract["summary_fresh_for_this_invocation"]
    assert not contract["target_extension_reached"]
