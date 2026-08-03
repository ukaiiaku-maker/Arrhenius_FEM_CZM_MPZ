from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from arrhenius_fracture import (
    mode_i_first_passage_v10_0_5_18_4_0_theta0_pf_parity as parity_entry,
)


def _load_compare_module():
    path = Path(__file__).resolve().parents[1] / "scripts" / "compare_v10051840_theta0_pf_parity.py"
    spec = importlib.util.spec_from_file_location("theta0_compare", path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_steps(case: Path, kj_values: list[float]) -> None:
    case.mkdir(parents=True, exist_ok=True)
    n = len(kj_values)
    frame = pd.DataFrame(
        {
            "step": np.arange(1, n + 1),
            "Uapp_m": np.arange(1, n + 1) * 2.0e-7,
            "Ftop_N": np.linspace(1.0e6, 2.0e6, n),
            "KJ_Pa_sqrtm": np.asarray(kj_values) * 1.0e6,
            "B": [0.2, 0.8, 0.0, 0.5, 0.0][:n],
            "N_em": np.linspace(10.0, 1.0, n),
            "crack_extension_m": np.asarray([0.0, 0.0, 5.0, 5.0, 10.0][:n]) * 1.0e-6,
            "n_fire": [0, 0, 1, 0, 1][:n],
        }
    )
    frame.to_csv(case / "steps_1000K.csv", index=False)


def test_entry_contract_records_matched_controls_and_bulk_mismatch(tmp_path: Path):
    kernel = tmp_path / "family.json"
    kernel.write_text("{}\n")
    contract = parity_entry._contract([], kernel, completed=False, error=None)
    assert contract["reference"]["case"] == "T1000K_th0_seed8666"
    assert contract["matched_controls"]["theta_deg"] == 0.0
    assert contract["matched_controls"]["n_stagger"] == 2
    assert contract["intentional_model_form_differences"]["bulk_plasticity_parity_complete"] is False
    assert contract["physics_guards"]["exact_stochastic_event_endpoint_preserved"] is True


def test_entry_rejects_nonzero_theta():
    args = ["--crystal-theta-deg", "30"]
    with pytest.raises(SystemExit, match="requires --crystal-theta-deg=0"):
        parity_entry._require_float(args, "--crystal-theta-deg", 0.0)


def test_comparison_aligns_by_extension_and_detects_slope_sign(tmp_path: Path):
    compare = _load_compare_module()
    pf = tmp_path / "pf"
    fem = tmp_path / "fem"
    _write_steps(pf, [20.0, 25.0, 30.0, 32.0, 40.0])
    _write_steps(fem, [20.0, 25.0, 30.0, 28.0, 24.0])

    pf_data = compare._trajectory(pf, "PF")
    fem_data = compare._trajectory(fem, "FEM")
    aligned = compare._aligned(pf_data, fem_data, 5.0)

    assert list(aligned["projected_extension_um"]) == [0.0, 5.0, 10.0]
    assert pf_data["metrics"]["event_R_curve_linear_slope_MPa_sqrt_m_per_um"] > 0.0
    assert fem_data["metrics"]["event_R_curve_linear_slope_MPa_sqrt_m_per_um"] < 0.0


def test_contract_comparison_reads_reference_and_fem_audits(tmp_path: Path):
    compare = _load_compare_module()
    pf = tmp_path / "pf"
    fem = tmp_path / "fem"
    pf.mkdir(); fem.mkdir()
    (pf / "run_args.json").write_text(
        json.dumps(
            {
                "crystal_theta_deg": 0.0,
                "dU": 2.0e-7,
                "dt": 8.4,
                "n_stagger": 2,
                "tip_h_fine": 1.0e-6,
                "tip_ratio": 1.2,
                "da_phys": 5.0e-6,
                "adaptive_event_target": 0.15,
                "nx": 36,
                "ny": 72,
                "max_fronts": 1,
            }
        )
    )
    (fem / parity_entry.AUDIT_FILE).write_text(
        json.dumps(
            {
                "matched_controls": {
                    "theta_deg": 0.0,
                    "dU_m": 2.0e-7,
                    "dt_s": 8.4,
                    "n_stagger": 2,
                    "tip_h_fine_m": 1.0e-6,
                    "tip_ratio": 1.2,
                    "da_phys_m": 5.0e-6,
                    "adaptive_event_target": 0.15,
                    "nx": 36,
                    "ny": 72,
                    "maximum_fronts": 1,
                },
                "intentional_model_form_differences": {
                    "bulk_plasticity_parity_complete": False
                },
            }
        )
    )
    audit = compare._contract_comparison(pf, fem)
    assert audit["all_listed_numerical_controls_match"] is True
    assert audit["intentional_model_form_differences"]["bulk_plasticity_parity_complete"] is False
