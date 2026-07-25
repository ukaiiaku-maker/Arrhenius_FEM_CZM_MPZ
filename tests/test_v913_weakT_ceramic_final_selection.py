from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd
import pytest

from arrhenius_fracture.emergent_gnd_contract_v913 import ACTIVE_CANDIDATE_PARAMETER_FIELDS


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "finalize_v913_weakT_ceramic_100um_selection.py"


def load_module():
    spec = importlib.util.spec_from_file_location("finalize_v913_weakT_ceramic", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


FINALIZER = load_module()


def active_values(seed: int) -> dict[str, float]:
    values = {}
    for index, field in enumerate(ACTIVE_CANDIDATE_PARAMETER_FIELDS, start=1):
        values[field] = float(seed + index) / 10.0
    values["Tref_K"] = 481.33
    values["peierls_nu0_s"] = 1.0e12
    values["taylor_nu0_s"] = 1.0e11
    values["rho_source0_m2"] = 1.0e14 + seed
    values["taylor_corr_rho_c_m2"] = 1.0e13 + seed
    return values


def row(candidate_id: str, target_class: str, score: float, strict: bool, seed: int):
    record = {
        "candidate_id": candidate_id,
        "target_class": target_class,
        "oneD_complete": True,
        "oneD_weakT_gate": strict if target_class == "weakT_FCC_like" else False,
        "oneD_ceramic_gate": strict if target_class == "ceramic_like" else False,
        "oneD_weakT_score": score if target_class == "weakT_FCC_like" else score + 100.0,
        "oneD_ceramic_score": score if target_class == "ceramic_like" else score + 100.0,
        "oneD_initial_span_MPa_sqrt_m": 0.2 + score,
        "oneD_developed_span_MPa_sqrt_m": 0.3 + score,
        "oneD_median_R_rise_MPa_sqrt_m": 2.0,
    }
    record.update(active_values(seed))
    return record


def write_rankings(root: Path, ceramic_strict_count: int = 3) -> None:
    weak = pd.DataFrame(
        [
            row("v913_zeroD_sobol_0000101", "weakT_FCC_like", 1.2, True, 1),
            row("v913_zeroD_sobol_0000102", "weakT_FCC_like", 0.4, True, 2),
            row("v913_zeroD_sobol_0000103", "weakT_FCC_like", 0.8, True, 3),
            row("v913_zeroD_sobol_0000104", "weakT_FCC_like", 0.1, False, 4),
        ]
    )
    ceramic_flags = [index < ceramic_strict_count for index in range(3)]
    ceramic = pd.DataFrame(
        [
            row("v913_zeroD_sobol_0000201", "ceramic_like", 0.9, ceramic_flags[0], 11),
            row("v913_zeroD_sobol_0000202", "ceramic_like", 0.3, ceramic_flags[1], 12),
            row("v913_zeroD_sobol_0000203", "ceramic_like", 0.6, ceramic_flags[2], 13),
        ]
    )
    weak.to_csv(root / "weakT_ranked.csv", index=False)
    ceramic.to_csv(root / "ceramic_ranked.csv", index=False)


def test_selects_two_strict_rows_and_promotes_lowest_score(tmp_path: Path):
    write_rankings(tmp_path)
    top2, primary, manifest = FINALIZER.build_outputs(tmp_path, 2)

    assert top2["candidate_id"].tolist() == [
        "v913_zeroD_sobol_0000102",
        "v913_zeroD_sobol_0000103",
        "v913_zeroD_sobol_0000202",
        "v913_zeroD_sobol_0000203",
    ]
    assert primary["candidate_id"].tolist() == [
        "v913_zeroD_sobol_0000102",
        "v913_zeroD_sobol_0000202",
    ]
    assert primary["option_key"].tolist() == [
        "v913_paper_weakT01_0000102_persistent_sites",
        "v913_paper_ceramic01_0000202_persistent_sites",
    ]
    assert top2["oneD_strict_gate_passed"].astype(bool).all()
    assert len(manifest["primary_candidates"]) == 2
    assert len(manifest["backup_candidates"]) == 2


def test_fails_closed_when_class_has_fewer_than_two_strict_rows(tmp_path: Path):
    write_rankings(tmp_path, ceramic_strict_count=1)
    with pytest.raises(RuntimeError, match="only 1 strict candidates"):
        FINALIZER.build_outputs(tmp_path, 2)
