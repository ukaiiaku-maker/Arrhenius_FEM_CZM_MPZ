from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pandas as pd

from arrhenius_fracture.emergent_gnd_contract_v913 import (
    ACTIVE_CANDIDATE_PARAMETER_FIELDS,
    candidate_parameter_fingerprint,
)


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "export_v913_weakT_ceramic_paper_handoff.py"
)
SPEC = importlib.util.spec_from_file_location("export_v913_weakT_ceramic_handoff", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


def _active(seed: float) -> dict[str, float]:
    values = {field: seed + index * 0.01 for index, field in enumerate(ACTIVE_CANDIDATE_PARAMETER_FIELDS)}
    values["Tref_K"] = 481.33
    values["peierls_nu0_s"] = 1.0e12
    values["taylor_nu0_s"] = 1.0e11
    return values


def _metrics(candidate_id: str, strict_field: str, score_field: str) -> dict[str, object]:
    return {
        "candidate_id": candidate_id,
        "complete_temperature_grid": True,
        "K50_temperature_span_MPa_sqrt_m": 1.0,
        "high_temperature_toughness_loss_MPa_sqrt_m": 1.0,
        "median_R_rise_first_to_50_MPa_sqrt_m": 1.0,
        "median_R_rise_25_to_50_MPa_sqrt_m": 0.0,
        "K50_peak_prominence_MPa_sqrt_m": 0.0,
        "weakT_gate": False,
        "weakT_score": 1.5,
        "ceramic_gate": False,
        "ceramic_score": 2.0,
        strict_field: True,
        score_field: 0.5,
    }


def test_build_handoff_preserves_exact_active_rows() -> None:
    weak_id = MODULE.EXPECTED["weakT_FCC_like"]["candidate_id"]
    ceramic_id = MODULE.EXPECTED["ceramic_like"]["candidate_id"]
    registry = pd.DataFrame(
        [
            {"candidate_id": weak_id, **_active(1.0)},
            {"candidate_id": ceramic_id, **_active(2.0)},
        ]
    )
    summary = {
        "schema": "synthetic",
        "selected": {
            "weakT_FCC_like": _metrics(weak_id, "weakT_gate", "weakT_score"),
            "ceramic_like": _metrics(ceramic_id, "ceramic_gate", "ceramic_score"),
        },
    }

    handoff, metadata = MODULE.build_handoff(registry, summary)
    assert handoff["candidate_id"].tolist() == [weak_id, ceramic_id]
    assert handoff["paper_material_class"].tolist() == [
        "weakT_FCC_like",
        "ceramic_like",
    ]
    assert len(metadata) == 2

    source_fingerprint = candidate_parameter_fingerprint(registry.to_dict(orient="records"))
    handoff_fingerprint = candidate_parameter_fingerprint(handoff.to_dict(orient="records"))
    assert handoff_fingerprint == source_fingerprint
    for field in ACTIVE_CANDIDATE_PARAMETER_FIELDS:
        assert handoff[field].tolist() == registry[field].tolist()


def test_expected_option_keys_are_stable() -> None:
    assert MODULE.EXPECTED["weakT_FCC_like"]["option_key"] == (
        "v913_paper_weakT01_0257068_persistent_sites"
    )
    assert MODULE.EXPECTED["ceramic_like"]["option_key"] == (
        "v913_paper_ceramic01_0189364_persistent_sites"
    )
