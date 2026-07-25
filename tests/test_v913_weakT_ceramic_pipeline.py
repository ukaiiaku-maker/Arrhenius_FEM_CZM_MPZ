from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
ANALYZER_PATH = ROOT / "scripts" / "analyze_v913_weakT_ceramic_1d.py"
SEARCH_PATH = ROOT / "scripts" / "run_v913_zero_d_weakT_ceramic_search.py"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


ANALYZER = load_module(ANALYZER_PATH, "analyze_v913_weakT_ceramic_1d")
SEARCH = load_module(SEARCH_PATH, "run_v913_zero_d_weakT_ceramic_search")


def synthetic_events():
    return [
        {
            "cumulative_projected_extension_m": 50.0e-6,
            "K_MPa_sqrt_m": 10.0,
        },
        {
            "cumulative_projected_extension_m": 250.0e-6,
            "K_MPa_sqrt_m": 14.0,
        },
        {
            "cumulative_projected_extension_m": 500.0e-6,
            "K_MPa_sqrt_m": 16.0,
        },
        {
            "cumulative_projected_extension_m": 1000.0e-6,
            "K_MPa_sqrt_m": 18.0,
        },
    ]


def test_checkpoint_and_extension_average_use_strict_target_coverage():
    events = synthetic_events()
    assert ANALYZER.checkpoint_K(events, 0.0) == pytest.approx(10.0)
    assert ANALYZER.checkpoint_K(events, 250.0) == pytest.approx(14.0)
    assert ANALYZER.checkpoint_K(events, 1000.0) == pytest.approx(18.0)
    assert ANALYZER.checkpoint_K(events, 1100.0) != ANALYZER.checkpoint_K(events, 1100.0)
    average = ANALYZER.extension_average_K(events, 1000.0)
    assert 10.0 < average < 18.0


def test_new_policy_expands_beyond_peak_oriented_temperature_and_source_bounds():
    policy = json.loads(
        (ROOT / "mpz_v9_13_zero_d_weakT_ceramic_search_policy.json").read_text()
    )
    dimensions = policy["search_dimensions"]
    assert float(dimensions["cleave_gT_eV_per_K"]["low"]) < 0.0
    assert float(dimensions["emit_gT_eV_per_K"]["low"]) < 0.0
    assert float(dimensions["rho_source0_m2"]["low"]) <= 1.0e11
    assert float(policy["local_anchor_fraction"]) < 0.5


def test_physical_surface_gate_rejects_negative_high_temperature_surface():
    frame = pd.DataFrame(
        [
            {
                "Tref_K": 481.33,
                "cleave_G00_eV": 2.0,
                "cleave_gT_eV_per_K": 0.0,
                "cleave_sigc0_GPa": 4.0,
                "cleave_sT_GPa_per_K": 0.0,
                "emit_G00_eV": 2.0,
                "emit_gT_eV_per_K": 0.0,
                "emit_sigc0_GPa": 4.0,
                "emit_sT_GPa_per_K": 0.0,
            },
            {
                "Tref_K": 481.33,
                "cleave_G00_eV": 1.0,
                "cleave_gT_eV_per_K": -0.01,
                "cleave_sigc0_GPa": 4.0,
                "cleave_sT_GPa_per_K": 0.0,
                "emit_G00_eV": 2.0,
                "emit_gT_eV_per_K": 0.0,
                "emit_sigc0_GPa": 4.0,
                "emit_sT_GPa_per_K": 0.0,
            },
        ]
    )
    gate = SEARCH._physical_surface_gate(frame, [300.0, 1300.0])
    assert gate.tolist() == [True, False]
