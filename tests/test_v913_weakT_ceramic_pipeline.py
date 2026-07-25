from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
ANALYZER_PATH = ROOT / "scripts" / "analyze_v913_weakT_ceramic_1d.py"
SEARCH_PATH = ROOT / "scripts" / "run_v913_zero_d_weakT_ceramic_search.py"
LAUNCHER_PATH = ROOT / "scripts" / "run_v913_weakT_ceramic_search.sh"


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
            "cumulative_projected_extension_m": 25.0e-6,
            "K_MPa_sqrt_m": 9.0,
        },
        {
            "cumulative_projected_extension_m": 50.0e-6,
            "K_MPa_sqrt_m": 10.0,
        },
        {
            "cumulative_projected_extension_m": 100.0e-6,
            "K_MPa_sqrt_m": 12.0,
        },
    ]


def test_checkpoint_and_extension_average_use_strict_target_coverage():
    events = synthetic_events()
    assert ANALYZER.checkpoint_K(events, 0.0) == pytest.approx(9.0)
    assert ANALYZER.checkpoint_K(events, 50.0) == pytest.approx(10.0)
    assert ANALYZER.checkpoint_K(events, 100.0) == pytest.approx(12.0)
    assert ANALYZER.checkpoint_K(events, 110.0) != ANALYZER.checkpoint_K(events, 110.0)
    average = ANALYZER.extension_average_K(events, 100.0)
    assert 9.0 < average < 12.0


def test_new_policy_expands_beyond_peak_oriented_temperature_and_source_bounds():
    policy = json.loads(
        (ROOT / "mpz_v9_13_zero_d_weakT_ceramic_search_policy.json").read_text()
    )
    dimensions = policy["search_dimensions"]
    assert float(dimensions["cleave_gT_eV_per_K"]["low"]) < 0.0
    assert float(dimensions["emit_gT_eV_per_K"]["low"]) < 0.0
    assert float(dimensions["rho_source0_m2"]["low"]) <= 1.0e11
    assert float(policy["local_anchor_fraction"]) < 0.5


def test_launcher_uses_reduced_five_temperature_100um_gate():
    text = LAUNCHER_PATH.read_text()
    assert "v9_13_long_map_exponential_110um_v2" in text
    assert 'TEMPERATURES_K="${TEMPERATURES_K:-300 600 900 1100 1200}"' in text
    assert 'TARGET_EXT_UM="${TARGET_EXT_UM:-100}"' in text
    assert 'SAMPLES="${SAMPLES:-131072}"' in text
    assert 'ZERO_D_EXACT_PER_CLASS="${ZERO_D_EXACT_PER_CLASS:-128}"' in text
    assert 'ZERO_D_PROMOTE_PER_CLASS="${ZERO_D_PROMOTE_PER_CLASS:-8}"' in text
    assert 'FINAL_PROMOTE_PER_CLASS="${FINAL_PROMOTE_PER_CLASS:-3}"' in text
    assert "one_d_${TARGET_TAG}um" in text
    assert "analysis_${TARGET_TAG}um" in text
    assert "1000um" not in text
    assert "LONG_LOADING_MAP" not in text
    assert text.count("scripts.run_v913_autonomous_dbtt_search") == 1
    assert "V913_WEAKT_CERAMIC_1D_START" in text


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
