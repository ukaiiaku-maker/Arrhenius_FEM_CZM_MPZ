from __future__ import annotations

from pathlib import Path
import sys

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import run_oneD_v2_natural_forward as natural
from reduced_fracture_v2.predictive import (
    FEMCZM_PREDICTIVE_LIFECYCLE_V2,
    PF_PREDICTIVE_LIFECYCLE_V2,
    ProviderMechanicsMap,
    SourceDriveMapBoundsError,
    provider_loading_map,
)


class _PFEngine:
    def __init__(self):
        self.call = None

    def predict_clock_increment(self, K, T, dt):
        self.call = (K, T, dt)
        return 0.125


class _FEMEngine:
    def __init__(self):
        self.call = None

    def sigma_tip(self, K):
        self.call = ("sigma_tip", K)
        return 2.0 * K

    def lambda_cleave(self, sigma, T):
        self.call = ("lambda_cleave", sigma, T)
        return 0.5, None, None


def test_natural_forward_uses_backend_owned_predictors():
    pf = _PFEngine()
    fem = _FEMEngine()
    assert natural._predict(pf, "pf", 1.0, 2.0, 1000.0, 3.0) == 0.125
    assert pf.call == (1.0, 1000.0, 3.0)
    assert natural._predict(fem, "fem", 1.0, 2.0, 1000.0, 3.0) == 1.5
    assert fem.call == ("lambda_cleave", 2.0, 1000.0)


def test_fem_natural_forward_excludes_rejected_aggregate_closure():
    source = (ROOT / "scripts/run_oneD_v2_natural_forward.py").read_text()
    assert "aggregate_emission=True" not in source
    assert "constructor=_fem_constructor" in source


def test_natural_forward_has_no_threshold_priming():
    source = (ROOT / "scripts/run_oneD_v2_natural_forward.py").read_text()
    assert "_threshold_prime" not in source
    assert '"threshold_priming":False' in source


def test_provider_lifecycle_reductions_are_distinct_and_candidate_independent():
    assert PF_PREDICTIVE_LIFECYCLE_V2.backend == "PF"
    assert FEMCZM_PREDICTIVE_LIFECYCLE_V2.backend == "FEMCZM"
    assert PF_PREDICTIVE_LIFECYCLE_V2.threshold_mode == "DETERMINISTIC_UNIT_ACTION"
    assert FEMCZM_PREDICTIVE_LIFECYCLE_V2.threshold_mode == "EXPONENTIAL_SOURCE_THRESHOLD"
    mechanics = ProviderMechanicsMap(
        "PF", np.array([0.0, 100.0]), np.ones(2), np.ones(2)
    )
    first = provider_loading_map(mechanics, seed=1, target_extension_m=20.0e-6)
    second = provider_loading_map(mechanics, seed=999, target_extension_m=20.0e-6)
    assert first.threshold_actions == second.threshold_actions == (1.0,) * 4
    assert first.path_advances_m == second.path_advances_m == (5.0e-6,) * 4


def test_terminal_native_baselines_match_topology_without_silent_bounds():
    path = ROOT / "analysis_outputs/oneD_v2_terminal_predictive_program/oneD_v2_native_baseline_results.csv"
    frame = pd.read_csv(path).set_index(["provider", "material_class"])
    assert frame.status.eq("TARGET_RIGHT_CENSORED").all()
    expected = {("PF", "Peak"): 1, ("PF", "DBTT"): 2,
                ("FEMCZM", "Peak"): 1, ("FEMCZM", "DBTT"): 3}
    assert frame.physical_avalanche_count.astype(int).to_dict() == expected
    assert frame.topology_match.astype(bool).all()


def test_source_map_out_of_domain_fails_closed():
    from scripts.run_oneD_v2_predictive_campaign import inputs

    _, _, providers = inputs()
    drive = providers["PF"][1]
    try:
        drive.evaluate(0.0, 1000.1e-6)
    except SourceDriveMapBoundsError:
        pass
    else:
        raise AssertionError("source-drive map silently extrapolated")
