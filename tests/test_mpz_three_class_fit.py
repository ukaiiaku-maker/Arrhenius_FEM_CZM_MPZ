from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pandas as pd

from arrhenius_fracture import sharp_front as sf
from arrhenius_fracture.config import make_emergent_config
from fit_mpz_three_classes import row_to_args, simulate, summarize_rcurve


def _row(name="weakT"):
    return pd.read_csv("mpz_three_class_initial_guesses.csv").set_index("target_class").loc[name]


def test_robust_rcurve_summary_recovers_plateau_and_slopes():
    rng = np.random.default_rng(7)
    a = np.arange(0.0, 1000.1, 10.0)
    K = 15.0 + 5.0 * (1.0 - np.exp(-a / 180.0)) + rng.normal(0.0, 0.12, a.size)
    q = summarize_rcurve(a, K, (20.0, 220.0), (700.0, 1000.0))
    assert 14.5 < q["K_init"] < 15.5
    assert 19.4 < q["K_plateau"] < 20.4
    assert 4.0 < q["delta_KR"] < 6.0
    assert q["early_rise_per_100um"] > 0.5
    assert abs(q["plateau_rise_per_100um"]) < 0.4
    assert 19.0 < q["sat_Kss"] < 21.0


def test_pre_renewal_spatial_state_is_reported():
    row = _row("DBTT")
    args = row_to_args(row, dK=0.1, Kdot=0.005, n_advances=1,
                       Kmax=65.0, da_um=5.0)
    eng = sf.build_engine(args, make_emergent_config().material)
    eng.mpz_state.retained[0, 0] = 2.0
    eng.B = 1.0
    info = eng._renew(1.0)
    assert info["n_fire"] == 1
    assert info["mpz_K_shield_pre_renewal_Pa_sqrt_m"] > 0.0
    assert info["mpz_retained_pre_renewal"] == 2.0
    assert 0.0 <= info["mpz_available_site_fraction_pre_renewal"] <= 1.0


def test_adaptive_first_passage_simulation_completes_without_event_staircase():
    row = _row("ceramic")
    opt = SimpleNamespace(
        dK=0.5, Kdot=0.005, n_advances=1, Kmax=40.0, da_um=5.0,
        early_window_um=(20.0, 220.0), plateau_window_um=(700.0, 1000.0),
        target_dB_substep=0.8, target_emission_hazard_substep=5.0,
        source_active_fraction_min=1.0e-3, min_substep_fraction=1.0e-7,
        max_substeps=100000,
    )
    q = simulate(row, 300.0, opt)
    assert q["n_events"] == 1
    assert np.isfinite(q["K_init"])
    assert not q["integration_stalled"]
    assert q["n_substeps"] > 1
