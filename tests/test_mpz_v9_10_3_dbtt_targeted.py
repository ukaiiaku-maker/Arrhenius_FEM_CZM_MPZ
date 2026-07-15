from __future__ import annotations

import pandas as pd

import optimize_mpz_v9_10_2_independent_shape_global as v102
import optimize_mpz_v9_10_3_dbtt_targeted_global as v103


def target_like_detail() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "T_K": 300.0,
                "completed": True,
                "K_init_proxy": 15.0,
                "K_plateau_proxy": 15.5,
                "delta_KR_proxy": 0.5,
                "early_rise_per_100um_proxy": 0.10,
                "plateau_rise_per_100um_proxy": 0.0,
                "K_init_target": 15.0,
                "K_init_scale": 2.0,
                "K_plateau_target": 15.5,
                "K_plateau_scale": 2.0,
                "early_rise_per_100um_target": 0.10,
                "early_rise_scale": 0.40,
                "plateau_rise_per_100um_target": 0.0,
                "plateau_rise_scale": 0.25,
                "delta_KR_min": 0.0,
                "delta_KR_max": 1.5,
                "weight": 1.0,
            },
            {
                "T_K": 700.0,
                "completed": True,
                "K_init_proxy": 19.0,
                "K_plateau_proxy": 22.0,
                "delta_KR_proxy": 3.0,
                "early_rise_per_100um_proxy": 0.50,
                "plateau_rise_per_100um_proxy": 0.0,
                "K_init_target": 19.0,
                "K_init_scale": 2.5,
                "K_plateau_target": 22.0,
                "K_plateau_scale": 3.0,
                "early_rise_per_100um_target": 0.50,
                "early_rise_scale": 0.60,
                "plateau_rise_per_100um_target": 0.0,
                "plateau_rise_scale": 0.35,
                "delta_KR_min": 1.0,
                "delta_KR_max": 5.0,
                "weight": 1.0,
            },
            {
                "T_K": 900.0,
                "completed": True,
                "K_init_proxy": 33.0,
                "K_plateau_proxy": 43.0,
                "delta_KR_proxy": 10.0,
                "early_rise_per_100um_proxy": 1.50,
                "plateau_rise_per_100um_proxy": 0.0,
                "K_init_target": 33.0,
                "K_init_scale": 3.0,
                "K_plateau_target": 43.0,
                "K_plateau_scale": 4.0,
                "early_rise_per_100um_target": 1.50,
                "early_rise_scale": 0.75,
                "plateau_rise_per_100um_target": 0.0,
                "plateau_rise_scale": 0.40,
                "delta_KR_min": 7.0,
                "delta_KR_max": 13.0,
                "weight": 1.0,
            },
            {
                "T_K": 1200.0,
                "completed": True,
                "K_init_proxy": 40.0,
                "K_plateau_proxy": 50.0,
                "delta_KR_proxy": 10.0,
                "early_rise_per_100um_proxy": 1.50,
                "plateau_rise_per_100um_proxy": 0.0,
                "K_init_target": 40.0,
                "K_init_scale": 3.0,
                "K_plateau_target": 50.0,
                "K_plateau_scale": 4.0,
                "early_rise_per_100um_target": 1.50,
                "early_rise_scale": 0.75,
                "plateau_rise_per_100um_target": 0.0,
                "plateau_rise_scale": 0.40,
                "delta_KR_min": 7.0,
                "delta_KR_max": 13.0,
                "weight": 1.0,
            },
        ]
    )


def test_target_aware_loss_removes_zero_shelf_shortcut() -> None:
    good = target_like_detail()
    bad = good.copy()
    bad.loc[bad.T_K == 300.0, ["K_init_proxy", "K_plateau_proxy"]] = 0.0
    good_loss = sum(v103.dbtt_target_components(good).values())
    bad_loss = sum(v103.dbtt_target_components(bad).values())
    assert good_loss == 0.0
    assert bad_loss > good_loss
    assert v103.dbtt_target_components(bad)["DBTT_low_shelf_guard_loss"] > 0.0


def test_dbtt_acceptance_requires_finite_low_shelf() -> None:
    good = target_like_detail()
    accepted, reason = v103.dbtt_acceptance("DBTT", good, {})
    assert accepted
    assert reason == "DBTT_target_aware_zeroD_gate_passed"

    bad = good.copy()
    bad.loc[bad.T_K == 300.0, ["K_init_proxy", "K_plateau_proxy"]] = 0.0
    accepted, reason = v103.dbtt_acceptance("DBTT", bad, {})
    assert not accepted
    assert reason == "DBTT_low_shelf_outside_window"


def test_dbtt_targeted_search_preserves_full_independent_shape_space() -> None:
    assert len(v102.PARAMETER_NAMES) == len(v102.BOUNDS) == 29
    assert "cleave_exp_a" in v102.PARAMETER_NAMES
    assert "emit_exp_a" in v102.PARAMETER_NAMES
    assert "peierls_exp_a" in v102.PARAMETER_NAMES
    assert "taylor_exp_a" in v102.PARAMETER_NAMES
