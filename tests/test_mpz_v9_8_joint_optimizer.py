import numpy as np
import pandas as pd

from optimize_mpz_v9_8_joint_response import (
    JointObjective,
    ObjectiveSettings,
    PARAMETER_NAMES,
    barrier_order_margin_eV,
    build_pt_model,
    load_targets,
    seed_vector,
    shape_from_seed,
    vector_to_parameters,
)


def synthetic_seed():
    return pd.Series(
        {
            "candidate_id": "synthetic",
            "cleave_G00_eV": 2.0,
            "cleave_gT_eV_per_K": 0.0,
            "cleave_sigc0_GPa": 4.0,
            "emit_G00_eV": 1.5,
            "emit_gT_eV_per_K": 0.0,
            "emit_sigc0_GPa": 2.5,
            "cleave_sT_GPa_per_K": 0.0,
            "cleave_exp_a": 0.2,
            "cleave_exp_n": 1.0,
            "cleave_floor_frac": 0.02,
            "emit_sT_GPa_per_K": 0.0,
            "emit_exp_a": 0.2,
            "emit_exp_n": 1.0,
            "emit_floor_frac": 0.02,
            "mpz_source_sites_per_system": 100.0,
            "mpz_source_refresh_length_m": 1.0e-6,
            "c_blunt": 0.5,
        }
    )


def test_ordered_absolute_barrier_parameterization():
    x = seed_vector(synthetic_seed())
    p = vector_to_parameters(x)
    assert p["taylor_H0_eV"] >= p["peierls_H0_eV"]
    assert p["source_sites_per_system"] > 0.0
    assert p["recovery_rate_s"] > 0.0


def test_barrier_order_margin_uses_same_resolved_stress():
    seed = synthetic_seed()
    x = seed_vector(seed)
    p = vector_to_parameters(x)
    p["peierls_H0_eV"] = 0.5
    p["delta_H_PT_eV"] = 1.0
    p["taylor_H0_eV"] = 1.5
    p["peierls_activation_entropy_kB"] = -10.0
    p["taylor_activation_entropy_kB"] = -10.0
    model = build_pt_model(p, shape_from_seed(seed), Tref_K=481.33)
    margin = barrier_order_margin_eV(
        model,
        np.array([300.0, 700.0, 1200.0]),
        np.linspace(0.0, 20.0e9, 21),
    )
    assert margin >= -1.0e-10


def test_target_interpolation(tmp_path):
    path = tmp_path / "targets.csv"
    pd.DataFrame(
        {
            "target_class": ["weakT", "weakT"],
            "T_K": [300.0, 1200.0],
            "K_init_target": [10.0, 20.0],
            "K_init_scale": [2.0, 2.0],
            "K_plateau_target": [15.0, 25.0],
            "K_plateau_scale": [2.0, 2.0],
            "early_rise_per_100um_target": [1.0, 1.0],
            "early_rise_scale": [0.5, 0.5],
            "plateau_rise_per_100um_target": [0.0, 0.0],
            "plateau_rise_scale": [0.3, 0.3],
            "delta_KR_min": [3.0, 3.0],
            "delta_KR_max": [7.0, 7.0],
            "weight": [1.0, 1.0],
        }
    ).to_csv(path, index=False)
    out = load_targets(path, "weakT", [300.0, 750.0, 1200.0])
    assert len(out) == 3
    assert np.isclose(out.loc[out.T_K == 750.0, "K_init_target"].iloc[0], 15.0)


def test_joint_objective_returns_finite_scalar(tmp_path):
    seed = synthetic_seed()
    targets = pd.DataFrame(
        {
            "target_class": ["weakT", "weakT"],
            "T_K": [300.0, 700.0],
            "K_init_target": [15.0, 15.0],
            "K_init_scale": [2.0, 2.0],
            "K_plateau_target": [20.0, 20.0],
            "K_plateau_scale": [3.0, 3.0],
            "early_rise_per_100um_target": [1.0, 1.0],
            "early_rise_scale": [0.5, 0.5],
            "plateau_rise_per_100um_target": [0.0, 0.0],
            "plateau_rise_scale": [0.3, 0.3],
            "delta_KR_min": [3.0, 3.0],
            "delta_KR_max": [7.0, 7.0],
            "weight": [1.0, 1.0],
        }
    )
    objective = JointObjective(
        ObjectiveSettings(
            target_class="weakT",
            temperatures=np.array([300.0, 700.0]),
            targets=targets,
            shape=shape_from_seed(seed),
            dK=2.0,
            Kmax=40.0,
        )
    )
    x = seed_vector(seed)
    value = objective(x)
    assert np.isfinite(value)
    assert value >= 0.0
    assert len(x) == len(PARAMETER_NAMES)
