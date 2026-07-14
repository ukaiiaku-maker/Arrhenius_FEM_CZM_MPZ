import numpy as np

import optimize_mpz_v9_10_1_shared_shape_global as shared
import optimize_mpz_v9_10_unified_global as base


def midpoint_vector():
    return np.asarray(
        [0.5 * (shared.BOUNDS[name][0] + shared.BOUNDS[name][1]) for name in shared.PARAMETER_NAMES],
        dtype=float,
    )


def test_shared_shape_parameterization_applies_to_all_four_barriers():
    p = shared.decode(midpoint_vector())
    assert p["shared_shape_all_four_active"] == 1.0
    assert p["cleave_exp_a"] == p["emit_exp_a"] == p["peierls_exp_a"] == p["taylor_exp_a"]
    assert p["cleave_exp_n"] == p["emit_exp_n"] == p["peierls_exp_n"] == p["taylor_exp_n"]
    assert p["taylor_H0_eV"] > p["peierls_H0_eV"]


def test_shared_search_is_full_space_and_lower_dimensional():
    assert len(shared.PARAMETER_NAMES) == 23
    assert "shared_exp_a" in shared.PARAMETER_NAMES
    assert "shared_exp_n" in shared.PARAMETER_NAMES
    assert "cleave_exp_a" not in shared.PARAMETER_NAMES
    assert "emit_exp_a" not in shared.PARAMETER_NAMES


def test_unified_model_uses_shared_shape_for_peierls_and_taylor():
    p = shared.decode(midpoint_vector())
    model = base.build_model(p, 481.33)
    assert model.cfg.parent.a == p["shared_exp_a"]
    assert model.cfg.parent.n == p["shared_exp_n"]

    # Peierls and Taylor inherit the common EXP-floor shape.  After normalizing
    # by their zero-stress and floor levels, their same-stress shape factors are
    # identical despite different barrier heights.
    T = 700.0
    sigma = 2.0e9
    gp0 = float(model.barrier_eV("peierls", 0.0, T))
    gt0 = float(model.barrier_eV("taylor", 0.0, T))
    gp = float(model.barrier_eV("peierls", sigma, T))
    gt = float(model.barrier_eV("taylor", sigma, T))
    fp = model.cfg.parent.floor_fraction * gp0
    ft = model.cfg.parent.floor_fraction * gt0
    shape_p = (gp - fp) / max(gp0 - fp, 1.0e-30)
    shape_t = (gt - ft) / max(gt0 - ft, 1.0e-30)
    assert np.isclose(shape_p, shape_t, rtol=1.0e-10, atol=1.0e-12)
