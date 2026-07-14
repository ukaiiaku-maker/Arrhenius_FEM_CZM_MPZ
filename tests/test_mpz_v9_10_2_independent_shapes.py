import numpy as np

import optimize_mpz_v9_10_2_independent_shape_global as opt
from arrhenius_fracture.emission_derived_plasticity import (
    EmissionDerivedPeierlsTaylorConfig,
    ExpFloorSurface,
)
from arrhenius_fracture.emission_derived_plasticity_v9102 import (
    EmissionDerivedPeierlsTaylorModel,
    IndependentShapeEntropyMechanismScale,
)
from arrhenius_fracture.moving_process_zone_v9102 import MovingProcessZoneState


def midpoint_vector():
    return np.asarray(
        [0.5 * (opt.BOUNDS[name][0] + opt.BOUNDS[name][1]) for name in opt.PARAMETER_NAMES],
        dtype=float,
    )


def test_independent_search_has_29_coordinates_and_strict_H_order():
    x = midpoint_vector()
    p = opt.decode(x)
    assert len(opt.PARAMETER_NAMES) == len(opt.bounds_list()) == 29
    assert p["taylor_H0_eV"] > p["peierls_H0_eV"]
    assert p["independent_shape_all_four_active"] == 1.0


def test_decode_preserves_four_distinct_shape_pairs():
    x = midpoint_vector()
    index = {name: i for i, name in enumerate(opt.PARAMETER_NAMES)}
    x[index["cleave_exp_a"]] = 0.11
    x[index["emit_exp_a"]] = 0.22
    x[index["peierls_exp_a"]] = 0.33
    x[index["taylor_exp_a"]] = 0.44
    x[index["cleave_exp_n"]] = 0.71
    x[index["emit_exp_n"]] = 0.82
    x[index["peierls_exp_n"]] = 0.93
    x[index["taylor_exp_n"]] = 1.04
    p = opt.decode(x)
    assert p["cleave_exp_a"] == 0.11
    assert p["emit_exp_a"] == 0.22
    assert p["peierls_exp_a"] == 0.33
    assert p["taylor_exp_a"] == 0.44
    assert p["cleave_exp_n"] == 0.71
    assert p["emit_exp_n"] == 0.82
    assert p["peierls_exp_n"] == 0.93
    assert p["taylor_exp_n"] == 1.04


def test_mechanism_shape_changes_stress_dependence_without_changing_H0():
    parent = ExpFloorSurface(
        G00_eV=2.0,
        gT_eV_per_K=0.0,
        sigc0_Pa=2.0e9,
        sT_Pa_per_K=0.0,
        a=0.2,
        n=1.0,
        floor_fraction=0.02,
    )
    model = EmissionDerivedPeierlsTaylorModel(
        EmissionDerivedPeierlsTaylorConfig(
            parent=parent,
            peierls=IndependentShapeEntropyMechanismScale(
                0.5, 0.0, 0.10, 0.75, rate_prefactor_s=1.0e12
            ),
            taylor=IndependentShapeEntropyMechanismScale(
                0.5, 0.0, 1.50, 2.00, rate_prefactor_s=1.0e11
            ),
        )
    )
    Gp0 = float(model.barrier_eV("peierls", 0.0, 700.0))
    Gt0 = float(model.barrier_eV("taylor", 0.0, 700.0))
    Gp = float(model.barrier_eV("peierls", 1.0e9, 700.0))
    Gt = float(model.barrier_eV("taylor", 1.0e9, 700.0))
    assert np.isclose(Gp0, Gt0)
    assert not np.isclose(Gp, Gt)


def test_v9102_spatial_state_has_no_encounters_when_peierls_is_frozen():
    rho = np.logspace(12.0, 16.0, 8)
    encounter = MovingProcessZoneState.encounter_rate_s(0.0, 1.0e-8, rho, 100.0)
    assert np.all(encounter == 0.0)
