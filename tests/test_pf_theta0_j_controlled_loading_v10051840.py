from __future__ import annotations

import pytest

from arrhenius_fracture import pf_theta0_j_controlled_loading_v10051840 as controller


def test_channel_power_matches_handoff_elastic_scaling():
    assert controller.channel_power_for("J") == 2.0
    assert controller.channel_power_for("KJ") == 1.0
    with pytest.raises(ValueError):
        controller.channel_power_for("bogus")


def test_elastic_predictor_reduces_to_handoff_sqrt_formula_for_J_channel():
    U_new = controller.elastic_predictor(
        1.0, 4.0, 9.0, power=2.0, achieved_min=1.0e-12
    )
    assert U_new == pytest.approx(1.0 * (9.0 / 4.0) ** 0.5)


def test_elastic_predictor_is_linear_for_KJ_channel():
    U_new = controller.elastic_predictor(
        2.0, 4.0, 8.0, power=1.0, achieved_min=1.0e-12
    )
    assert U_new == pytest.approx(2.0 * (8.0 / 4.0))


def test_perfectly_linear_KJ_response_converges_in_one_trial():
    k = 3.5e6

    def solve_trial(U: float) -> float:
        return k * U

    result = controller.solve_to_target(
        Uapp_initial=1.0e-7,
        achieved_initial=k * 1.0e-7,
        target=k * 5.0e-7,
        solve_trial=solve_trial,
        # 5x growth exceeds the default max_growth_ratio (4.0), which would
        # correctly safeguard-clip the first step; widen it here so this
        # test isolates the exact-linear-model convergence claim from the
        # separate safeguard-clamping behavior (covered by its own test).
        config=controller.JControllerConfig(channel_power=1.0, max_growth_ratio=10.0),
    )
    assert result.converged
    assert result.iterations == 1
    assert result.accepted_Uapp == pytest.approx(5.0e-7, rel=1.0e-9)
    assert result.achieved == pytest.approx(k * 5.0e-7, rel=1.0e-9)


def test_already_within_tolerance_takes_zero_trial_calls():
    calls = []

    def solve_trial(U: float) -> float:
        calls.append(U)
        return 10.0

    result = controller.solve_to_target(
        Uapp_initial=1.0,
        achieved_initial=10.0,
        target=10.0000001,
        solve_trial=solve_trial,
        config=controller.JControllerConfig(rel_tol=1.0e-3),
    )
    assert result.converged
    assert result.iterations == 0
    assert calls == []


def test_mildly_nonlinear_response_converges_via_secant_refinement():
    # KJ(U) = k*U - c*U**2 : softening response, not exactly linear or
    # quadratic, so the elastic predictor alone will under/overshoot and the
    # safeguarded secant correction must finish the job.
    k, c = 5.0e6, 1.0e12

    def solve_trial(U: float) -> float:
        return k * U - c * U * U

    target = solve_trial(4.0e-7) # pick an achievable target inside the safe branch
    result = controller.solve_to_target(
        Uapp_initial=1.0e-7,
        achieved_initial=solve_trial(1.0e-7),
        target=target,
        solve_trial=solve_trial,
        config=controller.JControllerConfig(rel_tol=1.0e-6, max_iterations=25),
    )
    assert result.converged
    assert result.achieved == pytest.approx(target, rel=1.0e-6)
    assert result.iterations >= 2 # nonlinear target should not resolve in one shot


def test_safeguard_bounds_growth_and_shrink_per_iteration():
    # A pathological predictor jump (huge target ratio) must be clamped to
    # the configured trust region rather than taking one huge unsafeguarded
    # step -- this is the "bounded secant/proportional correction" the
    # handoff requires.
    trial_Us = []

    def solve_trial(U: float) -> float:
        trial_Us.append(U)
        return 1.0e6 * U

    controller.solve_to_target(
        Uapp_initial=1.0e-7,
        achieved_initial=0.1,
        target=1.0e9, # would demand an enormous unsafeguarded jump
        solve_trial=solve_trial,
        config=controller.JControllerConfig(
            channel_power=1.0, max_growth_ratio=2.0, max_shrink_ratio=0.5, max_iterations=10
        ),
    )
    assert trial_Us[0] <= 1.0e-7 * 2.0 + 1.0e-30


def test_non_convergent_pathological_response_reports_converged_false_not_exception():
    def solve_trial(U: float) -> float:
        return 42.0 # constant: target is unreachable by construction

    result = controller.solve_to_target(
        Uapp_initial=1.0,
        achieved_initial=42.0 - 100.0, # force the initial-tolerance check to fail
        target=100.0,
        solve_trial=solve_trial,
        config=controller.JControllerConfig(max_iterations=5, rel_tol=1.0e-9),
    )
    assert result.converged is False
    assert result.iterations == 5
    assert len(result.trials) == 5


def test_achieved_min_floor_prevents_division_by_zero_from_rest_state():
    U_new = controller.elastic_predictor(
        0.0, 0.0, 5.0, power=1.0, achieved_min=1.0e-6
    )
    assert U_new == 0.0 # Uapp_old=0 dominates regardless of the floor

    U_new_from_nonzero_start = controller.elastic_predictor(
        1.0e-9, 0.0, 5.0, power=1.0, achieved_min=1.0e-6
    )
    assert U_new_from_nonzero_start == pytest.approx(1.0e-9 * (5.0 / 1.0e-6))


def test_rejected_trials_are_all_trials_before_the_accepted_one():
    k = 2.0e6

    def solve_trial(U: float) -> float:
        return k * U

    result = controller.solve_to_target(
        Uapp_initial=1.0e-7,
        achieved_initial=k * 1.0e-7,
        target=k * 3.0e-7,
        solve_trial=solve_trial,
    )
    assert result.converged
    accepted = [t for t in result.trials if t.accepted]
    assert len(accepted) == 1
    assert accepted[0] is result.trials[-1]
    assert result.rejected_trials == result.trials[:-1]


def test_config_rejects_invalid_bounds():
    with pytest.raises(ValueError):
        controller.JControllerConfig(rel_tol=0.0)
    with pytest.raises(ValueError):
        controller.JControllerConfig(max_iterations=0)
    with pytest.raises(ValueError):
        controller.JControllerConfig(channel_power=-1.0)
    with pytest.raises(ValueError):
        controller.JControllerConfig(max_growth_ratio=0.5, max_shrink_ratio=2.0)
