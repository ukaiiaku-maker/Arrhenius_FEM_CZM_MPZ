"""Transactional PF-reference driving-force (KJ/J target) loading controller.

Implements FEM_CZM_HANDOFF.md section 6 ("Immediate development priority:
driving-force-controlled loading") as pure, FEM-independent control logic: an
elastic square-root predictor followed by a safeguarded secant correction
that drives a scalar boundary-displacement trial toward a target KJ or J
within a tolerance, bounded by a per-iteration trust-region growth factor and
a hard iteration cap.

This module deliberately knows nothing about meshes, plasticity state, MPZ
state, hazard RNG, or cohesive geometry. It calls a caller-supplied
``solve_trial(Uapp) -> achieved`` callback and only ever reasons about the
scalar (Uapp, achieved) pairs it returns. Transactionality is therefore the
caller's responsibility -- and in this codebase that responsibility is
already discharged for free by the existing inner trial loop in
``sharp_front.run_2d`` (see the docstring of ``solve_to_target`` below for the
exact integration point), which snapshots and reverts ``u``, ``ep_gp``,
``rho_gp``, and ``Uapp`` around every trial before any front-engine hazard
state is touched. Hazard/RNG/MPZ state is committed only after a step is
accepted (outside the trial loop), so no additional snapshot machinery is
required to make target-seeking retries transactional in that loop.

No hazard, plasticity, MPZ, cohesive, or RNG state is read or written here.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

MODEL_ID = "pf_theta0_j_controlled_loading_v10_0_5_18_4_0"


@dataclass(frozen=True)
class JControllerConfig:
    rel_tol: float = 1.0e-3
    abs_floor: float = 1.0e-6
    max_iterations: int = 25
    achieved_min: float = 1.0e-12
    max_growth_ratio: float = 4.0
    max_shrink_ratio: float = 0.1
    channel_power: float = 1.0

    def __post_init__(self) -> None:
        if self.rel_tol <= 0.0:
            raise ValueError("rel_tol must be positive")
        if self.max_iterations < 1:
            raise ValueError("max_iterations must be >= 1")
        if self.channel_power <= 0.0:
            raise ValueError("channel_power must be positive")
        if not (0.0 < self.max_shrink_ratio <= 1.0 <= self.max_growth_ratio):
            raise ValueError("require max_shrink_ratio <= 1.0 <= max_growth_ratio")


@dataclass(frozen=True)
class TrialRecord:
    iteration: int
    Uapp: float
    achieved: float
    accepted: bool


@dataclass(frozen=True)
class JControlledStepResult:
    converged: bool
    accepted_Uapp: float
    achieved: float
    target: float
    iterations: int
    trials: tuple[TrialRecord, ...] = field(default_factory=tuple)

    @property
    def rejected_trials(self) -> tuple[TrialRecord, ...]:
        return tuple(t for t in self.trials if not t.accepted)


def channel_power_for(channel: str) -> float:
    """Exponent p such that achieved ~ Uapp**p in the elastic prefracture
    regime. FEM_CZM_HANDOFF.md section 6 states J ~ U**2; since KJ ~ sqrt(J),
    KJ ~ U**1.
    """
    if channel == "J":
        return 2.0
    if channel == "KJ":
        return 1.0
    raise ValueError(f"unknown channel {channel!r}; expected 'J' or 'KJ'")


def elastic_predictor(
    Uapp_old: float,
    achieved_old: float,
    target: float,
    *,
    power: float,
    achieved_min: float,
) -> float:
    """U_new = U_old * (target / max(achieved_old, achieved_min)) ** (1/power).

    Reduces to the handoff's stated square-root predictor
    (U_new = U_old * sqrt(J_target / max(J_old, J_min))) for power=2 (the J
    channel), and to the equivalent linear form for power=1 (the KJ channel).
    """
    denom = max(achieved_old, achieved_min)
    return Uapp_old * (target / denom) ** (1.0 / power)


def _safeguard(
    Uapp_candidate: float,
    Uapp_prev: float,
    *,
    max_growth_ratio: float,
    max_shrink_ratio: float,
) -> float:
    if Uapp_prev == 0.0:
        return Uapp_candidate
    lo = Uapp_prev * max_shrink_ratio
    hi = Uapp_prev * max_growth_ratio
    if Uapp_prev < 0.0:
        lo, hi = hi, lo
    return min(max(Uapp_candidate, lo), hi)


def _secant_update(
    U_a: float,
    f_a: float,
    U_b: float,
    f_b: float,
) -> float | None:
    """Root of the residual f(U) = achieved(U) - target via secant on the
    last two distinct trial points. Returns None if the points coincide or
    the residuals are degenerate (caller should fall back to the predictor).
    """
    if U_a == U_b:
        return None
    denom = f_b - f_a
    if denom == 0.0:
        return None
    return U_b - f_b * (U_b - U_a) / denom


def solve_to_target(
    Uapp_initial: float,
    achieved_initial: float,
    target: float,
    *,
    solve_trial: Callable[[float], float],
    config: JControllerConfig | None = None,
) -> JControlledStepResult:
    """Drive Uapp toward `target` by repeated calls to `solve_trial`.

    `solve_trial(Uapp)` must return the scalar J or KJ achieved at that trial
    displacement without committing any permanent state -- it is expected to
    be the caller's own re-solve of the same physical step at a new trial
    value, exactly as the inner `while True:` loop in
    `arrhenius_fracture/sharp_front.py::run_2d` (around lines 2456-2586)
    already does for its adaptive hazard-clock retry: `u`, `ep_gp`, `rho_gp`,
    and `Uapp` are saved before the trial and restored if the trial is
    rejected, and no front-engine hazard/RNG/MPZ state is touched until
    after the trial loop accepts a value. To wire this controller in, the
    intended integration replaces (or augments) that loop's
    `dU_step = cfg.loading.dU_top * trial_frac` / adaptive-clock accept
    condition with: predict `Uapp` from this module against the PF
    `KJ_target(t)` read by `pf_theta0_driving_force_trajectory_v10051840`,
    call `solve_trial` (the existing assemble/solve/update_plasticity/
    compute_J_integral block) inside the same revert-on-reject scaffolding,
    and only fall through to hazard/front-engine commit once this function
    reports `converged=True`.

    Returns a JControlledStepResult; on non-convergence within
    `max_iterations`, `converged` is False and `accepted_Uapp`/`achieved`
    report the last (best-effort) trial -- the caller decides whether to
    accept that as a right-censored step, per existing
    `control_state: right_censored_endpoint` conventions in this codebase,
    or fail the step.
    """
    cfg = config or JControllerConfig()
    trials: list[TrialRecord] = []

    def _within_tol(achieved: float) -> bool:
        return abs(achieved - target) <= max(cfg.rel_tol * abs(target), cfg.abs_floor)

    if _within_tol(achieved_initial):
        trials.append(TrialRecord(0, Uapp_initial, achieved_initial, True))
        return JControlledStepResult(
            converged=True,
            accepted_Uapp=Uapp_initial,
            achieved=achieved_initial,
            target=target,
            iterations=0,
            trials=tuple(trials),
        )

    U_prev, f_prev = Uapp_initial, achieved_initial
    U_curr = elastic_predictor(
        Uapp_initial,
        achieved_initial,
        target,
        power=cfg.channel_power,
        achieved_min=cfg.achieved_min,
    )
    U_curr = _safeguard(
        U_curr, Uapp_initial,
        max_growth_ratio=cfg.max_growth_ratio,
        max_shrink_ratio=cfg.max_shrink_ratio,
    )

    for iteration in range(1, cfg.max_iterations + 1):
        achieved = solve_trial(U_curr)
        converged = _within_tol(achieved)
        trials.append(TrialRecord(iteration, U_curr, achieved, converged))
        if converged:
            return JControlledStepResult(
                converged=True,
                accepted_Uapp=U_curr,
                achieved=achieved,
                target=target,
                iterations=iteration,
                trials=tuple(trials),
            )
        if iteration == cfg.max_iterations:
            break

        secant = _secant_update(U_prev, f_prev - target, U_curr, achieved - target)
        if secant is None:
            next_U = elastic_predictor(
                U_curr, achieved, target,
                power=cfg.channel_power, achieved_min=cfg.achieved_min,
            )
        else:
            next_U = secant
        next_U = _safeguard(
            next_U, U_curr,
            max_growth_ratio=cfg.max_growth_ratio,
            max_shrink_ratio=cfg.max_shrink_ratio,
        )
        U_prev, f_prev = U_curr, achieved
        U_curr = next_U

    last = trials[-1]
    return JControlledStepResult(
        converged=False,
        accepted_Uapp=last.Uapp,
        achieved=last.achieved,
        target=target,
        iterations=cfg.max_iterations,
        trials=tuple(trials),
    )


__all__ = [
    "MODEL_ID",
    "JControllerConfig",
    "TrialRecord",
    "JControlledStepResult",
    "channel_power_for",
    "elastic_predictor",
    "solve_to_target",
]
