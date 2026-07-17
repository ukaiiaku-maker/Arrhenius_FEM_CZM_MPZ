"""Exact state-change-aware timestep guard for the v9.10.4.4 reduced front.

The v9.10.4.3 count-aware guard estimated transport stiffness from
``active_count * (encounter + Taylor + escape + recovery rates)``.  That is
still overly restrictive when trapping/release is already close to its exact
equilibrium partition, or when most active content is retained and therefore
cannot escape.  In those states a very large microscopic rate can correspond
to almost no actual state change.

This guard derives a timestep from the exact exponential change used by each
operator: finite-source emission, retained/mobile equilibration, retained
recovery, and mobile escape.  A process imposes no restriction when its full
asymptotic state change is smaller than the allowed absolute count change.
Constitutive rates and state-update equations are unchanged.
"""
from __future__ import annotations

import math
from typing import Any

import numpy as np

from . import reduced_campaign_front_v9104 as _base
from .reduced_campaign_front_v9104 import ReducedCampaignFront


def _time_for_exact_change(
    asymptotic_change: np.ndarray | float,
    rate_s: np.ndarray | float,
    allowed_change: np.ndarray | float,
) -> float:
    """Return the minimum time at which an exact exponential change hits a bound.

    For ``change(t) = delta_inf * (1 - exp(-rate*t))``, no timestep limit is
    needed if ``delta_inf <= allowed_change``.  Arrays are handled componentwise
    and the smallest positive finite limit is returned.
    """
    delta, rate, allowed = np.broadcast_arrays(
        np.maximum(np.asarray(asymptotic_change, dtype=float), 0.0),
        np.maximum(np.asarray(rate_s, dtype=float), 0.0),
        np.maximum(np.asarray(allowed_change, dtype=float), 0.0),
    )
    active = (rate > 0.0) & (allowed > 0.0) & (delta > allowed)
    if not np.any(active):
        return math.inf
    fraction = np.clip(allowed[active] / delta[active], 0.0, 1.0 - 1.0e-15)
    limits = -np.log1p(-fraction) / rate[active]
    finite = limits[np.isfinite(limits) & (limits > 0.0)]
    return float(np.min(finite)) if finite.size else math.inf


def _limiter_diagnostics(self: ReducedCampaignFront, rates: dict[str, Any]) -> dict[str, float | str]:
    load_limit = self.s.max_dK_substep_MPa_sqrt_m / max(
        self.s.Kdot_MPa_sqrt_m_s, 1.0e-30
    )
    limits: dict[str, float] = {"loading": float(load_limit)}

    lam_c = max(float(rates["lambda_c_s"]), 0.0)
    if lam_c > 0.0:
        limits["cleavage_action"] = float(self.s.max_action_substep / lam_c)
        limits["cleavage_event"] = float(
            max(1.0 - self.B, self.s.event_tolerance) / lam_c
        )
    else:
        limits["cleavage_action"] = math.inf
        limits["cleavage_event"] = math.inf

    if self.plasticity_active:
        capacity = np.maximum(np.asarray(self.capacity, dtype=float), 0.0)
        available = np.maximum(np.asarray(self.available, dtype=float), 0.0)
        mobile = np.maximum(np.asarray(self.mobile, dtype=float), 0.0)
        retained = np.maximum(np.asarray(self.retained, dtype=float), 0.0)

        emit_allowed = max(float(self.s.max_emit_fraction_substep), 0.0) * capacity
        emit_rate = np.maximum(np.asarray(rates["lambda_e_s"], dtype=float), 0.0)
        limits["emission"] = _time_for_exact_change(
            available, emit_rate, emit_allowed
        )

        # Emission is itself bounded during the step.  Include that maximum
        # possible addition when estimating the subsequent exchange state.
        emitted_bound = np.minimum(available, emit_allowed)
        total_bound = mobile + retained + emitted_bound

        encounter = np.maximum(
            np.asarray(rates["encounter_rate_s"], dtype=float), 0.0
        )
        taylor = np.maximum(
            np.asarray(rates["taylor_rate_s"], dtype=float), 0.0
        )
        exchange = encounter + taylor
        retained_fraction_eq = np.divide(
            encounter,
            exchange,
            out=np.zeros_like(exchange),
            where=exchange > 0.0,
        )
        retained_eq_bound = retained_fraction_eq * total_bound
        exchange_delta = np.abs(retained_eq_bound - retained)
        exchange_allowed = (
            max(float(self.s.max_exchange_fraction_substep), 0.0) * capacity
        )
        limits["exchange"] = _time_for_exact_change(
            exchange_delta, exchange, exchange_allowed
        )

        # Recovery follows exchange in the exact plastic operator.  Bound it
        # using the largest retained population reachable in this step.
        retained_bound = np.maximum(retained, retained_eq_bound)
        recovery_rate = max(
            float(self.p.get("retained_recovery_rate_s", 0.0)), 0.0
        )
        limits["recovery"] = _time_for_exact_change(
            retained_bound, recovery_rate, exchange_allowed
        )

        # Escape acts only on mobile content, not on the full active count.
        mobile_eq_bound = np.maximum(total_bound - retained_eq_bound, 0.0)
        mobile_bound = np.maximum(mobile + emitted_bound, mobile_eq_bound)
        escape_rate = np.maximum(
            np.asarray(rates["velocity_m_s"], dtype=float), 0.0
        ) / max(self.s.L_pz_m, 1.0e-30)
        limits["escape"] = _time_for_exact_change(
            mobile_bound, escape_rate, exchange_allowed
        )

    finite_limits = {
        name: value for name, value in limits.items() if np.isfinite(value) and value > 0.0
    }
    if finite_limits:
        limiter = min(finite_limits, key=finite_limits.get)
        selected = float(finite_limits[limiter])
    else:
        limiter = "minimum_dt"
        selected = float(self.s.min_dt_s)

    return {
        **{f"dt_limit_{name}_s": float(value) for name, value in limits.items()},
        "dt_limiter": limiter,
        "dt_selected_s": max(selected, float(self.s.min_dt_s)),
        "available_total": float(np.sum(np.maximum(self.available, 0.0))),
        "mobile_total": float(np.sum(np.maximum(self.mobile, 0.0))),
        "retained_total": float(np.sum(np.maximum(self.retained, 0.0))),
    }


def _state_change_aware_choose_dt(
    self: ReducedCampaignFront, rates: dict[str, Any]
) -> float:
    diagnostic = _limiter_diagnostics(self, rates)
    h = float(diagnostic["dt_selected_s"])
    self._last_dt_diagnostics = diagnostic

    counts = getattr(self, "_dt_limiter_counts", None)
    if counts is None:
        counts = {}
        self._dt_limiter_counts = counts
    limiter = str(diagnostic["dt_limiter"])
    counts[limiter] = int(counts.get(limiter, 0)) + 1
    self._minimum_selected_dt_s = min(
        float(getattr(self, "_minimum_selected_dt_s", math.inf)), h
    )
    return max(h, float(self.s.min_dt_s))


_original_summarize_run = _base.summarize_run


def _diagnostic_summarize_run(
    front: ReducedCampaignFront,
    T_K: float,
    target_events: int,
    internal_steps: int,
) -> dict[str, Any]:
    result = dict(_original_summarize_run(front, T_K, target_events, internal_steps))
    rates = front.instantaneous_rates(front.K, float(T_K))
    diagnostics = dict(getattr(front, "_last_dt_diagnostics", {}))
    limiter_counts = dict(getattr(front, "_dt_limiter_counts", {}))
    if limiter_counts:
        dominant = max(limiter_counts, key=limiter_counts.get)
    else:
        dominant = "none"

    if len(front.events) >= target_events:
        termination = "target_events_reached"
    elif internal_steps >= front.s.max_internal_steps:
        termination = "max_internal_steps"
    elif front.K >= front.s.Kmax_MPa_sqrt_m:
        termination = "Kmax_reached"
    else:
        termination = "minimum_dt_or_unknown"

    result.update(
        {
            "termination_reason": termination,
            "terminal_K_MPa_sqrt_m": float(front.K),
            "terminal_time_s": float(front.time_s),
            "terminal_B": float(front.B),
            "terminal_a_um": float(front.a_um),
            "terminal_available_sites": float(np.sum(front.available)),
            "terminal_mobile_count": float(np.sum(front.mobile)),
            "terminal_retained_count": float(np.sum(front.retained)),
            "terminal_slip_count": float(np.sum(front.slip)),
            "terminal_K_shield_MPa_sqrt_m": float(rates["K_shield_MPa_sqrt_m"]),
            "terminal_sigma_back_max_Pa": float(np.max(rates["sigma_back_Pa"])),
            "terminal_lambda_c_s": float(rates["lambda_c_s"]),
            "terminal_lambda_e_max_s": float(
                np.max(np.asarray(rates["lambda_e_s"], dtype=float))
            ),
            "dominant_dt_limiter": dominant,
            "dt_limiter_counts": limiter_counts,
            "minimum_selected_dt_s": float(
                getattr(front, "_minimum_selected_dt_s", math.nan)
            ),
            **diagnostics,
        }
    )
    return result


ReducedCampaignFront._choose_dt = _state_change_aware_choose_dt
_base.summarize_run = _diagnostic_summarize_run


__all__ = [
    "_diagnostic_summarize_run",
    "_limiter_diagnostics",
    "_state_change_aware_choose_dt",
    "_time_for_exact_change",
]
