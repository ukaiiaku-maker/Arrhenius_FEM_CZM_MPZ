"""Asymptotic-preserving persistent-site transport for v10.0.5.14.5.

v10.0.5.14.4 removed the augmented diagnostic rows and made the high-temperature
transport path conservative.  The 300--800 K production campaign exposed a
second numerical regime: mobile Peierls glide and encounter storage are much
faster than Taylor release, so backward-Euler step doubling resolves the fast
initial layer down to nanoseconds even though the requested constitutive
macrostep is hundreds of seconds.

This release keeps every constitutive rate and finite-volume direction
unchanged.  When a measured timescale-separation gate is satisfied, it applies
the singular limit of the same equations:

* a mobile line undergoes the exact competing-risk sequence of storage versus
  forward transport/escape;
* Taylor release from each retained bin feeds that same mobile absorption map;
* the resulting retained-only substochastic generator is advanced exactly;
* state-dependent forest rates are closed by damped fixed-point iteration.

Outside that separated regime, the v10.0.5.14.4 physical backward-Euler solver
is used unchanged.  All paths remain nonnegative, conservative, and fail closed.
"""
from __future__ import annotations

import copy
from contextlib import contextmanager
from typing import Any, Iterator

import numpy as np
from scipy.linalg import expm

from . import persistent_site_transport_v1005144 as _v144
from .persistent_site_signed_transport_v100514 import (
    PersistentSiteSignedTransportMixin,
)

TRANSPORT_INTEGRATOR = (
    "hybrid_asymptotic_retained_exact__physical_BE_v10_0_5_14_5"
)
ASYMPTOTIC_MODEL = "competing_risk_mobile_elimination_retained_generator_v10_0_5_14_5"


def _stack_state(self, snapshot: dict[str, np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    mobile = np.vstack(
        [
            np.asarray(snapshot["mobile_positive"], dtype=float),
            np.asarray(snapshot["mobile_negative"], dtype=float),
        ]
    )
    retained = np.vstack(
        [
            np.asarray(snapshot["retained_positive"], dtype=float),
            np.asarray(snapshot["retained_negative"], dtype=float),
        ]
    )
    if mobile.shape != retained.shape or mobile.shape != (2 * self.n_systems, self.n_bins):
        raise RuntimeError("persistent-site stacked state has invalid shape")
    return np.maximum(mobile, 0.0), np.maximum(retained, 0.0)


def _unstack_state(
    self, mobile: np.ndarray, retained: np.ndarray
) -> dict[str, np.ndarray]:
    mobile = np.maximum(np.asarray(mobile, dtype=float), 0.0)
    retained = np.maximum(np.asarray(retained, dtype=float), 0.0)
    return {
        "mobile_positive": mobile[: self.n_systems].copy(),
        "mobile_negative": mobile[self.n_systems :].copy(),
        "retained_positive": retained[: self.n_systems].copy(),
        "retained_negative": retained[self.n_systems :].copy(),
    }


def _blend_snapshots(
    self,
    old: dict[str, np.ndarray],
    new: dict[str, np.ndarray],
    fraction: float,
) -> dict[str, np.ndarray]:
    q = float(np.clip(fraction, 0.0, 1.0))
    return {
        name: (1.0 - q) * np.asarray(old[name], float) + q * np.asarray(new[name], float)
        for name in self._TRANSPORT_ARRAY_NAMES
    }


def _rate_fields(
    self,
    snapshot: dict[str, np.ndarray],
    *,
    T_K: float,
    opening_stress_Pa: float,
) -> dict[str, np.ndarray]:
    forest = self._forest_from_snapshot(snapshot)
    radius = self.blunted_radius()
    stress = max(float(opening_stress_Pa), 0.0) * np.sqrt(
        radius / np.maximum(radius + self.x, radius)
    )
    rates = self._pt_model().rates(stress, forest, T_K, self.b_m)
    peierls = np.maximum(
        np.asarray(rates["peierls_rate_s"], dtype=float).reshape(-1), 0.0
    )
    release = np.maximum(
        np.asarray(rates["taylor_completion_rate_s"], dtype=float).reshape(-1),
        0.0,
    )
    jump = np.maximum(
        np.asarray(rates["jump_length_m"], dtype=float).reshape(-1), self.b_m
    )
    if not (peierls.shape == release.shape == jump.shape == (self.n_bins,)):
        raise RuntimeError("Peierls--Taylor transport rates do not match MPZ bins")
    if not (
        np.all(np.isfinite(peierls))
        and np.all(np.isfinite(release))
        and np.all(np.isfinite(jump))
    ):
        raise RuntimeError("nonfinite Peierls--Taylor transport coefficients")
    velocity = jump * peierls
    encounter = (
        float(self.candidate.encounter_efficiency)
        * velocity
        * np.sqrt(np.maximum(forest, 0.0))
    )
    faces = np.empty(self.n_bins + 1, dtype=float)
    faces[1:-1] = 0.5 * (velocity[:-1] + velocity[1:])
    faces[0], faces[-1] = velocity[0], velocity[-1]
    outflow = faces[1:] / max(self.dx, 1.0e-30)
    if not np.all(np.isfinite(encounter)) or not np.all(np.isfinite(outflow)):
        raise RuntimeError("nonfinite asymptotic transport coefficients")
    return {
        "forest": forest,
        "stress": stress,
        "peierls": peierls,
        "release": release,
        "jump": jump,
        "velocity": velocity,
        "encounter": encounter,
        "outflow": outflow,
    }


def _mobile_absorption_probabilities(
    outflow_rate_s: np.ndarray,
    encounter_rate_s: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return retained-bin and escape probabilities for each mobile start bin."""
    outflow = np.maximum(np.asarray(outflow_rate_s, dtype=float).reshape(-1), 0.0)
    encounter = np.maximum(
        np.asarray(encounter_rate_s, dtype=float).reshape(-1), 0.0
    )
    if outflow.shape != encounter.shape:
        raise ValueError("outflow and encounter arrays must have the same shape")
    n = int(outflow.size)
    retained = np.zeros((n, n), dtype=float)
    escaped = np.zeros(n, dtype=float)
    for start in range(n):
        survival = 1.0
        for j in range(start, n):
            total = outflow[j] + encounter[j]
            if total <= 0.0:
                # A genuinely immobile bin has no separated mobile limit.  The
                # selection gate rejects this path whenever it is reachable.
                break
            p_store = encounter[j] / total
            p_move = outflow[j] / total
            retained[start, j] = survival * p_store
            survival *= p_move
            if survival <= 1.0e-300:
                survival = 0.0
                break
        escaped[start] = survival
        row_sum = float(np.sum(retained[start]) + escaped[start])
        if row_sum > 0.0:
            retained[start] /= row_sum
            escaped[start] /= row_sum
    return retained, escaped


def _retained_generator(
    absorption: np.ndarray,
    escape_probability: np.ndarray,
    release_rate_s: np.ndarray,
) -> np.ndarray:
    """Build the retained-only column generator after fast mobile elimination."""
    P = np.maximum(np.asarray(absorption, dtype=float), 0.0)
    p_escape = np.maximum(np.asarray(escape_probability, dtype=float), 0.0)
    release = np.maximum(np.asarray(release_rate_s, dtype=float), 0.0)
    n = int(release.size)
    if P.shape != (n, n) or p_escape.shape != (n,):
        raise ValueError("retained-generator probability shapes do not match")
    Q = np.zeros((n, n), dtype=float)
    for source in range(n):
        for target in range(source + 1, n):
            Q[target, source] += release[source] * P[source, target]
        Q[source, source] -= release[source] * (
            p_escape[source] + float(np.sum(P[source, source + 1 :]))
        )
    return Q


def _selection_metrics(
    self,
    snapshot: dict[str, np.ndarray],
    fields: dict[str, np.ndarray],
    dt_s: float,
) -> dict[str, Any]:
    mobile, retained = _stack_state(self, snapshot)
    release = fields["release"]
    fast = fields["outflow"] + fields["encounter"]
    total_mass = max(float(np.sum(mobile) + np.sum(retained)), 1.0e-300)
    # Retained content only generates mobile content in proportion to its
    # release exposure over the macrostep.
    release_exposure = -np.expm1(-np.minimum(release * max(float(dt_s), 0.0), 700.0))
    source_weight = np.sum(mobile, axis=0) + np.sum(retained, axis=0) * release_exposure
    source_weight /= total_mass

    P, _ = _mobile_absorption_probabilities(fields["outflow"], fields["encounter"])
    reach = np.zeros(self.n_bins, dtype=float)
    for start in np.where(source_weight > 1.0e-16)[0]:
        probability_reaching = 1.0
        for j in range(int(start), self.n_bins):
            reach[j] = max(reach[j], float(source_weight[start]) * probability_reaching)
            total = fast[j]
            if total <= 0.0:
                break
            probability_reaching *= fields["outflow"][j] / total
            if probability_reaching * float(source_weight[start]) <= 1.0e-16:
                break

    reach_floor = max(
        float(getattr(self, "transport_qss_reach_fraction", 1.0e-12)), 1.0e-16
    )
    active = reach > reach_floor
    if not np.any(active):
        return {
            "selected": False,
            "reason": "no_reachable_mobile_source",
            "fast_courant_min": 0.0,
            "separation_min": 0.0,
            "reachable_bin_count": 0,
            "absorption": P,
        }
    fast_courant = fast[active] * max(float(dt_s), 0.0)
    separation = np.divide(
        fast[active],
        release[active],
        out=np.full(np.count_nonzero(active), np.inf),
        where=release[active] > 0.0,
    )
    min_courant = float(np.min(fast_courant))
    min_separation = float(np.min(separation))
    required_courant = max(
        float(getattr(self, "transport_qss_fast_courant_min", 50.0)), 1.0
    )
    required_separation = max(
        float(getattr(self, "transport_qss_separation_min", 20.0)), 1.0
    )
    selected = bool(
        np.all(fast[active] > 0.0)
        and min_courant >= required_courant
        and min_separation >= required_separation
    )
    return {
        "selected": selected,
        "reason": "timescale_separation" if selected else "separation_gate_not_met",
        "fast_courant_min": min_courant,
        "separation_min": min_separation,
        "required_fast_courant": required_courant,
        "required_separation": required_separation,
        "reachable_bin_count": int(np.count_nonzero(active)),
        "absorption": P,
    }


def _asymptotic_apply(
    self,
    initial_snapshot: dict[str, np.ndarray],
    coefficient_snapshot: dict[str, np.ndarray],
    *,
    dt_s: float,
    T_K: float,
    opening_stress_Pa: float,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    fields = _rate_fields(
        self,
        coefficient_snapshot,
        T_K=T_K,
        opening_stress_Pa=opening_stress_Pa,
    )
    selection = _selection_metrics(self, initial_snapshot, fields, dt_s)
    if not bool(selection["selected"]):
        raise RuntimeError(
            "asymptotic transport lost its timescale-separation gate during closure: "
            f"reason={selection['reason']} fast_courant_min="
            f"{selection['fast_courant_min']:.6e} separation_min="
            f"{selection['separation_min']:.6e}"
        )
    P = np.asarray(selection["absorption"], dtype=float)
    _, p_escape = _mobile_absorption_probabilities(
        fields["outflow"], fields["encounter"]
    )
    Q = _retained_generator(P, p_escape, fields["release"])

    mobile0, retained0 = _stack_state(self, initial_snapshot)
    retained_start = retained0 + mobile0 @ P
    initial_mobile_escape = float(np.sum(mobile0 @ p_escape))
    propagator = np.asarray(expm(max(float(dt_s), 0.0) * Q), dtype=float)
    retained_final = np.asarray((propagator @ retained_start.T).T, dtype=float)
    if not np.all(np.isfinite(retained_final)):
        raise RuntimeError("asymptotic retained propagator produced nonfinite state")
    magnitude = max(float(np.max(np.abs(retained_final))), float(np.max(retained_start)), 1.0e-300)
    negative_tolerance = 2.0e-11 * magnitude + 1.0e-300
    minimum = float(np.min(retained_final))
    if minimum < -negative_tolerance:
        raise RuntimeError(
            "asymptotic retained propagator violated nonnegative state: "
            f"minimum={minimum:.6e}, tolerance={negative_tolerance:.6e}"
        )
    retained_final = np.maximum(retained_final, 0.0)
    mobile_final = np.zeros_like(mobile0)
    result = _unstack_state(self, mobile_final, retained_final)

    initial_mass = self._snapshot_mass(initial_snapshot)
    final_mass = self._snapshot_mass(result)
    gain = final_mass - initial_mass
    gain_tolerance = 2.0e-10 * max(initial_mass, 1.0e-300) + 1.0e-300
    if gain > gain_tolerance:
        raise RuntimeError(
            "asymptotic persistent-site transport created line content: "
            f"gain={gain:.6e}, scale={initial_mass:.6e}"
        )
    escaped = max(initial_mass - final_mass, 0.0)
    initial_retained = float(np.sum(retained0))
    final_retained = float(np.sum(retained_final))
    signed_balance = initial_mass - final_mass - escaped
    diag: dict[str, Any] = {
        "dN_trapped": max(final_retained - initial_retained, 0.0),
        "dN_detrapped": max(initial_retained - final_retained, 0.0),
        "dN_escaped": escaped,
        "initial_mobile_direct_escape": initial_mobile_escape,
        "peierls_rate_min_s": float(np.min(fields["peierls"])),
        "peierls_rate_max_s": float(np.max(fields["peierls"])),
        "taylor_completion_rate_min_s": float(np.min(fields["release"])),
        "taylor_completion_rate_max_s": float(np.max(fields["release"])),
        "encounter_rate_min_s": float(np.min(fields["encounter"])),
        "encounter_rate_max_s": float(np.max(fields["encounter"])),
        "glide_velocity_max_m_s": float(np.max(fields["velocity"])),
        "rho_forest_min_m2": float(np.min(fields["forest"])),
        "rho_forest_max_m2": float(np.max(fields["forest"])),
        "max_frozen_courant": float(
            np.max(fields["velocity"]) * max(float(dt_s), 0.0) / max(self.dx, 1.0e-30)
        ),
        "line_content_conservation_error": abs(signed_balance),
        "line_content_conservation_signed_error": signed_balance,
        "line_content_conservation_scale": max(initial_mass, final_mass + escaped, 1.0e-300),
        "physical_generator_mass_gain": max(gain, 0.0),
        "physical_generator_action": ASYMPTOTIC_MODEL,
        "transport_asymptotic_fast_courant_min": selection["fast_courant_min"],
        "transport_asymptotic_separation_min": selection["separation_min"],
        "transport_asymptotic_reachable_bins": selection["reachable_bin_count"],
    }
    return result, diag


def _mechanics_state_error(
    self,
    first: dict[str, np.ndarray],
    second: dict[str, np.ndarray],
    scale: float,
) -> float:
    m1, r1 = _stack_state(self, first)
    m2, r2 = _stack_state(self, second)
    unsigned1 = np.sum(m1 + r1, axis=0)
    unsigned2 = np.sum(m2 + r2, axis=0)
    signed1 = np.sum(r1[: self.n_systems] - r1[self.n_systems :], axis=0)
    signed2 = np.sum(r2[: self.n_systems] - r2[self.n_systems :], axis=0)
    state_error = self._snapshot_difference(first, second) / scale
    unsigned_error = float(np.sum(np.abs(unsigned1 - unsigned2))) / scale
    signed_error = float(np.sum(np.abs(signed1 - signed2))) / scale
    return max(state_error, unsigned_error, signed_error)


def _transport_hybrid(
    self, *, dt_s: float, T_K: float, opening_stress_Pa: float
) -> dict[str, Any]:
    dt = max(float(dt_s), 0.0)
    initial = self._transport_snapshot()
    initial_mass = self._snapshot_mass(initial)
    if dt <= 0.0 or initial_mass <= 0.0:
        out = _v144._transport_physical(
            self, dt_s=dt, T_K=T_K, opening_stress_Pa=opening_stress_Pa
        )
        out["transport_integrator"] = TRANSPORT_INTEGRATOR
        out["transport_asymptotic_active"] = False
        out["transport_asymptotic_reason"] = "empty_or_zero_interval"
        return out

    fields = _rate_fields(
        self, initial, T_K=T_K, opening_stress_Pa=opening_stress_Pa
    )
    selection = _selection_metrics(self, initial, fields, dt)
    if not bool(selection["selected"]):
        out = _v144._transport_physical(
            self, dt_s=dt, T_K=T_K, opening_stress_Pa=opening_stress_Pa
        )
        out["transport_integrator"] = TRANSPORT_INTEGRATOR
        out["transport_asymptotic_active"] = False
        out["transport_asymptotic_reason"] = selection["reason"]
        out["transport_asymptotic_fast_courant_min"] = selection[
            "fast_courant_min"
        ]
        out["transport_asymptotic_separation_min"] = selection["separation_min"]
        return out

    rtol = max(float(getattr(self, "transport_nonlinear_rtol", 1.0e-3)), 1.0e-8)
    max_iterations = max(
        int(getattr(self, "transport_qss_max_iterations", 80)), 8
    )
    damping = float(
        np.clip(getattr(self, "transport_qss_damping", 0.5), 0.05, 1.0)
    )
    guess = copy.deepcopy(initial)
    maximum_residual = 0.0
    final_residual = float("inf")
    final = None
    final_diag: dict[str, Any] | None = None
    for iteration in range(1, max_iterations + 1):
        candidate, diag = _asymptotic_apply(
            self,
            initial,
            guess,
            dt_s=dt,
            T_K=T_K,
            opening_stress_Pa=opening_stress_Pa,
        )
        residual = _mechanics_state_error(
            self, candidate, guess, max(initial_mass, 1.0e-300)
        )
        maximum_residual = max(maximum_residual, residual)
        final_residual = residual
        if residual <= rtol:
            final = candidate
            final_diag = diag
            break
        guess = _blend_snapshots(self, guess, candidate, damping)
    if final is None or final_diag is None:
        raise RuntimeError(
            "persistent-site asymptotic closure failed to converge: "
            f"iterations={max_iterations}, residual={final_residual:.6e}, "
            f"rtol={rtol:.6e}, max_residual={maximum_residual:.6e}, "
            f"fast_courant_min={selection['fast_courant_min']:.6e}, "
            f"separation_min={selection['separation_min']:.6e}"
        )

    self._restore_transport_snapshot(final)
    escaped = float(final_diag["dN_escaped"])
    self.time_s += dt
    self.escaped_total += escaped
    out = {
        "dN_trapped": float(final_diag["dN_trapped"]),
        "dN_detrapped": float(final_diag["dN_detrapped"]),
        "dN_escaped": escaped,
        "dN_recovered": 0.0,
        "transport_substeps": 1,
        "transport_attempted_physical_solves": 0,
        "transport_attempted_exponentials": iteration,
        "transport_attempted_linear_solves": 0,
        "transport_rejected_intervals": 0,
        "transport_tail_accepts": 0,
        "transport_nonlinear_error_max": maximum_residual,
        "transport_nonlinear_error_final": final_residual,
        "transport_nonlinear_rtol": rtol,
        "transport_integrator": TRANSPORT_INTEGRATOR,
        "transport_asymptotic_active": True,
        "transport_asymptotic_model": ASYMPTOTIC_MODEL,
        "transport_asymptotic_reason": selection["reason"],
        "transport_asymptotic_iterations": iteration,
        "transport_asymptotic_damping": damping,
        "transport_cfl_limited": False,
        "explicit_recovery_active": False,
        **{
            key: value
            for key, value in final_diag.items()
            if key not in {"dN_trapped", "dN_detrapped", "dN_escaped"}
        },
    }
    self.last_transport = copy.deepcopy(out)
    return out


@contextmanager
def installed_asymptotic_transport_v1005145() -> Iterator[None]:
    with _v144.installed_split_transport_v1005144():
        old_transport = PersistentSiteSignedTransportMixin.transport
        PersistentSiteSignedTransportMixin.transport = _transport_hybrid
        try:
            yield
        finally:
            PersistentSiteSignedTransportMixin.transport = old_transport


__all__ = [
    "TRANSPORT_INTEGRATOR",
    "ASYMPTOTIC_MODEL",
    "_mobile_absorption_probabilities",
    "_retained_generator",
    "installed_asymptotic_transport_v1005145",
]
