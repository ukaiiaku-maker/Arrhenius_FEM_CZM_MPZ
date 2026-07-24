"""Exact PF v10.2.22 local signed-MPZ update map for FEM/CZM v10.0.5.15.

This module intentionally reproduces the discrete constitutive map used by the
active PF benchmark, rather than integrating a different coupled continuum
operator.  For one accepted local interval the order is

    persistent emission
    -> exact mobile/retained exchange
    -> zero explicit recovery
    -> population-weighted scalar Peierls advection
    -> absorbing escape
    -> wake exchange/advection.

The functions below are direct ports of ``UnifiedMPZState._exchange``,
``UnifiedMPZState._advect_forward``, and the signed scalar update installed by
``signed_burgers_shared_v1025.py`` at PF commit
198ece3aeb1d193a8c1c4857676fba720c088d27.
"""
from __future__ import annotations

import copy
import math
from typing import Any, Callable

import numpy as np

PF_REFERENCE_COMMIT = "198ece3aeb1d193a8c1c4857676fba720c088d27"
PF_UPDATE_MAP = "pf_v10_2_22_emit_exchange_zero_recovery_scalar_advection"


def exact_exchange(
    mobile: np.ndarray,
    retained: np.ndarray,
    encounter_rate_s: np.ndarray,
    release_rate_s: np.ndarray,
    dt_s: float,
) -> tuple[np.ndarray, np.ndarray, float, float]:
    """PF analytic two-state mobile/retained exchange, bin by bin."""
    m = np.maximum(np.asarray(mobile, dtype=float), 0.0)
    r = np.maximum(np.asarray(retained, dtype=float), 0.0)
    ke = np.maximum(np.asarray(encounter_rate_s, dtype=float), 0.0)[None, :]
    kr = np.maximum(np.asarray(release_rate_s, dtype=float), 0.0)[None, :]
    if m.shape != r.shape or m.ndim != 2:
        raise ValueError("mobile and retained arrays must have equal 2-D shapes")
    if ke.shape[1] != m.shape[1] or kr.shape[1] != m.shape[1]:
        raise ValueError("exchange-rate vectors must match the MPZ bin count")
    total = m + r
    rate = ke + kr
    equilibrium_retained = np.divide(
        ke,
        rate,
        out=np.zeros_like(rate),
        where=rate > 0.0,
    ) * total
    decay = np.exp(-np.minimum(rate * max(float(dt_s), 0.0), 700.0))
    new_retained = np.clip(
        equilibrium_retained + (r - equilibrium_retained) * decay,
        0.0,
        total,
    )
    new_mobile = total - new_retained
    trapped = float(np.sum(np.maximum(new_retained - r, 0.0)))
    released = float(np.sum(np.maximum(r - new_retained, 0.0)))
    return new_mobile, new_retained, trapped, released


def fractional_advect_forward(
    field: np.ndarray,
    distance_m: float,
    dx_m: float,
) -> tuple[np.ndarray, float]:
    """PF piecewise-linear cell-index remap away from the crack tip."""
    values = np.maximum(np.asarray(field, dtype=float), 0.0)
    distance = max(float(distance_m), 0.0)
    dx = max(float(dx_m), 1.0e-30)
    if values.ndim != 2:
        raise ValueError("advected field must be a 2-D system-by-bin array")
    if distance <= 0.0:
        return values.copy(), 0.0
    shift = distance / dx
    whole = int(math.floor(shift))
    fraction = shift - whole
    out = np.zeros_like(values)
    n_bins = values.shape[1]
    for index in range(n_bins):
        target = index + whole
        if target < n_bins:
            out[:, target] += values[:, index] * (1.0 - fraction)
        if fraction > 0.0 and target + 1 < n_bins:
            out[:, target + 1] += values[:, index] * fraction
    escaped = max(float(np.sum(values) - np.sum(out)), 0.0)
    return out, escaped


def _active_mass(state) -> float:
    return float(
        np.sum(state.mobile_positive)
        + np.sum(state.mobile_negative)
        + np.sum(state.retained_positive)
        + np.sum(state.retained_negative)
    )


def _wake_mass(state) -> float:
    return float(
        np.sum(state.wake_mobile_positive)
        + np.sum(state.wake_mobile_negative)
        + np.sum(state.wake_retained_positive)
        + np.sum(state.wake_retained_negative)
    )


def _retained_forest_density_m2(state, *, wake: bool) -> np.ndarray:
    """PF shared forest density: retained unsigned content only."""
    if wake:
        retained = state.wake_retained_positive + state.wake_retained_negative
        dx = float(state.wake_dx)
    else:
        retained = state.retained_positive + state.retained_negative
        dx = float(state.dx)
    width = max(float(state.blunting_length_m), dx, 1.0e-12)
    content = np.sum(np.maximum(retained, 0.0), axis=0)
    return np.maximum(
        float(state.candidate.rho_forest_floor_m2)
        + content / max(dx * width, 1.0e-30),
        1.0,
    )


def _stress_profile_Pa(state, opening_stress_Pa: float) -> np.ndarray:
    reference = max(float(state.blunting_length_m), float(state.dx), 1.0e-12)
    return max(float(opening_stress_Pa), 0.0) * np.sqrt(
        reference / np.maximum(reference + np.asarray(state.x, dtype=float), reference)
    )


def _rates(state, stress_profile_Pa: np.ndarray, rho_m2: np.ndarray, T_K: float):
    raw = state._pt_model().rates(
        np.asarray(stress_profile_Pa, dtype=float),
        np.asarray(rho_m2, dtype=float),
        float(T_K),
        float(state.b_m),
    )
    peierls = np.maximum(
        np.asarray(raw["peierls_rate_s"], dtype=float).reshape(-1), 0.0
    )
    release = np.maximum(
        np.asarray(raw["taylor_completion_rate_s"], dtype=float).reshape(-1),
        0.0,
    )
    jump = np.maximum(
        np.asarray(raw["jump_length_m"], dtype=float).reshape(-1),
        float(state.b_m),
    )
    expected = (int(state.n_bins),)
    if peierls.shape != expected or release.shape != expected or jump.shape != expected:
        raise RuntimeError("PF transport rates do not match the runtime MPZ grid")
    if not (
        np.all(np.isfinite(peierls))
        and np.all(np.isfinite(release))
        and np.all(np.isfinite(jump))
    ):
        raise RuntimeError("PF transport produced nonfinite rate coefficients")
    velocity = jump * peierls
    encounter = (
        max(float(state.candidate.encounter_efficiency), 0.0)
        * velocity
        * np.sqrt(np.maximum(rho_m2, 0.0))
    )
    return peierls, release, jump, velocity, encounter, raw


def _exchange_signed_species(
    state,
    encounter: np.ndarray,
    release: np.ndarray,
    dt_s: float,
    *,
    wake: bool,
) -> tuple[float, float]:
    prefix = "wake_" if wake else ""
    trapped = 0.0
    released = 0.0
    for sign in ("positive", "negative"):
        mobile_name = f"{prefix}mobile_{sign}"
        retained_name = f"{prefix}retained_{sign}"
        mobile, retained, dtrap, drelease = exact_exchange(
            getattr(state, mobile_name),
            getattr(state, retained_name),
            encounter,
            release,
            dt_s,
        )
        setattr(state, mobile_name, mobile)
        setattr(state, retained_name, retained)
        trapped += dtrap
        released += drelease
    return trapped, released


def _population_weighted_velocity(
    state,
    velocity_by_bin_m_s: np.ndarray,
    *,
    wake: bool,
) -> float:
    if wake:
        mobile = state.wake_mobile_positive + state.wake_mobile_negative
    else:
        mobile = state.mobile_positive + state.mobile_negative
    by_bin = np.sum(np.maximum(mobile, 0.0), axis=0)
    total = float(np.sum(by_bin))
    if total > 0.0:
        return max(float(np.sum(velocity_by_bin_m_s * by_bin) / total), 0.0)
    if wake:
        return 0.0
    nsrc = max(min(int(state._source_bin_count()), int(state.n_bins)), 1)
    return max(float(np.mean(velocity_by_bin_m_s[:nsrc])), 0.0)


def _evolve_wake_pf(state, *, dt_s: float, T_K: float) -> dict[str, float]:
    dt = max(float(dt_s), 0.0)
    if dt <= 0.0 or _wake_mass(state) <= 0.0:
        return {
            "wake_dN_trapped": 0.0,
            "wake_dN_released": 0.0,
            "wake_dN_recovered": 0.0,
            "wake_dN_mobile_transport_loss": 0.0,
        }
    rho = _retained_forest_density_m2(state, wake=True)
    zero_stress = np.zeros(int(state.wake_n_bins), dtype=float)
    peierls, release, jump, velocity, encounter, _ = _rates(
        state, zero_stress, rho, T_K
    )
    trapped, released = _exchange_signed_species(
        state, encounter, release, dt, wake=True
    )
    scalar_velocity = _population_weighted_velocity(state, velocity, wake=True)
    lost = 0.0
    for sign in ("positive", "negative"):
        name = f"wake_mobile_{sign}"
        moved, amount = fractional_advect_forward(
            getattr(state, name), scalar_velocity * dt, state.wake_dx
        )
        setattr(state, name, moved)
        lost += amount
    state.wake_discarded_total += lost
    return {
        "wake_dN_trapped": trapped,
        "wake_dN_released": released,
        "wake_dN_recovered": 0.0,
        "wake_dN_mobile_transport_loss": lost,
        "wake_glide_velocity_m_s": scalar_velocity,
        "wake_peierls_rate_max_s": float(np.max(peierls)),
    }


def evolve_pf_v10222(
    state,
    *,
    dt_s: float,
    T_K: float,
    opening_stress_Pa: float,
    drive_factors: np.ndarray,
    tau_signed_Pa: np.ndarray,
    emission_rate_function: Callable[[float, float], float],
) -> dict[str, Any]:
    """Apply one complete PF v10.2.22 persistent signed-state update."""
    dt = max(float(dt_s), 0.0)
    initial_inventory = _active_mass(state) + _wake_mass(state)
    escaped_before = float(state.escaped_total)
    discarded_before = float(state.wake_discarded_total)

    emission = state.emit_persistent(
        dt_s=dt,
        T_K=float(T_K),
        opening_stress_Pa=float(opening_stress_Pa),
        drive_factors=np.asarray(drive_factors, dtype=float),
        tau_signed_Pa=np.asarray(tau_signed_Pa, dtype=float),
        rate_function=emission_rate_function,
    )
    emitted = float(emission.get("dN_emit", 0.0))

    stress = _stress_profile_Pa(state, opening_stress_Pa)
    rho = _retained_forest_density_m2(state, wake=False)
    peierls, release, jump, velocity, encounter, raw = _rates(
        state, stress, rho, T_K
    )
    trapped, released = _exchange_signed_species(
        state, encounter, release, dt, wake=False
    )

    scalar_velocity = _population_weighted_velocity(state, velocity, wake=False)
    escaped = 0.0
    for sign in ("positive", "negative"):
        name = f"mobile_{sign}"
        moved, amount = fractional_advect_forward(
            getattr(state, name), scalar_velocity * dt, state.dx
        )
        setattr(state, name, moved)
        escaped += amount
    state.escaped_total += escaped
    state.time_s += dt
    wake = _evolve_wake_pf(state, dt_s=dt, T_K=T_K)

    final_inventory = _active_mass(state) + _wake_mass(state)
    escaped_increment = float(state.escaped_total) - escaped_before
    discarded_increment = float(state.wake_discarded_total) - discarded_before
    balance = initial_inventory + emitted - (
        final_inventory + escaped_increment + discarded_increment
    )
    scale = max(
        abs(initial_inventory + emitted),
        abs(final_inventory + escaped_increment + discarded_increment),
        1.0e-300,
    )
    relative_error = abs(balance) / scale
    if relative_error > 5.0e-11:
        raise RuntimeError(
            "PF update-map line-content conservation failed: "
            f"balance={balance:.9e}, relative_error={relative_error:.9e}"
        )

    diagnostics = {
        **emission,
        "dN_trapped": trapped,
        "dN_released": released,
        "dN_detrapped": released,
        "dN_recovered": 0.0,
        "dN_escaped": escaped,
        "peierls_rate_s": float(np.max(peierls)),
        "peierls_rate_min_s": float(np.min(peierls)),
        "taylor_completion_rate_s": float(np.max(release)),
        "taylor_completion_rate_min_s": float(np.min(release)),
        "encounter_rate_s": float(np.max(encounter)),
        "encounter_rate_min_s": float(np.min(encounter)),
        "glide_velocity_m_s": scalar_velocity,
        "glide_velocity_bin_max_m_s": float(np.max(velocity)),
        "jump_length_bin_max_m": float(np.max(jump)),
        "rho_forest_min_m2": float(np.min(rho)),
        "rho_forest_max_m2": float(np.max(rho)),
        "transport_integrator": PF_UPDATE_MAP,
        "transport_operator_order": (
            "emit_then_exact_exchange_then_zero_recovery_then_scalar_advection"
        ),
        "transport_substeps": 1,
        "transport_cfl_limited": False,
        "explicit_recovery_active": False,
        "pf_reference_commit": PF_REFERENCE_COMMIT,
        "line_content_balance_signed": balance,
        "line_content_balance_relative_error": relative_error,
        **wake,
    }
    state.last_transport = copy.deepcopy(diagnostics)
    return diagnostics


__all__ = [
    "PF_REFERENCE_COMMIT",
    "PF_UPDATE_MAP",
    "exact_exchange",
    "fractional_advect_forward",
    "evolve_pf_v10222",
]
