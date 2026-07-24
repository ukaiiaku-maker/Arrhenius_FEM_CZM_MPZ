"""Conservative stiff MPZ transport for v10.0.5.14.4.

The Peierls, encounter, Taylor-release, upwind transport, and absorbing-boundary
laws are unchanged.  Only the exact frozen physical state [mobile, retained] is
advanced.  Diagnostic accumulator rows are excluded.  The matrix exponential
is applied with sparse Krylov ``expm_multiply`` so strongly nonuniform real-rate
operators do not require formation of a dense propagator.
"""
from __future__ import annotations

import copy
from contextlib import contextmanager
from typing import Any, Iterator

import numpy as np
from scipy.linalg import expm
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import expm_multiply

from .persistent_site_signed_transport_v100514 import (
    PersistentSiteSignedTransportMixin,
)

TRANSPORT_INTEGRATOR = (
    "adaptive_physical_generator_krylov_exponential_v10_0_5_14_4"
)


def _exact_exchange_pair(
    mobile: np.ndarray,
    retained: np.ndarray,
    encounter_rate_s: np.ndarray,
    release_rate_s: np.ndarray,
    dt_s: float,
) -> tuple[np.ndarray, np.ndarray, float, float]:
    """Exact nonnegative local M <-> R solution, retained for audit tests."""
    m0 = np.maximum(np.asarray(mobile, dtype=float), 0.0)
    r0 = np.maximum(np.asarray(retained, dtype=float), 0.0)
    ke = np.maximum(np.asarray(encounter_rate_s, dtype=float), 0.0)[None, :]
    kr = np.maximum(np.asarray(release_rate_s, dtype=float), 0.0)[None, :]
    total = m0 + r0
    rate = ke + kr
    equilibrium = np.divide(ke * total, rate, out=np.zeros_like(total), where=rate > 0.0)
    decay = np.exp(-np.minimum(rate * max(float(dt_s), 0.0), 700.0))
    r1 = np.where(rate > 0.0, equilibrium + (r0 - equilibrium) * decay, r0)
    r1 = np.clip(r1, 0.0, total)
    m1 = total - r1
    return (
        m1,
        r1,
        float(np.sum(np.maximum(r1 - r0, 0.0))),
        float(np.sum(np.maximum(r0 - r1, 0.0))),
    )


def _advection_generator(velocity_m_s: np.ndarray, dx_m: float) -> np.ndarray:
    velocity = np.maximum(np.asarray(velocity_m_s, dtype=float).reshape(-1), 0.0)
    n = int(velocity.size)
    inv_dx = 1.0 / max(float(dx_m), 1.0e-30)
    faces = np.empty(n + 1, dtype=float)
    faces[1:-1] = 0.5 * (velocity[:-1] + velocity[1:])
    faces[0], faces[-1] = velocity[0], velocity[-1]
    A = np.zeros((n, n), dtype=float)
    for i in range(n):
        A[i, i] -= faces[i + 1] * inv_dx
        if i > 0:
            A[i, i - 1] += faces[i] * inv_dx
    return A


def _advect_mobile_exact(
    mobile: np.ndarray,
    velocity_m_s: np.ndarray,
    dx_m: float,
    dt_s: float,
) -> tuple[np.ndarray, float, float]:
    source = np.maximum(np.asarray(mobile, dtype=float), 0.0)
    dt = max(float(dt_s), 0.0)
    if dt <= 0.0 or float(np.sum(source)) <= 0.0:
        return source.copy(), 0.0, 0.0
    advanced = np.asarray(source @ expm(dt * _advection_generator(velocity_m_s, dx_m)).T)
    if not np.all(np.isfinite(advanced)):
        raise RuntimeError("physical exponential advection produced nonfinite state")
    magnitude = max(float(np.max(np.abs(advanced))), float(np.max(source)), 1.0e-300)
    tolerance = 2.0e-12 * magnitude + 1.0e-300
    if float(np.min(advanced)) < -tolerance:
        raise RuntimeError(
            "physical exponential advection violated nonnegative state: "
            f"minimum={float(np.min(advanced)):.6e}, tolerance={tolerance:.6e}"
        )
    advanced = np.maximum(advanced, 0.0)
    initial_mass = float(np.sum(source))
    final_mass = float(np.sum(advanced))
    gain = final_mass - initial_mass
    gain_tolerance = 2.0e-11 * max(initial_mass, 1.0e-300) + 1.0e-300
    if gain > gain_tolerance:
        raise RuntimeError(
            f"physical exponential advection created line content: gain={gain:.6e}"
        )
    return advanced, max(initial_mass - final_mass, 0.0), max(gain, 0.0)


def _physical_generator(
    velocity_m_s: np.ndarray,
    encounter_rate_s: np.ndarray,
    release_rate_s: np.ndarray,
    dx_m: float,
) -> csc_matrix:
    velocity = np.maximum(np.asarray(velocity_m_s, dtype=float).reshape(-1), 0.0)
    encounter = np.maximum(np.asarray(encounter_rate_s, dtype=float).reshape(-1), 0.0)
    release = np.maximum(np.asarray(release_rate_s, dtype=float).reshape(-1), 0.0)
    if not (velocity.shape == encounter.shape == release.shape):
        raise ValueError("physical transport coefficient shapes do not match")
    n = int(velocity.size)
    inv_dx = 1.0 / max(float(dx_m), 1.0e-30)
    faces = np.empty(n + 1, dtype=float)
    faces[1:-1] = 0.5 * (velocity[:-1] + velocity[1:])
    faces[0], faces[-1] = velocity[0], velocity[-1]
    A = np.zeros((2 * n, 2 * n), dtype=float)
    for i in range(n):
        A[i, i] -= faces[i + 1] * inv_dx + encounter[i]
        if i > 0:
            A[i, i - 1] += faces[i] * inv_dx
        A[i, n + i] += release[i]
        A[n + i, i] += encounter[i]
        A[n + i, n + i] -= release[i]
    return csc_matrix(A)


def _frozen_transport_step_physical(
    self,
    snapshot: dict[str, np.ndarray],
    *,
    dt_s: float,
    T_K: float,
    opening_stress_Pa: float,
) -> tuple[dict[str, np.ndarray], dict[str, float]]:
    dt = max(float(dt_s), 0.0)
    initial_mass = self._snapshot_mass(snapshot)
    if dt <= 0.0 or initial_mass <= 0.0:
        return copy.deepcopy(snapshot), {
            "dN_trapped": 0.0,
            "dN_detrapped": 0.0,
            "dN_escaped": 0.0,
            "max_frozen_courant": 0.0,
            "line_content_conservation_error": 0.0,
            "physical_generator_mass_gain": 0.0,
        }

    forest = self._forest_from_snapshot(snapshot)
    radius = self.blunted_radius()
    stress = max(float(opening_stress_Pa), 0.0) * np.sqrt(
        radius / np.maximum(radius + self.x, radius)
    )
    rates = self._pt_model().rates(stress, forest, T_K, self.b_m)
    peierls = np.maximum(np.asarray(rates["peierls_rate_s"], float).reshape(-1), 0.0)
    release = np.maximum(
        np.asarray(rates["taylor_completion_rate_s"], float).reshape(-1), 0.0
    )
    jump = np.maximum(np.asarray(rates["jump_length_m"], float).reshape(-1), self.b_m)
    if not (peierls.shape == release.shape == jump.shape == (self.n_bins,)):
        raise RuntimeError("Peierls--Taylor transport rates do not match MPZ bins")
    if not (np.all(np.isfinite(peierls)) and np.all(np.isfinite(release)) and np.all(np.isfinite(jump))):
        raise RuntimeError("nonfinite Peierls--Taylor transport coefficients")
    velocity = jump * peierls
    encounter = (
        float(self.candidate.encounter_efficiency)
        * velocity
        * np.sqrt(np.maximum(forest, 0.0))
    )
    if not np.all(np.isfinite(encounter)):
        raise RuntimeError("nonfinite encounter-storage coefficients")

    n = self.n_bins
    initial = np.zeros((2 * n, 2 * self.n_systems), dtype=float)
    initial_retained = 0.0
    for system in range(self.n_systems):
        initial[:n, system] = snapshot["mobile_positive"][system]
        initial[n:, system] = snapshot["retained_positive"][system]
        column = self.n_systems + system
        initial[:n, column] = snapshot["mobile_negative"][system]
        initial[n:, column] = snapshot["retained_negative"][system]
        initial_retained += float(np.sum(snapshot["retained_positive"][system]))
        initial_retained += float(np.sum(snapshot["retained_negative"][system]))

    operator = dt * _physical_generator(velocity, encounter, release, self.dx)
    advanced = np.asarray(expm_multiply(operator, initial), dtype=float)
    if not np.all(np.isfinite(advanced)):
        raise RuntimeError("physical Krylov transport produced nonfinite state")
    magnitude = max(float(np.max(np.abs(advanced))), float(np.max(initial)), 1.0e-300)
    negative_tolerance = 5.0e-10 * magnitude + 1.0e-300
    minimum = float(np.min(advanced))
    if minimum < -negative_tolerance:
        raise RuntimeError(
            "physical Krylov transport violated nonnegative state: "
            f"minimum={minimum:.6e}, tolerance={negative_tolerance:.6e}"
        )
    advanced = np.maximum(advanced, 0.0)

    result = {name: np.zeros_like(snapshot[name]) for name in self._TRANSPORT_ARRAY_NAMES}
    for system in range(self.n_systems):
        result["mobile_positive"][system] = advanced[:n, system]
        result["retained_positive"][system] = advanced[n:, system]
        column = self.n_systems + system
        result["mobile_negative"][system] = advanced[:n, column]
        result["retained_negative"][system] = advanced[n:, column]

    final_mass = self._snapshot_mass(result)
    final_retained = float(np.sum(result["retained_positive"]) + np.sum(result["retained_negative"]))
    gain = final_mass - initial_mass
    gain_tolerance = 2.0e-8 * max(initial_mass, 1.0e-300) + 1.0e-300
    if gain > gain_tolerance:
        raise RuntimeError(
            "physical Krylov transport created line content: "
            f"gain={gain:.6e}, scale={initial_mass:.6e}, tolerance={gain_tolerance:.6e}"
        )
    escaped = max(initial_mass - final_mass, 0.0)
    signed_balance = initial_mass - final_mass - escaped
    diagnostics = {
        "dN_trapped": max(final_retained - initial_retained, 0.0),
        "dN_detrapped": max(initial_retained - final_retained, 0.0),
        "dN_escaped": escaped,
        "peierls_rate_min_s": float(np.min(peierls)),
        "peierls_rate_max_s": float(np.max(peierls)),
        "taylor_completion_rate_min_s": float(np.min(release)),
        "taylor_completion_rate_max_s": float(np.max(release)),
        "encounter_rate_min_s": float(np.min(encounter)),
        "encounter_rate_max_s": float(np.max(encounter)),
        "glide_velocity_max_m_s": float(np.max(velocity)),
        "rho_forest_min_m2": float(np.min(forest)),
        "rho_forest_max_m2": float(np.max(forest)),
        "max_frozen_courant": float(np.max(velocity) * dt / max(self.dx, 1.0e-30)),
        "line_content_conservation_error": abs(signed_balance),
        "line_content_conservation_signed_error": signed_balance,
        "line_content_conservation_scale": max(initial_mass, final_mass + escaped, 1.0e-300),
        "physical_generator_mass_gain": max(gain, 0.0),
        "physical_generator_action": "sparse_krylov_expm_multiply",
    }
    return result, diagnostics


def _transport_physical(self, *, dt_s: float, T_K: float, opening_stress_Pa: float) -> dict[str, Any]:
    dt_total = max(float(dt_s), 0.0)
    initial = self._transport_snapshot()
    initial_global_mass = self._snapshot_mass(initial)
    if dt_total <= 0.0 or initial_global_mass <= 0.0:
        if dt_total > 0.0:
            self.time_s += dt_total
        out = {
            "dN_trapped": 0.0,
            "dN_detrapped": 0.0,
            "dN_escaped": 0.0,
            "dN_recovered": 0.0,
            "transport_substeps": 0,
            "transport_attempted_physical_exponentials": 0,
            "transport_attempted_exponentials": 0,
            "transport_rejected_intervals": 0,
            "transport_nonlinear_error_max": 0.0,
            "transport_integrator": TRANSPORT_INTEGRATOR,
            "transport_cfl_limited": False,
            "explicit_recovery_active": False,
        }
        self.last_transport = copy.deepcopy(out)
        return out

    rtol = max(float(getattr(self, "transport_nonlinear_rtol", 1.0e-3)), 1.0e-10)
    limit = max(int(self.max_transport_substeps), 12)
    min_interval = max(
        float(getattr(self, "transport_min_interval_s", 1.0e-12)),
        np.finfo(float).eps * max(dt_total, 1.0),
    )
    attempted = 0
    rejected = 0
    accepted_rows: list[dict[str, float]] = []
    max_error = 0.0

    def integrate(snapshot: dict[str, np.ndarray], interval: float) -> dict[str, np.ndarray]:
        nonlocal attempted, rejected, max_error
        current_mass = self._snapshot_mass(snapshot)
        if current_mass <= rtol * initial_global_mass * 1.0e-10:
            return copy.deepcopy(snapshot)
        if attempted + 3 > limit:
            raise RuntimeError(
                "persistent-site physical Krylov transport exceeded nonlinear solve budget: "
                f"attempted={attempted}, limit={limit}, interval_s={interval:.6e}, "
                f"max_error={max_error:.6e}"
            )
        full, _ = self._frozen_transport_step(
            snapshot, dt_s=interval, T_K=T_K, opening_stress_Pa=opening_stress_Pa
        )
        half, d1 = self._frozen_transport_step(
            snapshot, dt_s=0.5 * interval, T_K=T_K, opening_stress_Pa=opening_stress_Pa
        )
        two_half, d2 = self._frozen_transport_step(
            half, dt_s=0.5 * interval, T_K=T_K, opening_stress_Pa=opening_stress_Pa
        )
        attempted += 3
        scale = max(initial_global_mass, current_mass, self._snapshot_mass(two_half), 1.0e-300)
        error = self._snapshot_difference(full, two_half) / scale
        max_error = max(max_error, error)
        if error <= rtol or interval <= min_interval:
            accepted_rows.extend((d1, d2))
            return two_half
        rejected += 1
        middle = integrate(snapshot, 0.5 * interval)
        return integrate(middle, 0.5 * interval)

    final = integrate(initial, dt_total)
    self._restore_transport_snapshot(final)
    accepted = self._combine_transport_diagnostics(accepted_rows)
    escaped = float(accepted.get("dN_escaped", 0.0))
    self.time_s += dt_total
    self.escaped_total += escaped
    out = {
        "dN_trapped": float(accepted.get("dN_trapped", 0.0)),
        "dN_detrapped": float(accepted.get("dN_detrapped", 0.0)),
        "dN_escaped": escaped,
        "dN_recovered": 0.0,
        "transport_substeps": len(accepted_rows),
        "transport_attempted_physical_exponentials": attempted,
        "transport_attempted_exponentials": attempted,
        "transport_attempted_linear_solves": 0,
        "transport_rejected_intervals": rejected,
        "transport_nonlinear_error_max": max_error,
        "transport_nonlinear_rtol": rtol,
        "transport_integrator": TRANSPORT_INTEGRATOR,
        "transport_cfl_limited": False,
        "explicit_recovery_active": False,
        **{
            key: value
            for key, value in accepted.items()
            if key not in {"dN_trapped", "dN_detrapped", "dN_escaped"}
        },
    }
    self.last_transport = copy.deepcopy(out)
    return out


@contextmanager
def installed_split_transport_v1005144() -> Iterator[None]:
    old_frozen = PersistentSiteSignedTransportMixin._frozen_transport_step
    old_transport = PersistentSiteSignedTransportMixin.transport
    PersistentSiteSignedTransportMixin._frozen_transport_step = _frozen_transport_step_physical
    PersistentSiteSignedTransportMixin.transport = _transport_physical
    try:
        yield
    finally:
        PersistentSiteSignedTransportMixin._frozen_transport_step = old_frozen
        PersistentSiteSignedTransportMixin.transport = old_transport


__all__ = [
    "TRANSPORT_INTEGRATOR",
    "_exact_exchange_pair",
    "_advect_mobile_exact",
    "_physical_generator",
    "installed_split_transport_v1005144",
]
