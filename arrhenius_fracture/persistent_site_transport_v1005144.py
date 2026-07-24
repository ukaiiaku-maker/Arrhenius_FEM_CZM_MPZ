"""Conservative stiff MPZ transport for v10.0.5.14.4.

The PF-derived Peierls, encounter, Taylor-release, and absorbing-boundary laws
are unchanged.  v10.0.5.14.3 exponentiated the physical mobile/retained state
together with cumulative trapping, release, and escape counters.  At high
temperature those accumulator rows made the augmented operator badly scaled and
produced false conservation failures for small emitted populations.

v10.0.5.14.4 exponentiates only the exact frozen 2N physical generator
[mobile bins, retained bins].  The generator contains upwind transport,
encounter storage, Taylor release, and the absorbing outer-boundary sink.
Escape is evaluated as physical mass removed from that substochastic system.
Step doubling therefore measures only changes in the state-dependent rates.
"""
from __future__ import annotations

import copy
from contextlib import contextmanager
from typing import Any, Iterator

import numpy as np
from scipy.linalg import expm

from .persistent_site_signed_transport_v100514 import (
    PersistentSiteSignedTransportMixin,
)

TRANSPORT_INTEGRATOR = (
    "adaptive_physical_generator_exponential_v10_0_5_14_4"
)


def _exact_exchange_pair(
    mobile: np.ndarray,
    retained: np.ndarray,
    encounter_rate_s: np.ndarray,
    release_rate_s: np.ndarray,
    dt_s: float,
) -> tuple[np.ndarray, np.ndarray, float, float]:
    """Exact nonnegative local solution of M <-> R, retained for audit tests."""
    m0 = np.maximum(np.asarray(mobile, dtype=float), 0.0)
    r0 = np.maximum(np.asarray(retained, dtype=float), 0.0)
    ke = np.maximum(np.asarray(encounter_rate_s, dtype=float), 0.0)[None, :]
    kr = np.maximum(np.asarray(release_rate_s, dtype=float), 0.0)[None, :]
    dt = max(float(dt_s), 0.0)
    total = m0 + r0
    rate = ke + kr
    retained_equilibrium = np.divide(
        ke * total,
        rate,
        out=np.zeros_like(total),
        where=rate > 0.0,
    )
    decay = np.exp(-np.minimum(rate * dt, 700.0))
    r1 = np.where(
        rate > 0.0,
        retained_equilibrium + (r0 - retained_equilibrium) * decay,
        r0,
    )
    r1 = np.clip(r1, 0.0, total)
    m1 = total - r1
    trapped = float(np.sum(np.maximum(r1 - r0, 0.0)))
    released = float(np.sum(np.maximum(r0 - r1, 0.0)))
    return m1, r1, trapped, released


def _advection_generator(velocity_m_s: np.ndarray, dx_m: float) -> np.ndarray:
    """Return the positive upwind generator with an absorbing outer boundary."""
    velocity = np.maximum(np.asarray(velocity_m_s, dtype=float).reshape(-1), 0.0)
    n = int(velocity.size)
    inv_dx = 1.0 / max(float(dx_m), 1.0e-30)
    face_velocity = np.empty(n + 1, dtype=float)
    face_velocity[1:-1] = 0.5 * (velocity[:-1] + velocity[1:])
    face_velocity[0] = velocity[0]
    face_velocity[-1] = velocity[-1]
    generator = np.zeros((n, n), dtype=float)
    for i in range(n):
        generator[i, i] -= face_velocity[i + 1] * inv_dx
        if i > 0:
            generator[i, i - 1] += face_velocity[i] * inv_dx
    return generator


def _advect_mobile_exact(
    mobile: np.ndarray,
    velocity_m_s: np.ndarray,
    dx_m: float,
    dt_s: float,
) -> tuple[np.ndarray, float, float]:
    """Advance all mobile channels exactly and return absorbed line content."""
    source = np.maximum(np.asarray(mobile, dtype=float), 0.0)
    dt = max(float(dt_s), 0.0)
    if dt <= 0.0 or float(np.sum(source)) <= 0.0:
        return source.copy(), 0.0, 0.0
    propagator = expm(dt * _advection_generator(velocity_m_s, dx_m))
    advanced = np.asarray(source @ propagator.T, dtype=float)
    if not np.all(np.isfinite(advanced)):
        raise RuntimeError("physical exponential advection produced nonfinite state")
    magnitude = max(float(np.max(np.abs(advanced))), float(np.max(source)), 1.0e-300)
    negative_tolerance = 2.0e-12 * magnitude + 1.0e-300
    minimum = float(np.min(advanced))
    if minimum < -negative_tolerance:
        raise RuntimeError(
            "physical exponential advection violated nonnegative state: "
            f"minimum={minimum:.6e}, tolerance={negative_tolerance:.6e}"
        )
    advanced = np.maximum(advanced, 0.0)
    initial_mass = float(np.sum(source))
    final_mass = float(np.sum(advanced))
    mass_gain = final_mass - initial_mass
    gain_tolerance = 2.0e-11 * max(initial_mass, 1.0e-300) + 1.0e-300
    if mass_gain > gain_tolerance:
        raise RuntimeError(
            "physical exponential advection created line content: "
            f"gain={mass_gain:.6e}, tolerance={gain_tolerance:.6e}"
        )
    return advanced, max(initial_mass - final_mass, 0.0), max(mass_gain, 0.0)


def _physical_generator(
    velocity_m_s: np.ndarray,
    encounter_rate_s: np.ndarray,
    release_rate_s: np.ndarray,
    dx_m: float,
) -> np.ndarray:
    """Coupled mobile/retained generator with absorbing mobile escape."""
    velocity = np.maximum(np.asarray(velocity_m_s, dtype=float).reshape(-1), 0.0)
    encounter = np.maximum(
        np.asarray(encounter_rate_s, dtype=float).reshape(-1), 0.0
    )
    release = np.maximum(
        np.asarray(release_rate_s, dtype=float).reshape(-1), 0.0
    )
    if not (velocity.shape == encounter.shape == release.shape):
        raise ValueError("physical transport coefficient shapes do not match")
    n = int(velocity.size)
    inv_dx = 1.0 / max(float(dx_m), 1.0e-30)
    face_velocity = np.empty(n + 1, dtype=float)
    face_velocity[1:-1] = 0.5 * (velocity[:-1] + velocity[1:])
    face_velocity[0] = velocity[0]
    face_velocity[-1] = velocity[-1]
    generator = np.zeros((2 * n, 2 * n), dtype=float)
    for i in range(n):
        generator[i, i] -= face_velocity[i + 1] * inv_dx + encounter[i]
        if i > 0:
            generator[i, i - 1] += face_velocity[i] * inv_dx
        generator[i, n + i] += release[i]
        generator[n + i, i] += encounter[i]
        generator[n + i, n + i] -= release[i]
    return generator


def _frozen_transport_step_physical(
    self,
    snapshot: dict[str, np.ndarray],
    *,
    dt_s: float,
    T_K: float,
    opening_stress_Pa: float,
) -> tuple[dict[str, np.ndarray], dict[str, float]]:
    """Advance one frozen physical mobile/retained interval exactly."""
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
    peierls = np.maximum(
        np.asarray(rates["peierls_rate_s"], dtype=float).reshape(-1), 0.0
    )
    release_rate = np.maximum(
        np.asarray(rates["taylor_completion_rate_s"], dtype=float).reshape(-1),
        0.0,
    )
    jump = np.maximum(
        np.asarray(rates["jump_length_m"], dtype=float).reshape(-1), self.b_m
    )
    if not (
        peierls.shape == release_rate.shape == jump.shape == (self.n_bins,)
    ):
        raise RuntimeError("Peierls--Taylor transport rates do not match MPZ bins")
    if not (
        np.all(np.isfinite(peierls))
        and np.all(np.isfinite(release_rate))
        and np.all(np.isfinite(jump))
    ):
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
    n_columns = 2 * self.n_systems
    initial = np.zeros((2 * n, n_columns), dtype=float)
    initial_retained = 0.0
    for system in range(self.n_systems):
        initial[:n, system] = snapshot["mobile_positive"][system]
        initial[n:, system] = snapshot["retained_positive"][system]
        negative_column = self.n_systems + system
        initial[:n, negative_column] = snapshot["mobile_negative"][system]
        initial[n:, negative_column] = snapshot["retained_negative"][system]
        initial_retained += float(np.sum(snapshot["retained_positive"][system]))
        initial_retained += float(np.sum(snapshot["retained_negative"][system]))

    generator = _physical_generator(velocity, encounter, release_rate, self.dx)
    advanced = np.asarray(expm(dt * generator) @ initial, dtype=float)
    if not np.all(np.isfinite(advanced)):
        raise RuntimeError("physical transport exponential produced nonfinite state")
    magnitude = max(float(np.max(np.abs(advanced))), float(np.max(initial)), 1.0e-300)
    negative_tolerance = 2.0e-11 * magnitude + 1.0e-300
    minimum = float(np.min(advanced))
    if minimum < -negative_tolerance:
        raise RuntimeError(
            "physical transport exponential violated nonnegative state: "
            f"minimum={minimum:.6e}, tolerance={negative_tolerance:.6e}"
        )
    advanced = np.maximum(advanced, 0.0)

    result = {
        name: np.zeros_like(snapshot[name]) for name in self._TRANSPORT_ARRAY_NAMES
    }
    for system in range(self.n_systems):
        result["mobile_positive"][system] = advanced[:n, system]
        result["retained_positive"][system] = advanced[n:, system]
        negative_column = self.n_systems + system
        result["mobile_negative"][system] = advanced[:n, negative_column]
        result["retained_negative"][system] = advanced[n:, negative_column]

    final_mass = self._snapshot_mass(result)
    final_retained = float(
        np.sum(result["retained_positive"]) + np.sum(result["retained_negative"])
    )
    mass_gain = final_mass - initial_mass
    gain_tolerance = 2.0e-9 * max(initial_mass, 1.0e-300) + 1.0e-300
    if mass_gain > gain_tolerance:
        raise RuntimeError(
            "physical transport exponential created line content: "
            f"gain={mass_gain:.6e}, scale={initial_mass:.6e}, "
            f"tolerance={gain_tolerance:.6e}"
        )
    escaped = max(initial_mass - final_mass, 0.0)
    signed_balance = initial_mass - final_mass - escaped
    conservation_error = abs(signed_balance)
    conservation_scale = max(initial_mass, final_mass + escaped, 1.0e-300)

    diagnostics = {
        "dN_trapped": max(final_retained - initial_retained, 0.0),
        "dN_detrapped": max(initial_retained - final_retained, 0.0),
        "dN_escaped": escaped,
        "peierls_rate_min_s": float(np.min(peierls)),
        "peierls_rate_max_s": float(np.max(peierls)),
        "taylor_completion_rate_min_s": float(np.min(release_rate)),
        "taylor_completion_rate_max_s": float(np.max(release_rate)),
        "encounter_rate_min_s": float(np.min(encounter)),
        "encounter_rate_max_s": float(np.max(encounter)),
        "glide_velocity_max_m_s": float(np.max(velocity)),
        "rho_forest_min_m2": float(np.min(forest)),
        "rho_forest_max_m2": float(np.max(forest)),
        "max_frozen_courant": float(
            np.max(velocity) * dt / max(self.dx, 1.0e-30)
        ),
        "line_content_conservation_error": conservation_error,
        "line_content_conservation_signed_error": signed_balance,
        "line_content_conservation_scale": conservation_scale,
        "physical_generator_mass_gain": max(mass_gain, 0.0),
    }
    return result, diagnostics


def _transport_physical(
    self, *, dt_s: float, T_K: float, opening_stress_Pa: float
) -> dict[str, Any]:
    """Adaptive nonlinear transport using exact frozen physical generators."""
    dt_total = max(float(dt_s), 0.0)
    initial = self._transport_snapshot()
    initial_global_mass = self._snapshot_mass(initial)
    if dt_total <= 0.0:
        out = {
            "dN_trapped": 0.0,
            "dN_detrapped": 0.0,
            "dN_escaped": 0.0,
            "dN_recovered": 0.0,
            "transport_substeps": 0,
            "transport_attempted_physical_exponentials": 0,
            "transport_attempted_exponentials": 0,
            "transport_rejected_intervals": 0,
            "transport_integrator": TRANSPORT_INTEGRATOR,
            "transport_cfl_limited": False,
            "explicit_recovery_active": False,
        }
        self.last_transport = copy.deepcopy(out)
        return out
    if initial_global_mass <= 0.0:
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

    nonlinear_rtol = max(
        float(getattr(self, "transport_nonlinear_rtol", 1.0e-3)), 1.0e-10
    )
    max_exponentials = max(int(self.max_transport_substeps), 12)
    minimum_interval = max(
        float(getattr(self, "transport_min_interval_s", 1.0e-12)),
        np.finfo(float).eps * max(dt_total, 1.0),
    )
    attempted = 0
    rejected_intervals = 0
    accepted_diagnostics: list[dict[str, float]] = []
    maximum_error = 0.0

    def integrate_interval(
        snapshot: dict[str, np.ndarray], interval: float
    ) -> dict[str, np.ndarray]:
        nonlocal attempted, rejected_intervals, maximum_error
        current_mass = self._snapshot_mass(snapshot)
        if current_mass <= nonlinear_rtol * initial_global_mass * 1.0e-10:
            return copy.deepcopy(snapshot)
        if attempted + 3 > max_exponentials:
            raise RuntimeError(
                "persistent-site physical transport exceeded nonlinear solve budget: "
                f"attempted={attempted}, limit={max_exponentials}, "
                f"interval_s={interval:.6e}, max_error={maximum_error:.6e}"
            )
        full, _ = self._frozen_transport_step(
            snapshot,
            dt_s=interval,
            T_K=T_K,
            opening_stress_Pa=opening_stress_Pa,
        )
        half, first_diag = self._frozen_transport_step(
            snapshot,
            dt_s=0.5 * interval,
            T_K=T_K,
            opening_stress_Pa=opening_stress_Pa,
        )
        two_half, second_diag = self._frozen_transport_step(
            half,
            dt_s=0.5 * interval,
            T_K=T_K,
            opening_stress_Pa=opening_stress_Pa,
        )
        attempted += 3
        scale = max(
            initial_global_mass,
            current_mass,
            self._snapshot_mass(two_half),
            1.0e-300,
        )
        error = self._snapshot_difference(full, two_half) / scale
        maximum_error = max(maximum_error, error)
        if error <= nonlinear_rtol or interval <= minimum_interval:
            accepted_diagnostics.extend((first_diag, second_diag))
            return two_half
        rejected_intervals += 1
        midpoint = integrate_interval(snapshot, 0.5 * interval)
        return integrate_interval(midpoint, 0.5 * interval)

    final = integrate_interval(initial, dt_total)
    self._restore_transport_snapshot(final)
    accepted = self._combine_transport_diagnostics(accepted_diagnostics)
    escaped = float(accepted.get("dN_escaped", 0.0))
    self.time_s += dt_total
    self.escaped_total += escaped
    out = {
        "dN_trapped": float(accepted.get("dN_trapped", 0.0)),
        "dN_detrapped": float(accepted.get("dN_detrapped", 0.0)),
        "dN_escaped": escaped,
        "dN_recovered": 0.0,
        "transport_substeps": len(accepted_diagnostics),
        "transport_attempted_physical_exponentials": attempted,
        "transport_attempted_exponentials": attempted,
        "transport_attempted_linear_solves": 0,
        "transport_rejected_intervals": rejected_intervals,
        "transport_nonlinear_error_max": maximum_error,
        "transport_nonlinear_rtol": nonlinear_rtol,
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
    """Temporarily install the exact physical-state v10.0.5.14.4 transport."""
    old_frozen = PersistentSiteSignedTransportMixin._frozen_transport_step
    old_transport = PersistentSiteSignedTransportMixin.transport
    PersistentSiteSignedTransportMixin._frozen_transport_step = (
        _frozen_transport_step_physical
    )
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
