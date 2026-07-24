"""Stiff exponential transport correction for persistent signed MPZ v10.0.5.14.3.

The constitutive Peierls, encounter, Taylor-release, and absorbing-boundary
rates are unchanged.  For each nonlinear interval, the frozen finite-volume
mobile/retained generator is advanced exactly with a dense scaling-and-squaring
matrix exponential.  Step doubling therefore measures only the change of the
state-dependent coefficients, rather than the truncation error of backward
Euler in the infinitely stiff escape limit.
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

TRANSPORT_INTEGRATOR = "adaptive_frozen_generator_exponential_v10_0_5_14_3"


def _frozen_transport_step_exponential(
    self,
    snapshot: dict[str, np.ndarray],
    *,
    dt_s: float,
    T_K: float,
    opening_stress_Pa: float,
) -> tuple[dict[str, np.ndarray], dict[str, float]]:
    """Advance one frozen-coefficient finite-volume interval exactly."""
    dt = max(float(dt_s), 0.0)
    initial_mass = self._snapshot_mass(snapshot)
    if dt <= 0.0 or initial_mass <= 0.0:
        return copy.deepcopy(snapshot), {
            "dN_trapped": 0.0,
            "dN_detrapped": 0.0,
            "dN_escaped": 0.0,
            "max_frozen_courant": 0.0,
            "line_content_conservation_error": 0.0,
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
    inv_dx = 1.0 / max(self.dx, 1.0e-30)
    face_velocity = np.empty(n + 1, dtype=float)
    face_velocity[1:-1] = 0.5 * (velocity[:-1] + velocity[1:])
    face_velocity[0] = velocity[0]
    face_velocity[-1] = velocity[-1]

    # y = [mobile bins, retained bins, cumulative trapped, cumulative
    #      detrapped, cumulative escaped].  The generator is Metzler.  Including
    # the escape accumulator makes total line content exactly conservative.
    dim = 2 * n + 3
    generator = np.zeros((dim, dim), dtype=float)
    trap_row = 2 * n
    release_row = 2 * n + 1
    escape_row = 2 * n + 2
    for i in range(n):
        outflow = face_velocity[i + 1] * inv_dx
        generator[i, i] -= outflow + encounter[i]
        if i > 0:
            generator[i, i - 1] += face_velocity[i] * inv_dx
        generator[i, n + i] += release_rate[i]
        generator[n + i, i] += encounter[i]
        generator[n + i, n + i] -= release_rate[i]
        generator[trap_row, i] += encounter[i]
        generator[release_row, n + i] += release_rate[i]
    generator[escape_row, n - 1] += face_velocity[-1] * inv_dx

    n_columns = 2 * self.n_systems
    initial = np.zeros((dim, n_columns), dtype=float)
    for system in range(self.n_systems):
        initial[:n, system] = snapshot["mobile_positive"][system]
        initial[n : 2 * n, system] = snapshot["retained_positive"][system]
        negative_column = self.n_systems + system
        initial[:n, negative_column] = snapshot["mobile_negative"][system]
        initial[n : 2 * n, negative_column] = snapshot["retained_negative"][system]

    advanced = np.asarray(expm(dt * generator) @ initial, dtype=float)
    if not np.all(np.isfinite(advanced)):
        raise RuntimeError(
            "exponential persistent-site transport produced nonfinite state"
        )
    magnitude = max(float(np.max(np.abs(advanced))), 1.0)
    negative_tolerance = 5.0e-11 * magnitude
    minimum = float(np.min(advanced))
    if minimum < -negative_tolerance:
        raise RuntimeError(
            "exponential persistent-site transport violated nonnegative state: "
            f"minimum={minimum:.6e}"
        )
    advanced = np.maximum(advanced, 0.0)

    result = {
        name: np.zeros_like(snapshot[name]) for name in self._TRANSPORT_ARRAY_NAMES
    }
    for system in range(self.n_systems):
        result["mobile_positive"][system] = advanced[:n, system]
        result["retained_positive"][system] = advanced[n : 2 * n, system]
        negative_column = self.n_systems + system
        result["mobile_negative"][system] = advanced[:n, negative_column]
        result["retained_negative"][system] = advanced[
            n : 2 * n, negative_column
        ]

    trapped = float(np.sum(advanced[trap_row]))
    released = float(np.sum(advanced[release_row]))
    escaped = float(np.sum(advanced[escape_row]))
    final_mass = self._snapshot_mass(result)
    conservation_error = abs(initial_mass - final_mass - escaped)
    conservation_scale = max(initial_mass, final_mass + escaped, 1.0)
    if conservation_error > 2.0e-8 * conservation_scale:
        raise RuntimeError(
            "exponential persistent-site transport failed line-content conservation: "
            f"error={conservation_error:.6e}, scale={conservation_scale:.6e}"
        )

    diagnostics = {
        "dN_trapped": trapped,
        "dN_detrapped": released,
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
        "max_frozen_courant": float(np.max(velocity) * dt * inv_dx),
        "line_content_conservation_error": conservation_error,
    }
    return result, diagnostics


def _transport_exponential(
    self, *, dt_s: float, T_K: float, opening_stress_Pa: float
) -> dict[str, Any]:
    """Adaptive nonlinear transport with exact frozen-generator intervals."""
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
            "transport_attempted_exponentials": 0,
            "transport_attempted_linear_solves": 0,
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
            "transport_attempted_exponentials": 0,
            "transport_attempted_linear_solves": 0,
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
        if self._snapshot_mass(snapshot) <= nonlinear_rtol * initial_global_mass * 1e-10:
            return copy.deepcopy(snapshot)
        if attempted + 3 > max_exponentials:
            raise RuntimeError(
                "persistent-site exponential transport exceeded nonlinear solve budget: "
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
            self._snapshot_mass(snapshot),
            self._snapshot_mass(two_half),
            1.0e-30,
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
def installed_exponential_transport_v1005143() -> Iterator[None]:
    """Temporarily install the v10.0.5.14.3 transport on the shared state class."""
    old_frozen = PersistentSiteSignedTransportMixin._frozen_transport_step
    old_transport = PersistentSiteSignedTransportMixin.transport
    PersistentSiteSignedTransportMixin._frozen_transport_step = (
        _frozen_transport_step_exponential
    )
    PersistentSiteSignedTransportMixin.transport = _transport_exponential
    try:
        yield
    finally:
        PersistentSiteSignedTransportMixin._frozen_transport_step = old_frozen
        PersistentSiteSignedTransportMixin.transport = old_transport


__all__ = [
    "TRANSPORT_INTEGRATOR",
    "installed_exponential_transport_v1005143",
]
