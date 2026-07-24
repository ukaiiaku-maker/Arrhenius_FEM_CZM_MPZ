"""Peierls--Taylor transport mixin for the persistent signed MPZ state.

v10.0.5.14.2 keeps the same finite-volume upwind transport equations used by
v10.0.5.14.1 but replaces explicit CFL microstepping with an adaptive,
L-stable backward-Euler solve of the coupled mobile/retained system.  The
change is numerical only: Peierls glide, encounter storage, Taylor release,
and absorbing escape at the outer MPZ boundary are unchanged.
"""
from __future__ import annotations

import copy
from typing import Any, Callable

import numpy as np
from scipy.sparse import csc_matrix, eye, lil_matrix
from scipy.sparse.linalg import spsolve

from .emission_derived_plasticity import (
    CorrelatedTaylorConfig,
    EmissionDerivedPeierlsTaylorConfig,
    ExpFloorSurface,
)
from .emission_derived_plasticity_v9102 import (
    EmissionDerivedPeierlsTaylorModel,
    IndependentShapeEntropyMechanismScale,
)


class PersistentSiteSignedTransportMixin:
    _TRANSPORT_ARRAY_NAMES = (
        "mobile_positive",
        "mobile_negative",
        "retained_positive",
        "retained_negative",
    )

    def _pt_model(self) -> EmissionDerivedPeierlsTaylorModel:
        c = self.candidate
        parent = ExpFloorSurface(
            G00_eV=c.emit_G00_eV,
            gT_eV_per_K=c.emit_gT_eV_per_K,
            sigc0_Pa=c.emit_sigc0_GPa * 1.0e9,
            sT_Pa_per_K=c.emit_sT_GPa_per_K * 1.0e9,
            Tref_K=c.Tref_K,
            a=c.emit_exp_a,
            n=c.emit_exp_n,
            floor_fraction=c.emit_floor_frac,
            floor_min_eV=1.0e-4,
            floor_max_fraction=0.95,
        )
        return EmissionDerivedPeierlsTaylorModel(
            EmissionDerivedPeierlsTaylorConfig(
                parent=parent,
                peierls=IndependentShapeEntropyMechanismScale(
                    energy_ratio=c.peierls_H0_eV / max(c.emit_G00_eV, 1.0e-30),
                    activation_entropy_kB=c.peierls_activation_entropy_kB,
                    exp_a=c.peierls_exp_a,
                    exp_n=c.peierls_exp_n,
                    stress_ratio=1.0,
                    rate_prefactor_s=c.peierls_nu0_s,
                ),
                taylor=IndependentShapeEntropyMechanismScale(
                    energy_ratio=c.taylor_H0_eV / max(c.emit_G00_eV, 1.0e-30),
                    activation_entropy_kB=c.taylor_activation_entropy_kB,
                    exp_a=c.taylor_exp_a,
                    exp_n=c.taylor_exp_n,
                    stress_ratio=1.0,
                    rate_prefactor_s=c.taylor_nu0_s,
                ),
                correlated_taylor=CorrelatedTaylorConfig(
                    rho_c_m2=c.taylor_corr_rho_c_m2,
                    renewal_time_s=1.0,
                    m_exponent=1.0,
                    m_scale=c.taylor_corr_scale,
                    m_cap=float("inf"),
                ),
                peierls_stress_fraction=c.peierls_stress_fraction,
                taylor_stress_fraction=c.taylor_stress_fraction,
                taylor_phi_max=float("inf"),
                mobile_fraction_low_density=0.01,
                mobile_saturation_density_m2=float("inf"),
                mobile_density_floor_m2=0.0,
                jump_fraction_of_forest_spacing=1.0,
                jump_length_min_m=self.b_m,
                rate_cap_s=float("inf"),
            )
        )

    def _transport_snapshot(self) -> dict[str, np.ndarray]:
        return {
            name: np.asarray(getattr(self, name), dtype=float).copy()
            for name in self._TRANSPORT_ARRAY_NAMES
        }

    def _restore_transport_snapshot(self, snapshot: dict[str, np.ndarray]) -> None:
        for name in self._TRANSPORT_ARRAY_NAMES:
            setattr(self, name, np.asarray(snapshot[name], dtype=float).copy())

    def _forest_from_snapshot(self, snapshot: dict[str, np.ndarray]) -> np.ndarray:
        width = max(self.blunting_length_m, self.dx)
        unsigned = np.sum(
            snapshot["mobile_positive"]
            + snapshot["mobile_negative"]
            + snapshot["retained_positive"]
            + snapshot["retained_negative"],
            axis=0,
        )
        return np.maximum(
            float(self.candidate.rho_forest_floor_m2)
            + unsigned / max(self.dx * width, 1.0e-30),
            1.0,
        )

    @staticmethod
    def _snapshot_mass(snapshot: dict[str, np.ndarray]) -> float:
        return float(sum(np.sum(snapshot[name]) for name in snapshot))

    @staticmethod
    def _snapshot_difference(
        first: dict[str, np.ndarray], second: dict[str, np.ndarray]
    ) -> float:
        return float(
            sum(
                np.sum(np.abs(first[name] - second[name]))
                for name in PersistentSiteSignedTransportMixin._TRANSPORT_ARRAY_NAMES
            )
        )

    def _frozen_transport_step(
        self,
        snapshot: dict[str, np.ndarray],
        *,
        dt_s: float,
        T_K: float,
        opening_stress_Pa: float,
    ) -> tuple[dict[str, np.ndarray], dict[str, float]]:
        """Advance one frozen-coefficient coupled transport interval implicitly."""
        dt = max(float(dt_s), 0.0)
        if dt <= 0.0:
            return copy.deepcopy(snapshot), {
                "dN_trapped": 0.0,
                "dN_detrapped": 0.0,
                "dN_escaped": 0.0,
                "max_frozen_courant": 0.0,
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
        n = self.n_bins
        face_velocity = np.empty(n + 1, dtype=float)
        face_velocity[1:-1] = 0.5 * (velocity[:-1] + velocity[1:])
        face_velocity[0] = velocity[0]
        face_velocity[-1] = velocity[-1]

        # y = [mobile bins, retained bins, cumulative trapped, cumulative
        #      detrapped, cumulative escaped].  The generator is Metzler and
        # conserves line content when the absorbing escape accumulator is included.
        dim = 2 * n + 3
        generator = lil_matrix((dim, dim), dtype=float)
        trap_row = 2 * n
        release_row = 2 * n + 1
        escape_row = 2 * n + 2
        for i in range(n):
            outflow = face_velocity[i + 1] / max(self.dx, 1.0e-30)
            generator[i, i] -= outflow + encounter[i]
            if i > 0:
                generator[i, i - 1] += face_velocity[i] / max(self.dx, 1.0e-30)
            generator[i, n + i] += release_rate[i]
            generator[n + i, i] += encounter[i]
            generator[n + i, n + i] -= release_rate[i]
            generator[trap_row, i] += encounter[i]
            generator[release_row, n + i] += release_rate[i]
        generator[escape_row, n - 1] += face_velocity[-1] / max(self.dx, 1.0e-30)

        n_columns = 2 * self.n_systems
        initial = np.zeros((dim, n_columns), dtype=float)
        for system in range(self.n_systems):
            initial[:n, system] = snapshot["mobile_positive"][system]
            initial[n : 2 * n, system] = snapshot["retained_positive"][system]
            negative_column = self.n_systems + system
            initial[:n, negative_column] = snapshot["mobile_negative"][system]
            initial[n : 2 * n, negative_column] = snapshot["retained_negative"][system]

        system_matrix = eye(dim, format="csc") - dt * csc_matrix(generator)
        advanced = np.asarray(spsolve(system_matrix, initial), dtype=float)
        if advanced.ndim == 1:
            advanced = advanced[:, None]
        if not np.all(np.isfinite(advanced)):
            raise RuntimeError("implicit persistent-site transport produced nonfinite state")
        magnitude = max(float(np.max(np.abs(advanced))), 1.0)
        negative_tolerance = 1.0e-11 * magnitude
        if float(np.min(advanced)) < -negative_tolerance:
            raise RuntimeError(
                "implicit persistent-site transport violated nonnegative state: "
                f"minimum={float(np.min(advanced)):.6e}"
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
        initial_mass = self._snapshot_mass(snapshot)
        final_mass = self._snapshot_mass(result)
        conservation_error = abs(initial_mass - final_mass - escaped)
        conservation_scale = max(initial_mass, final_mass + escaped, 1.0)
        if conservation_error > 1.0e-8 * conservation_scale:
            raise RuntimeError(
                "implicit persistent-site transport failed line-content conservation: "
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
            "max_frozen_courant": float(
                np.max(velocity) * dt / max(self.dx, 1.0e-30)
            ),
            "line_content_conservation_error": conservation_error,
        }
        return result, diagnostics

    @staticmethod
    def _combine_transport_diagnostics(
        diagnostics: list[dict[str, float]],
    ) -> dict[str, float]:
        if not diagnostics:
            return {}
        summed = {
            key: float(sum(item.get(key, 0.0) for item in diagnostics))
            for key in (
                "dN_trapped",
                "dN_detrapped",
                "dN_escaped",
            )
        }
        minima = {
            key: float(min(item[key] for item in diagnostics if key in item))
            for key in (
                "peierls_rate_min_s",
                "taylor_completion_rate_min_s",
                "encounter_rate_min_s",
                "rho_forest_min_m2",
            )
        }
        maxima = {
            key: float(max(item[key] for item in diagnostics if key in item))
            for key in (
                "peierls_rate_max_s",
                "taylor_completion_rate_max_s",
                "encounter_rate_max_s",
                "glide_velocity_max_m_s",
                "rho_forest_max_m2",
                "max_frozen_courant",
                "line_content_conservation_error",
            )
        }
        return {**summed, **minima, **maxima}

    def transport(
        self, *, dt_s: float, T_K: float, opening_stress_Pa: float
    ) -> dict[str, Any]:
        dt_total = max(float(dt_s), 0.0)
        initial = self._transport_snapshot()
        if dt_total <= 0.0 or self._snapshot_mass(initial) <= 0.0:
            out = {
                "dN_trapped": 0.0,
                "dN_detrapped": 0.0,
                "dN_escaped": 0.0,
                "dN_recovered": 0.0,
                "transport_substeps": 0,
                "transport_attempted_linear_solves": 0,
                "transport_rejected_intervals": 0,
                "transport_integrator": "adaptive_backward_euler_upwind_v10_0_5_14_2",
                "transport_cfl_limited": False,
                "explicit_recovery_active": False,
            }
            self.last_transport = copy.deepcopy(out)
            return out

        nonlinear_rtol = max(
            float(getattr(self, "transport_nonlinear_rtol", 1.0e-3)), 1.0e-10
        )
        max_solves = max(int(self.max_transport_substeps), 12)
        minimum_interval = max(
            float(getattr(self, "transport_min_interval_s", 1.0e-12)),
            np.finfo(float).eps * max(dt_total, 1.0),
        )
        attempted_solves = 0
        rejected_intervals = 0
        accepted_diagnostics: list[dict[str, float]] = []
        maximum_error = 0.0

        def integrate_interval(
            snapshot: dict[str, np.ndarray], interval: float
        ) -> dict[str, np.ndarray]:
            nonlocal attempted_solves, rejected_intervals, maximum_error
            if attempted_solves + 3 > max_solves:
                raise RuntimeError(
                    "persistent-site implicit transport exceeded nonlinear solve budget: "
                    f"attempted={attempted_solves}, limit={max_solves}, "
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
            attempted_solves += 3
            scale = max(
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
            "transport_attempted_linear_solves": attempted_solves,
            "transport_rejected_intervals": rejected_intervals,
            "transport_nonlinear_error_max": maximum_error,
            "transport_nonlinear_rtol": nonlinear_rtol,
            "transport_integrator": "adaptive_backward_euler_upwind_v10_0_5_14_2",
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

    def evolve(
        self,
        *,
        dt_s: float,
        T_K: float,
        opening_stress_Pa: float,
        drive_factors: np.ndarray,
        tau_signed_Pa: np.ndarray,
        emission_rate_function: Callable[[float, float], float],
    ) -> dict[str, Any]:
        emission = self.emit_persistent(
            dt_s=dt_s,
            T_K=T_K,
            opening_stress_Pa=opening_stress_Pa,
            drive_factors=drive_factors,
            tau_signed_Pa=tau_signed_Pa,
            rate_function=emission_rate_function,
        )
        transport = self.transport(
            dt_s=dt_s,
            T_K=T_K,
            opening_stress_Pa=opening_stress_Pa,
        )
        return {**emission, **transport}


__all__ = ["PersistentSiteSignedTransportMixin"]
