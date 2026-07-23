"""Peierls--Taylor transport mixin for persistent signed MPZ state."""
from __future__ import annotations

import copy
from typing import Any, Callable

import numpy as np

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

    @staticmethod
    def _exchange(
        mobile: np.ndarray,
        retained: np.ndarray,
        encounter: np.ndarray,
        release: np.ndarray,
        dt: float,
    ) -> tuple[np.ndarray, np.ndarray, float, float]:
        ke = np.maximum(np.asarray(encounter, dtype=float), 0.0)[None, :]
        kt = np.maximum(np.asarray(release, dtype=float), 0.0)[None, :]
        total = np.maximum(mobile, 0.0) + np.maximum(retained, 0.0)
        rate = ke + kt
        frac_r = np.divide(ke, rate, out=np.zeros_like(rate), where=rate > 0.0)
        equilibrium = frac_r * total
        decay = np.exp(-np.minimum(rate * max(float(dt), 0.0), 700.0))
        new_r = np.clip(
            equilibrium + (retained - equilibrium) * decay, 0.0, total
        )
        new_m = total - new_r
        trapped = float(np.sum(np.maximum(new_r - retained, 0.0)))
        released = float(np.sum(np.maximum(retained - new_r, 0.0)))
        return new_m, new_r, trapped, released

    @staticmethod
    def _advect_forward(
        field: np.ndarray,
        velocity_m_s: np.ndarray,
        dx: float,
        dt: float,
    ) -> tuple[np.ndarray, float]:
        source = np.asarray(field, dtype=float)
        velocity = np.asarray(velocity_m_s, dtype=float).reshape(-1)
        out = source.copy()
        escaped = 0.0
        for system in range(source.shape[0]):
            f = source[system]
            v = np.maximum(velocity, 0.0)
            face_v = np.empty(f.size + 1)
            face_v[1:-1] = 0.5 * (v[:-1] + v[1:])
            face_v[0], face_v[-1] = v[0], v[-1]
            flux = np.zeros(f.size + 1)
            for j in range(1, f.size):
                flux[j] = face_v[j] * f[j - 1]
            flux[-1] = face_v[-1] * f[-1]
            out[system] = np.maximum(
                f - float(dt) * (flux[1:] - flux[:-1]) / float(dx), 0.0
            )
            escaped += max(float(dt) * flux[-1] / float(dx), 0.0)
        return out, escaped

    def transport(
        self, *, dt_s: float, T_K: float, opening_stress_Pa: float
    ) -> dict[str, Any]:
        dt_total = max(float(dt_s), 0.0)
        remaining = dt_total
        trapped = released = escaped = 0.0
        substeps = 0
        last: dict[str, Any] = {}
        while remaining > 0.0:
            substeps += 1
            if substeps > self.max_transport_substeps:
                raise RuntimeError("persistent-site transport exceeded max substeps")
            forest = self.forest_density_profile_m2()
            radius = self.blunted_radius()
            stress = max(float(opening_stress_Pa), 0.0) * np.sqrt(
                radius / np.maximum(radius + self.x, radius)
            )
            rates = self._pt_model().rates(stress, forest, T_K, self.b_m)
            peierls = np.maximum(
                np.asarray(rates["peierls_rate_s"], dtype=float).reshape(-1),
                0.0,
            )
            release_rate = np.maximum(
                np.asarray(
                    rates["taylor_completion_rate_s"], dtype=float
                ).reshape(-1),
                0.0,
            )
            jump = np.maximum(
                np.asarray(rates["jump_length_m"], dtype=float).reshape(-1),
                self.b_m,
            )
            velocity = jump * peierls
            encounter = (
                float(self.candidate.encounter_efficiency)
                * velocity
                * np.sqrt(np.maximum(forest, 0.0))
            )
            max_rate = max(
                float(np.max(encounter)),
                float(np.max(release_rate)),
                float(np.max(velocity)) / max(self.dx, 1.0e-30),
            )
            dt = (
                remaining
                if max_rate <= 0.0
                else min(remaining, self.max_transport_cfl / max_rate)
            )
            for m_name, r_name in (
                ("mobile_positive", "retained_positive"),
                ("mobile_negative", "retained_negative"),
            ):
                m, r, t, rel = self._exchange(
                    getattr(self, m_name),
                    getattr(self, r_name),
                    encounter,
                    release_rate,
                    dt,
                )
                setattr(self, m_name, m)
                setattr(self, r_name, r)
                trapped += t
                released += rel
                m, esc = self._advect_forward(
                    getattr(self, m_name), velocity, self.dx, dt
                )
                setattr(self, m_name, m)
                escaped += esc
            remaining -= dt
            last = {
                "peierls_rate_min_s": float(np.min(peierls)),
                "peierls_rate_max_s": float(np.max(peierls)),
                "taylor_completion_rate_min_s": float(np.min(release_rate)),
                "taylor_completion_rate_max_s": float(np.max(release_rate)),
                "encounter_rate_min_s": float(np.min(encounter)),
                "encounter_rate_max_s": float(np.max(encounter)),
                "glide_velocity_max_m_s": float(np.max(velocity)),
                "rho_forest_min_m2": float(np.min(forest)),
                "rho_forest_max_m2": float(np.max(forest)),
            }
        self.time_s += dt_total
        self.escaped_total += escaped
        out = {
            "dN_trapped": trapped,
            "dN_detrapped": released,
            "dN_escaped": escaped,
            "dN_recovered": 0.0,
            "transport_substeps": substeps,
            "explicit_recovery_active": False,
            **last,
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
