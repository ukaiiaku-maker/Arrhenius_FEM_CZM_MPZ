"""State-coupled exact stochastic emission for FEM/CZM v10.0.5.18.2.

This release keeps the v10.0.5.18 event law: one independent exponential
integrated-hazard clock per signed channel and one mechanically normalized line
packet per completed activation.  It adds a non-mutating constitutive predictor
for the outer adaptive load controller.  The predictor rejects and subdivides a
load increment when the endpoint approximation would change the emission rate,
Taylor backstress, tip radius, physical front width, source area, or persistent
site multiplicity too strongly.

No emission count is capped and no post-onset mean-field burst is introduced.
The accepted update remains the exact event-localized v10.0.5.18 process, which
recomputes geometry and backstress after every activation.
"""
from __future__ import annotations

import copy
import math
import os
from typing import Any

import numpy as np

from .persistent_site_stochastic_emission_v100518 import (
    EMISSION_HAZARD_SCHEMA,
    EMISSION_TRANSPORT_SCHEMA,
    PersistentSiteStochasticEmissionMovingTipFrontEngineV100518,
    _emission_hazards,
)

MODEL_ID = "FEM_CZM_state_coupled_exact_stochastic_emission_v10_0_5_18_2"
ADAPTIVE_SCHEMA = (
    "v10.0.5.18.2_outer_load_subdivision_from_emission_state_change"
)


def two_channel_absolute_opening_drives(self, fallback_K: float):
    """Return cleavage and emission K drives without double anisotropic scaling.

    The inherited mixed-mode wrapper supplies a scalar emission multiplier.  The
    persistent-site model separately resolves two signed channel factors from FEM
    tensor probes.  Applying both factors multiplies the same anisotropic emission
    drive twice.  For the persistent two-channel engine the opening stress is
    therefore formed from the absolute J-derived K; the signed channel factors are
    applied exactly once inside the source hazard.
    """

    KJ = max(float(fallback_K), 0.0)
    latest = dict(getattr(getattr(self, "_mm", None), "latest", {}) or {})
    fc = float(latest.get("cleavage_factor", 1.0))
    fe = float(latest.get("emission_factor", 1.0))
    resolved = bool(latest.get("two_channel_drive_reliable", False))
    use_absolute = bool(
        getattr(self, "persistent_site_source_active", False) and resolved
    )
    K_cleave = KJ * fc
    K_emit = KJ if use_absolute else KJ * fe
    metadata = {
        "KJ": KJ,
        "fc": fc,
        "fe": fe,
        **latest,
        "persistent_two_channel_absolute_opening_drive": use_absolute,
        "scalar_emission_factor_applied_to_Kemit": not use_absolute,
        "two_channel_factor_applied_inside_emission_hazard": use_absolute,
    }
    return K_cleave, K_emit, metadata


def _relative_change(initial: float, final: float, floor: float) -> float:
    return abs(float(final) - float(initial)) / max(abs(float(initial)), float(floor))


def _maximum_log_change(initial: np.ndarray, final: np.ndarray) -> float:
    a = np.asarray(initial, dtype=float).reshape(-1)
    b = np.asarray(final, dtype=float).reshape(-1)
    if a.shape != b.shape:
        raise ValueError("hazard arrays must have the same shape")
    active = np.maximum(a, b) > 1.0e-300
    if not np.any(active):
        return 0.0
    la = np.log(np.maximum(a[active], 1.0e-300))
    lb = np.log(np.maximum(b[active], 1.0e-300))
    return float(np.max(np.abs(lb - la)))


class PersistentSiteStateCoupledStochasticEmissionFrontEngineV1005182(
    PersistentSiteStochasticEmissionMovingTipFrontEngineV100518
):
    """Exact stochastic events with adaptive outer load-ramp resolution."""

    state_coupled_emission_adaptivity_active = True
    persistent_two_channel_absolute_opening_drive = True

    _adaptive_reference_target_default = 0.05
    _adaptive_max_log_hazard_change_default = 0.25
    _adaptive_max_backstress_fraction_default = 0.05
    _adaptive_max_geometry_fraction_default = 0.05
    _adaptive_stress_floor_Pa_default = 1.0e6

    @classmethod
    def configure_stochastic(cls, **kwargs) -> None:
        super().configure_stochastic(**kwargs)
        cls._adaptive_reference_target_default = max(
            float(os.environ.get("EMISSION_ADAPT_REFERENCE_TARGET", "0.05")),
            1.0e-12,
        )
        cls._adaptive_max_log_hazard_change_default = max(
            float(os.environ.get("EMISSION_ADAPT_MAX_LOG_HAZARD_CHANGE", "0.25")),
            1.0e-8,
        )
        cls._adaptive_max_backstress_fraction_default = max(
            float(os.environ.get("EMISSION_ADAPT_MAX_BACKSTRESS_FRACTION", "0.05")),
            1.0e-8,
        )
        cls._adaptive_max_geometry_fraction_default = max(
            float(os.environ.get("EMISSION_ADAPT_MAX_GEOMETRY_FRACTION", "0.05")),
            1.0e-8,
        )
        cls._adaptive_stress_floor_Pa_default = max(
            float(os.environ.get("EMISSION_ADAPT_STRESS_FLOOR_PA", "1e6")),
            1.0,
        )

    def reset(self):
        super().reset()
        self.emission_K_prev_Pa_sqrt_m = 0.0
        self.emission_adaptive_reference_target = float(
            type(self)._adaptive_reference_target_default
        )
        self.emission_adaptive_max_log_hazard_change = float(
            type(self)._adaptive_max_log_hazard_change_default
        )
        self.emission_adaptive_max_backstress_fraction = float(
            type(self)._adaptive_max_backstress_fraction_default
        )
        self.emission_adaptive_max_geometry_fraction = float(
            type(self)._adaptive_max_geometry_fraction_default
        )
        self.emission_adaptive_stress_floor_Pa = float(
            type(self)._adaptive_stress_floor_Pa_default
        )
        self.emission_adaptive_prediction_count = 0
        self.emission_adaptive_rejection_requests = 0
        self.emission_last_adaptive_prediction: dict[str, Any] = {}

    def _capture_state(self) -> dict[str, Any]:
        snapshot = super()._capture_state()
        if hasattr(self, "emission_K_prev_Pa_sqrt_m"):
            snapshot.update(
                {
                    "emission_K_prev_Pa_sqrt_m": float(
                        self.emission_K_prev_Pa_sqrt_m
                    ),
                    "emission_adaptive_prediction_count": int(
                        self.emission_adaptive_prediction_count
                    ),
                    "emission_adaptive_rejection_requests": int(
                        self.emission_adaptive_rejection_requests
                    ),
                    "emission_last_adaptive_prediction": copy.deepcopy(
                        self.emission_last_adaptive_prediction
                    ),
                }
            )
        return snapshot

    def _restore_state(self, snapshot: dict[str, Any]) -> None:
        super()._restore_state(snapshot)
        if "emission_K_prev_Pa_sqrt_m" not in snapshot:
            return
        self.emission_K_prev_Pa_sqrt_m = float(
            snapshot["emission_K_prev_Pa_sqrt_m"]
        )
        self.emission_adaptive_prediction_count = int(
            snapshot["emission_adaptive_prediction_count"]
        )
        self.emission_adaptive_rejection_requests = int(
            snapshot["emission_adaptive_rejection_requests"]
        )
        self.emission_last_adaptive_prediction = copy.deepcopy(
            snapshot["emission_last_adaptive_prediction"]
        )

    def _state_change_prediction(
        self,
        *,
        K_emit: float,
        T_K: float,
        dt_s: float,
        drive_factors: np.ndarray,
        tau_signed_Pa: np.ndarray,
    ) -> dict[str, Any]:
        dt = max(float(dt_s), 0.0)
        factors = np.asarray(drive_factors, dtype=float).reshape(-1)
        tau = np.asarray(tau_signed_Pa, dtype=float).reshape(-1)

        geometry0 = self.mpz_state.source_geometry()
        rho0, back0 = self.mpz_state.backstress()
        opening0 = self._opening_stress(self.emission_K_prev_Pa_sqrt_m)[0]
        opening1 = self._opening_stress(float(K_emit))[0]

        hazard_start = _emission_hazards(
            self.mpz_state,
            T_K=float(T_K),
            opening_stress_Pa=opening0,
            drive_factors=factors,
            tau_signed_Pa=tau,
            rate_function=self._emission_rate_per_site,
        )
        hazard_endpoint = _emission_hazards(
            self.mpz_state,
            T_K=float(T_K),
            opening_stress_Pa=opening1,
            drive_factors=factors,
            tau_signed_Pa=tau,
            rate_function=self._emission_rate_per_site,
        )

        trial_state = self.mpz_state.copy()
        deterministic = trial_state.emit_persistent(
            dt_s=dt,
            T_K=float(T_K),
            opening_stress_Pa=opening1,
            drive_factors=factors,
            tau_signed_Pa=tau,
            rate_function=self._emission_rate_per_site,
        )
        geometry1 = trial_state.source_geometry()
        rho1, back1 = trial_state.backstress()
        hazard_feedback = _emission_hazards(
            trial_state,
            T_K=float(T_K),
            opening_stress_Pa=opening1,
            drive_factors=factors,
            tau_signed_Pa=tau,
            rate_function=self._emission_rate_per_site,
        )

        drive1 = np.asarray(hazard_endpoint["drive_Pa"], dtype=float)
        backstress_fraction = float(
            np.max(
                np.abs(np.asarray(back1) - np.asarray(back0))
                / np.maximum(np.abs(drive1), self.emission_adaptive_stress_floor_Pa)
            )
        )
        radius_fraction = _relative_change(
            geometry0["tip_radius_m"],
            geometry1["tip_radius_m"],
            self.mpz_state.b_m,
        )
        width_fraction = _relative_change(
            geometry0["front_width_m"],
            geometry1["front_width_m"],
            self.mpz_state.b_m,
        )
        area_fraction = _relative_change(
            geometry0["source_area_m2"],
            geometry1["source_area_m2"],
            self.mpz_state.b_m ** 2,
        )
        multiplicity_fraction = _relative_change(
            geometry0["multiplicity_per_system"],
            geometry1["multiplicity_per_system"],
            1.0,
        )
        geometry_fraction = max(
            radius_fraction,
            width_fraction,
            area_fraction,
            multiplicity_fraction,
        )
        load_log_change = _maximum_log_change(
            hazard_start["aggregate_hazard_s"],
            hazard_endpoint["aggregate_hazard_s"],
        )
        feedback_log_change = _maximum_log_change(
            hazard_endpoint["aggregate_hazard_s"],
            hazard_feedback["aggregate_hazard_s"],
        )
        log_hazard_change = max(load_log_change, feedback_log_change)

        ratios = {
            "log_hazard": (
                log_hazard_change / self.emission_adaptive_max_log_hazard_change
            ),
            "backstress": (
                backstress_fraction
                / self.emission_adaptive_max_backstress_fraction
            ),
            "geometry": (
                geometry_fraction / self.emission_adaptive_max_geometry_fraction
            ),
        }
        controlling = max(ratios, key=ratios.get)
        normalized = self.emission_adaptive_reference_target * max(
            max(ratios.values()), 0.0
        )
        return {
            "adaptive_increment": float(normalized),
            "controlling_component": controlling,
            "component_ratios": ratios,
            "load_log_hazard_change": float(load_log_change),
            "feedback_log_hazard_change": float(feedback_log_change),
            "backstress_fraction": float(backstress_fraction),
            "geometry_fraction": float(geometry_fraction),
            "tip_radius_fraction": float(radius_fraction),
            "front_width_fraction": float(width_fraction),
            "source_area_fraction": float(area_fraction),
            "multiplicity_fraction": float(multiplicity_fraction),
            "opening_stress_start_Pa": float(opening0),
            "opening_stress_endpoint_Pa": float(opening1),
            "K_emit_start_Pa_sqrt_m": float(self.emission_K_prev_Pa_sqrt_m),
            "K_emit_endpoint_Pa_sqrt_m": float(K_emit),
            "predicted_source_activations": float(
                deterministic.get("source_activations", 0.0)
            ),
            "predicted_activations_by_system": np.asarray(
                deterministic.get("activations_by_system", np.zeros(self.mpz_state.n_systems)),
                dtype=float,
            ),
            "geometry_initial": geometry0,
            "geometry_predicted": geometry1,
            "rho_back_initial_by_system_m2": np.asarray(rho0, dtype=float),
            "rho_back_predicted_by_system_m2": np.asarray(rho1, dtype=float),
            "sigma_back_initial_by_system_Pa": np.asarray(back0, dtype=float),
            "sigma_back_predicted_by_system_Pa": np.asarray(back1, dtype=float),
            "aggregate_hazard_start_by_system_s": np.asarray(
                hazard_start["aggregate_hazard_s"], dtype=float
            ),
            "aggregate_hazard_endpoint_by_system_s": np.asarray(
                hazard_endpoint["aggregate_hazard_s"], dtype=float
            ),
            "aggregate_hazard_feedback_by_system_s": np.asarray(
                hazard_feedback["aggregate_hazard_s"], dtype=float
            ),
        }

    def predict_clock_increment_drives(self, K_cleave, K_emit, T, dt):
        factors, tau_signed, _ = self._two_channel_drive()
        prediction = self._state_change_prediction(
            K_emit=float(K_emit),
            T_K=float(T),
            dt_s=max(float(dt), 0.0),
            drive_factors=factors,
            tau_signed_Pa=tau_signed,
        )
        self.emission_adaptive_prediction_count += 1
        self.emission_last_adaptive_prediction = copy.deepcopy(prediction)
        adaptive_increment = float(prediction["adaptive_increment"])
        if adaptive_increment > self.emission_adaptive_reference_target:
            self.emission_adaptive_rejection_requests += 1
            return adaptive_increment
        cleavage_increment = super().predict_clock_increment_drives(
            K_cleave,
            K_emit,
            T,
            dt,
        )
        return float(max(cleavage_increment, adaptive_increment))

    def step_drives(self, K_cleave, K_emit, T, dt, metadata=None):
        out = super().step_drives(
            K_cleave,
            K_emit,
            T,
            dt,
            metadata=metadata,
        )
        self.emission_K_prev_Pa_sqrt_m = float(K_emit)
        out.update(
            {
                "emission_state_coupled_adaptivity_active": True,
                "emission_adaptive_schema": ADAPTIVE_SCHEMA,
                "emission_adaptive_prediction_count": int(
                    self.emission_adaptive_prediction_count
                ),
                "emission_adaptive_rejection_requests": int(
                    self.emission_adaptive_rejection_requests
                ),
                "emission_last_adaptive_prediction": copy.deepcopy(
                    self.emission_last_adaptive_prediction
                ),
                "emission_scalar_directional_factor_applied_to_opening_K": False,
                "emission_two_channel_factor_applied_once": True,
                "emission_accepted_update_exact_event_localized": True,
                "emission_post_onset_mean_field_burst": False,
            }
        )
        return out

    @classmethod
    def audit_payload(cls) -> dict[str, Any]:
        payload = super().audit_payload()
        payload["model_id"] = MODEL_ID
        emission = dict(payload.get("stochastic_emission", {}) or {})
        emission.update(
            {
                "adaptive_schema": ADAPTIVE_SCHEMA,
                "accepted_update": "exact_event_localized_one_activation_packets",
                "post_onset_mean_field_burst": False,
                "emission_count_cap_applied_to_physics": False,
                "outer_load_increment_subdivided_from_state_change": True,
                "adaptive_state_variables": [
                    "aggregate_hazard",
                    "Taylor_backstress",
                    "tip_radius",
                    "physical_front_width",
                    "source_area",
                    "persistent_site_multiplicity",
                ],
                "scalar_emission_factor_applied_to_opening_K": False,
                "signed_channel_factors_applied_once": True,
                "backstress_recomputed_after_each_event": True,
                "tip_radius_recomputed_after_each_event": True,
                "front_width_recomputed_after_each_event": True,
                "site_multiplicity_recomputed_after_each_event": True,
            }
        )
        payload["stochastic_emission"] = emission
        return payload


__all__ = [
    "ADAPTIVE_SCHEMA",
    "EMISSION_HAZARD_SCHEMA",
    "EMISSION_TRANSPORT_SCHEMA",
    "MODEL_ID",
    "PersistentSiteStateCoupledStochasticEmissionFrontEngineV1005182",
    "two_channel_absolute_opening_drives",
]
