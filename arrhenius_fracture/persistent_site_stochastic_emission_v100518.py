"""Stochastic signed-emission clocks for FEM/CZM v10.0.5.18.

Cleavage retains the v10.0.5.16 exponential first-passage threshold and
threshold-correlated crack-advance reward.  This release replaces the
continuous expected emission increment with one exponential integrated-hazard
clock per signed slip channel.  Emission events are localized inside each
plastic half-step and interleaved with the existing PF transport map.
"""
from __future__ import annotations

import copy
import hashlib
import math
import os
from typing import Any, Callable

import numpy as np

from .persistent_site_pf_update_v100515 import (
    PF_REFERENCE_COMMIT,
    PF_UPDATE_MAP,
    _active_mass,
    _evolve_wake_pf,
    _exchange_signed_species,
    _population_weighted_velocity,
    _rates,
    _retained_forest_density_m2,
    _stress_profile_Pa,
    _wake_mass,
    fractional_advect_forward,
)
from .persistent_site_stochastic_tip_v100516 import (
    PersistentSitePFStochasticMovingTipFrontEngineV100516,
    draw_hazard_threshold,
)

MODEL_ID = "FEM_CZM_stochastic_cleavage_and_signed_emission_v10_0_5_18"
EMISSION_HAZARD_SCHEMA = (
    "v10.0.5.18_independent_signed_channel_exponential_integrated_hazard"
)
EMISSION_TRANSPORT_SCHEMA = (
    "v10.0.5.18_event_localized_emission_interleaved_PF_v10_2_22_transport"
)


def derive_stream_seed(case_seed: int, stream_name: str, channel: int) -> int:
    """Derive a stable domain-separated 64-bit RNG seed."""
    token = (
        f"v10.0.5.18|case_seed={int(case_seed)}|stream={stream_name}|"
        f"channel={int(channel)}"
    ).encode("utf-8")
    return int.from_bytes(hashlib.sha256(token).digest()[:8], "big")


def _draw_exponential(rng: np.random.Generator, floor: float) -> float:
    return draw_hazard_threshold("exponential", rng, floor)


def _sum_numeric(target: dict[str, float], source: dict[str, Any]) -> None:
    for key, value in source.items():
        if isinstance(value, (bool, np.bool_)):
            continue
        if isinstance(value, (int, float, np.integer, np.floating)):
            target[key] = target.get(key, 0.0) + float(value)


def advance_pf_transport_only_v10222(
    state,
    *,
    dt_s: float,
    T_K: float,
    opening_stress_Pa: float,
) -> dict[str, Any]:
    """Advance the unchanged PF transport map without creating an emission event."""
    dt = max(float(dt_s), 0.0)
    if dt <= 0.0:
        return {
            "dN_emit": 0.0,
            "dN_trapped": 0.0,
            "dN_released": 0.0,
            "dN_detrapped": 0.0,
            "dN_recovered": 0.0,
            "dN_escaped": 0.0,
            "transport_substeps": 0,
            "transport_integrator": PF_UPDATE_MAP,
            "explicit_recovery_active": False,
            "line_content_balance_signed": 0.0,
            "line_content_balance_relative_error": 0.0,
        }

    initial_inventory = _active_mass(state) + _wake_mass(state)
    escaped_before = float(state.escaped_total)
    discarded_before = float(state.wake_discarded_total)

    stress = _stress_profile_Pa(state, opening_stress_Pa)
    rho = _retained_forest_density_m2(state, wake=False)
    peierls, release, jump, velocity, encounter, _ = _rates(
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
    balance = initial_inventory - (
        final_inventory + escaped_increment + discarded_increment
    )
    scale = max(
        abs(initial_inventory),
        abs(final_inventory + escaped_increment + discarded_increment),
        1.0e-300,
    )
    relative_error = abs(balance) / scale
    if relative_error > 5.0e-11:
        raise RuntimeError(
            "stochastic-emission transport-only line-content conservation failed: "
            f"balance={balance:.9e}, relative_error={relative_error:.9e}"
        )

    diagnostics = {
        "dN_emit": 0.0,
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
            "event_localized_emission_interleaved_with_exact_exchange_"
            "zero_recovery_scalar_advection"
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


def _emission_hazards(
    state,
    *,
    T_K: float,
    opening_stress_Pa: float,
    drive_factors: np.ndarray,
    tau_signed_Pa: np.ndarray,
    rate_function: Callable[[float, float], float],
) -> dict[str, Any]:
    factors = np.asarray(drive_factors, dtype=float).reshape(-1)
    tau = np.asarray(tau_signed_Pa, dtype=float).reshape(-1)
    if factors.shape != (state.n_systems,) or tau.shape != (state.n_systems,):
        raise ValueError("one stochastic emission clock is required per signed channel")
    if not np.all(np.isfinite(factors)) or not np.all(np.isfinite(tau)):
        raise ValueError("signed emission drive channels must be finite")

    geometry = state.source_geometry()
    multiplicity = float(geometry["multiplicity_per_system"])
    rho, sigma_back = state.backstress()
    drive = np.maximum(factors * max(float(opening_stress_Pa), 0.0), 0.0)
    signs = np.sign(tau)
    sigma = np.maximum(drive - sigma_back, 0.0)
    rates = np.zeros(state.n_systems, dtype=float)
    for system in range(state.n_systems):
        if signs[system] != 0.0 and sigma[system] > 0.0:
            rates[system] = max(float(rate_function(sigma[system], T_K)), 0.0)
    hazards = multiplicity * rates
    if not np.all(np.isfinite(hazards)) or np.any(hazards < 0.0):
        raise RuntimeError("stochastic emission produced an invalid aggregate hazard")
    return {
        "geometry": geometry,
        "multiplicity": multiplicity,
        "rho_back_m2": rho,
        "sigma_back_Pa": sigma_back,
        "drive_Pa": drive,
        "sigma_effective_Pa": sigma,
        "rate_per_site_s": rates,
        "aggregate_hazard_s": hazards,
        "signs": signs,
    }


def _deposit_one_emission_event(state, system: int, sign: float) -> float:
    if sign == 0.0:
        raise RuntimeError("cannot deposit an emission event with zero Burgers sign")
    channel = int(system)
    line_content = float(state.activation_to_line_content_by_system[channel])
    if not math.isfinite(line_content) or line_content <= 0.0:
        raise RuntimeError("activation-to-line conversion must be positive and finite")
    nsrc = max(min(int(state._source_bin_count()), int(state.n_bins)), 1)
    amount = line_content / float(nsrc)
    if sign > 0.0:
        state.mobile_positive[channel, :nsrc] += amount
        state.accumulated_slip_positive[channel, :nsrc] += amount
    else:
        state.mobile_negative[channel, :nsrc] += amount
        state.accumulated_slip_negative[channel, :nsrc] += amount
    state.emitted_total += line_content
    return line_content


def _linear_hazard_crossing_time(
    rate0: float,
    rate1: float,
    interval_s: float,
    required_action: float,
) -> float | None:
    """Crossing time for a linearly varying nonnegative aggregate hazard."""
    h = max(float(interval_s), 0.0)
    need = max(float(required_action), 0.0)
    l0 = max(float(rate0), 0.0)
    l1 = max(float(rate1), 0.0)
    if need <= 0.0:
        return 0.0
    total = 0.5 * (l0 + l1) * h
    if h <= 0.0 or total + 1.0e-14 * max(total, need, 1.0) < need:
        return None
    slope = (l1 - l0) / h
    if abs(slope) <= 1.0e-14 * max(l0, l1, 1.0) / max(h, 1.0e-300):
        if l0 <= 0.0:
            return h
        return min(max(need / l0, 0.0), h)
    discriminant = l0 * l0 + 2.0 * slope * need
    if discriminant < 0.0 and discriminant > -1.0e-12 * max(l0 * l0, 1.0):
        discriminant = 0.0
    if discriminant < 0.0:
        return h
    denominator = l0 + math.sqrt(discriminant)
    if denominator <= 0.0:
        return h
    return min(max(2.0 * need / denominator, 0.0), h)


class PersistentSiteStochasticEmissionMovingTipFrontEngineV100518(
    PersistentSitePFStochasticMovingTipFrontEngineV100516
):
    """Stochastic cleavage plus independent stochastic signed-emission clocks."""

    stochastic_emission_active = True
    _emission_case_seed_default = 0
    _emission_minimum_threshold_default = 1.0e-12
    _emission_max_action_substep_default = 0.05
    _emission_max_events_per_half_step_default = 10000

    @classmethod
    def configure_stochastic(cls, **kwargs) -> None:
        super().configure_stochastic(**kwargs)
        cleavage_seed = int(kwargs.get("hazard_seed", 0))
        raw = os.environ.get("EMISSION_HAZARD_SEED", "").strip()
        case_seed = cleavage_seed if not raw else int(raw)
        if case_seed < 0:
            raise ValueError("EMISSION_HAZARD_SEED must be nonnegative")
        cls._emission_case_seed_default = case_seed
        cls._emission_minimum_threshold_default = max(
            float(os.environ.get("EMISSION_HAZARD_MIN_THRESHOLD", "1e-12")),
            1.0e-300,
        )
        cls._emission_max_action_substep_default = max(
            float(os.environ.get("EMISSION_MAX_ACTION_SUBSTEP", "0.05")),
            1.0e-6,
        )
        cls._emission_max_events_per_half_step_default = max(
            int(os.environ.get("EMISSION_MAX_EVENTS_PER_HALF_STEP", "10000")),
            1,
        )

    def reset(self):
        super().reset()
        self.emission_case_seed = int(type(self)._emission_case_seed_default)
        self.emission_minimum_threshold = float(
            type(self)._emission_minimum_threshold_default
        )
        self.emission_max_action_substep = float(
            type(self)._emission_max_action_substep_default
        )
        self.emission_max_events_per_half_step = int(
            type(self)._emission_max_events_per_half_step_default
        )
        self.emission_stream_seeds = np.asarray(
            [
                derive_stream_seed(self.emission_case_seed, "signed_emission", channel)
                for channel in range(self.mpz_state.n_systems)
            ],
            dtype=np.uint64,
        )
        self._emission_rngs = [
            np.random.default_rng(int(seed)) for seed in self.emission_stream_seeds
        ]
        self.emission_threshold_action = np.asarray(
            [
                _draw_exponential(rng, self.emission_minimum_threshold)
                for rng in self._emission_rngs
            ],
            dtype=float,
        )
        self.emission_action_current = np.zeros(self.mpz_state.n_systems, dtype=float)
        self.emission_event_index = np.zeros(self.mpz_state.n_systems, dtype=np.int64)
        self.emission_event_count_total = 0
        self.emission_event_history: list[dict[str, Any]] = []
        self.emission_internal_substeps_total = 0
        self.emission_transport_time_total_s = 0.0

    def _capture_state(self) -> dict[str, Any]:
        snapshot = super()._capture_state()
        if hasattr(self, "_emission_rngs"):
            snapshot.update(
                {
                    "emission_rng_states": [
                        copy.deepcopy(rng.bit_generator.state)
                        for rng in self._emission_rngs
                    ],
                    "emission_threshold_action": self.emission_threshold_action.copy(),
                    "emission_action_current": self.emission_action_current.copy(),
                    "emission_event_index": self.emission_event_index.copy(),
                    "emission_event_count_total": int(self.emission_event_count_total),
                    "emission_event_history": copy.deepcopy(
                        self.emission_event_history
                    ),
                    "emission_internal_substeps_total": int(
                        self.emission_internal_substeps_total
                    ),
                    "emission_transport_time_total_s": float(
                        self.emission_transport_time_total_s
                    ),
                }
            )
        return snapshot

    def _restore_state(self, snapshot: dict[str, Any]) -> None:
        super()._restore_state(snapshot)
        if "emission_rng_states" not in snapshot:
            return
        rngs: list[np.random.Generator] = []
        for state in snapshot["emission_rng_states"]:
            rng = np.random.default_rng()
            rng.bit_generator.state = copy.deepcopy(state)
            rngs.append(rng)
        self._emission_rngs = rngs
        self.emission_threshold_action = np.asarray(
            snapshot["emission_threshold_action"], dtype=float
        ).copy()
        self.emission_action_current = np.asarray(
            snapshot["emission_action_current"], dtype=float
        ).copy()
        self.emission_event_index = np.asarray(
            snapshot["emission_event_index"], dtype=np.int64
        ).copy()
        self.emission_event_count_total = int(
            snapshot["emission_event_count_total"]
        )
        self.emission_event_history = copy.deepcopy(
            snapshot["emission_event_history"]
        )
        self.emission_internal_substeps_total = int(
            snapshot["emission_internal_substeps_total"]
        )
        self.emission_transport_time_total_s = float(
            snapshot["emission_transport_time_total_s"]
        )

    def _stochastic_emission_transport(
        self,
        *,
        dt_s: float,
        T_K: float,
        opening_stress_Pa: float,
        drive_factors: np.ndarray,
        tau_signed_Pa: np.ndarray,
    ) -> dict[str, Any]:
        remaining = max(float(dt_s), 0.0)
        diagnostics: dict[str, float] = {}
        events_this_step: list[dict[str, Any]] = []
        emitted_lines = 0.0
        substeps = 0
        elapsed = 0.0
        last_hazard = None

        while remaining > 0.0:
            substeps += 1
            if substeps > 20 * self.emission_max_events_per_half_step + 1000:
                raise RuntimeError(
                    "stochastic emission exceeded the internal substep budget; "
                    "reduce the outer timestep without changing the event law"
                )
            start = _emission_hazards(
                self.mpz_state,
                T_K=T_K,
                opening_stress_Pa=opening_stress_Pa,
                drive_factors=drive_factors,
                tau_signed_Pa=tau_signed_Pa,
                rate_function=self._emission_rate_per_site,
            )
            hazard0 = np.asarray(start["aggregate_hazard_s"], dtype=float)
            residual = np.maximum(
                self.emission_threshold_action - self.emission_action_current,
                0.0,
            )
            positive = hazard0 > 0.0
            if not np.any(positive):
                transport = advance_pf_transport_only_v10222(
                    self.mpz_state,
                    dt_s=remaining,
                    T_K=T_K,
                    opening_stress_Pa=opening_stress_Pa,
                )
                _sum_numeric(diagnostics, transport)
                elapsed += remaining
                remaining = 0.0
                last_hazard = start
                break

            max_hazard = float(np.max(hazard0))
            h = min(
                remaining,
                self.emission_max_action_substep / max(max_hazard, 1.0e-300),
            )
            predicted = np.divide(
                residual,
                hazard0,
                out=np.full_like(residual, np.inf),
                where=positive,
            )
            next_predicted = float(np.min(predicted))
            if math.isfinite(next_predicted):
                h = min(h, max(1.25 * next_predicted, 1.0e-15))
            h = min(max(h, min(remaining, 1.0e-15)), remaining)

            state_start = self.mpz_state.copy()
            proposal = advance_pf_transport_only_v10222(
                self.mpz_state,
                dt_s=h,
                T_K=T_K,
                opening_stress_Pa=opening_stress_Pa,
            )
            end = _emission_hazards(
                self.mpz_state,
                T_K=T_K,
                opening_stress_Pa=opening_stress_Pa,
                drive_factors=drive_factors,
                tau_signed_Pa=tau_signed_Pa,
                rate_function=self._emission_rate_per_site,
            )
            hazard1 = np.asarray(end["aggregate_hazard_s"], dtype=float)
            crossing_times = np.full(self.mpz_state.n_systems, np.inf)
            for system in range(self.mpz_state.n_systems):
                crossing = _linear_hazard_crossing_time(
                    hazard0[system],
                    hazard1[system],
                    h,
                    residual[system],
                )
                if crossing is not None:
                    crossing_times[system] = crossing

            if not np.any(np.isfinite(crossing_times)):
                self.emission_action_current += 0.5 * (hazard0 + hazard1) * h
                _sum_numeric(diagnostics, proposal)
                elapsed += h
                remaining = max(remaining - h, 0.0)
                last_hazard = end
                continue

            system = int(np.argmin(crossing_times))
            tau = min(max(float(crossing_times[system]), 0.0), h)
            self.mpz_state = state_start
            partial = advance_pf_transport_only_v10222(
                self.mpz_state,
                dt_s=tau,
                T_K=T_K,
                opening_stress_Pa=opening_stress_Pa,
            )
            at_event = _emission_hazards(
                self.mpz_state,
                T_K=T_K,
                opening_stress_Pa=opening_stress_Pa,
                drive_factors=drive_factors,
                tau_signed_Pa=tau_signed_Pa,
                rate_function=self._emission_rate_per_site,
            )
            hazard_event = np.asarray(
                at_event["aggregate_hazard_s"], dtype=float
            )
            self.emission_action_current += 0.5 * (hazard0 + hazard_event) * tau
            self.emission_action_current[system] = self.emission_threshold_action[
                system
            ]
            threshold_completed = float(self.emission_threshold_action[system])
            sign = float(np.asarray(at_event["signs"])[system])
            line_content = _deposit_one_emission_event(
                self.mpz_state,
                system,
                sign,
            )
            emitted_lines += line_content
            event_index = int(self.emission_event_index[system])
            event_record = {
                "system": system,
                "sign": sign,
                "event_index": event_index,
                "time_within_half_step_s": elapsed + tau,
                "threshold_action": threshold_completed,
                "aggregate_hazard_at_event_s": float(hazard_event[system]),
                "rate_per_site_at_event_s": float(
                    np.asarray(at_event["rate_per_site_s"])[system]
                ),
                "multiplicity_at_event": float(at_event["multiplicity"]),
                "line_content_added": line_content,
                "tip_radius_before_event_m": float(
                    start["geometry"]["tip_radius_m"]
                ),
                "tip_radius_after_event_m": float(
                    self.mpz_state.blunted_radius()
                ),
                "front_width_at_event_m": float(
                    at_event["geometry"]["front_width_m"]
                ),
            }
            events_this_step.append(event_record)
            self.emission_event_history.append(copy.deepcopy(event_record))
            self.emission_event_index[system] += 1
            self.emission_event_count_total += 1
            self.emission_action_current[system] = 0.0
            self.emission_threshold_action[system] = _draw_exponential(
                self._emission_rngs[system],
                self.emission_minimum_threshold,
            )
            _sum_numeric(diagnostics, partial)
            elapsed += tau
            remaining = max(remaining - tau, 0.0)
            last_hazard = at_event
            if len(events_this_step) > self.emission_max_events_per_half_step:
                raise RuntimeError(
                    "stochastic emission exceeded EMISSION_MAX_EVENTS_PER_HALF_STEP; "
                    "reduce the outer timestep rather than batching state-changing events"
                )
            if tau <= 0.0 and remaining > 0.0:
                remaining = max(remaining - min(remaining, 1.0e-15), 0.0)

        self.emission_internal_substeps_total += substeps
        self.emission_transport_time_total_s += elapsed
        geometry = self.mpz_state.source_geometry()
        if last_hazard is None:
            last_hazard = _emission_hazards(
                self.mpz_state,
                T_K=T_K,
                opening_stress_Pa=opening_stress_Pa,
                drive_factors=drive_factors,
                tau_signed_Pa=tau_signed_Pa,
                rate_function=self._emission_rate_per_site,
            )
        result: dict[str, Any] = {
            **diagnostics,
            "dN_emit": emitted_lines,
            "source_activations": float(len(events_this_step)),
            "stochastic_emission_events": events_this_step,
            "stochastic_emission_event_count": len(events_this_step),
            "stochastic_emission_event_count_total": int(
                self.emission_event_count_total
            ),
            "stochastic_emission_threshold_action": (
                self.emission_threshold_action.copy()
            ),
            "stochastic_emission_action_current": (
                self.emission_action_current.copy()
            ),
            "stochastic_emission_event_index": self.emission_event_index.copy(),
            "stochastic_emission_stream_seeds": self.emission_stream_seeds.copy(),
            "aggregate_hazard_final_by_system_s": np.asarray(
                last_hazard["aggregate_hazard_s"], dtype=float
            ),
            "rate_final_by_system_s": np.asarray(
                last_hazard["rate_per_site_s"], dtype=float
            ),
            "persistent_site_geometry": geometry,
            "persistent_site_multiplicity_per_system": float(
                geometry["multiplicity_per_system"]
            ),
            "stochastic_emission_schema": EMISSION_HAZARD_SCHEMA,
            "stochastic_emission_transport_schema": EMISSION_TRANSPORT_SCHEMA,
            "emission_internal_substeps": substeps,
            "emission_transport_time_s": elapsed,
            "finite_source_inventory_active": False,
            "source_depletion_active": False,
            "source_refresh_active": False,
            "explicit_recovery_active": False,
        }
        self.mpz_state.last_emission = copy.deepcopy(result)
        return result

    def _plastic_half_step(
        self,
        *,
        dt_s: float,
        T_K: float,
        K_emit: float,
        drive_factors: np.ndarray,
        tau_signed_Pa: np.ndarray,
    ) -> tuple[dict[str, Any], float]:
        dt = max(float(dt_s), 0.0)
        if dt <= 0.0:
            return {
                "dN_emit": 0.0,
                "dN_trapped": 0.0,
                "dN_released": 0.0,
                "dN_escaped": 0.0,
                "dN_recovered": 0.0,
                "stochastic_emission_event_count": 0,
                "transport_integrator": PF_UPDATE_MAP,
            }, 0.0
        opening, _, _, _ = self._opening_stress(K_emit)
        result = self._stochastic_emission_transport(
            dt_s=dt,
            T_K=float(T_K),
            opening_stress_Pa=opening,
            drive_factors=drive_factors,
            tau_signed_Pa=tau_signed_Pa,
        )
        self.W_emit += (
            opening
            * float(self.b)
            * float(self.f.L_pz)
            * max(float(result.get("dN_emit", 0.0)), 0.0)
        )
        self._sync_compat()
        return result, opening

    def predict_clock_increment_drives(self, K_cleave, K_emit, T, dt):
        factors, tau_signed, _ = self._two_channel_drive()
        trial = copy.copy(self)
        trial.mpz_state = self.mpz_state.copy()
        trial._checkpoint_origin_snapshot = None
        trial._geometry_veto_snapshot = None
        trial._checkpoint_fired_last_step = False
        trial.hazard_threshold_history = list(self.hazard_threshold_history)
        trial.stochastic_event_length_history = list(
            self.stochastic_event_length_history
        )
        hazard_rng = np.random.default_rng()
        hazard_rng.bit_generator.state = copy.deepcopy(
            self._hazard_rng.bit_generator.state
        )
        trial._hazard_rng = hazard_rng
        trial._emission_rngs = []
        for source in self._emission_rngs:
            rng = np.random.default_rng()
            rng.bit_generator.state = copy.deepcopy(source.bit_generator.state)
            trial._emission_rngs.append(rng)
        trial.emission_threshold_action = self.emission_threshold_action.copy()
        trial.emission_action_current = self.emission_action_current.copy()
        trial.emission_event_index = self.emission_event_index.copy()
        trial.emission_event_history = copy.deepcopy(self.emission_event_history)
        predicted = trial._integrate_coupled(
            K_cleave=float(K_cleave),
            K_emit=float(K_emit),
            T_K=float(T),
            dt_s=max(float(dt), 0.0),
            drive_factors=factors,
            tau_signed_Pa=tau_signed,
        )
        self.kinetic_prediction_calls += 1
        return float(max(predicted["dB"], 0.0))

    def step_drives(self, K_cleave, K_emit, T, dt, metadata=None):
        out = super().step_drives(K_cleave, K_emit, T, dt, metadata=metadata)
        out.update(
            {
                "stochastic_emission_enabled": True,
                "stochastic_emission_schema": EMISSION_HAZARD_SCHEMA,
                "stochastic_emission_transport_schema": EMISSION_TRANSPORT_SCHEMA,
                "emission_case_seed": int(self.emission_case_seed),
                "emission_stream_seeds": self.emission_stream_seeds.tolist(),
                "emission_threshold_action": (
                    self.emission_threshold_action.tolist()
                ),
                "emission_action_current": self.emission_action_current.tolist(),
                "emission_event_index": self.emission_event_index.tolist(),
                "emission_event_count_total": int(
                    self.emission_event_count_total
                ),
                "emission_rng": "numpy_default_rng_domain_separated_seed",
                "emission_event_packet_semantics": (
                    "one_activation_times_mechanical_line_conversion"
                ),
                "emission_noise_added_to_stress_or_barrier": False,
            }
        )
        return out

    @classmethod
    def audit_payload(cls) -> dict[str, Any]:
        payload = super().audit_payload()
        case_seed = int(cls._emission_case_seed_default)
        payload.update(
            {
                "model_id": MODEL_ID,
                "stochastic_emission": {
                    "schema": EMISSION_HAZARD_SCHEMA,
                    "transport_schema": EMISSION_TRANSPORT_SCHEMA,
                    "case_seed": case_seed,
                    "stream_seeds": [
                        derive_stream_seed(case_seed, "signed_emission", channel)
                        for channel in range(2)
                    ],
                    "distribution": "exponential_unit_mean",
                    "random_transform": "minus_log_uniform",
                    "independent_signed_channel_streams": True,
                    "event_packet": (
                        "one_activation_times_mechanical_line_conversion"
                    ),
                    "backstress_recomputed_after_each_event": True,
                    "tip_radius_recomputed_after_each_event": True,
                    "site_multiplicity_recomputed_after_each_event": True,
                    "transport_interleaved_with_event_localization": True,
                    "trial_rejection_restores_rng_and_clock_state": True,
                    "finite_source_inventory": False,
                    "source_refresh": False,
                    "explicit_recovery": False,
                },
            }
        )
        return payload


__all__ = [
    "EMISSION_HAZARD_SCHEMA",
    "EMISSION_TRANSPORT_SCHEMA",
    "MODEL_ID",
    "PersistentSiteStochasticEmissionMovingTipFrontEngineV100518",
    "advance_pf_transport_only_v10222",
    "derive_stream_seed",
]
