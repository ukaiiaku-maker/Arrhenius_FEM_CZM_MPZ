"""Cycle-block loading adapter for the v10.0.5.2 kinetic-CZM front engine.

This module changes only the applied loading history. Every phase is advanced by
``super().integrate_kinetics`` so the certified MPZ transport, source depletion,
parallel opening/emission channels, active shielding, cleavage clock, and source
refresh laws remain authoritative.
"""
from __future__ import annotations

import copy
import math
from dataclasses import dataclass, asdict
from typing import Any, Mapping

import numpy as np


@dataclass(frozen=True)
class FatigueLoadingConfigV10053:
    R: float = 0.1
    frequency_Hz: float = 1.0e3
    n_phase: int = 96
    closure_clip: bool = True
    quadrature: str = "midpoint"

    def validate(self) -> "FatigueLoadingConfigV10053":
        if not math.isfinite(self.R):
            raise ValueError("R must be finite")
        if not math.isfinite(self.frequency_Hz) or self.frequency_Hz <= 0.0:
            raise ValueError("frequency_Hz must be finite and positive")
        if int(self.n_phase) < 8:
            raise ValueError("n_phase must be at least 8")
        if self.quadrature != "midpoint":
            raise ValueError("only midpoint cycle quadrature is supported")
        return self

    @property
    def period_s(self) -> float:
        return 1.0 / float(self.frequency_Hz)

    def phase_factors(self) -> np.ndarray:
        phase = (np.arange(int(self.n_phase), dtype=float) + 0.5) * (
            2.0 * np.pi / int(self.n_phase)
        )
        mean = 0.5 * (1.0 + float(self.R))
        amplitude = 0.5 * (1.0 - float(self.R))
        factors = mean + amplitude * np.cos(phase)
        if self.closure_clip:
            factors = np.maximum(factors, 0.0)
        return factors


class FatigueKineticsMixinV10053:
    """Interpret an outer physical-time interval as a cyclic K history."""

    supports_progressive_fatigue_v10053 = True
    fatigue_loading_changes_constitutive_physics = False

    def configure_fatigue_v10053(
        self, config: FatigueLoadingConfigV10053 | None
    ) -> None:
        self._fatigue_v10053_config = (
            None if config is None else copy.deepcopy(config).validate()
        )
        self._fatigue_v10053_active = config is not None
        self._fatigue_v10053_last = {}

    @property
    def fatigue_config_v10053(self) -> FatigueLoadingConfigV10053 | None:
        return getattr(self, "_fatigue_v10053_config", None)

    @staticmethod
    def _sum_mapping(target: dict[str, Any], source: Mapping[str, Any]) -> None:
        """Fallback diagnostics-aware aggregator."""
        for key, value in source.items():
            if isinstance(value, (bool, np.bool_)):
                continue
            if isinstance(value, (int, float, np.integer, np.floating)):
                target[key] = float(target.get(key, 0.0)) + float(value)
                continue
            if key == "dN_emit_per_system":
                current = np.asarray(value, dtype=float)
                previous = np.asarray(target.get(key, np.zeros_like(current)), dtype=float)
                if previous.shape != current.shape:
                    raise RuntimeError(f"fatigue channel shape changed for {key}")
                target[key] = (previous + current).tolist()
            elif key in {
                "lambda_emit_per_system_s-1",
                "emission_probability_per_system",
                "slip_system_drive_factors",
                "sigma_emission_backstress_per_system_Pa",
                "sigma_emission_effective_per_system_Pa",
                "emission_backstress_density_per_system_m2",
            }:
                target[key] = np.asarray(value, dtype=float).tolist()

    def _aggregate_mapping(
        self, target: dict[str, Any], source: Mapping[str, Any]
    ) -> None:
        aggregator = getattr(self, "_sum_numeric", None)
        if callable(aggregator):
            aggregator(target, source)
        else:  # pragma: no cover
            self._sum_mapping(target, source)

    def _integrate_fatigue_interval_v10053(
        self,
        K_open_max: float,
        K_cleave_max: float,
        T_K: float,
        dt_s: float,
        *,
        system_weights: np.ndarray | None = None,
    ) -> dict[str, Any]:
        cfg = self.fatigue_config_v10053
        if cfg is None:
            return super().integrate_kinetics(
                K_open_max,
                K_cleave_max,
                T_K,
                dt_s,
                system_weights=system_weights,
            )

        requested = max(float(dt_s), 0.0)
        factors = cfg.phase_factors()
        phase_dt = requested / float(factors.size) if factors.size else requested
        plastic: dict[str, Any] = {}
        advance: dict[str, Any] = {}
        dB = 0.0
        micro_advance = 0.0
        consumed = 0.0
        internal = 0
        fired = False
        last: dict[str, Any] | None = None
        sigma_cleave_values: list[float] = []
        sigma_emit_values: list[float] = []

        for phase_index, factor in enumerate(factors):
            if phase_dt <= 0.0:
                break
            Ko = max(float(K_open_max) * float(factor), 0.0)
            Kc = max(float(K_cleave_max) * float(factor), 0.0)
            row = super().integrate_kinetics(
                Ko,
                Kc,
                T_K,
                phase_dt,
                system_weights=system_weights,
            )
            last = row
            channels_phase = dict(row.get("channels", {}))
            sigma_cleave_values.append(
                float(channels_phase.get("sigma_cleave_eff_Pa", 0.0))
            )
            sigma_emit_values.append(
                float(channels_phase.get("sigma_emission_effective_Pa", 0.0))
            )
            self._aggregate_mapping(plastic, row.get("plastic", {}))
            self._aggregate_mapping(advance, row.get("advance", {}))
            dB += float(row.get("dB", 0.0))
            micro_advance += float(row.get("micro_advance_step_m", 0.0))
            consumed += float(row.get("dt_consumed_s", phase_dt))
            internal += int(row.get("internal_substeps", 0))
            if bool(row.get("fired", False)):
                fired = True
                break

        if last is None:
            channels = self.stress_channels(
                float(K_open_max), float(K_cleave_max), system_weights
            )
            last = {
                "lambda_c_effective_s-1": 0.0,
                "lambda_c_raw_s-1": 0.0,
                "G_cleave_eff_eV": 0.0,
                "channels": channels,
            }

        cycles_requested = requested * float(cfg.frequency_Hz)
        cycles_consumed = consumed * float(cfg.frequency_Hz)
        Kmax_exact = max(float(K_cleave_max), 0.0)
        Kmin_exact = float(cfg.R) * Kmax_exact
        if cfg.closure_clip:
            Kmin_exact = max(Kmin_exact, 0.0)
        result = {
            "fired": fired,
            "n_fire": 1 if fired else 0,
            "dB": dB,
            "micro_advance_step_m": micro_advance,
            "micro_advance_total_m": float(self.micro_advance_total_m),
            "checkpoint_committed_total_m": float(self.checkpoint_advance_total_m),
            "dt_consumed_s": consumed,
            "dt_unused_s": max(requested - consumed, 0.0),
            "internal_substeps": internal,
            "lambda_c_effective_s-1": float(
                last.get("lambda_c_effective_s-1", 0.0)
            ),
            "lambda_c_raw_s-1": float(last.get("lambda_c_raw_s-1", 0.0)),
            "G_cleave_eff_eV": float(last.get("G_cleave_eff_eV", 0.0)),
            "plastic": plastic,
            "advance": advance,
            "channels": dict(last.get("channels", {})),
            "fatigue_loading_v10053": True,
            "fatigue_cycles_requested": cycles_requested,
            "fatigue_cycles_consumed": cycles_consumed,
            "fatigue_R": float(cfg.R),
            "fatigue_frequency_Hz": float(cfg.frequency_Hz),
            "fatigue_n_phase": int(cfg.n_phase),
            "fatigue_Kmax_Pa_sqrt_m": Kmax_exact,
            "fatigue_Kmin_Pa_sqrt_m": Kmin_exact,
            "fatigue_DeltaK_Pa_sqrt_m": Kmax_exact - Kmin_exact,
            "fatigue_avg_sigma_cleave_eff_Pa": (
                float(np.mean(sigma_cleave_values)) if sigma_cleave_values else 0.0
            ),
            "fatigue_max_sigma_cleave_eff_Pa": (
                float(np.max(sigma_cleave_values)) if sigma_cleave_values else 0.0
            ),
            "fatigue_avg_sigma_emission_effective_Pa": (
                float(np.mean(sigma_emit_values)) if sigma_emit_values else 0.0
            ),
            "fatigue_phase_blocks_completed": phase_index + 1 if factors.size else 0,
            "fatigue_cycle_jump_operator_split": True,
        }
        self._fatigue_v10053_last = copy.deepcopy(result)
        return result

    def integrate_kinetics(
        self,
        K_open: float,
        K_cleave: float,
        T_K: float,
        dt_s: float,
        *,
        system_weights: np.ndarray | None = None,
    ) -> dict[str, Any]:
        if not bool(getattr(self, "_fatigue_v10053_active", False)):
            return super().integrate_kinetics(
                K_open,
                K_cleave,
                T_K,
                dt_s,
                system_weights=system_weights,
            )
        return self._integrate_fatigue_interval_v10053(
            K_open,
            K_cleave,
            T_K,
            dt_s,
            system_weights=system_weights,
        )

    def predict_fatigue_cycle(
        self,
        waveform: Any,
        T_K: float,
        n_phase: int,
        *,
        system_weights: np.ndarray | None = None,
    ) -> dict[str, Any]:
        """Transactional one-cycle predictor used only for adaptive block sizing."""
        configured = self.fatigue_config_v10053
        cfg = FatigueLoadingConfigV10053(
            R=float(getattr(waveform, "R", configured.R if configured else 0.1)),
            frequency_Hz=float(
                getattr(
                    waveform,
                    "frequency_Hz",
                    configured.frequency_Hz if configured else 1.0e3,
                )
            ),
            n_phase=int(n_phase),
            closure_clip=bool(
                getattr(
                    waveform,
                    "closure_clip",
                    configured.closure_clip if configured else True,
                )
            ),
        ).validate()
        snapshot = self.snapshot_kinetic_state()
        active0 = bool(getattr(self, "_fatigue_v10053_active", False))
        cfg0 = configured
        last0 = copy.deepcopy(getattr(self, "_fatigue_v10053_last", {}))
        try:
            self.configure_fatigue_v10053(cfg)
            result = self._integrate_fatigue_interval_v10053(
                float(getattr(waveform, "Kmax")),
                float(getattr(waveform, "Kmax")),
                float(T_K),
                cfg.period_s,
                system_weights=system_weights,
            )
            plastic = dict(result.get("plastic", {}))
            channels = dict(result.get("channels", {}))
            return {
                "dN_emit_per_cycle": float(plastic.get("dN_emit", 0.0)),
                "dN_store_per_cycle": float(plastic.get("dN_trapped", 0.0)),
                "dN_mobile_per_cycle": float(plastic.get("dN_emit", 0.0)),
                "dN_escape_per_cycle": float(plastic.get("dN_escaped", 0.0)),
                "dN_peierls_per_cycle": float(plastic.get("peierls_events", 0.0)),
                "dN_taylor_per_cycle": float(plastic.get("taylor_completions", 0.0)),
                "mu_cleave_per_cycle": float(result.get("dB", 0.0)),
                "avg_sigma_tip": float(
                    result.get("fatigue_avg_sigma_cleave_eff_Pa", 0.0)
                ),
                "max_sigma_tip": float(
                    result.get("fatigue_max_sigma_cleave_eff_Pa", 0.0)
                ),
                "avg_sigma_emit_eff": float(
                    result.get(
                        "fatigue_avg_sigma_emission_effective_Pa",
                        channels.get("sigma_emission_effective_Pa", 0.0),
                    )
                ),
                "fatigue_predictor_uses_v10052_state_updates": True,
            }
        finally:
            self.restore_kinetic_state(snapshot)
            self._fatigue_v10053_active = active0
            self._fatigue_v10053_config = cfg0
            self._fatigue_v10053_last = last0

    def audit_payload(self) -> dict[str, Any]:
        payload = super().audit_payload()
        cfg = self.fatigue_config_v10053
        payload.update(
            {
                "progressive_fatigue_v10053_supported": True,
                "fatigue_loading_config_v10053": None if cfg is None else asdict(cfg),
                "fatigue_loading_changes_constitutive_physics": False,
                "fatigue_cycle_jump_operator_split": True,
            }
        )
        return payload


__all__ = ["FatigueLoadingConfigV10053", "FatigueKineticsMixinV10053"]
