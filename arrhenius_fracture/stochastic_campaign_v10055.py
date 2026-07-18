"""Stochastic VHCF campaign state and hybrid cycle-block scheduler.

v10.0.5.5 preserves the v10.0.5.4 cycle-integrated FEM/CZM/MPZ model but
realizes the finite source-emission history stochastically.  The authoritative
one-cycle predictor remains mean-field so block selection is stable and does
not consume random numbers.  Accepted blocks use bounded binomial source
sampling; transport, trapping, release and recovery then evolve conditional on
the realized source history.
"""
from __future__ import annotations

import copy
from dataclasses import asdict, dataclass
import math
import os
from typing import Any, Callable

import numpy as np

from . import mixed_mode_first_passage_v9_11 as v911
from .kinetic_campaign_czm import (
    CampaignKineticMPZState,
    KineticCampaignCZMConfig,
    apply_pf_manifest_to_mpz_config,
)
from .kinetic_campaign_czm_v1003 import CampaignAwareV1003TipEngineMixin
from .kinetic_campaign_czm_v1005 import STATE_MODEL
from .kinetic_campaign_czm_v10052 import (
    ChannelDiagnosticParallelOpeningEmissionCZMFrontEngine,
)
from .kinetic_fatigue_v10053 import (
    FatigueKineticsMixinV10053,
    FatigueLoadingConfigV10053,
)
from .stochastic_kinetics_v911 import (
    normalize_event_statistics,
    sample_effective_binomial,
)

POINT_RELEASE = "10.0.5.5"
MODEL_ID = "FEM_CZM_parallel_opening_stochastic_fatigue_v10_0_5_5"
ENGINE_REGISTRY_V10055: list[Any] = []


@dataclass(frozen=True)
class HybridSchedulerConfigV10055:
    """Numerical event-count controls for stochastic cycle blocks.

    These values control temporal resolution, not material kinetics.  The
    Arrhenius hazards and finite source capacities remain unchanged.
    """

    enabled: bool = True
    rare_event_target: float = 0.25
    tau_leap_target: float = 3.0
    tau_switch_expected_events: float = 10.0

    def validate(self) -> "HybridSchedulerConfigV10055":
        if self.rare_event_target <= 0.0:
            raise ValueError("rare_event_target must be positive")
        if self.tau_leap_target < self.rare_event_target:
            raise ValueError("tau_leap_target must be >= rare_event_target")
        if self.tau_switch_expected_events < self.tau_leap_target:
            raise ValueError(
                "tau_switch_expected_events must be >= tau_leap_target"
            )
        return self

    @classmethod
    def from_environment(cls) -> "HybridSchedulerConfigV10055":
        return cls(
            enabled=os.environ.get("ARRHENIUS_STOCHASTIC_BLOCKS", "1") != "0",
            rare_event_target=float(
                os.environ.get("ARRHENIUS_RARE_EVENT_TARGET", "0.25")
            ),
            tau_leap_target=float(
                os.environ.get("ARRHENIUS_TAU_LEAP_TARGET", "3.0")
            ),
            tau_switch_expected_events=float(
                os.environ.get("ARRHENIUS_TAU_SWITCH_EXPECTED_EVENTS", "10.0")
            ),
        ).validate()


class StochasticCampaignKineticMPZStateV10055(CampaignKineticMPZState):
    """Campaign MPZ state with stochastic finite-source emission.

    The inherited v9.11 state already owns a reproducible RNG stream and restart
    serialization.  v10.0.5.4 bypassed that sampler in ``emit_exact``; this
    override restores bounded Bernoulli/binomial realizations while preserving
    the deterministic model's mean exactly.
    """

    state_model = "kinetic_campaign_czm_stochastic_v10055"

    def emit_exact(
        self,
        dt_s: float,
        sigma_opening_Pa: float,
        T_K: float,
        system_weights: np.ndarray | None = None,
    ) -> dict[str, Any]:
        mean_field = bool(getattr(self, "_v10055_force_mean_field", False))
        stochastic = (
            normalize_event_statistics(getattr(self, "event_statistics", "deterministic"))
            == "stochastic"
            and bool(getattr(self, "stochastic_emission", False))
            and not mean_field
        )
        if not stochastic:
            out = super().emit_exact(
                dt_s, sigma_opening_Pa, T_K, system_weights
            )
            out["stochastic_source_commit"] = 0.0
            out["stochastic_source_mean_field_predictor"] = float(mean_field)
            return out

        dt = max(float(dt_s), 0.0)
        weights = self._weights(system_weights, self.n_systems)
        sigma_back = self.taylor_backstress_Pa()
        sigma_emit = np.maximum(
            weights * max(float(sigma_opening_Pa), 0.0) - sigma_back, 0.0
        )
        rates = np.maximum(self.manifest.emission.rate(sigma_emit, T_K), 0.0)
        available0 = np.maximum(self.available_sites.copy(), 0.0)
        probability = 1.0 - np.exp(-np.minimum(rates * dt, 700.0))

        emitted_system = np.zeros(self.n_systems, dtype=float)
        for index in range(self.n_systems):
            # Deterministic v10 mean: available * probability * directional weight.
            p_effective = float(
                np.clip(probability[index] * weights[index], 0.0, 1.0)
            )
            emitted_system[index] = sample_effective_binomial(
                self._emission_rng,
                float(available0[index]),
                p_effective,
            )

        self.available_sites = np.clip(
            available0 - emitted_system, 0.0, self.site_capacity
        )
        nsrc = max(min(int(self.cfg.source_bin_count), self.n_bins), 1)
        self.mobile[:, :nsrc] += emitted_system[:, None] / nsrc
        self.accumulated_slip[:, :nsrc] += emitted_system[:, None] / nsrc
        emitted = float(np.sum(emitted_system))
        self.emitted_total += emitted
        self.cumulative_emitted += emitted
        self.last_emitted_per_system = emitted_system.copy()
        self.last_emission_rate_per_system_s = rates.copy()
        self.last_sigma_emit_per_system_Pa = sigma_emit.copy()
        self.stochastic_emission_events += int(np.count_nonzero(emitted_system > 0.0))

        return {
            "dN_emit": emitted,
            "dN_emit_per_system": emitted_system.tolist(),
            "lambda_emit_per_system_s-1": rates.tolist(),
            "sigma_emission_backstress_per_system_Pa": sigma_back.tolist(),
            "sigma_emission_effective_per_system_Pa": sigma_emit.tolist(),
            "emission_backstress_density_per_system_m2": (
                self.last_local_density_per_system_m2.tolist()
            ),
            "stochastic_source_commit": 1.0,
            "stochastic_source_mean_field_predictor": 0.0,
            "stochastic_source_events_block": float(
                np.count_nonzero(emitted_system > 0.0)
            ),
            "stochastic_source_count_block": emitted,
        }

    def diagnostics_campaign(self) -> dict[str, Any]:
        out = super().diagnostics_campaign()
        out.update(
            {
                "event_statistics": str(
                    getattr(self, "event_statistics", "deterministic")
                ),
                "stochastic_emission_active": float(
                    normalize_event_statistics(
                        getattr(self, "event_statistics", "deterministic")
                    )
                    == "stochastic"
                    and bool(getattr(self, "stochastic_emission", False))
                ),
                "stochastic_emission_events": float(
                    getattr(self, "stochastic_emission_events", 0)
                ),
                "stochastic_seed": float(getattr(self, "stochastic_seed", 1)),
                "stochastic_stream": float(getattr(self, "stochastic_stream", 17011)),
            }
        )
        return out


class StochasticPredictorMeanFieldMixinV10055:
    """Keep adaptive prediction deterministic while accepted commits are random."""

    supports_stochastic_fatigue_v10055 = True

    def predict_fatigue_cycle(self, *args, **kwargs):
        state = self.mpz_state
        old = bool(getattr(state, "_v10055_force_mean_field", False))
        state._v10055_force_mean_field = True
        try:
            result = super().predict_fatigue_cycle(*args, **kwargs)
            self._v10055_predictor_mean_field_calls = int(
                getattr(self, "_v10055_predictor_mean_field_calls", 0)
            ) + 1
            result["v10055_predictor_mean_field"] = True
            return result
        finally:
            # The transactional predictor may restore a deep-copied MPZ object.
            self.mpz_state._v10055_force_mean_field = old


def engine_factory_v10055(
    original_build,
    context,
    mm,
    row,
    manifest,
    kinetic_cfg: KineticCampaignCZMConfig,
    registry,
):
    """Build v10.0.5.5 using unchanged v10.0.5.2 constitutive surfaces."""

    def build(parsed_args, material):
        base = original_build(parsed_args, material)
        cfg = v911.build_mpz_config(mm, row)
        apply_pf_manifest_to_mpz_config(cfg, manifest, kinetic_cfg)
        base.f.r0 = 1.0e-6
        base.f.L_pz = cfg.length_m
        base.f.c_blunt = float(manifest.c_blunt)
        base.f.nu0_c = float(manifest.cleavage.attempt_frequency_s)
        base.f.nu0_e = float(manifest.emission.attempt_frequency_s)
        base.f.m_hits = 3.0
        base.f.tau_c = 1.0e-6
        base.f.sigma_cap = 0.0
        base.f.dN_cap = math.inf
        base.f.N_sat = math.inf
        base.f.recover_k = 0.0
        base.f.k_shield = 0.0
        base.f.chi_shield = 0.0
        base.f.v_emb_b3 = 0.0

        Engine = type(
            "StochasticParallelOpeningEmissionKineticCZMEngineV10055",
            (
                StochasticPredictorMeanFieldMixinV10055,
                FatigueKineticsMixinV10053,
                CampaignAwareV1003TipEngineMixin,
                ChannelDiagnosticParallelOpeningEmissionCZMFrontEngine,
            ),
            {
                "supports_progressive_kinetic_czm": True,
                "supports_tensor_resolved_parallel_coupling": True,
                "supports_complete_per_channel_diagnostics": True,
                "supports_progressive_fatigue_v10053": True,
                "supports_stochastic_fatigue_v10055": True,
                "integration_point_release": POINT_RELEASE,
            },
        )
        eng = Engine(
            base.f,
            base.cb,
            base.eb,
            base.G,
            base.nu,
            base.b,
            cfg,
            manifest,
            kinetic_cfg,
        )
        # Replace the expected-value v10 campaign state with the stochastic
        # campaign subclass.  No material or transport parameter is altered.
        eng.mpz_state = StochasticCampaignKineticMPZStateV10055(
            cfg,
            manifest,
            b=base.b,
            G_Pa=base.G,
            kinetic_cfg=kinetic_cfg,
        )
        eng._sync_compat()
        eng._mm_init(context)
        fatigue_active = bool(getattr(parsed_args, "fatigue_cycles", False))
        if fatigue_active:
            eng.configure_fatigue_v10053(
                FatigueLoadingConfigV10053(
                    R=float(getattr(parsed_args, "R", 0.1) or 0.0),
                    frequency_Hz=float(
                        getattr(parsed_args, "frequency_Hz", 1.0e3) or 1.0e3
                    ),
                    n_phase=int(getattr(parsed_args, "n_phase", 96) or 96),
                    closure_clip=not bool(
                        getattr(parsed_args, "no_closure_clip", False)
                    ),
                )
            )
        else:
            eng.configure_fatigue_v10053(None)

        if getattr(eng, "state_model", None) != STATE_MODEL:
            raise RuntimeError(
                "v10.0.5.5 factory constructed the wrong state model: "
                f"{getattr(eng, 'state_model', None)!r}"
            )
        required = (
            "tensor_resolved_parallel_coupling",
            "per_channel_strang_diagnostics_complete",
            "supports_progressive_fatigue_v10053",
            "supports_stochastic_fatigue_v10055",
        )
        missing = [name for name in required if not bool(getattr(eng, name, False))]
        if missing:
            raise RuntimeError(
                "v10.0.5.5 engine lacks required capabilities: " + ", ".join(missing)
            )
        registry.append(eng)
        ENGINE_REGISTRY_V10055.append(eng)
        return eng

    return build


def hybrid_choose_block_factory_v10055(
    original: Callable[..., dict[str, Any]],
    config: HybridSchedulerConfigV10055,
    audit: dict[str, Any],
) -> Callable[..., dict[str, Any]]:
    """Wrap the existing state/hazard controller with event-count resolution."""

    config = copy.deepcopy(config).validate()
    audit.clear()
    audit.update(
        {
            "config": asdict(config),
            "selection_calls": 0,
            "rare_event_calls": 0,
            "tau_leap_calls": 0,
            "quiet_calls": 0,
            "records": [],
        }
    )

    def choose(controller, pred, user_block_cycles=None):
        base = dict(original(controller, pred, user_block_cycles))
        audit["selection_calls"] += 1
        if not config.enabled:
            return base

        # Emission creates new mobile content. Escape is a second state-changing
        # count. Peierls and Taylor are subchannels of escape and are not added
        # again, avoiding event-rate double counting.
        event_rate = max(float(pred.mu_emit), 0.0) + max(
            float(pred.escape_per_cycle), 0.0
        )
        base_cycles = max(float(base.get("cycles", 0.0)), 0.0)
        expected_base = event_rate * base_cycles
        mode = "quiet"
        target = math.inf

        if event_rate > 0.0 and math.isfinite(event_rate):
            if expected_base > config.tau_switch_expected_events:
                mode = "tau_leap"
                target = config.tau_leap_target
            elif expected_base > config.rare_event_target:
                mode = "rare_event"
                target = config.rare_event_target

        if math.isfinite(target):
            event_limit = target / event_rate
            if event_limit < base_cycles:
                min_cycles = max(float(controller.cfg.min_block_cycles), 0.0)
                cycles = max(event_limit, min_cycles)
                base["cycles"] = cycles
                base["limiter"] = f"stochastic_{mode}"
                base["unlimited_cycles"] = event_limit
                limits = dict(base.get("candidate_limits", {}))
                limits[f"stochastic_{mode}"] = event_limit
                base["candidate_limits"] = limits
                base_cycles = cycles
                expected_base = event_rate * cycles

        if mode == "rare_event":
            audit["rare_event_calls"] += 1
        elif mode == "tau_leap":
            audit["tau_leap_calls"] += 1
        else:
            audit["quiet_calls"] += 1

        record = {
            "mode": mode,
            "event_rate_per_cycle": event_rate,
            "cycles_selected": base_cycles,
            "expected_state_events": expected_base,
            "limiter": str(base.get("limiter", "unknown")),
        }
        audit["records"].append(record)
        base["v10055_stochastic_mode"] = mode
        base["v10055_event_rate_per_cycle"] = event_rate
        base["v10055_expected_state_events"] = expected_base
        return base

    choose._v10055_hybrid_stochastic_scheduler = True
    return choose


__all__ = [
    "POINT_RELEASE",
    "MODEL_ID",
    "ENGINE_REGISTRY_V10055",
    "HybridSchedulerConfigV10055",
    "StochasticCampaignKineticMPZStateV10055",
    "StochasticPredictorMeanFieldMixinV10055",
    "engine_factory_v10055",
    "hybrid_choose_block_factory_v10055",
]
