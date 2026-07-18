"""Diagnostics-complete v10.0.5.2 parallel opening/emission engine.

This point release changes no constitutive equation or state update.  It only
preserves vector-valued per-channel diagnostics while the inherited Strang
integrator aggregates its two plastic half-steps and any internal substeps.
"""
from __future__ import annotations

import math
from typing import Any, Mapping

import numpy as np

from . import mixed_mode_first_passage_v9_11 as v911
from .kinetic_campaign_czm import (
    KineticCampaignCZMConfig,
    apply_pf_manifest_to_mpz_config,
)
from .kinetic_campaign_czm_v1003 import CampaignAwareV1003TipEngineMixin
from .kinetic_campaign_czm_v1005 import (
    STATE_MODEL,
    ParallelOpeningEmissionCZMFrontEngine,
)

MODEL_ID = "FEM_CZM_parallel_opening_tensor_emission_v10_0_5_2"

_VECTOR_SUM_KEYS = frozenset({
    "dN_emit_per_system",
})
_VECTOR_LAST_KEYS = frozenset({
    "lambda_emit_per_system_s-1",
    "emission_probability_per_system",
    "slip_system_drive_factors",
    "sigma_emission_backstress_per_system_Pa",
    "sigma_emission_effective_per_system_Pa",
    "emission_backstress_density_per_system_m2",
})


def _finite_vector(value: Any, key: str) -> np.ndarray:
    array = np.asarray(value, dtype=float).reshape(-1)
    if array.size < 1 or not np.all(np.isfinite(array)):
        raise RuntimeError(f"nonfinite or empty per-channel diagnostic {key}")
    return array


class ChannelDiagnosticParallelOpeningEmissionCZMFrontEngine(
    ParallelOpeningEmissionCZMFrontEngine
):
    """v10.0.5 mechanics/kinetics with lossless channel diagnostics."""

    state_model_detail = (
        "pf_v10_1_7_1_parallel_opening_tensor_emission_reset_safe_v10052"
    )
    per_channel_strang_diagnostics_complete = True
    per_channel_emitted_increment_semantics = "sum_over_all_strang_half_steps"
    per_channel_hazard_semantics = "last_rate_after_final_strang_half_step"

    @staticmethod
    def _sum_numeric(target: dict[str, Any], source: Mapping[str, Any]) -> None:
        """Aggregate scalars as before and retain selected finite vectors.

        The inherited v10.0.5 state update already applies each half-step.  This
        method only controls the returned diagnostic dictionary.  Emitted counts
        are additive.  Rates, stresses, probabilities, and drive factors are
        instantaneous diagnostics and therefore retain the final half-step value.
        """

        for key, value in source.items():
            if isinstance(value, (bool, np.bool_)):
                continue
            if isinstance(value, (int, float, np.integer, np.floating)):
                target[key] = float(target.get(key, 0.0)) + float(value)
                continue
            if key in _VECTOR_SUM_KEYS:
                current = _finite_vector(value, key)
                previous = target.get(key)
                if previous is None:
                    target[key] = current.tolist()
                else:
                    old = _finite_vector(previous, key)
                    if old.shape != current.shape:
                        raise RuntimeError(
                            f"per-channel diagnostic shape changed for {key}: "
                            f"{old.shape} -> {current.shape}"
                        )
                    target[key] = (old + current).tolist()
                continue
            if key in _VECTOR_LAST_KEYS:
                target[key] = _finite_vector(value, key).tolist()

    def audit_payload(self) -> dict[str, Any]:
        payload = super().audit_payload()
        payload.update({
            "model": MODEL_ID,
            "state_model_detail": self.state_model_detail,
            "per_channel_strang_diagnostics_complete": True,
            "per_channel_emitted_increment_semantics": (
                self.per_channel_emitted_increment_semantics
            ),
            "per_channel_hazard_semantics": self.per_channel_hazard_semantics,
            "constitutive_physics_changed_in_v10052": False,
        })
        return payload


def engine_factory_v10052(
    original_build,
    context,
    mm,
    row,
    manifest,
    kinetic_cfg: KineticCampaignCZMConfig,
    registry,
):
    """Build the v10.0.5 engine with diagnostics-complete aggregation."""

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
            "ParallelOpeningEmissionKineticCZMEngineV10052",
            (
                CampaignAwareV1003TipEngineMixin,
                ChannelDiagnosticParallelOpeningEmissionCZMFrontEngine,
            ),
            {
                "supports_progressive_kinetic_czm": True,
                "supports_tensor_resolved_parallel_coupling": True,
                "supports_complete_per_channel_diagnostics": True,
                "integration_point_release": "10.0.5.2",
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
        eng._mm_init(context)
        if getattr(eng, "state_model", None) != STATE_MODEL:
            raise RuntimeError(
                "v10.0.5.2 factory constructed the wrong state model: "
                f"{getattr(eng, 'state_model', None)!r}"
            )
        if not bool(getattr(eng, "tensor_resolved_parallel_coupling", False)):
            raise RuntimeError(
                "v10.0.5.2 factory constructed an engine without parallel coupling"
            )
        if not bool(getattr(eng, "per_channel_strang_diagnostics_complete", False)):
            raise RuntimeError(
                "v10.0.5.2 factory constructed an engine without complete diagnostics"
            )
        registry.append(eng)
        return eng

    return build


__all__ = [
    "MODEL_ID",
    "ChannelDiagnosticParallelOpeningEmissionCZMFrontEngine",
    "engine_factory_v10052",
]
