"""v10.0.5.3 engine factory: v10.0.5.2 physics under cyclic loading."""
from __future__ import annotations

import math

from . import mixed_mode_first_passage_v9_11 as v911
from .kinetic_campaign_czm import (
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

MODEL_ID = "FEM_CZM_parallel_opening_tensor_emission_fatigue_v10_0_5_3"


def engine_factory_v10053(
    original_build,
    context,
    mm,
    row,
    manifest,
    kinetic_cfg: KineticCampaignCZMConfig,
    registry,
):
    """Build the diagnostics-complete engine with a loading-only fatigue mixin."""

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
            "ParallelOpeningEmissionKineticCZMEngineV10053",
            (
                FatigueKineticsMixinV10053,
                CampaignAwareV1003TipEngineMixin,
                ChannelDiagnosticParallelOpeningEmissionCZMFrontEngine,
            ),
            {
                "supports_progressive_kinetic_czm": True,
                "supports_tensor_resolved_parallel_coupling": True,
                "supports_complete_per_channel_diagnostics": True,
                "supports_progressive_fatigue_v10053": True,
                "integration_point_release": "10.0.5.3",
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
                "v10.0.5.3 factory constructed the wrong state model: "
                f"{getattr(eng, 'state_model', None)!r}"
            )
        required = (
            "tensor_resolved_parallel_coupling",
            "per_channel_strang_diagnostics_complete",
            "supports_progressive_fatigue_v10053",
        )
        missing = [name for name in required if not bool(getattr(eng, name, False))]
        if missing:
            raise RuntimeError(
                "v10.0.5.3 engine lacks required capabilities: " + ", ".join(missing)
            )
        registry.append(eng)
        return eng

    return build


__all__ = ["MODEL_ID", "engine_factory_v10053"]
