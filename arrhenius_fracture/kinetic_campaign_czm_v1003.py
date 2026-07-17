"""v10.0.3 campaign-aware mixed-mode dispatch and engine factory.

The inherited v8/v9.11 mixed-mode mixin recognizes only the compatibility token
``moving_pz``.  A v10 campaign engine has ``state_model=kinetic_campaign_czm``;
without an explicit dispatch override, adaptive clock prediction and accepted
steps fall through to the legacy scalar N_em pathway.
"""
from __future__ import annotations

from typing import Any

from . import mixed_mode_first_passage_v9_11 as v911
from .kinetic_campaign_czm import apply_pf_manifest_to_mpz_config
from .kinetic_campaign_czm_v1001 import (
    ResetSafeDevelopedStateDiagnosticCZMFrontEngine,
)

STATE_MODEL = "kinetic_campaign_czm"


class CampaignAwareV1003TipEngineMixin(v911.CalibratedV911TipEngineMixin):
    """Dispatch v10 campaign clocks through the separated-drive API."""

    supports_progressive_kinetic_czm = True

    @staticmethod
    def _metadata(mm: dict[str, Any]) -> dict[str, Any]:
        return {
            "anisotropic_KJ_Pa_sqrt_m": mm.get("KJ"),
            "anisotropic_Kcleave_Pa_sqrt_m": mm.get("KJ", 0.0) * mm.get("fc", 1.0),
            "anisotropic_Kemit_Pa_sqrt_m": mm.get("KJ", 0.0) * mm.get("fe", 1.0),
            "anisotropic_cleavage_factor": mm.get("fc"),
            "anisotropic_emission_factor": mm.get("fe"),
            "anisotropic_reference_phase_deg": mm.get("reference_traction_phase_deg"),
            "anisotropic_candidate_angle_deg": mm.get("candidate_angle_deg"),
            "anisotropic_candidate_sigma_nn_Pa": mm.get("candidate_sigma_nn_Pa"),
            "anisotropic_probe_sigma1_Pa": mm.get("probe_sigma1_Pa"),
            "anisotropic_probe_tau_slip_abs_Pa": mm.get("slip_tau_abs_Pa"),
            "anisotropic_probe_sigma_cleave_overdrive_Pa": mm.get("candidate_overdrive_Pa"),
            "anisotropic_slip_system": mm.get("slip_system_name"),
            "anisotropic_directional_factor_cap_active": mm.get(
                "directional_factor_cap_active", False
            ),
        }

    def predict_clock_increment(self, K, T, dt):
        if getattr(self, "state_model", None) == STATE_MODEL:
            Kc, Ke, _ = self._mm_drives(K)
            return self.predict_clock_increment_drives(Kc, Ke, T, dt)
        return super().predict_clock_increment(K, T, dt)

    def step(self, K, T, dt):
        if getattr(self, "state_model", None) == STATE_MODEL:
            Kc, Ke, mm = self._mm_drives(K)
            self._mm_prev_Kcleave = Kc
            return self.step_drives(
                Kc,
                Ke,
                T,
                dt,
                metadata=self._metadata(mm),
            )
        return super().step(K, T, dt)


def engine_factory_v1003(
    original_build,
    context,
    mm,
    row,
    manifest,
    kinetic_cfg,
    registry,
):
    """Return a live v10 campaign engine builder for the v9.11 wrapper."""

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
        base.f.dN_cap = float("inf")
        base.f.N_sat = float("inf")
        base.f.recover_k = 0.0
        base.f.k_shield = 0.0
        base.f.chi_shield = 0.0
        base.f.v_emb_b3 = 0.0

        Engine = type(
            "PFEquivalentModeIKineticCZMEngineV1003",
            (
                CampaignAwareV1003TipEngineMixin,
                ResetSafeDevelopedStateDiagnosticCZMFrontEngine,
            ),
            {
                "supports_progressive_kinetic_czm": True,
                "integration_point_release": "10.0.3",
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
                "v10.0.3 factory constructed the wrong state model: "
                f"{getattr(eng, 'state_model', None)!r}"
            )
        registry.append(eng)
        return eng

    return build


__all__ = [
    "STATE_MODEL",
    "CampaignAwareV1003TipEngineMixin",
    "engine_factory_v1003",
]
