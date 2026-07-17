"""Mode-I entry point for the PF-equivalent kinetic CZM branch.

Stage-B ``abrupt`` runs use the mature 2-D adaptive-CZM loop with the new
front-local kinetic state and exact PF v10.1.7.1 material rows.  The progressive
trial-interface mechanics are intentionally exposed through
:mod:`kinetic_cohesive_stepper`; they are not silently routed through the older
v9.17 renew-then-open controller.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys
from typing import Any

from . import mixed_mode_first_passage_v9_11 as v911
from . import mode_i_first_passage_v9_11 as modei911
from .kinetic_campaign_czm import (
    DevelopedStateDiagnosticCZMFrontEngine,
    KineticCampaignCZMConfig,
    apply_pf_manifest_to_mpz_config,
)
from .pf_equivalent_material_manifest import (
    PF_SOURCE,
    load_material_manifest,
    pf_manifest_path,
)

MODEL_ID = "FEM_CZM_Mode_I_kinetic_campaign_czm_v10_0"


def _option_value(argv: list[str], name: str) -> str | None:
    for i, token in enumerate(argv):
        if token == name and i + 1 < len(argv):
            return argv[i + 1]
        if token.startswith(name + "="):
            return token.split("=", 1)[1]
    return None


def _replace_option(argv: list[str], name: str, value: str) -> None:
    i = 0
    while i < len(argv):
        token = argv[i]
        if token == name:
            del argv[i:i + 2]
            continue
        if token.startswith(name + "="):
            del argv[i]
            continue
        i += 1
    argv.extend([name, str(value)])


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
    p.add_argument("--v10-material-source", default=PF_SOURCE)
    p.add_argument("--v10-material-class", default=None)
    p.add_argument(
        "--czm-opening-coupling",
        default="abrupt",
        choices=("abrupt", "clock_linear"),
    )
    p.add_argument("--max-action-substep", type=float, default=0.02)
    p.add_argument("--max-translation-substep-m", type=float, default=1.0e-7)
    p.add_argument("--max-internal-steps", type=int, default=20000)
    p.add_argument("--min-kinetic-substep-s", type=float, default=1.0e-15)
    return p


def _engine_factory_v10(original_build, context, mm, row, manifest, kinetic_cfg, registry):
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
            "PFEquivalentModeIKineticCZMEngine",
            (v911.CalibratedV911TipEngineMixin, DevelopedStateDiagnosticCZMFrontEngine),
            {},
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
        registry.append(eng)
        return eng

    return build


def main(argv: list[str] | None = None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    opts, remaining = parser().parse_known_args(user_args)
    material_class = (
        opts.v10_material_class
        or _option_value(remaining, "--mpz-material-class")
        or "ceramic"
    )
    if opts.v10_material_source != PF_SOURCE:
        raise SystemExit(
            "v10.0 equivalence validation currently requires --v10-material-source "
            f"{PF_SOURCE}; legacy rows remain available through the v9 entry points"
        )
    if opts.czm_opening_coupling == "clock_linear":
        raise SystemExit(
            "clock_linear is not routed through the old renew-then-open v9.17 loop. "
            "Use the kinetic_cohesive_stepper Stage-C harness until the dedicated "
            "v10 progressive 2-D loop passes its transaction tests."
        )

    manifest = load_material_manifest(material_class, parameter_source=PF_SOURCE)
    manifest_path = pf_manifest_path(material_class)
    _replace_option(remaining, "--mpz-material-manifest", str(manifest_path))
    _replace_option(remaining, "--mpz-material-class", manifest.name)
    if _option_value(remaining, "--mpz-length-um") is None:
        remaining.extend(["--mpz-length-um", "100"])
    if _option_value(remaining, "--mpz-n-bins") is None:
        remaining.extend(["--mpz-n-bins", "200"])

    kinetic_cfg = KineticCampaignCZMConfig(
        max_action_substep=opts.max_action_substep,
        max_translation_substep_m=opts.max_translation_substep_m,
        min_substep_s=opts.min_kinetic_substep_s,
        max_internal_steps=opts.max_internal_steps,
        coupling_scheme="strang",
        wake_shielding=False,
        active_shielding=True,
        signed_active_shielding=True,
        mobile_shield_fraction=1.0,
        backstress_scale=1.0,
        source_refresh_scale=1.0,
    ).validate()

    original_factory = v911._engine_factory
    engines: list[Any] = []

    def patched_factory(original_build, context, mm, row):
        return _engine_factory_v10(
            original_build, context, mm, row, manifest, kinetic_cfg, engines
        )

    v911._engine_factory = patched_factory
    try:
        results = modei911.main(remaining)
    finally:
        v911._engine_factory = original_factory

    out_value = _option_value(remaining, "--out")
    if out_value is not None:
        out = Path(out_value)
        out.mkdir(parents=True, exist_ok=True)
        payload = {
            "model": MODEL_ID,
            "front_state_model": "kinetic_campaign_czm",
            "material_parameter_source": PF_SOURCE,
            "material": manifest.as_dict(),
            "kinetic_config": vars(kinetic_cfg),
            "opening_coupling": "abrupt",
            "wake_shielding_active": False,
            "stress_channels_separated": True,
            "continuous_mpz_translation_active": True,
            "source_refresh_from_advance_only": True,
            "full_progressive_trial_loop_active": False,
            "stage": "B_abrupt_regression",
            "engine_audits": [eng.audit_payload() for eng in engines],
        }
        (out / "kinetic_campaign_czm_v10_0_audit.json").write_text(
            json.dumps(payload, indent=2, default=str)
        )
    return results


if __name__ == "__main__":
    main()
