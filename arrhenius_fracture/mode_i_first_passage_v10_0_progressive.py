"""Progressive Mode-I PF-equivalent kinetic CZM entry point.

This module combines the v10 campaign front state and guarded progressive
single-front loop with the established production geometry safeguards.  It does
not import the v9.17/v9.18 renew-then-open hazard controller.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
import sys
from typing import Any

from . import crack_backend as _cb
from . import fem as _fem
from . import mesh as _mesh
from . import mixed_mode_first_passage_v9_11 as v911
from . import mode_i_first_passage_v9_11 as modei911
from . import mode_i_first_passage_v9_18_3 as _v9183
from . import mode_i_first_passage_v9_18_5 as _v9185
from . import mode_i_first_passage_v9_18_5_6 as _v91856
from . import sharp_front as _sharp_front
from .kinetic_campaign_czm import KineticCampaignCZMConfig
from .kinetic_progressive_2d_v10 import (
    build_progressive_run_2d,
    progressive_runtime_payload,
    reset_progressive_runtime,
    write_progressive_runtime_audit,
)
from .mode_i_first_passage_v10_0 import (
    MODEL_ID,
    _engine_factory_v10,
    _option_value,
    _replace_option,
    parser,
)
from .pf_equivalent_material_manifest import (
    PF_SOURCE,
    load_material_manifest,
    pf_manifest_path,
)

PROGRESSIVE_MODEL_ID = MODEL_ID + "_progressive_clock_linear"


def _set_env(name: str, value: str, saved: dict[str, str | None]) -> None:
    if name not in saved:
        saved[name] = os.environ.get(name)
    os.environ[name] = str(value)


def _restore_env(saved: dict[str, str | None]) -> None:
    for name, value in saved.items():
        if value is None:
            os.environ.pop(name, None)
        else:
            os.environ[name] = value


def main(argv: list[str] | None = None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    opts, remaining = parser().parse_known_args(user_args)
    if opts.czm_opening_coupling != "clock_linear":
        raise SystemExit(
            "mode_i_first_passage_v10_0_progressive requires "
            "--czm-opening-coupling clock_linear"
        )
    if opts.v10_material_source != PF_SOURCE:
        raise SystemExit(
            "progressive equivalence validation requires PF v10.1.7.1 material rows"
        )

    material_class = (
        opts.v10_material_class
        or _option_value(remaining, "--mpz-material-class")
        or "ceramic"
    )
    manifest = load_material_manifest(material_class, parameter_source=PF_SOURCE)
    _replace_option(
        remaining,
        "--mpz-material-manifest",
        str(pf_manifest_path(material_class)),
    )
    _replace_option(remaining, "--mpz-material-class", manifest.name)
    _replace_option(remaining, "--crack-backend", "adaptive_czm")
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
    original_run_2d = _sharp_front.run_2d
    original_make_mesh = _mesh.make_tri_mesh
    original_assemble = _fem.assemble_mechanics
    original_solve = _fem.solve_dirichlet
    original_insert = _cb.AdaptiveCZMBackend._insert_target_in_incident_triangle
    original_advance = _cb.AdaptiveCZMBackend.advance
    engines: list[Any] = []
    env_saved: dict[str, str | None] = {}
    error: BaseException | None = None

    def patched_factory(original_build, context, mm, row):
        return _engine_factory_v10(
            original_build,
            context,
            mm,
            row,
            manifest,
            kinetic_cfg,
            engines,
        )

    reset_progressive_runtime()
    _v9185._RUNTIME["mesh"] = None
    _v9185._RUNTIME["corridor_centers"] = []
    _v9185._RUNTIME["component_anchor_history"] = []
    _v9185._RUNTIME["quality_vetoes"] = []
    _v91856._AUDIT["accepted_events"] = []
    _v91856._AUDIT["resolution_warnings"] = []
    _v91856._AUDIT["quality_vetoes"] = []
    _v91856._AUDIT["consecutive_veto_abort"] = None

    _v9185._make_corridor_mesh._original = original_make_mesh
    _v9185._assemble_with_mesh_capture._original = original_assemble
    _v91856._strict_quality_advance_v91856._original = original_advance

    _set_env("ARRHENIUS_CZM_OPENING_COUPLING", "clock_linear", env_saved)
    _set_env(
        "ARRHENIUS_MAX_TRIAL_DAMAGE_CHANGE",
        os.environ.get("ARRHENIUS_MAX_TRIAL_DAMAGE_CHANGE", "0.05"),
        env_saved,
    )
    target = _option_value(remaining, "--target-crack-extension-um")
    if target is not None and "ARRHENIUS_COMMITTED_TARGET_EXTENSION_UM" not in os.environ:
        _set_env("ARRHENIUS_COMMITTED_TARGET_EXTENSION_UM", target, env_saved)

    v911._engine_factory = patched_factory
    _mesh.make_tri_mesh = _v9185._make_corridor_mesh
    _fem.assemble_mechanics = _v9185._assemble_with_mesh_capture
    _fem.solve_dirichlet = _v9185._component_anchored_solve
    _cb.AdaptiveCZMBackend._insert_target_in_incident_triangle = (
        _v9183._edge_aware_insert_target_in_incident_triangle
    )
    _cb.AdaptiveCZMBackend.advance = _v91856._strict_quality_advance_v91856
    _sharp_front.run_2d = build_progressive_run_2d(original_run_2d)

    results = None
    try:
        results = modei911.main(remaining)
    except BaseException as exc:
        error = exc
        raise
    finally:
        v911._engine_factory = original_factory
        _sharp_front.run_2d = original_run_2d
        _mesh.make_tri_mesh = original_make_mesh
        _fem.assemble_mechanics = original_assemble
        _fem.solve_dirichlet = original_solve
        _cb.AdaptiveCZMBackend._insert_target_in_incident_triangle = original_insert
        _cb.AdaptiveCZMBackend.advance = original_advance
        _restore_env(env_saved)
        _v91856._write_audit(remaining, error)

    out_value = _option_value(remaining, "--out")
    if out_value is not None:
        out = Path(out_value)
        out.mkdir(parents=True, exist_ok=True)
        progressive_path = write_progressive_runtime_audit(out)
        progressive = progressive_runtime_payload()
        payload = {
            "model": PROGRESSIVE_MODEL_ID,
            "front_state_model": "kinetic_campaign_czm",
            "material_parameter_source": PF_SOURCE,
            "material": manifest.as_dict(),
            "kinetic_config": vars(kinetic_cfg),
            "opening_coupling": "clock_linear",
            "wake_shielding_active": False,
            "stress_channels_separated": True,
            "continuous_mpz_translation_active": True,
            "source_refresh_from_advance_only": True,
            "full_progressive_trial_loop_active": bool(
                progressive["full_progressive_trial_loop_active"]
            ),
            "progressive_runtime_audit": progressive_path.name,
            "trial_insertions": progressive["trial_insertions"],
            "committed_events": progressive["committed_events"],
            "damage_rejections": progressive["damage_rejections"],
            "full_rollbacks": progressive["full_rollbacks"],
            "geometry_quality_wrapper": "v9.18.5.6_explicit_quality_chain",
            "prefined_mode_i_corridor_enabled": True,
            "component_wise_incremental_x_anchor_enabled": True,
            "engine_audits": [eng.audit_payload() for eng in engines],
        }
        (out / "kinetic_campaign_czm_v10_0_audit.json").write_text(
            json.dumps(payload, indent=2, default=str)
        )
        if not payload["full_progressive_trial_loop_active"]:
            raise RuntimeError(
                "progressive run returned without activating the dedicated trial loop"
            )
    return results


if __name__ == "__main__":
    main()
