"""v10.0.5.18.3.3 four-class production entry with target-aware corridor mesh."""
from __future__ import annotations

import json
import os
from pathlib import Path
import sys

from . import mode_i_first_passage_v9_18_5_3 as _v91853
from . import mode_i_first_passage_v10_0_5_18_3_2_four_class_joint_K_single_trial_stochastic_emission as _base
from .long_growth_corridor_v10051833 import (
    CORRIDOR_SCHEMA,
    MODEL_ID as CORRIDOR_MODEL,
    target_aware_long_growth_corridor_mesh,
)

POINT_RELEASE = "10.0.5.18.3.3"
MODEL_ID = "FEM_CZM_four_class_joint_K_single_trial_target_aware_corridor_v10_0_5_18_3_3"
PRODUCTION_MANIFEST = "persistent_site_production_manifest_v10_0_5_18_3_3.json"
SELECTION_MANIFEST = "persistent_site_parameter_selection_v10_0_5_18_3_3.json"
TRANSFER_MANIFEST = "four_class_parameter_transfer_v10_0_5_18_3_3.json"
SEED_MANIFEST = "stochastic_seed_manifest_v10_0_5_18_3_3.json"

_OLD_FILES = {
    "persistent_site_production_manifest_v10_0_5_18_3_2.json": PRODUCTION_MANIFEST,
    "persistent_site_parameter_selection_v10_0_5_18_3_2.json": SELECTION_MANIFEST,
    "four_class_parameter_transfer_v10_0_5_18_3_2.json": TRANSFER_MANIFEST,
    "stochastic_seed_manifest_v10_0_5_18_3_2.json": SEED_MANIFEST,
}


def _option_value(argv: list[str], name: str) -> str | None:
    for index, token in enumerate(argv):
        if token == name and index + 1 < len(argv):
            return argv[index + 1]
        if token.startswith(name + "="):
            return token.split("=", 1)[1]
    return None


def _required_positive_option(argv: list[str], name: str) -> float:
    raw = _option_value(argv, name)
    if raw is None:
        raise SystemExit(f"v10.0.5.18.3.3 requires {name}")
    value = float(raw)
    if value <= 0.0:
        raise SystemExit(f"v10.0.5.18.3.3 requires positive {name}")
    return value


def _out_path(argv: list[str]) -> Path | None:
    raw = _option_value(argv, "--out")
    return None if raw is None else Path(raw).expanduser().resolve()


def _rewrite_outputs(out: Path | None, target_um: float, da_um: float, lpz_um: float) -> None:
    if out is None:
        return
    corridor_path = out / "compact_corridor_mesh_v91852.json"
    corridor = {}
    if corridor_path.is_file():
        try:
            corridor = json.loads(corridor_path.read_text())
        except Exception:
            corridor = {}

    fields = {
        "target_aware_long_growth_corridor": True,
        "target_aware_corridor_schema": CORRIDOR_SCHEMA,
        "target_aware_corridor_model": CORRIDOR_MODEL,
        "committed_target_extension_um_propagated_before_mesh": target_um,
        "physical_event_length_um_propagated_before_mesh": da_um,
        "process_zone_length_um_propagated_before_mesh": lpz_um,
        "corridor_acceptance_metric": "triangle_quality_and_h_tip_over_L_pz",
        "tip_h_over_da_role": "audit_warning_only",
        "child_area_ratio_floor_relaxed": False,
        "constitutive_physics_changed": False,
    }

    for old_name, new_name in _OLD_FILES.items():
        old = out / old_name
        if not old.is_file():
            continue
        payload = json.loads(old.read_text())
        payload["point_release"] = POINT_RELEASE
        payload["model"] = MODEL_ID
        if old_name.startswith("persistent_site_production_manifest"):
            payload["schema"] = "persistent_site_production_manifest_v10_0_5_18_3_3"
            physics = dict(payload.get("physics_contract", {}) or {})
            physics.update(fields)
            payload["physics_contract"] = physics
            payload["long_growth_corridor"] = corridor
        elif old_name.startswith("persistent_site_parameter_selection"):
            payload["schema"] = MODEL_ID
            policy = dict(payload.get("policy", {}) or {})
            policy.update(fields)
            payload["policy"] = policy
        else:
            payload.update(fields)
        new = out / new_name
        new.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n")
        old.unlink()


def main(argv: list[str] | None = None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    target_um = _required_positive_option(user_args, "--target-crack-extension-um")
    da_um = _required_positive_option(user_args, "--da-phys")
    lpz_raw = _option_value(user_args, "--mpz-length-um")
    lpz_um = 50.0 if lpz_raw is None else float(lpz_raw)
    if lpz_um <= 0.0:
        raise SystemExit("v10.0.5.18.3.3 requires positive --mpz-length-um")
    out = _out_path(user_args)

    requested_env = {
        "ARRHENIUS_COMMITTED_TARGET_EXTENSION_UM": str(target_um),
        "ARRHENIUS_PHYSICAL_DA_UM": str(da_um),
        "ARRHENIUS_CORRIDOR_PROCESS_ZONE_UM": str(lpz_um),
        "ARRHENIUS_CORRIDOR_MAX_CENTER_GAP_UM": os.environ.get(
            "ARRHENIUS_CORRIDOR_MAX_CENTER_GAP_UM", "100"
        ),
        "ARRHENIUS_MAX_CORRIDOR_H_OVER_LPZ": os.environ.get(
            "ARRHENIUS_MAX_CORRIDOR_H_OVER_LPZ", "0.25"
        ),
    }
    saved_env = {key: os.environ.get(key) for key in requested_env}
    saved_corridor = _v91853._quality_selected_corridor_mesh
    target_aware_long_growth_corridor_mesh._original = (
        saved_corridor._original
        if hasattr(saved_corridor, "_original")
        else None
    )
    # v9.18.5.3 assigns this function to the lower mesh slot and then attaches
    # the raw mesh constructor as ``_original`` at runtime.
    _v91853._quality_selected_corridor_mesh = target_aware_long_growth_corridor_mesh
    for key, value in requested_env.items():
        os.environ[key] = value

    try:
        return _base.main(user_args)
    finally:
        _v91853._quality_selected_corridor_mesh = saved_corridor
        for key, value in saved_env.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value
        _rewrite_outputs(out, target_um, da_um, lpz_um)


if __name__ == "__main__":
    main()


__all__ = [
    "MODEL_ID",
    "POINT_RELEASE",
    "PRODUCTION_MANIFEST",
    "SEED_MANIFEST",
    "SELECTION_MANIFEST",
    "TRANSFER_MANIFEST",
    "main",
]
