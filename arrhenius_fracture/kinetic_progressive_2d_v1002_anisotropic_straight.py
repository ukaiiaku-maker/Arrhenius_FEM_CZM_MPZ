"""Anisotropic-elastic, straight-path adapter for the v10.0.2 lifecycle.

The inherited v9.11 direct Mode-I entry point requires ``--crystal-aniso`` and
``--crystal-compete`` because its calibrated integration keeps cubic anisotropic
elasticity active in FEM equilibrium and J evaluation.  In the legacy
``sharp_front.run_2d`` implementation, however, ``crystal_aniso`` also sets the
unrelated ``deflect`` path-selection switch.

This adapter changes only that source-level coupling for the guarded v10.0.2
progressive function: anisotropic elasticity has already been constructed from
``args.crystal_aniso`` before the ``deflect`` assignment, while the crack path is
forced to the prescribed straight single-front Mode-I checkpoint.  No material
parameter, stress channel, J definition, or topology-quality rule is changed.
"""
from __future__ import annotations

import copy
import inspect
import json
from pathlib import Path
import textwrap
from typing import Any

from . import kinetic_progressive_2d_v1002 as _base

SCHEMA = "kinetic_campaign_czm_progressive_2d_v10_0_2_anisotropic_straight"

_DEFLECT_ANCHOR = "        deflect = bool(getattr(args, 'crystal_aniso', False))\n"
_DEFLECT_REPLACEMENT = (
    "        # v10.0.2 progressive gate: keep cubic anisotropic FEM/J active,\n"
    "        # but prescribe one straight Mode-I checkpoint path.\n"
    "        deflect = False\n"
)


def build_progressive_run_2d_v1002_anisotropic_straight(original_run_2d):
    source = textwrap.dedent(inspect.getsource(original_run_2d))
    count = source.count(_DEFLECT_ANCHOR)
    if count != 1:
        raise RuntimeError(
            "v10.0.2 anisotropic-straight adapter requires exactly one "
            f"deflect anchor; found {count}"
        )
    modified_source = source.replace(_DEFLECT_ANCHOR, _DEFLECT_REPLACEMENT)

    original_getsource = _base.inspect.getsource

    def patched_getsource(obj):
        if obj is original_run_2d:
            return modified_source
        return original_getsource(obj)

    _base.inspect.getsource = patched_getsource
    try:
        transformed = _base.build_progressive_run_2d_v1002(original_run_2d)
    finally:
        _base.inspect.getsource = original_getsource

    transformed._v1002_anisotropic_elasticity_preserved = True
    transformed._v1002_path_deflection_forced_off = True
    transformed._v1002_crystal_compete_used_for_v911_validation_only = True
    return transformed


def reset_progressive_runtime_v1002_anisotropic_straight() -> None:
    _base.reset_progressive_runtime_v1002()


def progressive_runtime_payload_v1002_anisotropic_straight() -> dict[str, Any]:
    payload = copy.deepcopy(_base.progressive_runtime_payload_v1002())
    payload.update({
        "schema": SCHEMA,
        "anisotropic_elasticity_preserved": True,
        "anisotropic_J_preserved": True,
        "path_deflection_forced_off": True,
        "straight_single_front_mode_I_checkpoint": True,
        "crystal_compete_role": "v9.11_direct_mode_validation_contract_only",
        "anisotropic_path_selection_active": False,
    })
    return payload


def write_progressive_runtime_audit_v1002_anisotropic_straight(
    out: str | Path,
) -> Path:
    root = Path(out)
    root.mkdir(parents=True, exist_ok=True)
    path = root / "kinetic_campaign_czm_progressive_2d_v10_0_2.json"
    path.write_text(json.dumps(
        progressive_runtime_payload_v1002_anisotropic_straight(),
        indent=2,
        default=str,
    ))
    return path


__all__ = [
    "SCHEMA",
    "build_progressive_run_2d_v1002_anisotropic_straight",
    "reset_progressive_runtime_v1002_anisotropic_straight",
    "progressive_runtime_payload_v1002_anisotropic_straight",
    "write_progressive_runtime_audit_v1002_anisotropic_straight",
]
