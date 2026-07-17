"""Consolidated, fail-closed run_2d source adapter for v10.0.3.

This adapter is intentionally limited to compatibility/accounting changes around
the established v10.0.2 lifecycle transform.  It does not alter material
parameters, rates, J evaluation, mechanics, or topology-quality rules.
"""
from __future__ import annotations

import inspect
import textwrap

from . import kinetic_progressive_2d_v1002 as _v1002


_DEFLECT = "        deflect = bool(getattr(args, 'crystal_aniso', False))\n"
_DEFLECT_NEW = (
    "        # v10.0.3: anisotropic FEM/J remains active, while this first\n"
    "        # lifecycle gate prescribes one straight Mode-I checkpoint path.\n"
    "        deflect = False\n"
)

_SUMMARY = (
    "        defl_deg = 0.0; path_dy_mm = 0.0; branched = False\n"
    "        n_primary_adv = 0; n_branch_adv = 0; n_fronts_final = 1\n"
)
_SUMMARY_NEW = (
    "        defl_deg = 0.0; path_dy_mm = 0.0; branched = False\n"
    "        # Non-deflecting runs still own a real front engine and must report\n"
    "        # committed checkpoints rather than the legacy zero initializer.\n"
    "        n_primary_adv = int(round(getattr(eng, 'n_adv', 0)))\n"
    "        n_branch_adv = 0; n_fronts_final = 1\n"
)


def _replace_exact(source: str, old: str, new: str, name: str, expected: int = 1) -> str:
    count = source.count(old)
    if count != expected:
        raise RuntimeError(
            f"v10.0.3 source adapter expected {expected} {name} anchor(s); found {count}"
        )
    return source.replace(old, new)


def _campaign_state_compatibility(source: str) -> str:
    eq = "== 'moving_pz'"
    ne = "!= 'moving_pz'"
    eq_count = source.count(eq)
    ne_count = source.count(ne)
    if eq_count < 2 or ne_count < 1:
        raise RuntimeError(
            "v10.0.3 moving-state compatibility anchors changed: "
            f"eq={eq_count}, ne={ne_count}"
        )
    source = source.replace(
        eq,
        "in ('moving_pz', 'kinetic_campaign_czm')",
    )
    source = source.replace(
        ne,
        "not in ('moving_pz', 'kinetic_campaign_czm')",
    )
    return source


def build_progressive_run_2d_v1003_source(original_run_2d):
    source = textwrap.dedent(inspect.getsource(original_run_2d))
    source = _replace_exact(source, _DEFLECT, _DEFLECT_NEW, "deflect")
    source = _replace_exact(source, _SUMMARY, _SUMMARY_NEW, "summary accounting")
    source = _campaign_state_compatibility(source)

    original_getsource = _v1002.inspect.getsource

    def patched_getsource(obj):
        if obj is original_run_2d:
            return source
        return original_getsource(obj)

    _v1002.inspect.getsource = patched_getsource
    try:
        transformed = _v1002.build_progressive_run_2d_v1002(original_run_2d)
    finally:
        _v1002.inspect.getsource = original_getsource

    transformed._v1003_source_adapter = True
    transformed._v1003_anisotropic_elasticity_preserved = True
    transformed._v1003_path_deflection_forced_off = True
    transformed._v1003_campaign_state_compatibility = True
    transformed._v1003_nondeflect_summary_accounting = True
    return transformed


__all__ = ["build_progressive_run_2d_v1003_source"]
