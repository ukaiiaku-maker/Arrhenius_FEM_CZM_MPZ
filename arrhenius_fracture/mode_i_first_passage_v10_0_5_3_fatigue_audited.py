"""Audited v10.0.5.3 progressive-fatigue entry point.

This module replaces the original v10.0.5.3 implementation's brittle attempt to
rewrite the Python source of ``build_progressive_run_2d_v1002`` itself.  Instead,
it applies four fail-closed loading/accounting edits to the actual ``run_2d``
source and then delegates to the unchanged, established v10.0.2 progressive
trial-CZM lifecycle builder.

No material parameter, barrier, FEM constitutive law, tensor projection,
cohesive update, MPZ transport law, source-depletion law, or topology-quality
rule is changed here.
"""
from __future__ import annotations

import inspect
from typing import Any

from . import kinetic_progressive_2d_v1002 as _v1002
from . import kinetic_progressive_2d_v1003_source as _v1003_source
from . import mode_i_first_passage_v10_0_5_3_fatigue as _original_entry
from . import sharp_front as _sharp_front

POINT_RELEASE = "10.0.5.3"
MODEL_ID = "FEM_CZM_Mode_I_progressive_cycle_block_fatigue_v10_0_5_3_audited"


def _replace_unique(source: str, old: str, new: str, name: str) -> str:
    count = source.count(old)
    if count != 1:
        raise RuntimeError(
            f"v10.0.5.3 audited adapter expected exactly one {name} anchor; "
            f"found {count}"
        )
    return source.replace(old, new)


_BACKEND_ANCHOR = (
    "        crack_backend = build_crack_backend(args, cfg.geometry)\n"
    "        cohesive_network = crack_backend.cohesive_network\n"
)

_ACCEPTED_CYCLES_ANCHOR = (
    "            fatigue_cycles_accepted = float(locals().get('fatigue_cycles_trial', 0.0)) "
    "if fatigue_mode else 0.0\n"
)

_SINGLE_FRONT_FATIGUE_ANCHOR = (
    "            else:\n"
    "                if fatigue_mode:\n"
    "                    wave_acc = FatigueWaveform(\n"
)

_DENSITY_TRANSPORT_ANCHOR = (
    "            # ---- density transport (common to both paths) ----\n"
)


def patch_run_2d_source_v10053(source: str) -> str:
    """Apply only the cyclic-loading and accepted-time plumbing edits.

    The established v10.0.2 builder still owns trial insertion, rollback,
    cohesive damage correction, checkpoint commit, remeshing, and time carry.
    """

    # The v10.0.2 builder deliberately rejects fatigue at its construction gate.
    # Temporarily hide the fatigue flag only while that gate constructs the same
    # progressive backend; restore it immediately afterward for all load/cycle logic.
    backend_wrapped = (
        "        progressive_fatigue_v10053 = bool(fatigue_mode)\n"
        "        if progressive_fatigue_v10053:\n"
        "            fatigue_mode = False\n"
        + _BACKEND_ANCHOR
        + "        if progressive_fatigue_v10053:\n"
        "            fatigue_mode = True\n"
    )
    source = _replace_unique(
        source,
        _BACKEND_ANCHOR,
        backend_wrapped,
        "backend construction",
    )

    # Convert the adaptive cycle block selected by the existing fatigue controller
    # to the physical interval consumed by the progressive kinetic-CZM lifecycle.
    accepted_replacement = (
        _ACCEPTED_CYCLES_ANCHOR
        + "            if kinetic_progressive and fatigue_mode:\n"
        "                frequency_v10053 = max(\n"
        "                    float(getattr(args, 'frequency_Hz', 1.0e3) or 1.0e3),\n"
        "                    1.0e-300,\n"
        "                )\n"
        "                dt_cur = fatigue_cycles_accepted / frequency_v10053\n"
    )
    source = _replace_unique(
        source,
        _ACCEPTED_CYCLES_ANCHOR,
        accepted_replacement,
        "accepted cycle-block conversion",
    )

    # Bypass the retired cycle_step_front state commit for the progressive engine.
    # The dynamic lifecycle target permits at most one checkpoint per equilibrium
    # state; the outer loop then re-equilibrates the changed crack geometry.
    dispatch_replacement = (
        "            else:\n"
        "                if kinetic_progressive and fatigue_mode:\n"
        "                    checkpoint_now_v10053 = float(\n"
        "                        getattr(eng, 'checkpoint_advance_total_m', 0.0)\n"
        "                    )\n"
        "                    os.environ['ARRHENIUS_COMMITTED_TARGET_EXTENSION_UM'] = (\n"
        "                        f\"{(checkpoint_now_v10053 + float(da_phys)) * 1.0e6:.17g}\"\n"
        "                    )\n"
        "                if fatigue_mode and not kinetic_progressive:\n"
        "                    wave_acc = FatigueWaveform(\n"
    )
    source = _replace_unique(
        source,
        _SINGLE_FRONT_FATIGUE_ANCHOR,
        dispatch_replacement,
        "single-front fatigue dispatch",
    )

    # If a checkpoint occurs before the requested cycle block ends, count only the
    # physical time consumed by the accepted lifecycle.  The unused interval is not
    # silently credited as fatigue cycles on the old geometry.
    accounting_replacement = (
        "            if kinetic_progressive and fatigue_mode:\n"
        "                frequency_v10053 = max(\n"
        "                    float(getattr(args, 'frequency_Hz', 1.0e3) or 1.0e3),\n"
        "                    1.0e-300,\n"
        "                )\n"
        "                requested_cycles_v10053 = float(fatigue_cycles_accepted)\n"
        "                consumed_dt_v10053 = max(\n"
        "                    float(info.get('dt_consumed_s', dt_cur)), 0.0\n"
        "                )\n"
        "                dt_cur = consumed_dt_v10053\n"
        "                fatigue_cycles_accepted = consumed_dt_v10053 * frequency_v10053\n"
        "                info['cycles_requested_v10053'] = requested_cycles_v10053\n"
        "                info['cycles'] = fatigue_cycles_accepted\n"
        "                info['cycle_limiter'] = str(\n"
        "                    locals().get('fatigue_cycle_limiter_trial', 'progressive_lifecycle')\n"
        "                )\n"
        "                info['cycle_unlimited'] = float(\n"
        "                    locals().get(\n"
        "                        'fatigue_cycle_unlimited_trial', requested_cycles_v10053\n"
        "                    )\n"
        "                )\n"
        + _DENSITY_TRANSPORT_ANCHOR
    )
    source = _replace_unique(
        source,
        _DENSITY_TRANSPORT_ANCHOR,
        accounting_replacement,
        "consumed cycle accounting",
    )
    return source


def build_progressive_run_2d_v10053_audited(original_run_2d):
    """Delegate to the unchanged v10.0.2 builder using audited run_2d source."""

    original_getsource = _v1002.inspect.getsource

    def patched_getsource(obj):
        source = original_getsource(obj)
        if obj is original_run_2d:
            source = patch_run_2d_source_v10053(source)
        return source

    _v1002.inspect.getsource = patched_getsource
    try:
        transformed = _original_entry._ORIGINAL_V1002_BUILDER(original_run_2d)
    finally:
        _v1002.inspect.getsource = original_getsource

    transformed._v10053_progressive_fatigue = True
    transformed._v10053_audited_run_2d_adapter = True
    transformed._v10053_legacy_fatigue_commit_bypassed = True
    transformed._v10053_one_checkpoint_per_outer_block = True
    transformed._v10053_consumed_cycle_accounting = True
    transformed._v10053_constitutive_physics_changed = False
    return transformed


def validate_source_transform_v10053() -> dict[str, Any]:
    """Compile the exact current run_2d/v10.0.3/v10.0.2 transform chain.

    This is a construction preflight only; it does not solve FEM equations or
    mutate a material state.  It exists specifically to prevent campaign launch
    when any source anchor or wrapper composition has drifted.
    """

    saved_builder = _v1002.build_progressive_run_2d_v1002
    _v1002.reset_progressive_runtime_v1002()
    _v1002.build_progressive_run_2d_v1002 = build_progressive_run_2d_v10053_audited
    try:
        transformed = _v1003_source.build_progressive_run_2d_v1003_source(
            _sharp_front.run_2d
        )
    finally:
        _v1002.build_progressive_run_2d_v1002 = saved_builder

    required = {
        "v1003_source_adapter": bool(
            getattr(transformed, "_v1003_source_adapter", False)
        ),
        "v1002_event_lifecycle": bool(
            getattr(transformed, "_v1002_event_lifecycle", False)
        ),
        "v10053_audited_adapter": bool(
            getattr(transformed, "_v10053_audited_run_2d_adapter", False)
        ),
        "legacy_fatigue_commit_bypassed": bool(
            getattr(transformed, "_v10053_legacy_fatigue_commit_bypassed", False)
        ),
        "consumed_cycle_accounting": bool(
            getattr(transformed, "_v10053_consumed_cycle_accounting", False)
        ),
        "constitutive_physics_changed": bool(
            getattr(transformed, "_v10053_constitutive_physics_changed", True)
        ),
    }
    failed = [
        name
        for name, value in required.items()
        if (name == "constitutive_physics_changed" and value)
        or (name != "constitutive_physics_changed" and not value)
    ]
    if failed:
        raise RuntimeError(
            "v10.0.5.3 source-transform preflight failed: " + ", ".join(failed)
        )
    return {
        "point_release": POINT_RELEASE,
        "model": MODEL_ID,
        "source_transform_preflight_passed": True,
        **required,
    }


def main(argv: list[str] | None = None):
    validate_source_transform_v10053()
    _original_entry._build_progressive_run_2d_v10053 = (
        build_progressive_run_2d_v10053_audited
    )
    return _original_entry.main(argv)


if __name__ == "__main__":
    main()


__all__ = [
    "POINT_RELEASE",
    "MODEL_ID",
    "patch_run_2d_source_v10053",
    "build_progressive_run_2d_v10053_audited",
    "validate_source_transform_v10053",
    "main",
]
