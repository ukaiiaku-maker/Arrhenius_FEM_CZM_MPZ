"""v10.0.5.5 stochastic VHCF first-passage entry point.

This point release retains the v10.0.5.4 hazard-integrated cycle horizon and
adds three numerical capabilities:

1. reproducible stochastic finite-source emission and stochastic cleavage
   thresholds using the existing Arrhenius hazards;
2. hybrid rare-event / bounded tau-leap cycle blocks selected from expected
   state-changing event counts;
3. conservative reuse of an unchanged maximum-load FEM state in tip-only,
   load-hold fatigue calculations.

The cycle horizon is a user-selected experimental cutoff.  1e14 cycles is
supported but is not imposed as a physical target or default stopping event.
"""
from __future__ import annotations

from datetime import datetime, timezone
import inspect
import json
import math
import os
from pathlib import Path
import sys
from typing import Any

from . import fatigue_v1 as _fatigue_v1
from . import kinetic_progressive_2d_v1002 as _v1002
from . import mode_i_first_passage_v10_0_5_3_fatigue as _v10053_original
from . import mode_i_first_passage_v10_0_5_3_fatigue_audited as _v10053_audited
from . import mode_i_first_passage_v10_0_5_4_vhcf as _v10054
from .stochastic_campaign_v10055 import (
    ENGINE_REGISTRY_V10055,
    HybridSchedulerConfigV10055,
    engine_factory_v10055,
    hybrid_choose_block_factory_v10055,
)

POINT_RELEASE = "10.0.5.5"
MODEL_ID = "FEM_CZM_Mode_I_stochastic_VHCF_hybrid_v10_0_5_5"
COMPLETION_MANIFEST = "run_completion_v10_0_5_5_stochastic_vhcf.json"
STOCHASTIC_AUDIT = "stochastic_vhcf_v10_0_5_5.json"
FEM_CACHE_AUDIT = "vhcf_fem_cache_v10_0_5_5.json"


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _option_value(args: list[str], option: str, default: str | None = None):
    try:
        index = args.index(option)
    except ValueError:
        return default
    return args[index + 1] if index + 1 < len(args) else default


def _replace_unique(source: str, old: str, new: str, name: str) -> str:
    count = source.count(old)
    if count != 1:
        raise RuntimeError(
            f"v10.0.5.5 expected exactly one {name} anchor; found {count}"
        )
    return source.replace(old, new)


_CACHE_INIT_ANCHOR = """        prev_a_tip_for_block = float(a_tip)

        while step < args.steps:
"""

_CACHE_MECHANICS_ANCHOR = """                sigma_gp = np.zeros((3, mesh.ne)); psi_gp = np.zeros(mesh.ne); Ftop = 0.0
                for it in range(args.n_stagger):
                    Kmat, Rint, sigma_gp, seq_gp, s1_gp, psi_gp = assemble_mechanics(
                        mesh, u, ep_gp, rho_gp, d, D, mat, cohesive_network=cohesive_network)
                    u, Ftop = solve_dirichlet(Kmat, Rint, u, bnd, Uy_top, Uy_bot)
                    Kmat, Rint, sigma_gp, seq_gp, s1_gp, psi_gp = assemble_mechanics(
                        mesh, u, ep_gp, rho_gp, d, D, mat, cohesive_network=cohesive_network)
                    if fatigue_mode and cyclic_mechanics_enabled:
                        # In fatigue mode the plastic strain/dislocation state is
                        # advanced by the explicit cyclic mechanics block after the
                        # Kmax/J predictor chooses an accepted cycle count.  Do not
                        # also apply a monotonic max-load plastic increment here.
                        dot_ep = np.zeros(mesh.ne)
                    else:
                        ep_gp, rho_gp, dot_ep = update_plasticity(
                            ep_gp, rho_gp, sigma_gp, mat, T, dt_cur,
                            plast_model, cfg.dislocations)
"""

_CACHE_AUDIT_ANCHOR = """        tag = f"{int(T):04d}K"
"""


def patch_run_2d_source_v10055(source: str) -> str:
    """Cache only mechanically identical tip-only load-hold equilibrium states."""

    cache_init = """        prev_a_tip_for_block = float(a_tip)
        vhcf_cache_requested_v10055 = os.environ.get('ARRHENIUS_VHCF_FEM_CACHE', '1') != '0'
        vhcf_cache_enabled_v10055 = (
            vhcf_cache_requested_v10055
            and fatigue_mode
            and bool(getattr(args, 'fatigue_hold_load', False))
            and not cyclic_mechanics_enabled
            and str(getattr(args, 'bulk_plasticity_mode_v911', 'tip_only')) == 'tip_only'
        )
        vhcf_cache_v10055 = None
        vhcf_fem_solves_v10055 = 0
        vhcf_fem_reuses_v10055 = 0

        def _vhcf_cohesive_signature_v10055():
            if cohesive_network is None:
                return ()
            values = []
            for name in ('damage', 'd', 'state', 'normal_damage', 'tangential_damage'):
                value = getattr(cohesive_network, name, None)
                try:
                    arr = np.asarray(value, dtype=float)
                except Exception:
                    continue
                if arr.size > 0 and np.all(np.isfinite(arr)):
                    values.append((name, int(arr.size), float(np.sum(arr)), float(np.max(arr))))
            return tuple(values)

        def _vhcf_mechanics_key_v10055(U_value):
            return (
                float(U_value), int(mesh.nn), int(mesh.ne), float(np.sum(d)),
                float(a_tip), _vhcf_cohesive_signature_v10055(),
            )

        while step < args.steps:
"""
    source = _replace_unique(
        source, _CACHE_INIT_ANCHOR, cache_init, "FEM cache initialization"
    )

    cache_mechanics = """                cache_key_v10055 = _vhcf_mechanics_key_v10055(Uapp)
                reuse_mechanics_v10055 = (
                    vhcf_cache_enabled_v10055
                    and vhcf_cache_v10055 is not None
                    and vhcf_cache_v10055['key'] == cache_key_v10055
                )
                if reuse_mechanics_v10055:
                    u = vhcf_cache_v10055['u'].copy()
                    sigma_gp = vhcf_cache_v10055['sigma_gp'].copy()
                    seq_gp = vhcf_cache_v10055['seq_gp'].copy()
                    s1_gp = vhcf_cache_v10055['s1_gp'].copy()
                    psi_gp = vhcf_cache_v10055['psi_gp'].copy()
                    Ftop = float(vhcf_cache_v10055['Ftop'])
                    dot_ep = np.zeros(mesh.ne)
                    vhcf_fem_reuses_v10055 += 1
                else:
                    sigma_gp = np.zeros((3, mesh.ne)); psi_gp = np.zeros(mesh.ne); Ftop = 0.0
                    for it in range(args.n_stagger):
                        Kmat, Rint, sigma_gp, seq_gp, s1_gp, psi_gp = assemble_mechanics(
                            mesh, u, ep_gp, rho_gp, d, D, mat, cohesive_network=cohesive_network)
                        u, Ftop = solve_dirichlet(Kmat, Rint, u, bnd, Uy_top, Uy_bot)
                        Kmat, Rint, sigma_gp, seq_gp, s1_gp, psi_gp = assemble_mechanics(
                            mesh, u, ep_gp, rho_gp, d, D, mat, cohesive_network=cohesive_network)
                        if fatigue_mode and cyclic_mechanics_enabled:
                            dot_ep = np.zeros(mesh.ne)
                        else:
                            ep_gp, rho_gp, dot_ep = update_plasticity(
                                ep_gp, rho_gp, sigma_gp, mat, T, dt_cur,
                                plast_model, cfg.dislocations)
                    vhcf_fem_solves_v10055 += 1
                    if vhcf_cache_enabled_v10055:
                        vhcf_cache_v10055 = {
                            'key': cache_key_v10055,
                            'u': u.copy(),
                            'sigma_gp': sigma_gp.copy(),
                            'seq_gp': seq_gp.copy(),
                            's1_gp': s1_gp.copy(),
                            'psi_gp': psi_gp.copy(),
                            'Ftop': float(Ftop),
                        }
"""
    source = _replace_unique(
        source, _CACHE_MECHANICS_ANCHOR, cache_mechanics, "maximum-load mechanics"
    )

    cache_audit = """        with open(os.path.join(args.out, 'vhcf_fem_cache_v10_0_5_5.json'), 'w') as fp:
            json.dump({
                'schema': 'vhcf_fem_cache_v10_0_5_5',
                'cache_requested': bool(vhcf_cache_requested_v10055),
                'cache_enabled': bool(vhcf_cache_enabled_v10055),
                'fem_equilibrium_solves': int(vhcf_fem_solves_v10055),
                'fem_equilibrium_reuses': int(vhcf_fem_reuses_v10055),
                'bulk_mode': str(getattr(args, 'bulk_plasticity_mode_v911', 'tip_only')),
                'cyclic_mechanics_enabled': bool(cyclic_mechanics_enabled),
                'fatigue_hold_load': bool(getattr(args, 'fatigue_hold_load', False)),
                'cache_invalidates_on': [
                    'remote_displacement', 'mesh_size', 'damage_sum',
                    'crack_tip_position', 'cohesive_state_signature'
                ],
            }, fp, indent=2)

        tag = f"{int(T):04d}K"
"""
    source = _replace_unique(
        source, _CACHE_AUDIT_ANCHOR, cache_audit, "FEM cache audit output"
    )
    return source


def build_progressive_run_2d_v10055(original_run_2d):
    """Compose the cache patch ahead of the audited v10.0.5.3 lifecycle patch."""

    original_getsource = _v1002.inspect.getsource

    def patched_getsource(obj):
        source = original_getsource(obj)
        if obj is original_run_2d:
            source = patch_run_2d_source_v10055(source)
        return source

    _v1002.inspect.getsource = patched_getsource
    try:
        transformed = _ORIGINAL_AUDITED_BUILDER(original_run_2d)
    finally:
        _v1002.inspect.getsource = original_getsource

    transformed._v10055_stochastic_vhcf = True
    transformed._v10055_fem_cache = True
    transformed._v10055_constitutive_physics_changed = False
    return transformed


_ORIGINAL_AUDITED_BUILDER = _v10053_audited.build_progressive_run_2d_v10053_audited


def validate_source_transform_v10055() -> dict[str, Any]:
    saved = _v10053_audited.build_progressive_run_2d_v10053_audited
    _v10053_audited.build_progressive_run_2d_v10053_audited = (
        build_progressive_run_2d_v10055
    )
    try:
        result = _v10053_audited.validate_source_transform_v10053()
    finally:
        _v10053_audited.build_progressive_run_2d_v10053_audited = saved
    result = dict(result)
    result.update(
        {
            "point_release": POINT_RELEASE,
            "v10055_source_transform_preflight_passed": True,
            "stochastic_vhcf_adapter": True,
            "fem_cache_adapter": True,
        }
    )
    return result


def _engine_audit() -> list[dict[str, Any]]:
    rows = []
    for index, engine in enumerate(ENGINE_REGISTRY_V10055):
        state = getattr(engine, "mpz_state", None)
        diagnostics = (
            state.diagnostics_campaign()
            if state is not None and hasattr(state, "diagnostics_campaign")
            else {}
        )
        rows.append(
            {
                "engine_index": index,
                "event_statistics": str(getattr(engine, "event_statistics", "unknown")),
                "cleavage_threshold": float(getattr(engine, "B_target", 1.0)),
                "cleavage_event_index": int(
                    getattr(getattr(engine, "_threshold_stream", None), "event_index", 0)
                ),
                "predictor_mean_field_calls": int(
                    getattr(engine, "_v10055_predictor_mean_field_calls", 0)
                ),
                **diagnostics,
            }
        )
    return rows


def main(argv: list[str] | None = None):
    args = list(sys.argv[1:] if argv is None else argv)
    out_value = _option_value(args, "--out")
    if out_value is None:
        raise SystemExit("v10.0.5.5 stochastic VHCF requires --out")
    out = Path(out_value)
    out.mkdir(parents=True, exist_ok=True)
    status_path = out / COMPLETION_MANIFEST
    status = {
        "schema": "authoritative_run_completion_v10_0_5_5_stochastic_vhcf",
        "point_release": POINT_RELEASE,
        "model": MODEL_ID,
        "started_utc": _utc_now(),
        "completed_utc": None,
        "status": "running",
        "run_completed_without_exception": False,
        "constitutive_surfaces_changed": False,
    }
    status_path.write_text(json.dumps(status, indent=2))

    preflight = validate_source_transform_v10055()
    scheduler_cfg = HybridSchedulerConfigV10055.from_environment()
    scheduler_audit: dict[str, Any] = {}
    ENGINE_REGISTRY_V10055.clear()

    saved_factory = _v10053_original.engine_factory_v10053
    saved_choose = _fatigue_v1.FatigueCycleHazardController.choose_block_cycles_diagnostic
    saved_builder = _v10053_audited.build_progressive_run_2d_v10053_audited
    _v10053_original.engine_factory_v10053 = engine_factory_v10055
    _fatigue_v1.FatigueCycleHazardController.choose_block_cycles_diagnostic = (
        hybrid_choose_block_factory_v10055(
            saved_choose, scheduler_cfg, scheduler_audit
        )
    )
    _v10053_audited.build_progressive_run_2d_v10053_audited = (
        build_progressive_run_2d_v10055
    )

    try:
        results = _v10054.main(args)
        engines = _engine_audit()
        requested_statistics = os.environ.get(
            "ARRHENIUS_EVENT_STATISTICS", "stochastic"
        )
        if not engines:
            raise RuntimeError("v10.0.5.5 completed without a registered engine")
        if not any(row["predictor_mean_field_calls"] > 0 for row in engines):
            raise RuntimeError(
                "v10.0.5.5 completed without mean-field stochastic block prediction"
            )
        if requested_statistics.strip().lower() == "stochastic":
            bad = [
                row for row in engines
                if row.get("event_statistics") != "stochastic"
                or not bool(row.get("stochastic_emission_active", 0.0))
            ]
            if bad:
                raise RuntimeError(
                    "v10.0.5.5 stochastic mode was requested but not active"
                )

        cache_path = out / FEM_CACHE_AUDIT
        audit = {
            "schema": "stochastic_vhcf_v10_0_5_5",
            "point_release": POINT_RELEASE,
            "model": MODEL_ID,
            "event_statistics_requested": requested_statistics,
            "stochastic_emission_requested": os.environ.get(
                "ARRHENIUS_STOCHASTIC_EMISSION", "1"
            ) != "0",
            "stochastic_seed": int(
                os.environ.get("ARRHENIUS_STOCHASTIC_SEED", "1")
            ),
            "scheduler": scheduler_audit,
            "engines": engines,
            "source_transform": preflight,
            "fem_cache_audit": FEM_CACHE_AUDIT,
            "fem_cache_audit_exists": cache_path.exists(),
            "stochastic_scope": (
                "finite source emission and cleavage renewal; Peierls/Taylor "
                "transport remains deterministic conditional on realized source history"
            ),
            "cycle_horizon_role": "user-selected experimental censoring horizon",
            "constitutive_surfaces_changed": False,
        }
        (out / STOCHASTIC_AUDIT).write_text(
            json.dumps(audit, indent=2, default=str)
        )
        source_completion = out / _v10054.COMPLETION_MANIFEST
        source_status = (
            json.loads(source_completion.read_text())
            if source_completion.exists()
            else {}
        )
        status.update(
            {
                "completed_utc": _utc_now(),
                "status": source_status.get("status", "complete"),
                "termination": source_status.get("termination"),
                "right_censored": source_status.get("right_censored", False),
                "run_completed_without_exception": True,
                "stochastic_audit": STOCHASTIC_AUDIT,
                "source_v10054_completion": _v10054.COMPLETION_MANIFEST,
            }
        )
        status_path.write_text(json.dumps(status, indent=2, default=str))
        return results
    except BaseException as exc:
        status.update(
            {
                "completed_utc": _utc_now(),
                "status": "failed",
                "runtime_error_type": type(exc).__name__,
                "runtime_error": str(exc),
            }
        )
        status_path.write_text(json.dumps(status, indent=2, default=str))
        raise
    finally:
        _v10053_original.engine_factory_v10053 = saved_factory
        _fatigue_v1.FatigueCycleHazardController.choose_block_cycles_diagnostic = saved_choose
        _v10053_audited.build_progressive_run_2d_v10053_audited = saved_builder


if __name__ == "__main__":
    main()


__all__ = [
    "POINT_RELEASE",
    "MODEL_ID",
    "COMPLETION_MANIFEST",
    "STOCHASTIC_AUDIT",
    "FEM_CACHE_AUDIT",
    "patch_run_2d_source_v10055",
    "build_progressive_run_2d_v10055",
    "validate_source_transform_v10055",
    "main",
]
