"""v10.0.2 progressive 2-D event lifecycle with retry and time carry.

This guarded source transform retains the v10.0 trial-cohesive topology and
v9.18.5.6 geometry-quality chain while closing two lifecycle gaps:

* a rejected trial-damage increment is retried transactionally at reduced dt;
* unused physical time after a committed checkpoint is carried at the same
  applied load into a newly inserted and re-equilibrated trial segment.
"""
from __future__ import annotations

import copy
import inspect
import json
import os
from pathlib import Path
import textwrap
from typing import Any

import numpy as np

from .cohesive_trial_state import KineticTrialAdaptiveCZMBackend
from .kinetic_cohesive_stepper import (
    KineticCohesiveStepper,
    KineticCohesiveStepperConfig,
)
from .kinetic_event_lifecycle_v1002 import (
    EventLifecycleConfig,
    KineticEventLifecycleController,
)
from .kinetic_progressive_2d_v10 import _v10_format_progressive_info

SCHEMA = "kinetic_campaign_czm_progressive_2d_v10_0_2"

_RUNTIME: dict[str, Any] = {
    "active": False,
    "source_transform_installed": False,
    "records": [],
    "trial_insertions": 0,
    "committed_events": 0,
    "damage_rejections": 0,
    "accepted_substeps": 0,
    "full_rollbacks": 0,
    "carried_time_s": 0.0,
    "max_commits_in_outer_interval": 0,
    "anchor_counts": {},
}


def reset_progressive_runtime_v1002() -> None:
    _RUNTIME.update({
        "active": False,
        "source_transform_installed": False,
        "records": [],
        "trial_insertions": 0,
        "committed_events": 0,
        "damage_rejections": 0,
        "accepted_substeps": 0,
        "full_rollbacks": 0,
        "carried_time_s": 0.0,
        "max_commits_in_outer_interval": 0,
        "anchor_counts": {},
    })


def _require_unique(source: str, anchor: str, name: str) -> None:
    count = source.count(anchor)
    _RUNTIME["anchor_counts"][name] = count
    if count != 1:
        raise RuntimeError(
            f"v10.0.2 progressive run_2d transform requires exactly one "
            f"{name} anchor; found {count}"
        )


_INCREMENT_KEYS = (
    "dN_emit_block",
    "dN_store_block",
    "dN_mobile_block",
    "dN_escape_block",
    "dN_peierls_block",
    "dN_taylor_block",
    "dB_block",
    "mu_emit",
    "mu_escape",
    "mu_cleave_pred",
    "micro_advance_step_m",
)


def _merge_progressive_infos(
    infos: list[dict[str, Any]],
    lifecycle: Any,
) -> dict[str, Any]:
    if not infos:
        raise RuntimeError(
            "progressive event lifecycle returned no accepted kinetic substep"
        )
    out = copy.deepcopy(infos[-1])
    for key in _INCREMENT_KEYS:
        out[key] = float(sum(float(row.get(key, 0.0)) for row in infos))
    out["fired"] = int(lifecycle.committed_events) > 0
    out["n_fire"] = int(lifecycle.committed_events)
    out["n_fire_available"] = int(lifecycle.committed_events)
    out["dt_consumed_s"] = float(lifecycle.consumed_dt_s)
    out["dt_unused_s"] = float(lifecycle.unused_dt_s)
    out["event_lifecycle_rejected_attempts"] = int(
        lifecycle.rejected_attempts
    )
    out["event_lifecycle_accepted_substeps"] = len(infos)
    out["event_lifecycle_committed_events"] = int(
        lifecycle.committed_events
    )
    out["event_lifecycle_stopped_at_target"] = bool(
        lifecycle.stopped_at_target
    )
    out["v_crack"] = (
        float(out["micro_advance_step_m"])
        / max(float(lifecycle.consumed_dt_s), 1.0e-300)
    )
    return out


def build_progressive_run_2d_v1002(original_run_2d):
    source = textwrap.dedent(inspect.getsource(original_run_2d))

    backend_anchor = """        crack_backend = build_crack_backend(args, cfg.geometry)\n        cohesive_network = crack_backend.cohesive_network\n"""
    adaptive_anchor = """            adaptive_target = min(adaptive_target, 0.8)\n"""
    step_anchor = """                else:\n                    info = eng.step(KJ, T, dt_cur)\n                if info['fired']:\n                    if Kc_first is None:\n"""
    _require_unique(source, backend_anchor, "backend_construction")
    _require_unique(source, adaptive_anchor, "adaptive_target")
    _require_unique(source, step_anchor, "single_front_step")

    backend_replacement = """        kinetic_progressive = (\n            getattr(eng, 'state_model', '') == 'kinetic_campaign_czm'\n            and os.environ.get('ARRHENIUS_CZM_OPENING_COUPLING', 'abrupt').strip().lower() == 'clock_linear'\n        )\n        if kinetic_progressive:\n            if deflect or fatigue_mode or max_fronts != 1:\n                raise RuntimeError(\n                    'v10.0.2 progressive kinetic CZM requires single-front monotonic Mode I with branching disabled'\n                )\n            if str(getattr(args, 'crack_backend', 'adaptive_czm')).lower() != 'adaptive_czm':\n                raise RuntimeError('v10.0.2 progressive kinetic CZM requires --crack-backend adaptive_czm')\n            crack_backend = KineticTrialAdaptiveCZMBackend(\n                geom=cfg.geometry,\n                penalty_normal_Pa_per_m=float(getattr(args, 'czm_penalty_normal', 1.0e18)),\n                penalty_tangent_Pa_per_m=float(getattr(args, 'czm_penalty_tangent', 1.0e18)),\n                max_angle_error_deg=float(getattr(args, 'czm_max_angle_error_deg', 35.0)),\n                min_area_ratio=float(getattr(args, 'czm_min_area_ratio', 0.08)),\n                min_triangle_quality=float(getattr(args, 'czm_min_triangle_quality', 0.035)),\n                max_node_move_factor=float(getattr(args, 'czm_max_node_move_factor', 1.75)),\n                max_hrefine_subsegments=int(getattr(args, 'czm_max_hrefine_subsegments', 512)),\n                opening_coupling='clock_linear',\n            )\n            kinetic_stepper = KineticCohesiveStepper(KineticCohesiveStepperConfig(\n                opening_coupling='clock_linear',\n                maximum_damage_change=float(os.environ.get('ARRHENIUS_MAX_TRIAL_DAMAGE_CHANGE', '0.05')),\n                correction_skip_threshold=float(os.environ.get('ARRHENIUS_TRIAL_CORRECTION_SKIP', '1e-4')),\n            ))\n            kinetic_lifecycle = KineticEventLifecycleController(EventLifecycleConfig(\n                min_retry_dt_s=float(os.environ.get('ARRHENIUS_MIN_TRIAL_RETRY_DT_S', '1e-18')),\n                max_retries_per_substep=int(os.environ.get('ARRHENIUS_MAX_TRIAL_RETRIES', '64')),\n                max_accepted_substeps_per_interval=int(os.environ.get('ARRHENIUS_MAX_ACCEPTED_SUBSTEPS_PER_INTERVAL', '10000')),\n            ))\n            _V1002_RUNTIME['active'] = True\n            _V1002_RUNTIME['source_transform_installed'] = True\n        else:\n            crack_backend = build_crack_backend(args, cfg.geometry)\n            kinetic_stepper = None\n            kinetic_lifecycle = None\n        cohesive_network = crack_backend.cohesive_network\n"""
    source = source.replace(backend_anchor, backend_replacement)

    adaptive_replacement = """            adaptive_target = min(adaptive_target, 0.8)\n            if kinetic_progressive:\n                adaptive_target = min(\n                    adaptive_target,\n                    0.8 * float(os.environ.get('ARRHENIUS_MAX_TRIAL_DAMAGE_CHANGE', '0.05')),\n                )\n"""
    source = source.replace(adaptive_anchor, adaptive_replacement)

    step_replacement = """                else:\n                    if kinetic_progressive:\n                        interval_infos_v1002 = []\n                        target_um_v1002 = float(os.environ.get(\n                            'ARRHENIUS_COMMITTED_TARGET_EXTENSION_UM', 'inf'))\n                        target_m_v1002 = target_um_v1002 * 1.0e-6\n\n                        def target_reached_v1002():\n                            return (\n                                np.isfinite(target_m_v1002)\n                                and eng.checkpoint_advance_total_m\n                                >= target_m_v1002 - max(1.0e-12 * da_phys, 1.0e-15)\n                            )\n\n                        def ensure_trial_v1002():\n                            nonlocal mesh, bnd, d, u, ep_gp, rho_gp, dot_ep\n                            nonlocal pz_store_gp, pz_mobile_gp, pz_escape_gp, pz_emit_gp\n                            nonlocal x, y, cx_e, cy_e, adj, cohesive_network\n                            if crack_backend.active_trial(0) is None:\n                                p0_v1002 = np.array([a_tip, 0.0])\n                                p1_v1002 = np.array([\n                                    min(a_tip + da_phys, cfg.geometry.Lx - 2e-5),\n                                    0.0,\n                                ])\n                                if float(np.linalg.norm(p1_v1002 - p0_v1002)) <= 1.0e-15:\n                                    raise RuntimeError(\n                                        'v10.0.2 cannot insert another physical checkpoint before the domain boundary'\n                                    )\n                                rr_v1002 = crack_backend.begin_trial_segment(\n                                    mesh=mesh, boundary=bnd, damage=d, displacement=u,\n                                    p0=p0_v1002, p1=p1_v1002,\n                                    direction=np.array([1.0, 0.0]),\n                                    front_id=0, front_engine=eng,\n                                    bulk_history={\n                                        'ep_gp': ep_gp, 'rho_gp': rho_gp,\n                                        'dot_ep': dot_ep,\n                                        'pz_store_gp': pz_store_gp,\n                                        'pz_mobile_gp': pz_mobile_gp,\n                                        'pz_escape_gp': pz_escape_gp,\n                                        'pz_emit_gp': pz_emit_gp,\n                                    },\n                                    front_position=p0_v1002,\n                                    front_path=[p0_v1002.copy()],\n                                    kill_r=max(mesh.hbar_tip, 0.5e-6),\n                                )\n                                if not rr_v1002.inserted:\n                                    raise RuntimeError(\n                                        'v10.0.2 trial insertion failed before consuming cleavage action: '\n                                        + str(rr_v1002.reason)\n                                    )\n                                parent_v1002 = rr_v1002.elem_parent_map\n                                mesh, bnd, d, u = (\n                                    rr_v1002.mesh, rr_v1002.boundary,\n                                    rr_v1002.damage, rr_v1002.displacement,\n                                )\n                                if parent_v1002 is not None:\n                                    pm_v1002 = np.asarray(parent_v1002, dtype=int)\n                                    ep_gp = np.ascontiguousarray(ep_gp[:, pm_v1002])\n                                    rho_gp = np.ascontiguousarray(rho_gp[pm_v1002])\n                                    dot_ep = np.ascontiguousarray(dot_ep[pm_v1002])\n                                    pz_store_gp = np.ascontiguousarray(pz_store_gp[pm_v1002])\n                                    pz_mobile_gp = np.ascontiguousarray(pz_mobile_gp[pm_v1002])\n                                    pz_escape_gp = np.ascontiguousarray(pz_escape_gp[pm_v1002])\n                                    pz_emit_gp = np.ascontiguousarray(pz_emit_gp[pm_v1002])\n                                elif mesh.ne != rho_gp.size:\n                                    raise RuntimeError(\n                                        'v10.0.2 trial insertion changed bulk element count without a parent map'\n                                    )\n                                x = mesh.nodes[:, 0]; y = mesh.nodes[:, 1]\n                                cx_e = mesh.nodes[mesh.elems].mean(axis=1)[:, 0]\n                                cy_e = mesh.nodes[mesh.elems].mean(axis=1)[:, 1]\n                                adj = build_elem_adjacency(mesh) if rho_transport_c > 0.0 else None\n                                cohesive_network = crack_backend.cohesive_network\n                                _V1002_RUNTIME['trial_insertions'] += 1\n\n                            trial_v1002 = crack_backend.active_trial(0)\n                            return {\n                                'trial_logs': tuple(trial_v1002.log_indices),\n                                'trial_event_id': int(trial_v1002.event_index),\n                                'N_em_pre': float(eng.N_em),\n                                'a_tip_m_before': float(a_tip),\n                            }\n\n                        def mechanics_v1002():\n                            nonlocal u, Ftop, sigma_gp, seq_gp, s1_gp, psi_gp, KJ\n                            for _it_v1002 in range(max(int(args.n_stagger), 1)):\n                                Kmat_v1002, Rint_v1002, sigma_gp, seq_gp, s1_gp, psi_gp = assemble_mechanics(\n                                    mesh, u, ep_gp, rho_gp, d, D, mat,\n                                    cohesive_network=cohesive_network)\n                                u, Ftop = solve_dirichlet(\n                                    Kmat_v1002, Rint_v1002, u, bnd, Uy_top, Uy_bot)\n                                Kmat_v1002, Rint_v1002, sigma_gp, seq_gp, s1_gp, psi_gp = assemble_mechanics(\n                                    mesh, u, ep_gp, rho_gp, d, D, mat,\n                                    cohesive_network=cohesive_network)\n                            h_v1002 = mesh.hbar_tip if mesh.hbar_tip > 0 else mesh.hbar\n                            _J_v1002, K_v1002, _Jinfo_v1002 = compute_J_integral(\n                                mesh, u, sigma_gp, psi_gp, d,\n                                np.array([a_tip, 0.0]), np.array([1.0, 0.0]),\n                                mat, ell=max(r_J_cluster_ell, 3.0 * h_v1002),\n                                crack_segments=_backend_crack_segments())\n                            KJ = max(float(K_v1002), 0.0)\n                            return {\n                                'K_open_Pa_sqrt_m': KJ,\n                                'K_cleave_input_Pa_sqrt_m': KJ,\n                                'slip_system_weights': np.ones(2),\n                            }\n\n                        def mechanics_snapshot_v1002():\n                            return {\n                                'u': u.copy(), 'Ftop': float(Ftop),\n                                'sigma_gp': sigma_gp.copy(),\n                                'seq_gp': seq_gp.copy(),\n                                's1_gp': s1_gp.copy(),\n                                'psi_gp': psi_gp.copy(),\n                                'KJ': float(KJ),\n                            }\n\n                        def mechanics_restore_v1002(state_v1002):\n                            nonlocal u, Ftop, sigma_gp, seq_gp, s1_gp, psi_gp, KJ\n                            u = state_v1002['u'].copy()\n                            Ftop = float(state_v1002['Ftop'])\n                            sigma_gp = state_v1002['sigma_gp'].copy()\n                            seq_gp = state_v1002['seq_gp'].copy()\n                            s1_gp = state_v1002['s1_gp'].copy()\n                            psi_gp = state_v1002['psi_gp'].copy()\n                            KJ = float(state_v1002['KJ'])\n\n                        def full_rollback_v1002(payload_v1002):\n                            nonlocal mesh, bnd, d, u, ep_gp, rho_gp, dot_ep\n                            nonlocal pz_store_gp, pz_mobile_gp, pz_escape_gp, pz_emit_gp\n                            nonlocal x, y, cx_e, cy_e, adj, cohesive_network\n                            mesh = payload_v1002['mesh']\n                            bnd = payload_v1002['boundary']\n                            d = payload_v1002['damage']\n                            u = payload_v1002['displacement']\n                            bulk_v1002 = payload_v1002['bulk_history']\n                            ep_gp = bulk_v1002['ep_gp'].copy()\n                            rho_gp = bulk_v1002['rho_gp'].copy()\n                            dot_ep = bulk_v1002['dot_ep'].copy()\n                            pz_store_gp = bulk_v1002['pz_store_gp'].copy()\n                            pz_mobile_gp = bulk_v1002['pz_mobile_gp'].copy()\n                            pz_escape_gp = bulk_v1002['pz_escape_gp'].copy()\n                            pz_emit_gp = bulk_v1002['pz_emit_gp'].copy()\n                            x = mesh.nodes[:, 0]; y = mesh.nodes[:, 1]\n                            cx_e = mesh.nodes[mesh.elems].mean(axis=1)[:, 0]\n                            cy_e = mesh.nodes[mesh.elems].mean(axis=1)[:, 1]\n                            adj = build_elem_adjacency(mesh) if rho_transport_c > 0.0 else None\n                            cohesive_network = crack_backend.cohesive_network\n                            _V1002_RUNTIME['full_rollbacks'] += 1\n\n                        def advance_trial_v1002(dt_trial_v1002):\n                            return kinetic_stepper.advance(\n                                backend=crack_backend, front_engine=eng, front_id=0,\n                                T_K=T, dt_s=dt_trial_v1002,\n                                solve_mechanics=mechanics_v1002,\n                                external_snapshot=mechanics_snapshot_v1002,\n                                external_restore=mechanics_restore_v1002,\n                                on_full_rollback=full_rollback_v1002,\n                            )\n\n                        def on_accepted_v1002(context_v1002, result_v1002,\n                                              requested_v1002, retries_v1002):\n                            info_v1002 = _v10_format_progressive_info(\n                                eng, result_v1002, KJ, context_v1002['N_em_pre'])\n                            interval_infos_v1002.append(info_v1002)\n                            _V1002_RUNTIME['accepted_substeps'] += 1\n                            _V1002_RUNTIME['records'].append({\n                                'temperature_K': float(T),\n                                'step': int(step),\n                                'a_tip_m_before': context_v1002['a_tip_m_before'],\n                                'trial_requested_dt_s': float(requested_v1002),\n                                'retry_count': int(retries_v1002),\n                                'carry_sequence_index': len(interval_infos_v1002) - 1,\n                                **copy.deepcopy(info_v1002),\n                            })\n\n                        def on_commit_v1002(context_v1002, result_v1002):\n                            nonlocal a_tip, Kc_first, Kc_first_step\n                            nonlocal x, y, cx_e, cy_e, adj\n                            _V1002_RUNTIME['committed_events'] += 1\n                            if Kc_first is None:\n                                Kc_first = KJ\n                                Kc_first_step = len(hist['sigma_back'])\n                            trial_logs_v1002 = context_v1002['trial_logs']\n                            if trial_logs_v1002:\n                                a_tip = max(\n                                    float(crack_backend.advance_log[i]['x1'])\n                                    for i in trial_logs_v1002\n                                )\n                            x = mesh.nodes[:, 0]; y = mesh.nodes[:, 1]\n                            cx_e = mesh.nodes[mesh.elems].mean(axis=1)[:, 0]\n                            cy_e = mesh.nodes[mesh.elems].mean(axis=1)[:, 1]\n                            adj = build_elem_adjacency(mesh) if rho_transport_c > 0.0 else None\n\n                        lifecycle_v1002 = kinetic_lifecycle.consume_interval(\n                            total_dt_s=dt_cur,\n                            ensure_trial=ensure_trial_v1002,\n                            advance_trial=advance_trial_v1002,\n                            on_accepted=on_accepted_v1002,\n                            on_commit=on_commit_v1002,\n                            target_reached=target_reached_v1002,\n                        )\n                        _V1002_RUNTIME['damage_rejections'] += int(\n                            lifecycle_v1002.rejected_attempts)\n                        _V1002_RUNTIME['max_commits_in_outer_interval'] = max(\n                            int(_V1002_RUNTIME['max_commits_in_outer_interval']),\n                            int(lifecycle_v1002.committed_events),\n                        )\n                        for accepted_v1002 in lifecycle_v1002.accepted_steps:\n                            if bool(accepted_v1002.result.committed):\n                                _V1002_RUNTIME['carried_time_s'] += max(\n                                    float(accepted_v1002.requested_dt_s)\n                                    - float(accepted_v1002.result.dt_consumed_s),\n                                    0.0,\n                                )\n                        info = _merge_progressive_infos(\n                            interval_infos_v1002, lifecycle_v1002)\n                    else:\n                        info = eng.step(KJ, T, dt_cur)\n                if info['fired'] and not kinetic_progressive:\n                    if Kc_first is None:\n"""
    source = source.replace(step_anchor, step_replacement)

    namespace = dict(original_run_2d.__globals__)
    namespace.update({
        "KineticTrialAdaptiveCZMBackend": KineticTrialAdaptiveCZMBackend,
        "KineticCohesiveStepper": KineticCohesiveStepper,
        "KineticCohesiveStepperConfig": KineticCohesiveStepperConfig,
        "EventLifecycleConfig": EventLifecycleConfig,
        "KineticEventLifecycleController": KineticEventLifecycleController,
        "_v10_format_progressive_info": _v10_format_progressive_info,
        "_merge_progressive_infos": _merge_progressive_infos,
        "_V1002_RUNTIME": _RUNTIME,
    })
    code = compile(source, "<kinetic_progressive_2d_v1002>", "exec")
    exec(code, namespace)
    transformed = namespace[original_run_2d.__name__]
    transformed.__name__ = original_run_2d.__name__
    transformed.__qualname__ = original_run_2d.__qualname__
    transformed.__doc__ = original_run_2d.__doc__
    transformed._v10_progressive_source_transform = True
    transformed._v1002_event_lifecycle = True
    transformed._v10_original = original_run_2d
    return transformed


def progressive_runtime_payload_v1002() -> dict[str, Any]:
    return {
        "schema": SCHEMA,
        "active": bool(_RUNTIME["active"]),
        "source_transform_installed": bool(
            _RUNTIME["source_transform_installed"]
        ),
        "anchor_counts": dict(_RUNTIME["anchor_counts"]),
        "trial_insertions": int(_RUNTIME["trial_insertions"]),
        "committed_events": int(_RUNTIME["committed_events"]),
        "damage_rejections": int(_RUNTIME["damage_rejections"]),
        "accepted_substeps": int(_RUNTIME["accepted_substeps"]),
        "full_rollbacks": int(_RUNTIME["full_rollbacks"]),
        "carried_time_s": float(_RUNTIME["carried_time_s"]),
        "max_commits_in_outer_interval": int(
            _RUNTIME["max_commits_in_outer_interval"]
        ),
        "records": copy.deepcopy(_RUNTIME["records"]),
        "full_progressive_trial_loop_active": bool(
            _RUNTIME["active"] and _RUNTIME["source_transform_installed"]
        ),
        "rejected_step_retry_active": True,
        "unused_time_carry_active": True,
        "same_load_re_equilibration_after_commit": True,
        "one_topology_event_per_equilibrium_state": True,
        "continuous_mpz_translation_before_commit": True,
        "mpz_advance_on_commit_m": 0.0,
        "wake_shielding_active": False,
        "independent_cohesive_failure_criterion_added": False,
    }


def write_progressive_runtime_audit_v1002(out: str | Path) -> Path:
    root = Path(out)
    root.mkdir(parents=True, exist_ok=True)
    path = root / "kinetic_campaign_czm_progressive_2d_v10_0_2.json"
    path.write_text(
        json.dumps(progressive_runtime_payload_v1002(), indent=2, default=str)
    )
    return path


__all__ = [
    "SCHEMA",
    "reset_progressive_runtime_v1002",
    "build_progressive_run_2d_v1002",
    "progressive_runtime_payload_v1002",
    "write_progressive_runtime_audit_v1002",
]
