"""Guarded progressive single-front 2-D loop for v10 kinetic CZM.

The mature ``sharp_front.run_2d`` function is large and shared by legacy,
fatigue, branching, and mixed-mode campaigns.  Rather than copying that driver
or changing legacy behavior, this module compiles a versioned transformation of
three exact source anchors:

1. adaptive-CZM backend construction;
2. the accepted cleavage-action increment limit;
3. the single-front monotonic engine-step/post-fire geometry block.

The transformation aborts if any anchor is absent or non-unique.  Progressive
mode is restricted to one front, branching off, monotonic Mode I, and
``adaptive_czm``.  The old path is byte-for-byte preserved when progressive mode
is inactive.
"""
from __future__ import annotations

import copy
import inspect
import json
import math
import os
from pathlib import Path
import textwrap
from typing import Any, Mapping

import numpy as np

from .cohesive_trial_state import KineticTrialAdaptiveCZMBackend
from .kinetic_cohesive_stepper import (
    KineticCohesiveStepper,
    KineticCohesiveStepperConfig,
)

SCHEMA = "kinetic_campaign_czm_progressive_2d_v10_0"

_V10_PROGRESSIVE_RUNTIME: dict[str, Any] = {
    "active": False,
    "source_transform_installed": False,
    "records": [],
    "trial_insertions": 0,
    "committed_events": 0,
    "damage_rejections": 0,
    "full_rollbacks": 0,
    "anchor_counts": {},
}


def reset_progressive_runtime() -> None:
    _V10_PROGRESSIVE_RUNTIME.update({
        "active": False,
        "source_transform_installed": False,
        "records": [],
        "trial_insertions": 0,
        "committed_events": 0,
        "damage_rejections": 0,
        "full_rollbacks": 0,
        "anchor_counts": {},
    })


def _mean_or_zero(value: Any) -> float:
    try:
        array = np.asarray(value, dtype=float)
    except Exception:
        return 0.0
    return float(np.mean(array)) if array.size else 0.0


def _v10_format_progressive_info(
    engine: Any,
    step_result: Any,
    KJ: float,
    N_em_pre: float,
) -> dict[str, Any]:
    kinetics = dict(step_result.kinetics)
    channels = dict(kinetics.get("channels", {}))
    plastic = dict(kinetics.get("plastic", {}))
    advance = dict(kinetics.get("advance", {}))
    state = engine.mpz_state.diagnostics_campaign()
    generic = engine.mpz_state.diagnostics(
        engine.G,
        engine.nu,
        engine.b,
        engine.f.r0,
        float(engine.manifest.c_blunt),
    )
    rates = np.asarray(
        getattr(engine.mpz_state, "last_emission_rate_per_system_s", []),
        dtype=float,
    )
    site_capacity = np.asarray(engine.mpz_state.site_capacity, dtype=float)
    available = np.asarray(engine.mpz_state.available_sites, dtype=float)
    site_fraction = float(
        np.mean(np.divide(
            available,
            site_capacity,
            out=np.ones_like(available),
            where=site_capacity > 0.0,
        ))
    ) if site_capacity.size else 1.0
    Gc = float(kinetics.get("G_cleave_eff_eV", 0.0))
    out = {
        "fired": bool(step_result.committed),
        "n_fire": 1 if step_result.committed else 0,
        "n_fire_available": 1 if step_result.committed else 0,
        "sigma_tip": float(channels.get("sigma_cleave_eff_Pa", 0.0)),
        "sigma_emit_tip": float(channels.get("sigma_opening_tip_Pa", 0.0)),
        "sigma_back": float(channels.get("sigma_emission_backstress_Pa", 0.0)),
        "sigma_cleave_eff_Pa": float(channels.get("sigma_cleave_eff_Pa", 0.0)),
        "sigma_opening_tip_Pa": float(channels.get("sigma_opening_tip_Pa", 0.0)),
        "sigma_emission_backstress_Pa": float(
            channels.get("sigma_emission_backstress_Pa", 0.0)
        ),
        "sigma_emission_effective_Pa": float(
            channels.get("sigma_emission_effective_Pa", 0.0)
        ),
        "lambda_c": float(kinetics.get("lambda_c_effective_s-1", 0.0)),
        "lambda_c_raw": float(kinetics.get("lambda_c_raw_s-1", 0.0)),
        "lambda_e": float(np.mean(rates)) if rates.size else 0.0,
        "lambda_emit_per_system_s-1": rates.tolist(),
        "B": float(engine.B),
        "cleavage_clock_B": float(engine.B),
        "N_em": float(engine.N_em),
        "N_em_pre_renewal": float(N_em_pre),
        "N_em_shed_to_wake": 0.0,
        "r_eff": float(engine.r_eff()),
        "r_eff_m": float(engine.r_eff()),
        "W_emit": float(engine.W_emit),
        "v_crack": (
            float(kinetics.get("micro_advance_step_m", 0.0))
            / max(float(kinetics.get("dt_consumed_s", 0.0)), 1.0e-300)
        ),
        "G_cleave_raw_eV": Gc,
        "G_cleave_eff_eV": Gc,
        "S_cleave_kB": 0.0,
        "dGcleave_dsigma_eV_per_GPa": 0.0,
        "vstar_cleave_b3": 0.0,
        "cleave_barrier_kind_code": 1.0,
        "front_state_model_code": 2.0,
        "front_state_model": "kinetic_campaign_czm",
        "mpz_K_shield_Pa_sqrt_m": float(
            channels.get("K_shield_effective_Pa_sqrt_m", 0.0)
        ),
        "mpz_active_K_shield_Pa_sqrt_m": float(
            channels.get("K_shield_effective_Pa_sqrt_m", 0.0)
        ),
        "mpz_wake_K_shield_Pa_sqrt_m": 0.0,
        "mpz_mobile_count": float(state.get("mobile_count", 0.0)),
        "mpz_retained_count": float(state.get("retained_count", 0.0)),
        "mpz_available_site_fraction": site_fraction,
        "mpz_local_slip_count": float(engine.mpz_state.local_slip_count()),
        "mpz_escaped_total": float(state.get("cumulative_escaped", 0.0)),
        "mpz_recovered_total": float(state.get("cumulative_recovered", 0.0)),
        "mpz_wake_retained_total": 0.0,
        "dN_emit_block": float(plastic.get("dN_emit", 0.0)),
        "dN_store_block": float(plastic.get("dN_trapped", 0.0)),
        "dN_mobile_block": float(plastic.get("dN_emit", 0.0)),
        "dN_escape_block": float(plastic.get("dN_escaped", 0.0)),
        "dN_peierls_block": float(plastic.get("peierls_events", 0.0)),
        "dN_taylor_block": float(plastic.get("taylor_completions", 0.0)),
        "dB_block": float(kinetics.get("dB", 0.0)),
        "storage_fraction": (
            float(state.get("retained_fraction", 0.0))
        ),
        "mu_emit": float(plastic.get("dN_emit", 0.0)),
        "mu_escape": float(plastic.get("dN_escaped", 0.0)),
        "mu_cleave_pred": float(kinetics.get("dB", 0.0)),
        "anisotropic_KJ_Pa_sqrt_m": float(KJ),
        "K_open_Pa_sqrt_m": float(channels.get("K_open_Pa_sqrt_m", KJ)),
        "K_cleave_input_Pa_sqrt_m": float(
            channels.get("K_cleave_input_Pa_sqrt_m", KJ)
        ),
        "K_shield_raw_Pa_sqrt_m": float(
            channels.get("K_shield_raw_Pa_sqrt_m", 0.0)
        ),
        "K_shield_effective_Pa_sqrt_m": float(
            channels.get("K_shield_effective_Pa_sqrt_m", 0.0)
        ),
        "micro_advance_step_m": float(
            kinetics.get("micro_advance_step_m", 0.0)
        ),
        "micro_advance_total_m": float(
            kinetics.get("micro_advance_total_m", engine.micro_advance_total_m)
        ),
        "checkpoint_committed_total_m": float(
            kinetics.get(
                "checkpoint_committed_total_m",
                engine.checkpoint_advance_total_m,
            )
        ),
        "dt_consumed_s": float(step_result.dt_consumed_s),
        "dt_unused_s": float(step_result.dt_unused_s),
        "internal_substeps": int(kinetics.get("internal_substeps", 0)),
        "trial_event_id": int(step_result.trial_event_id),
        "trial_cohesive_damage": float(step_result.damage_after),
        "trial_status": "committed" if step_result.committed else "trial",
        "progressive_opening_active": True,
        "wake_shielding_mechanically_active": False,
        **state,
        **generic,
        **advance,
    }
    return out


def _require_unique(source: str, anchor: str, name: str) -> None:
    count = source.count(anchor)
    _V10_PROGRESSIVE_RUNTIME["anchor_counts"][name] = count
    if count != 1:
        raise RuntimeError(
            f"v10 progressive run_2d transform requires exactly one {name} anchor; "
            f"found {count}"
        )


def build_progressive_run_2d(original_run_2d):
    source = textwrap.dedent(inspect.getsource(original_run_2d))

    backend_anchor = """        crack_backend = build_crack_backend(args, cfg.geometry)\n        cohesive_network = crack_backend.cohesive_network\n"""
    adaptive_anchor = """            adaptive_target = min(adaptive_target, 0.8)\n"""
    step_anchor = """                else:\n                    info = eng.step(KJ, T, dt_cur)\n                if info['fired']:\n                    if Kc_first is None:\n"""
    _require_unique(source, backend_anchor, "backend_construction")
    _require_unique(source, adaptive_anchor, "adaptive_target")
    _require_unique(source, step_anchor, "single_front_step")

    backend_replacement = """        kinetic_progressive = (\n            getattr(eng, 'state_model', '') == 'kinetic_campaign_czm'\n            and os.environ.get('ARRHENIUS_CZM_OPENING_COUPLING', 'abrupt').strip().lower() == 'clock_linear'\n        )\n        if kinetic_progressive:\n            if deflect or fatigue_mode or max_fronts != 1:\n                raise RuntimeError(\n                    'v10 progressive kinetic CZM currently requires single-front monotonic Mode I with branching disabled'\n                )\n            if str(getattr(args, 'crack_backend', 'adaptive_czm')).lower() != 'adaptive_czm':\n                raise RuntimeError('v10 progressive kinetic CZM requires --crack-backend adaptive_czm')\n            crack_backend = KineticTrialAdaptiveCZMBackend(\n                geom=cfg.geometry,\n                penalty_normal_Pa_per_m=float(getattr(args, 'czm_penalty_normal', 1.0e18)),\n                penalty_tangent_Pa_per_m=float(getattr(args, 'czm_penalty_tangent', 1.0e18)),\n                max_angle_error_deg=float(getattr(args, 'czm_max_angle_error_deg', 35.0)),\n                min_area_ratio=float(getattr(args, 'czm_min_area_ratio', 0.08)),\n                min_triangle_quality=float(getattr(args, 'czm_min_triangle_quality', 0.035)),\n                max_node_move_factor=float(getattr(args, 'czm_max_node_move_factor', 1.75)),\n                max_hrefine_subsegments=int(getattr(args, 'czm_max_hrefine_subsegments', 512)),\n                opening_coupling='clock_linear',\n            )\n            kinetic_stepper = KineticCohesiveStepper(KineticCohesiveStepperConfig(\n                opening_coupling='clock_linear',\n                maximum_damage_change=float(os.environ.get('ARRHENIUS_MAX_TRIAL_DAMAGE_CHANGE', '0.05')),\n                correction_skip_threshold=float(os.environ.get('ARRHENIUS_TRIAL_CORRECTION_SKIP', '1e-4')),\n            ))\n            _V10_PROGRESSIVE_RUNTIME['active'] = True\n            _V10_PROGRESSIVE_RUNTIME['source_transform_installed'] = True\n        else:\n            crack_backend = build_crack_backend(args, cfg.geometry)\n            kinetic_stepper = None\n        cohesive_network = crack_backend.cohesive_network\n"""
    source = source.replace(backend_anchor, backend_replacement)

    adaptive_replacement = """            adaptive_target = min(adaptive_target, 0.8)\n            if kinetic_progressive:\n                adaptive_target = min(\n                    adaptive_target,\n                    0.8 * float(os.environ.get('ARRHENIUS_MAX_TRIAL_DAMAGE_CHANGE', '0.05')),\n                )\n"""
    source = source.replace(adaptive_anchor, adaptive_replacement)

    step_replacement = """                else:\n                    if kinetic_progressive:\n                        if crack_backend.active_trial(0) is None:\n                            p0_v10 = np.array([a_tip, 0.0])\n                            p1_v10 = np.array([\n                                min(a_tip + da_phys, cfg.geometry.Lx - 2e-5),\n                                0.0,\n                            ])\n                            rr_v10 = crack_backend.begin_trial_segment(\n                                mesh=mesh, boundary=bnd, damage=d, displacement=u,\n                                p0=p0_v10, p1=p1_v10, direction=np.array([1.0, 0.0]),\n                                front_id=0, front_engine=eng,\n                                bulk_history={\n                                    'ep_gp': ep_gp, 'rho_gp': rho_gp,\n                                    'pz_store_gp': pz_store_gp,\n                                    'pz_mobile_gp': pz_mobile_gp,\n                                    'pz_escape_gp': pz_escape_gp,\n                                    'pz_emit_gp': pz_emit_gp,\n                                },\n                                front_position=p0_v10,\n                                front_path=[p0_v10.copy()],\n                                kill_r=max(mesh.hbar_tip, 0.5e-6),\n                            )\n                            if not rr_v10.inserted:\n                                raise RuntimeError(\n                                    'v10 progressive trial insertion failed before consuming cleavage action: '\n                                    + str(rr_v10.reason)\n                                )\n                            parent_v10 = rr_v10.elem_parent_map\n                            mesh, bnd, d, u = (\n                                rr_v10.mesh, rr_v10.boundary,\n                                rr_v10.damage, rr_v10.displacement,\n                            )\n                            if parent_v10 is not None:\n                                pm_v10 = np.asarray(parent_v10, dtype=int)\n                                ep_gp = np.ascontiguousarray(ep_gp[:, pm_v10])\n                                rho_gp = np.ascontiguousarray(rho_gp[pm_v10])\n                                dot_ep = np.ascontiguousarray(dot_ep[pm_v10])\n                                pz_store_gp = np.ascontiguousarray(pz_store_gp[pm_v10])\n                                pz_mobile_gp = np.ascontiguousarray(pz_mobile_gp[pm_v10])\n                                pz_escape_gp = np.ascontiguousarray(pz_escape_gp[pm_v10])\n                                pz_emit_gp = np.ascontiguousarray(pz_emit_gp[pm_v10])\n                            elif mesh.ne != rho_gp.size:\n                                raise RuntimeError(\n                                    'v10 progressive trial insertion changed bulk element count without a parent map'\n                                )\n                            x = mesh.nodes[:, 0]; y = mesh.nodes[:, 1]\n                            cx_e = mesh.nodes[mesh.elems].mean(axis=1)[:, 0]\n                            cy_e = mesh.nodes[mesh.elems].mean(axis=1)[:, 1]\n                            adj = build_elem_adjacency(mesh) if rho_transport_c > 0.0 else None\n                            cohesive_network = crack_backend.cohesive_network\n                            _V10_PROGRESSIVE_RUNTIME['trial_insertions'] += 1\n\n                        trial_v10 = crack_backend.active_trial(0)\n                        trial_logs_v10 = tuple(trial_v10.log_indices)\n                        N_pre_v10 = float(eng.N_em)\n\n                        def mechanics_v10():\n                            nonlocal u, Ftop, sigma_gp, seq_gp, s1_gp, psi_gp, KJ\n                            for _it_v10 in range(max(int(args.n_stagger), 1)):\n                                Kmat_v10, Rint_v10, sigma_gp, seq_gp, s1_gp, psi_gp = assemble_mechanics(\n                                    mesh, u, ep_gp, rho_gp, d, D, mat,\n                                    cohesive_network=cohesive_network)\n                                u, Ftop = solve_dirichlet(\n                                    Kmat_v10, Rint_v10, u, bnd, Uy_top, Uy_bot)\n                                Kmat_v10, Rint_v10, sigma_gp, seq_gp, s1_gp, psi_gp = assemble_mechanics(\n                                    mesh, u, ep_gp, rho_gp, d, D, mat,\n                                    cohesive_network=cohesive_network)\n                            h_v10 = mesh.hbar_tip if mesh.hbar_tip > 0 else mesh.hbar\n                            _J_v10, K_v10, _Jinfo_v10 = compute_J_integral(\n                                mesh, u, sigma_gp, psi_gp, d,\n                                np.array([a_tip, 0.0]), np.array([1.0, 0.0]),\n                                mat, ell=max(r_J_cluster_ell, 3.0 * h_v10),\n                                crack_segments=_backend_crack_segments())\n                            KJ = max(float(K_v10), 0.0)\n                            return {\n                                'K_open_Pa_sqrt_m': KJ,\n                                'K_cleave_input_Pa_sqrt_m': KJ,\n                                'slip_system_weights': np.ones(2),\n                            }\n\n                        def mechanics_snapshot_v10():\n                            return {\n                                'u': u.copy(), 'Ftop': float(Ftop),\n                                'sigma_gp': sigma_gp.copy(),\n                                'seq_gp': seq_gp.copy(),\n                                's1_gp': s1_gp.copy(),\n                                'psi_gp': psi_gp.copy(),\n                                'KJ': float(KJ),\n                            }\n\n                        def mechanics_restore_v10(state_v10):\n                            nonlocal u, Ftop, sigma_gp, seq_gp, s1_gp, psi_gp, KJ\n                            u = state_v10['u'].copy()\n                            Ftop = float(state_v10['Ftop'])\n                            sigma_gp = state_v10['sigma_gp'].copy()\n                            seq_gp = state_v10['seq_gp'].copy()\n                            s1_gp = state_v10['s1_gp'].copy()\n                            psi_gp = state_v10['psi_gp'].copy()\n                            KJ = float(state_v10['KJ'])\n\n                        def full_rollback_v10(payload_v10):\n                            nonlocal mesh, bnd, d, u, ep_gp, rho_gp\n                            nonlocal pz_store_gp, pz_mobile_gp, pz_escape_gp, pz_emit_gp\n                            nonlocal x, y, cx_e, cy_e, adj, cohesive_network\n                            mesh = payload_v10['mesh']\n                            bnd = payload_v10['boundary']\n                            d = payload_v10['damage']\n                            u = payload_v10['displacement']\n                            bulk_v10 = payload_v10['bulk_history']\n                            ep_gp = bulk_v10['ep_gp'].copy()\n                            rho_gp = bulk_v10['rho_gp'].copy()\n                            pz_store_gp = bulk_v10['pz_store_gp'].copy()\n                            pz_mobile_gp = bulk_v10['pz_mobile_gp'].copy()\n                            pz_escape_gp = bulk_v10['pz_escape_gp'].copy()\n                            pz_emit_gp = bulk_v10['pz_emit_gp'].copy()\n                            x = mesh.nodes[:, 0]; y = mesh.nodes[:, 1]\n                            cx_e = mesh.nodes[mesh.elems].mean(axis=1)[:, 0]\n                            cy_e = mesh.nodes[mesh.elems].mean(axis=1)[:, 1]\n                            adj = build_elem_adjacency(mesh) if rho_transport_c > 0.0 else None\n                            cohesive_network = crack_backend.cohesive_network\n                            _V10_PROGRESSIVE_RUNTIME['full_rollbacks'] += 1\n\n                        result_v10 = kinetic_stepper.advance(\n                            backend=crack_backend, front_engine=eng, front_id=0,\n                            T_K=T, dt_s=dt_cur, solve_mechanics=mechanics_v10,\n                            external_snapshot=mechanics_snapshot_v10,\n                            external_restore=mechanics_restore_v10,\n                            on_full_rollback=full_rollback_v10,\n                        )\n                        if not result_v10.accepted:\n                            _V10_PROGRESSIVE_RUNTIME['damage_rejections'] += 1\n                            raise RuntimeError(\n                                'v10 progressive outer step exceeded the accepted damage increment; '\n                                'reduce the loading/time increment to approximately '\n                                + str(result_v10.recommended_dt_s) + ' s'\n                            )\n                        info = _v10_format_progressive_info(\n                            eng, result_v10, KJ, N_pre_v10)\n                        _V10_PROGRESSIVE_RUNTIME['records'].append({\n                            'temperature_K': float(T),\n                            'step': int(step),\n                            'a_tip_m_before': float(a_tip),\n                            **copy.deepcopy(info),\n                        })\n                        if result_v10.committed:\n                            _V10_PROGRESSIVE_RUNTIME['committed_events'] += 1\n                            if Kc_first is None:\n                                Kc_first = KJ; Kc_first_step = len(hist['sigma_back'])\n                            if trial_logs_v10:\n                                a_tip = max(\n                                    float(crack_backend.advance_log[i]['x1'])\n                                    for i in trial_logs_v10\n                                )\n                            x = mesh.nodes[:, 0]; y = mesh.nodes[:, 1]\n                            cx_e = mesh.nodes[mesh.elems].mean(axis=1)[:, 0]\n                            cy_e = mesh.nodes[mesh.elems].mean(axis=1)[:, 1]\n                            adj = build_elem_adjacency(mesh) if rho_transport_c > 0.0 else None\n                    else:\n                        info = eng.step(KJ, T, dt_cur)\n                if info['fired'] and not kinetic_progressive:\n                    if Kc_first is None:\n"""
    source = source.replace(step_anchor, step_replacement)

    namespace = dict(original_run_2d.__globals__)
    namespace.update({
        "KineticTrialAdaptiveCZMBackend": KineticTrialAdaptiveCZMBackend,
        "KineticCohesiveStepper": KineticCohesiveStepper,
        "KineticCohesiveStepperConfig": KineticCohesiveStepperConfig,
        "_v10_format_progressive_info": _v10_format_progressive_info,
        "_V10_PROGRESSIVE_RUNTIME": _V10_PROGRESSIVE_RUNTIME,
    })
    code = compile(source, "<kinetic_progressive_2d_v10>", "exec")
    exec(code, namespace)
    transformed = namespace[original_run_2d.__name__]
    transformed.__name__ = original_run_2d.__name__
    transformed.__qualname__ = original_run_2d.__qualname__
    transformed.__doc__ = original_run_2d.__doc__
    transformed._v10_progressive_source_transform = True
    transformed._v10_original = original_run_2d
    return transformed


def progressive_runtime_payload() -> dict[str, Any]:
    return {
        "schema": SCHEMA,
        "active": bool(_V10_PROGRESSIVE_RUNTIME["active"]),
        "source_transform_installed": bool(
            _V10_PROGRESSIVE_RUNTIME["source_transform_installed"]
        ),
        "anchor_counts": dict(_V10_PROGRESSIVE_RUNTIME["anchor_counts"]),
        "trial_insertions": int(_V10_PROGRESSIVE_RUNTIME["trial_insertions"]),
        "committed_events": int(_V10_PROGRESSIVE_RUNTIME["committed_events"]),
        "damage_rejections": int(_V10_PROGRESSIVE_RUNTIME["damage_rejections"]),
        "full_rollbacks": int(_V10_PROGRESSIVE_RUNTIME["full_rollbacks"]),
        "records": copy.deepcopy(_V10_PROGRESSIVE_RUNTIME["records"]),
        "full_progressive_trial_loop_active": bool(
            _V10_PROGRESSIVE_RUNTIME["active"]
            and _V10_PROGRESSIVE_RUNTIME["source_transform_installed"]
        ),
        "one_topology_event_per_geometry_solve": True,
        "continuous_mpz_translation_before_commit": True,
        "mpz_advance_on_commit_m": 0.0,
        "wake_shielding_active": False,
        "independent_cohesive_failure_criterion_added": False,
    }


def write_progressive_runtime_audit(out: str | Path) -> Path:
    root = Path(out)
    root.mkdir(parents=True, exist_ok=True)
    path = root / "kinetic_campaign_czm_progressive_2d_v10_0.json"
    path.write_text(json.dumps(progressive_runtime_payload(), indent=2, default=str))
    return path


__all__ = [
    "SCHEMA",
    "reset_progressive_runtime",
    "build_progressive_run_2d",
    "progressive_runtime_payload",
    "write_progressive_runtime_audit",
]
