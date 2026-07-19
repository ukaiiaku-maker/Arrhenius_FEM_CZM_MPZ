"""One-step elastic probe through the actual v10.0.5.5 production run_2d path.

This entry composes on top of the audited stochastic-VHCF source-transform stack.
It inserts a read-only recorder immediately after the production maximum-load FEM
equilibrium solve.  Plastic evolution is replaced by an elastic no-op only for
this audit entry.  The production mesh, notch stamp, boundary constraints, crack
backend, contour selection, line-of-sight policy and exclusion radius remain the
ones constructed by ``sharp_front.run_2d``.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
import sys
from typing import Any

import numpy as np

from . import mode_i_first_passage_v10_0_5_5_stochastic_vhcf_audited as _v10055
from . import plasticity as _plasticity
from .production_j_parity_v10059 import PROBE_JSON

POINT_RELEASE = "10.0.5.9"
MODEL_ID = "FEM_CZM_production_initialization_J_probe_v10_0_5_9"


_AUDIT_INSERT_ANCHOR = """                    if vhcf_cache_enabled_v10055:
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


_ORIGINAL_AUDITED_PATCH = _v10055.patch_run_2d_source_v10055_audited


def _replace_unique(source: str, old: str, new: str, name: str) -> str:
    count = source.count(old)
    if count != 1:
        raise RuntimeError(
            f"v10.0.5.9 expected exactly one {name} anchor; found {count}"
        )
    return source.replace(old, new)


def patch_run_2d_source_v10059(source: str) -> str:
    """Add a read-only post-equilibrium recorder to the transformed production source."""
    patched = _ORIGINAL_AUDITED_PATCH(source)
    injection = _AUDIT_INSERT_ANCHOR + """
                probe_path_v10059 = os.environ.get('ARRHENIUS_V10059_PROBE_PATH', '')
                if probe_path_v10059 and not os.path.exists(probe_path_v10059):
                    from arrhenius_fracture.production_j_parity_v10059 import record_production_j_probe_v10059
                    contour_text_v10059 = os.environ.get(
                        'ARRHENIUS_V10059_CONTOURS_UM', '180 240 300'
                    )
                    contours_v10059 = [
                        float(token) * 1.0e-6
                        for token in contour_text_v10059.replace(',', ' ').split()
                        if token
                    ]
                    h_local_v10059 = mesh.hbar_tip if mesh.hbar_tip > 0 else mesh.hbar
                    if deflect and fronts is not None:
                        root_v10059 = fronts[0]
                        tip_v10059 = np.asarray(root_v10059['xy'], dtype=float)
                        direction_v10059 = np.asarray(
                            root_v10059.get('fwd', np.array([1.0, 0.0])), dtype=float
                        )
                        production_ell_v10059 = max(r_J_cluster_ell, 3.0 * h_local_v10059)
                        production_segments_v10059 = []
                        production_exclude_v10059 = 2.0 * kill_r
                        production_path_v10059 = 'anisotropic_root_cluster_with_2killr_exclusion'
                    else:
                        tip_v10059 = np.array([a_tip, 0.0], dtype=float)
                        direction_v10059 = np.array([1.0, 0.0], dtype=float)
                        production_ell_v10059 = max(r_J_cluster_ell, 3.0 * h_local_v10059)
                        production_segments_v10059 = _backend_crack_segments()
                        production_exclude_v10059 = 0.0
                        production_path_v10059 = 'straight_cluster_no_exclusion'
                    record_production_j_probe_v10059(
                        path=probe_path_v10059,
                        mesh=mesh,
                        u=u,
                        ep_gp=ep_gp,
                        sigma_gp=sigma_gp,
                        psi_tension_gp=psi_gp,
                        d=d,
                        D=D,
                        Kmat=Kmat,
                        mat=mat,
                        Ftop_N_per_thickness=Ftop,
                        Uapp_m=Uapp,
                        tip_xy=tip_v10059,
                        direction=direction_v10059,
                        half_thickness_m=half_h,
                        kill_r_m=kill_r,
                        production_ell_m=production_ell_v10059,
                        production_segments=production_segments_v10059,
                        production_exclude_radius_m=production_exclude_v10059,
                        production_path=production_path_v10059,
                        contour_outer_m=contours_v10059,
                        specimen_width_m=cfg.geometry.Lx,
                        specimen_height_m=cfg.geometry.Ly,
                        requested_tip_h_m=float(getattr(args, 'tip_h_fine', 0.0) or 0.0),
                        crack_backend_name=crack_backend.name,
                        crystal_anisotropic=bool(getattr(args, 'crystal_aniso', False)),
                        crystal_theta_deg=float(getattr(args, 'crystal_theta_deg', 0.0) or 0.0),
                    )
"""
    return _replace_unique(
        patched,
        _AUDIT_INSERT_ANCHOR,
        injection,
        "post-equilibrium production recorder",
    )


def validate_source_transform_v10059() -> dict[str, Any]:
    source = __import__("inspect").getsource(
        __import__("arrhenius_fracture.sharp_front", fromlist=["run_2d"]).run_2d
    )
    patched = patch_run_2d_source_v10059(source)
    compile(patched, "<v10.0.5.9-production-j-probe>", "exec")
    required = {
        "production_recorder": "record_production_j_probe_v10059(" in patched,
        "production_exclusion": "production_exclude_v10059 = 2.0 * kill_r" in patched,
        "full_audited_v10055_stack": "cohesive_elements" in patched,
    }
    failed = [name for name, value in required.items() if not value]
    if failed:
        raise RuntimeError(
            "v10.0.5.9 source-transform preflight failed: " + ", ".join(failed)
        )
    return {
        "point_release": POINT_RELEASE,
        "model": MODEL_ID,
        "source_transform_preflight_passed": True,
        "constitutive_physics_changed": False,
        **required,
    }


def _elastic_update_plasticity(
    ep_gp,
    rho_gp,
    sigma_gp,
    mat,
    T,
    dt,
    plast_model,
    dislocation_cfg,
    *args,
    **kwargs,
):
    """Audit-only elastic no-op with the production update function's return contract."""
    del sigma_gp, mat, T, dt, plast_model, dislocation_cfg, args
    ep_out = np.asarray(ep_gp, dtype=float).copy()
    rho_out = np.asarray(rho_gp, dtype=float).copy()
    dot = np.zeros_like(rho_out)
    if bool(kwargs.get("return_info", False)):
        return ep_out, rho_out, dot, {
            "dWp_accepted_gp": np.zeros_like(rho_out),
            "audit_elastic_noop": True,
        }
    return ep_out, rho_out, dot


def _option_value(args: list[str], option: str, default: str | None = None):
    try:
        index = args.index(option)
    except ValueError:
        return default
    return args[index + 1] if index + 1 < len(args) else default


def _has_option(args: list[str], option: str) -> bool:
    return any(token == option or token.startswith(option + "=") for token in args)


def _ensure_v911_probe_contract(args: list[str]) -> list[str]:
    """Satisfy the complete v9.11 production-entry preflight without altering physics.

    The v9.11 Mode-I path requires 2-D anisotropic crystal competition with one
    non-branching front.  The campaign runner already supplies mode=2d and
    max-fronts=1.  This audit entry adds the required competition switch so the
    recorded state follows the same anisotropic root-front production semantics.
    """
    resolved = list(args)
    if _has_option(resolved, "--no-crystal-aniso"):
        raise SystemExit(
            "v10.0.5.9 production J parity requires --crystal-aniso to match "
            "the v10.0.5.8 anisotropic fixed-grip reference"
        )
    if not _has_option(resolved, "--crystal-aniso"):
        resolved.append("--crystal-aniso")
    if not _has_option(resolved, "--crystal-compete"):
        resolved.append("--crystal-compete")
    if _has_option(resolved, "--crystal-branch"):
        raise SystemExit(
            "v10.0.5.9 production J parity requires branching disabled"
        )
    max_fronts = _option_value(resolved, "--max-fronts", "1")
    if int(max_fronts or "1") != 1:
        raise SystemExit(
            "v10.0.5.9 production J parity requires --max-fronts 1"
        )
    return resolved


def main(argv: list[str] | None = None):
    args = _ensure_v911_probe_contract(
        list(sys.argv[1:] if argv is None else argv)
    )
    out_value = _option_value(args, "--out")
    if out_value is None:
        raise SystemExit("v10.0.5.9 production J probe requires --out")
    out = Path(out_value).resolve()
    out.mkdir(parents=True, exist_ok=True)
    probe_path = out / PROBE_JSON
    if probe_path.exists():
        probe_path.unlink()

    validate_source_transform_v10059()
    os.environ["ARRHENIUS_V10059_PROBE_PATH"] = str(probe_path)
    os.environ.setdefault("ARRHENIUS_V10059_CONTOURS_UM", "180 240 300")
    os.environ.setdefault("ARRHENIUS_EVENT_STATISTICS", "mean_field")
    os.environ.setdefault("ARRHENIUS_STOCHASTIC_EMISSION", "0")
    os.environ.setdefault("ARRHENIUS_VHCF_FEM_CACHE", "0")

    saved_patch = _v10055.patch_run_2d_source_v10055_audited
    saved_update = _plasticity.update_plasticity
    _v10055.patch_run_2d_source_v10055_audited = patch_run_2d_source_v10059
    _plasticity.update_plasticity = _elastic_update_plasticity
    try:
        result = _v10055.main(args)
    finally:
        _v10055.patch_run_2d_source_v10055_audited = saved_patch
        _plasticity.update_plasticity = saved_update
        os.environ.pop("ARRHENIUS_V10059_PROBE_PATH", None)

    if not probe_path.exists():
        raise RuntimeError("production path completed without writing the v10.0.5.9 probe")
    payload = json.loads(probe_path.read_text())
    payload["source_transform"] = validate_source_transform_v10059()
    payload["v911_probe_contract"] = {
        "mode_2d": _option_value(args, "--mode") == "2d",
        "crystal_anisotropic": _has_option(args, "--crystal-aniso"),
        "crystal_competition": _has_option(args, "--crystal-compete"),
        "crystal_branching": _has_option(args, "--crystal-branch"),
        "max_fronts": int(_option_value(args, "--max-fronts", "1") or "1"),
    }
    payload["base_run_returned"] = True
    probe_path.write_text(json.dumps(payload, indent=2, default=str))
    print(probe_path)
    return result


if __name__ == "__main__":
    main()


__all__ = [
    "POINT_RELEASE",
    "MODEL_ID",
    "patch_run_2d_source_v10059",
    "validate_source_transform_v10059",
    "_ensure_v911_probe_contract",
    "main",
]
