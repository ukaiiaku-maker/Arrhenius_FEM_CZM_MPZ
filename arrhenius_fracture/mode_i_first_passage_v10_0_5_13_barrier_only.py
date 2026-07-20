"""v10.0.5.13 full 2-D FEM/CZM barrier-only monotonic entry.

This entry preserves the existing v9.18.5.6 adaptive-CZM/quality lifecycle,
the state-coupled ``bulk_same_pt_km`` continuum update, and the productionized
330 um refinement support.  Only the opening, emission, Peierls, and Taylor
barrier surfaces are selected from the v9.11.1 response registry.

Candidate-specific source inventory, source refresh, encounter/retention,
recovery, shielding, blunting, MPZ-grid, and developed-state fields are ignored.
Those controls are one common 2-D solver policy for all four options.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import sys
from typing import Any, Mapping

from . import bulk_state_v911 as _bulk
from . import crack_backend as _crack_backend
from . import mesh as _mesh
from . import mixed_mode_first_passage_v9_11 as _v911
from . import mode_i_first_passage_v9_18_5 as _v9185
from . import mode_i_first_passage_v9_18_5_6 as _solver
from .barrier_only_response_registry_v100513 import (
    PARAMETER_SOURCE,
    TWO_D_STATE_POLICY,
    load_barrier_option,
)
from .mode_i_first_passage_v10_0 import _option_value, _replace_option
from .mode_i_first_passage_v10_0_5_12_phase_c import (
    CLUSTER_J_RADIUS_TO_LEGACY_ELL,
)
from .mode_i_first_passage_v10_0_5_12_3_phase_c import _annotate_mesh
from .moving_process_zone import build_mpz_config_from_namespace
from .physical_refinement_mesh_v100510 import (
    clear_physical_refinement_v100510,
    configure_physical_refinement_v100510,
    make_physical_refinement_mesh_v100510,
)

POINT_RELEASE = "10.0.5.13"
MODEL_ID = "FEM_CZM_full_2D_barrier_only_state_coupled_v10_0_5_13"
PRODUCTION_MANIFEST = "barrier_only_production_manifest_v10_0_5_13.json"


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(dict(payload), indent=2, default=str))


def _wrapper_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
    parser.add_argument("--barrier-option", required=True)
    parser.add_argument("--barrier-registry", type=Path, default=None)
    parser.add_argument("--tip-refinement-radius-um", type=float, default=330.0)
    parser.add_argument("--selected-cluster-J-outer-um", type=float, default=240.0)
    parser.add_argument("--local-J-outer-um", type=float, default=100.0)
    return parser


def _float_option(argv: list[str], name: str) -> float | None:
    value = _option_value(argv, name)
    return None if value is None else float(value)


def _require_close_or_set(
    argv: list[str],
    name: str,
    expected: float,
    *,
    rel_tol: float = 1.0e-12,
    abs_tol: float = 1.0e-15,
) -> None:
    existing = _float_option(argv, name)
    if existing is not None and not math.isclose(
        existing, expected, rel_tol=rel_tol, abs_tol=abs_tol
    ):
        raise SystemExit(
            f"{name}={existing:.16g} conflicts with v10.0.5.13 value {expected:.16g}"
        )
    _replace_option(argv, name, f"{expected:.16g}")


def _require_text_or_set(argv: list[str], name: str, expected: str) -> None:
    existing = _option_value(argv, name)
    if existing is not None and str(existing) != str(expected):
        raise SystemExit(
            f"{name}={existing!r} conflicts with v10.0.5.13 value {expected!r}"
        )
    _replace_option(argv, name, str(expected))


def _require_int_or_set(argv: list[str], name: str, expected: int) -> None:
    existing = _option_value(argv, name)
    if existing is not None and int(existing) != int(expected):
        raise SystemExit(
            f"{name}={existing!r} conflicts with v10.0.5.13 value {expected}"
        )
    _replace_option(argv, name, str(int(expected)))


def cluster_j_legacy_length_m(physical_outer_radius_um: float) -> float:
    radius_um = float(physical_outer_radius_um)
    if not math.isfinite(radius_um) or radius_um <= 0.0:
        raise ValueError("cluster-J physical outer radius must be finite and positive")
    return radius_um * 1.0e-6 / CLUSTER_J_RADIUS_TO_LEGACY_ELL


def _apply_barrier_pt_config(cfg: Any, row: Mapping[str, Any]) -> None:
    """Install only Peierls/Taylor barrier kinetics plus common 2-D controls."""
    emit0 = max(float(row["emit_G00_eV"]), 1.0e-30)
    values = {
        "bulk_kinetics_model": "emission_derived_peierls_taylor_multihit",
        "bulk_kinetics_model_detail": "v100513_barriers_only_common_2d_state",
        "pt_emit_G00_eV": float(row["emit_G00_eV"]),
        "pt_emit_gT_eV_per_K": float(row["emit_gT_eV_per_K"]),
        "pt_emit_sigc0_Pa": float(row["emit_sigc0_GPa"]) * 1.0e9,
        "pt_emit_sT_Pa_per_K": float(row["emit_sT_GPa_per_K"]) * 1.0e9,
        "pt_emit_Tref_K": float(row["Tref_K"]),
        "pt_emit_exp_a": float(row["emit_exp_a"]),
        "pt_emit_exp_n": float(row["emit_exp_n"]),
        "pt_emit_floor_frac": float(row["emit_floor_frac"]),
        "pt_emit_floor_min_eV": 1.0e-4,
        "pt_emit_floor_max_frac": 0.95,
        "pt_peierls_energy_ratio": float(row["peierls_H0_eV"]) / emit0,
        "pt_peierls_entropy_ratio": float(row["peierls_activation_entropy_kB"]),
        "pt_peierls_activation_entropy_kB": float(
            row["peierls_activation_entropy_kB"]
        ),
        "pt_peierls_exp_a": float(row["peierls_exp_a"]),
        "pt_peierls_exp_n": float(row["peierls_exp_n"]),
        "pt_peierls_stress_ratio": 1.0,
        "pt_peierls_stress_fraction": float(
            TWO_D_STATE_POLICY["peierls_stress_fraction"]
        ),
        "pt_peierls_nu0_s": float(row["peierls_nu0_s"]),
        "pt_taylor_energy_ratio": float(row["taylor_H0_eV"]) / emit0,
        "pt_taylor_entropy_ratio": float(row["taylor_activation_entropy_kB"]),
        "pt_taylor_activation_entropy_kB": float(
            row["taylor_activation_entropy_kB"]
        ),
        "pt_taylor_exp_a": float(row["taylor_exp_a"]),
        "pt_taylor_exp_n": float(row["taylor_exp_n"]),
        "pt_taylor_stress_ratio": 1.0,
        "pt_taylor_stress_fraction": float(
            TWO_D_STATE_POLICY["taylor_stress_fraction"]
        ),
        "pt_taylor_nu0_s": float(row["taylor_nu0_s"]),
        "pt_taylor_corr_rho_c": float(
            TWO_D_STATE_POLICY["taylor_corr_rho_c_m2"]
        ),
        "pt_taylor_renewal_time_s": float(
            TWO_D_STATE_POLICY["taylor_renewal_time_s"]
        ),
        "pt_taylor_m_exponent": float(TWO_D_STATE_POLICY["taylor_m_exponent"]),
        "pt_taylor_m_scale": float(TWO_D_STATE_POLICY["taylor_m_scale"]),
        "pt_taylor_m_cap": float(TWO_D_STATE_POLICY["taylor_m_cap"]),
        "pt_encounter_efficiency": float(
            TWO_D_STATE_POLICY["encounter_efficiency"]
        ),
        "pt_forest_density_floor_m2": float(
            TWO_D_STATE_POLICY["forest_density_floor_m2"]
        ),
    }
    for name, value in values.items():
        setattr(cfg, name, value)


def _build_barrier_only_mpz_config(args: Any, row: Mapping[str, Any]):
    length_m = float(TWO_D_STATE_POLICY["mpz_length_um"]) * 1.0e-6
    cfg = build_mpz_config_from_namespace(args, default_length_m=length_m)
    cfg.length_m = length_m
    cfg.n_bins = int(TWO_D_STATE_POLICY["mpz_n_bins"])
    cfg.n_systems = 2

    # Common solver-native state policy; no candidate-specific 1-D state field.
    cfg.source_sites_per_system = float(
        TWO_D_STATE_POLICY["source_sites_per_system"]
    )
    cfg.source_recovery_rate_s = float(
        TWO_D_STATE_POLICY["source_recovery_rate_s"]
    )
    cfg.source_refresh_length_m = float(
        TWO_D_STATE_POLICY["source_refresh_length_um"]
    ) * 1.0e-6
    cfg.source_bin_count = int(TWO_D_STATE_POLICY["source_bin_count"])
    cfg.shielding_orientation_factors = tuple(
        TWO_D_STATE_POLICY["shielding_orientation_factors"]
    )
    cfg.mobile_shield_fraction = float(
        TWO_D_STATE_POLICY["mobile_shield_fraction"]
    )
    cfg.shielding_core_m = float(TWO_D_STATE_POLICY["shielding_core_m"])
    cfg.retained_recovery_nu0_s = float(
        TWO_D_STATE_POLICY["retained_recovery_nu0_s"]
    )
    cfg.retained_recovery_barrier_eV = float(
        TWO_D_STATE_POLICY["retained_recovery_barrier_eV"]
    )
    cfg.retained_recovery_activation_volume_b3 = float(
        TWO_D_STATE_POLICY["retained_recovery_activation_volume_b3"]
    )
    cfg.mobile_recovery_rate_s = float(
        TWO_D_STATE_POLICY["mobile_recovery_rate_s"]
    )
    cfg.pair_annihilation_rate_per_count_s = float(
        TWO_D_STATE_POLICY["pair_annihilation_rate_per_count_s"]
    )
    cfg.blunting_length_m = float(TWO_D_STATE_POLICY["blunting_length_um"]) * 1.0e-6
    cfg.blunting_slip_fraction = float(
        TWO_D_STATE_POLICY["blunting_slip_fraction"]
    )
    _apply_barrier_pt_config(cfg, row)
    return cfg


def _set_bulk_barrier_defaults(args: Any, row: Mapping[str, Any], bulk_mode: str) -> None:
    _apply_barrier_pt_config(args, row)
    values = {
        "front_state_model": "moving_pz",
        "pz_store_to_rho_scale": 0.0,
        "tip_source_rho_per_emit": 0.0,
        "bulk_mult_frac": 0.0,
        "exhaustion": False,
        "bulk_plasticity_mode_v911": str(bulk_mode),
    }
    for name, value in values.items():
        setattr(args, name, value)


def _compact_audit(row: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "candidate_id": row["candidate_id"],
        "target_class": row["target_class"],
        "parameter_source": PARAMETER_SOURCE,
        "barrier_fingerprint_sha256": row["barrier_fingerprint_sha256"],
        "parameter_fingerprint_sha256": row["parameter_fingerprint_sha256"],
        "barrier_only_transfer": True,
        "barrier_fields_transferred": row["barrier_fields_transferred"],
        "candidate_state_fields_ignored": row["candidate_state_fields_ignored"],
        "two_d_state_policy": TWO_D_STATE_POLICY,
        "bulk_state_evolves_in_fem": True,
        "shielding_derived_from_evolving_state": True,
    }


def _mesh_payload(mesh: Any, requested_radius_m: float) -> dict[str, Any]:
    actual = getattr(mesh, "production_refinement_radius_m", None) if mesh is not None else None
    centers = getattr(mesh, "production_refinement_centers_m", None) if mesh is not None else None
    policy = getattr(mesh, "production_refinement_policy", None) if mesh is not None else None
    actual_value = None if actual is None else float(actual)
    verified = actual_value is not None and math.isclose(
        actual_value, requested_radius_m, rel_tol=1.0e-12, abs_tol=1.0e-15
    )
    return {
        "requested_radius_m": requested_radius_m,
        "requested_radius_um": requested_radius_m * 1.0e6,
        "actual_radius_m": actual_value,
        "actual_radius_um": None if actual_value is None else actual_value * 1.0e6,
        "policy": policy,
        "centers_m": centers,
        "actual_radius_verified": verified,
    }


def _tip_centers(args: tuple[Any, ...], kwargs: dict[str, Any]) -> Any:
    if "tip_centers" in kwargs:
        return kwargs["tip_centers"]
    return args[2] if len(args) >= 3 else None


def main(argv: list[str] | None = None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    wrapper, remaining = _wrapper_parser().parse_known_args(user_args)
    out_value = _option_value(remaining, "--out")
    if out_value is None:
        raise SystemExit("v10.0.5.13 requires --out")
    out = Path(out_value).resolve()
    out.mkdir(parents=True, exist_ok=True)

    option = load_barrier_option(wrapper.barrier_option, wrapper.barrier_registry)
    barrier_csv = option.write_barrier_csv(
        out / "barrier_only_inputs" / "selected_barrier_manifest.csv"
    )
    row = option.legacy_row(str(barrier_csv.resolve()))
    row.update(
        {
            "barrier_fields_transferred": option.audit_payload()[
                "barrier_fields_transferred"
            ],
            "candidate_state_fields_ignored": option.ignored_candidate_state,
        }
    )

    radius_um = float(wrapper.tip_refinement_radius_um)
    cluster_um = float(wrapper.selected_cluster_J_outer_um)
    local_um = float(wrapper.local_J_outer_um)
    if not all(math.isfinite(x) and x > 0.0 for x in (radius_um, cluster_um, local_um)):
        raise SystemExit("refinement and J radii must be finite and positive")
    if max(cluster_um, local_um) > radius_um:
        raise SystemExit(
            "J contour is not supported by physical refinement: "
            f"radius={radius_um:g} um cluster={cluster_um:g} um local={local_um:g} um"
        )

    _require_text_or_set(remaining, "--mpz-material-manifest", str(barrier_csv.resolve()))
    _require_text_or_set(remaining, "--mpz-material-class", option.canonical_class)
    _require_text_or_set(
        remaining, "--bulk-plasticity-mode", TWO_D_STATE_POLICY["bulk_plasticity_mode"]
    )
    _require_close_or_set(
        remaining, "--mpz-length-um", float(TWO_D_STATE_POLICY["mpz_length_um"])
    )
    _require_int_or_set(
        remaining, "--mpz-n-bins", int(TWO_D_STATE_POLICY["mpz_n_bins"])
    )
    _require_close_or_set(
        remaining, "--rJ-cluster", cluster_j_legacy_length_m(cluster_um)
    )
    _require_close_or_set(remaining, "--rJ-outer", local_um * 1.0e-6)
    _require_text_or_set(remaining, "--crack-backend", "adaptive_czm")
    _require_int_or_set(remaining, "--max-fronts", 1)

    radius_m = radius_um * 1.0e-6
    configure_physical_refinement_v100510(radius_m)

    status = {
        "schema": "barrier_only_production_manifest_v10_0_5_13",
        "model": MODEL_ID,
        "point_release": POINT_RELEASE,
        "started_utc": _utc_now(),
        "completed_utc": None,
        "status": "running",
        "run_completed_without_exception": False,
        "runtime_error_type": None,
        "runtime_error": None,
        "barrier_option": option.audit_payload(),
        "selected_barrier_csv": str(barrier_csv.resolve()),
        "selected_barrier_csv_sha256": _sha256(barrier_csv),
        "two_d_state_policy": TWO_D_STATE_POLICY,
        "candidate_state_fields_applied": False,
        "bulk_plasticity_mode": TWO_D_STATE_POLICY["bulk_plasticity_mode"],
        "selected_cluster_J_outer_um": cluster_um,
        "cluster_J_legacy_length_m": cluster_j_legacy_length_m(cluster_um),
        "local_J_outer_um": local_um,
        "refinement": {
            "requested_radius_um": radius_um,
            "ordinary_corridor_and_remesh_path": True,
        },
        "underlying_solver": "mode_i_first_passage_v9_18_5_6",
        "branching_enabled": False,
        "max_fronts": 1,
        "argv": remaining,
    }
    manifest_path = out / PRODUCTION_MANIFEST
    _write(manifest_path, status)

    saved = {
        "load_selected_row": _v911.load_selected_row,
        "apply_exact_barrier_args": _v911.apply_exact_barrier_args,
        "build_mpz_config": _v911.build_mpz_config,
        "set_bulk_defaults": _v911._set_bulk_pt_namespace_defaults,
        "compact_audit": _v911.compact_audit,
        "bulk_apply_pt": _bulk.apply_pt_dislocation_config,
        "make_mesh": _mesh.make_tri_mesh,
        "rebuild_mesh": _crack_backend.rebuild_tri_mesh,
    }
    env_saved = {
        name: os.environ.get(name)
        for name in (
            "ARRHENIUS_EVENT_STATISTICS",
            "ARRHENIUS_STOCHASTIC_EMISSION",
            "ARRHENIUS_PREFINED_MODE_I_CORRIDOR",
        )
    }

    def load_selected(_path, material_class=None):
        if material_class is not None and str(material_class) != option.canonical_class:
            raise ValueError(
                f"material class {material_class!r} conflicts with {option.canonical_class!r}"
            )
        return dict(row)

    original_rebuild = _crack_backend.rebuild_tri_mesh

    def rebuild_with_refinement_metadata(*args, **kwargs):
        mesh = original_rebuild(*args, **kwargs)
        return _annotate_mesh(mesh, radius_m, _tip_centers(args, kwargs))

    _v911.load_selected_row = load_selected
    # Opening/emission barriers use the existing exact v9.11 function; its
    # non-barrier numerical controls are the already validated v9.11 controls.
    _v911.apply_exact_barrier_args = saved["apply_exact_barrier_args"]
    _v911.build_mpz_config = _build_barrier_only_mpz_config
    _v911._set_bulk_pt_namespace_defaults = _set_bulk_barrier_defaults
    _v911.compact_audit = _compact_audit
    _bulk.apply_pt_dislocation_config = _apply_barrier_pt_config
    _mesh.make_tri_mesh = make_physical_refinement_mesh_v100510
    _crack_backend.rebuild_tri_mesh = rebuild_with_refinement_metadata
    os.environ["ARRHENIUS_EVENT_STATISTICS"] = "deterministic"
    os.environ["ARRHENIUS_STOCHASTIC_EMISSION"] = "0"
    os.environ["ARRHENIUS_PREFINED_MODE_I_CORRIDOR"] = "1"

    try:
        result = _solver.main(remaining)
        mesh_data = _mesh_payload(_v9185._RUNTIME.get("mesh"), radius_m)
        if not mesh_data["actual_radius_verified"]:
            raise RuntimeError(
                "run returned without verified physical refinement metadata: "
                + json.dumps(mesh_data, default=str)
            )
        status.update(
            {
                "completed_utc": _utc_now(),
                "status": "complete",
                "run_completed_without_exception": True,
                "mesh_refinement_runtime": mesh_data,
            }
        )
        _write(manifest_path, status)
        return result
    except BaseException as exc:
        status.update(
            {
                "completed_utc": _utc_now(),
                "status": "failed",
                "run_completed_without_exception": False,
                "runtime_error_type": type(exc).__name__,
                "runtime_error": str(exc),
                "mesh_refinement_runtime": _mesh_payload(
                    _v9185._RUNTIME.get("mesh"), radius_m
                ),
            }
        )
        _write(manifest_path, status)
        raise
    finally:
        _v911.load_selected_row = saved["load_selected_row"]
        _v911.apply_exact_barrier_args = saved["apply_exact_barrier_args"]
        _v911.build_mpz_config = saved["build_mpz_config"]
        _v911._set_bulk_pt_namespace_defaults = saved["set_bulk_defaults"]
        _v911.compact_audit = saved["compact_audit"]
        _bulk.apply_pt_dislocation_config = saved["bulk_apply_pt"]
        _mesh.make_tri_mesh = saved["make_mesh"]
        _crack_backend.rebuild_tri_mesh = saved["rebuild_mesh"]
        clear_physical_refinement_v100510()
        for name, value in env_saved.items():
            if value is None:
                os.environ.pop(name, None)
            else:
                os.environ[name] = value


if __name__ == "__main__":
    main()


__all__ = [
    "POINT_RELEASE",
    "MODEL_ID",
    "PRODUCTION_MANIFEST",
    "cluster_j_legacy_length_m",
    "_apply_barrier_pt_config",
    "_build_barrier_only_mpz_config",
    "_set_bulk_barrier_defaults",
    "main",
]
