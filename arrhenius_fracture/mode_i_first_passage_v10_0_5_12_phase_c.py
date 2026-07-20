"""Phase-C monotonic entry with production refinement and v9.11.1 options.

The wrapper composes the diagnostics-complete v10.0.5.2 progressive monotonic
solver with two production integrations that were audit-only at v10.0.5.11:

* a fixed physical refinement radius propagated through the ordinary corridor and
  remesh path; and
* an exact four-option v9.11.1 response registry, including the peak option.

No barrier is fitted or mutated here.  The selected registry row is fingerprinted,
materialized as a one-row v9.11 CSV, and consumed by the existing constitutive
parser.  Mobile shielding is forced to the registry value zero; retained state is
the only explicit unresolved shielding contribution.

``sharp_front --rJ-cluster`` accepts a legacy domain length ``ell`` whose actual
outer contour radius is approximately ``8*ell``.  Phase C exposes the requested
physical outer radius and converts it exactly once before validating the solver
argument.
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
from typing import Any

from . import mesh as _mesh
from . import mode_i_first_passage_v10_0_3_progressive as _v1003
from . import mode_i_first_passage_v10_0_5_2_parallel as _v10052
from . import mode_i_first_passage_v9_18_5 as _v9185
from .kinetic_campaign_czm import KineticCampaignCZMConfig as _BaseKineticConfig
from .mode_i_first_passage_v10_0 import _option_value, _replace_option
from .mpz_response_registry_v100512 import (
    PARAMETER_SOURCE,
    load_option,
)
from .physical_refinement_mesh_v100510 import (
    clear_physical_refinement_v100510,
    configure_physical_refinement_v100510,
    make_physical_refinement_mesh_v100510,
)

POINT_RELEASE = "10.0.5.12.1"
MODEL_ID = "FEM_CZM_Phase_C_four_option_monotonic_v10_0_5_12_1"
PRODUCTION_MANIFEST = "phase_c_production_manifest_v10_0_5_12.json"
CLUSTER_J_RADIUS_TO_LEGACY_ELL = 8.0


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


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
    if existing is not None and not math.isclose(existing, expected, rel_tol=rel_tol, abs_tol=abs_tol):
        raise SystemExit(f"{name}={existing:.16g} conflicts with Phase-C value {expected:.16g}")
    _replace_option(argv, name, f"{expected:.16g}")


def _require_int_or_set(argv: list[str], name: str, expected: int) -> None:
    value = _option_value(argv, name)
    if value is not None and int(value) != int(expected):
        raise SystemExit(f"{name}={value} conflicts with Phase-C value {expected}")
    _replace_option(argv, name, str(int(expected)))


def cluster_j_legacy_length_m(physical_outer_radius_um: float) -> float:
    radius_um = float(physical_outer_radius_um)
    if not math.isfinite(radius_um) or radius_um <= 0.0:
        raise ValueError("cluster-J physical outer radius must be finite and positive")
    return radius_um * 1.0e-6 / CLUSTER_J_RADIUS_TO_LEGACY_ELL


class _RegistryKineticConfig(_BaseKineticConfig):
    """Preserve the registry contract that mobile content does not shield."""

    def __init__(self, *args, **kwargs):
        kwargs["mobile_shield_fraction"] = 0.0
        super().__init__(*args, **kwargs)


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


def _write(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, default=str))


def _wrapper_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
    parser.add_argument("--phase-c-option", required=True)
    parser.add_argument("--phase-c-registry", type=Path, default=None)
    parser.add_argument("--tip-refinement-radius-um", type=float, default=330.0)
    parser.add_argument("--selected-cluster-J-outer-um", type=float, default=240.0)
    parser.add_argument("--local-J-outer-um", type=float, default=100.0)
    return parser


def main(argv: list[str] | None = None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    wrapper, remaining = _wrapper_parser().parse_known_args(user_args)
    out_value = _option_value(remaining, "--out")
    if out_value is None:
        raise SystemExit("Phase-C production entry requires --out")
    out = Path(out_value).resolve()
    out.mkdir(parents=True, exist_ok=True)
    audit_path = out / PRODUCTION_MANIFEST

    option = load_option(wrapper.phase_c_option, wrapper.phase_c_registry)
    radius_um = float(wrapper.tip_refinement_radius_um)
    cluster_um = float(wrapper.selected_cluster_J_outer_um)
    local_um = float(wrapper.local_J_outer_um)
    if not all(math.isfinite(x) and x > 0.0 for x in (radius_um, cluster_um, local_um)):
        raise SystemExit("Phase-C refinement and J radii must be finite and positive")
    if max(cluster_um, local_um) > radius_um:
        raise SystemExit(
            "Phase-C J contour is not supported by the physical refinement radius: "
            f"radius={radius_um:g} um, cluster={cluster_um:g} um, local={local_um:g} um"
        )

    input_dir = out / "phase_c_inputs"
    manifest_csv = option.write_selected_csv(input_dir / "selected_parameter_manifest.csv")
    manifest = option.material_manifest(source_path=str(manifest_csv.resolve()))

    cluster_ell_m = cluster_j_legacy_length_m(cluster_um)
    _replace_option(remaining, "--v10-material-source", PARAMETER_SOURCE)
    _replace_option(remaining, "--v10-material-class", option.option_key)
    _replace_option(remaining, "--mpz-material-manifest", str(manifest_csv.resolve()))
    _replace_option(remaining, "--mpz-material-class", option.canonical_class)
    _require_close_or_set(remaining, "--mpz-length-um", option.mpz_length_um)
    _require_int_or_set(remaining, "--mpz-n-bins", option.mpz_n_bins)
    _require_close_or_set(remaining, "--rJ-cluster", cluster_ell_m)
    _require_close_or_set(remaining, "--rJ-outer", local_um * 1.0e-6)
    _replace_option(remaining, "--crack-backend", "adaptive_czm")
    _replace_option(remaining, "--max-fronts", "1")

    radius_m = radius_um * 1.0e-6
    configure_physical_refinement_v100510(radius_m)

    saved = {
        "source": _v1003.PF_SOURCE,
        "loader": _v1003.load_material_manifest,
        "path": _v1003.pf_manifest_path,
        "config": _v1003.KineticCampaignCZMConfig,
        "mesh": _mesh.make_tri_mesh,
    }
    env_saved: dict[str, str | None] = {}

    def load_selected(material_class: str, *, parameter_source: str = PARAMETER_SOURCE, **_):
        if parameter_source != PARAMETER_SOURCE:
            raise ValueError(
                f"Phase-C loader requires parameter_source={PARAMETER_SOURCE!r}, got {parameter_source!r}"
            )
        if str(material_class) != option.option_key:
            raise ValueError(
                f"Phase-C invocation selected {material_class!r}, expected {option.option_key!r}"
            )
        return manifest

    def selected_path(material_class: str) -> Path:
        if str(material_class) != option.option_key:
            raise ValueError(
                f"Phase-C manifest path requested for {material_class!r}, expected {option.option_key!r}"
            )
        return manifest_csv.resolve()

    status: dict[str, Any] = {
        "schema": "phase_c_production_manifest_v10_0_5_12_1",
        "model": MODEL_ID,
        "point_release": POINT_RELEASE,
        "started_utc": _utc_now(),
        "completed_utc": None,
        "status": "running",
        "run_completed_without_exception": False,
        "runtime_error_type": None,
        "runtime_error": None,
        "option": option.audit_payload(),
        "selected_manifest_csv": str(manifest_csv.resolve()),
        "selected_manifest_sha256": _sha256(manifest_csv),
        "selected_cluster_J_outer_um": cluster_um,
        "cluster_J_legacy_length_m": cluster_ell_m,
        "cluster_J_radius_to_legacy_ell": CLUSTER_J_RADIUS_TO_LEGACY_ELL,
        "local_J_outer_um": local_um,
        "refinement": {
            "requested_radius_um": radius_um,
            "audit_only": False,
            "ordinary_corridor_and_remesh_path": True,
        },
        "event_statistics": "mean_field",
        "stochastic_emission": False,
        "mobile_shield_fraction": 0.0,
        "branching_enabled": False,
        "max_fronts": 1,
        "underlying_solver": _v10052.MODEL_ID,
        "underlying_point_release": _v10052.POINT_RELEASE,
        "argv": remaining,
    }
    _write(audit_path, status)

    _v1003.PF_SOURCE = PARAMETER_SOURCE
    _v1003.load_material_manifest = load_selected
    _v1003.pf_manifest_path = selected_path
    _v1003.KineticCampaignCZMConfig = _RegistryKineticConfig
    _mesh.make_tri_mesh = make_physical_refinement_mesh_v100510
    _set_env("ARRHENIUS_EVENT_STATISTICS", "mean_field", env_saved)
    _set_env("ARRHENIUS_STOCHASTIC_EMISSION", "0", env_saved)
    _set_env("ARRHENIUS_VHCF_FEM_CACHE", "0", env_saved)
    _set_env("ARRHENIUS_PREFINED_MODE_I_CORRIDOR", "1", env_saved)

    try:
        result = _v10052.main(remaining)
        mesh_data = _mesh_payload(_v9185._RUNTIME.get("mesh"), radius_m)
        if not mesh_data["actual_radius_verified"]:
            raise RuntimeError(
                "Phase-C solver returned without verified physical refinement metadata: "
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
        _write(audit_path, status)
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
        _write(audit_path, status)
        raise
    finally:
        _v1003.PF_SOURCE = saved["source"]
        _v1003.load_material_manifest = saved["loader"]
        _v1003.pf_manifest_path = saved["path"]
        _v1003.KineticCampaignCZMConfig = saved["config"]
        _mesh.make_tri_mesh = saved["mesh"]
        clear_physical_refinement_v100510()
        _restore_env(env_saved)


if __name__ == "__main__":
    main()


__all__ = [
    "POINT_RELEASE",
    "MODEL_ID",
    "PRODUCTION_MANIFEST",
    "CLUSTER_J_RADIUS_TO_LEGACY_ELL",
    "cluster_j_legacy_length_m",
    "main",
]
