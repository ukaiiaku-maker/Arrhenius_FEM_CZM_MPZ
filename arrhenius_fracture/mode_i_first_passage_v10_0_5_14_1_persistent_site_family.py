"""v10.0.5.14.1: PF v10.2.22 persistent sites with the real v10.2.14 kernel family."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys
from typing import Any

import numpy as np

from . import barrier_only_response_registry_v100513 as _legacy_registry
from . import mixed_mode_first_passage_v8 as _mm
from . import mixed_mode_first_passage_v9_11 as _v911
from . import mode_i_first_passage_v10_0_5_13_barrier_only as _core
from . import mode_i_first_passage_v10_0_5_13_4_barrier_only as _policy_entry
from . import mode_i_first_passage_v10_0_5_13_5_barrier_only as _base
from . import mode_i_first_passage_v10_0_5_14_persistent_site as _v100514
from .anisotropic_two_channel_drive_v100514 import augmented_j_wrapper_factory
from .persistent_site_front_engine_v100514 import (
    PersistentSiteMovingProcessZoneFrontEngineV100514,
)
from .persistent_site_registry_v100514 import ROWS, select_persistent_site_row
from .persistent_site_signed_support_v100514 import SignedShieldingKernelV100514
from .signed_kernel_family_v1005141 import (
    FAMILY_SCHEMA,
    SignedShieldingKernelFamilyV1005141,
    load_signed_shielding_artifact_v1005141,
)

POINT_RELEASE = "10.0.5.14.1"
MODEL_ID = "FEM_CZM_full_2D_PF_v10_2_22_persistent_site_kernel_family_v10_0_5_14_1"
PRODUCTION_MANIFEST = "persistent_site_production_manifest_v10_0_5_14_1.json"


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
    parser.add_argument(
        "--persistent-site-option", required=True, choices=tuple(ROWS)
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--signed-kernel-family", type=Path)
    group.add_argument("--signed-shielding-kernel", type=Path)
    return parser


def _option_audit(option) -> dict[str, Any]:
    payload = option.audit_payload()
    payload["point_release"] = POINT_RELEASE
    return payload


def _artifact_audit(artifact) -> dict[str, Any]:
    if isinstance(artifact, SignedShieldingKernelFamilyV1005141):
        return artifact.audit_payload()
    return {
        "schema": artifact.metadata.get("schema", "fixed_signed_kernel"),
        "point_release": POINT_RELEASE,
        "source_path": artifact.source_path,
        "artifact_kind": "fixed_kernel_compatibility",
        "active_shape": list(
            artifact.active_kernel_Pa_sqrt_m_per_signed_line.shape
        ),
        "wake_shape": list(
            artifact.wake_kernel_Pa_sqrt_m_per_signed_line.shape
        ),
        "activation_to_line_content_by_system": (
            artifact.activation_to_line_content_by_system.tolist()
        ),
        "fixed_kernel_approximation": True,
    }


def _compact_audit_factory(option, policy, artifact):
    artifact_audit = _artifact_audit(artifact)

    def compact(row):
        return {
            "model": MODEL_ID,
            "point_release": POINT_RELEASE,
            "candidate_id": option.candidate_id,
            "target_class": option.canonical_class,
            "parameter_source": _v100514.PARAMETER_SOURCE,
            "candidate_state_fields_applied": True,
            "persistent_site_option": _option_audit(option),
            "two_d_state_policy": policy,
            "signed_kernel_artifact": artifact_audit,
            "bulk_state_evolves_in_fem": False,
            "moving_crack_tip_mpz_active": True,
            "shielding_derived_from_signed_retained_state": True,
            "finite_source_inventory": False,
            "source_depletion": False,
            "source_refresh": False,
            "kernel_interpolated_by_cumulative_crack_path_extension": (
                artifact_audit.get("artifact_kind")
                == "crack_extension_kernel_family"
            ),
        }

    return compact


def _validate_artifact(artifact, candidate) -> None:
    n_systems = int(candidate.n_slip_channels)
    n_bins = int(candidate.n_bins_recommended)
    length_m = float(candidate.L_pz_um_recommended) * 1.0e-6
    dx = length_m / n_bins
    active_x = (np.arange(n_bins, dtype=float) + 0.5) * dx
    wake_x = active_x.copy()
    if isinstance(artifact, SignedShieldingKernelFamilyV1005141):
        artifact.validate(
            n_systems,
            n_bins,
            active_x_m=active_x,
            wake_x_m=wake_x,
        )
        # Force exact initial-state construction before the FEM solve starts.
        artifact.snapshot(0.0, active_x, wake_x)
    else:
        artifact.validate(n_systems, n_bins)


def _update_manifest(out, option, policy, artifact, completed: bool) -> None:
    inherited = out / _base.PRODUCTION_MANIFEST
    inherited_payload: dict[str, Any] = {}
    if inherited.is_file():
        try:
            inherited_payload = json.loads(inherited.read_text())
        except Exception:
            inherited_payload = {}
    artifact_audit = _artifact_audit(artifact)
    payload = {
        **inherited_payload,
        "schema": "persistent_site_production_manifest_v10_0_5_14_1",
        "model": MODEL_ID,
        "point_release": POINT_RELEASE,
        "status": (
            "complete" if completed else inherited_payload.get("status", "failed")
        ),
        "completed_utc": _v100514._utc_now(),
        "run_completed_without_exception": bool(completed),
        "persistent_site_option": _option_audit(option),
        "two_d_state_policy": policy,
        "candidate_state_fields_applied": True,
        "signed_kernel_artifact": artifact_audit,
        "physics_contract": {
            "reference": "PF v10.2.22 commit 198ece3aeb1d193a8c1c4857676fba720c088d27",
            "persistent_source_sites": True,
            "finite_source_inventory": False,
            "source_depletion": False,
            "source_refresh": False,
            "available_site_fraction_required": 1.0,
            "source_sites_refreshed_required": 0.0,
            "signed_mobile_retained_state": True,
            "unsigned_backstress_signed_shielding_split": True,
            "mobile_shield_fraction": 0.0,
            "front_width_grid_independent": True,
            "implicit_backstress_complementarity": True,
            "fractional_moving_frame": True,
            "trial_commit_state": True,
            "wake_shielding": False,
            "bulk_plasticity_mode": "tip_only",
            "kernel_family_schema": (
                FAMILY_SCHEMA
                if isinstance(artifact, SignedShieldingKernelFamilyV1005141)
                else None
            ),
            "kernel_interpolation_coordinate": (
                "cumulative_crack_path_extension_m"
                if isinstance(artifact, SignedShieldingKernelFamilyV1005141)
                else "fixed_kernel"
            ),
            "kernel_extrapolation_allowed": False,
            "kernel_spatial_projection": (
                "piecewise_linear_with_endpoint_hold"
                if isinstance(artifact, SignedShieldingKernelFamilyV1005141)
                else "preprojected_fixed_grid"
            ),
            "constitutive_K_shield_cap": False,
        },
    }
    (out / PRODUCTION_MANIFEST).write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n"
    )


def main(argv: list[str] | None = None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    wrapper, remaining = _parser().parse_known_args(user_args)
    candidate = select_persistent_site_row(wrapper.persistent_site_option)
    option = _v100514.PersistentSiteOptionAdapterV100514(candidate)
    artifact_path = (
        wrapper.signed_kernel_family
        if wrapper.signed_kernel_family is not None
        else wrapper.signed_shielding_kernel
    )
    artifact = load_signed_shielding_artifact_v1005141(artifact_path)
    if wrapper.signed_kernel_family is not None and not isinstance(
        artifact, SignedShieldingKernelFamilyV1005141
    ):
        raise SystemExit(
            "--signed-kernel-family requires schema " + FAMILY_SCHEMA
        )
    _validate_artifact(artifact, candidate)
    policy = _v100514.persistent_site_policy(candidate)
    policy = {
        **policy,
        "policy_id": "PF_v10_2_22_persistent_sites_signed_kernel_family_v1005141",
        "kernel_artifact_kind": (
            "crack_extension_family"
            if isinstance(artifact, SignedShieldingKernelFamilyV1005141)
            else "fixed_kernel_compatibility"
        ),
        "kernel_interpolation_coordinate": (
            "cumulative_crack_path_extension_m"
            if isinstance(artifact, SignedShieldingKernelFamilyV1005141)
            else "fixed"
        ),
        "kernel_extrapolation_allowed": False,
        "wake_shielding_active": False,
        "constitutive_K_shield_cap": False,
    }

    out_value = None
    if "--out" in remaining:
        index = remaining.index("--out")
        if index + 1 < len(remaining):
            out_value = remaining[index + 1]
    if out_value is None:
        raise SystemExit("v10.0.5.14.1 requires --out")
    out = Path(out_value).expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)

    _v100514._replace_option(remaining, "--barrier-option", candidate.option_key)
    _v100514._replace_option(remaining, "--bulk-plasticity-mode", "tip_only")
    _v100514._replace_option(
        remaining, "--mpz-length-um", str(candidate.L_pz_um_recommended)
    )
    _v100514._replace_option(
        remaining, "--mpz-n-bins", str(candidate.n_bins_recommended)
    )
    _v100514._replace_option(remaining, "--c-blunt", str(candidate.c_blunt))

    saved = {
        "core_load": _core.load_barrier_option,
        "core_policy": _core.TWO_D_STATE_POLICY,
        "core_audit": _core._compact_audit,
        "registry_policy": _legacy_registry.TWO_D_STATE_POLICY,
        "entry_policy": _policy_entry.TIP_ONLY_POLICY,
        "v911_engine": _v911.MovingProcessZone2DFrontEngine,
        "j_wrapper_factory": _mm._j_wrapper_factory,
    }
    PersistentSiteMovingProcessZoneFrontEngineV100514.configure(
        candidate, artifact
    )
    _core.load_barrier_option = lambda value, path=None: option
    _core.TWO_D_STATE_POLICY = policy
    _core._compact_audit = _compact_audit_factory(option, policy, artifact)
    _legacy_registry.TWO_D_STATE_POLICY = policy
    _policy_entry.TIP_ONLY_POLICY = policy
    _v911.MovingProcessZone2DFrontEngine = (
        PersistentSiteMovingProcessZoneFrontEngineV100514
    )
    _mm._j_wrapper_factory = augmented_j_wrapper_factory(
        saved["j_wrapper_factory"]
    )

    selection = {
        "schema": MODEL_ID,
        "point_release": POINT_RELEASE,
        "created_utc": _v100514._utc_now(),
        "option": _option_audit(option),
        "policy": policy,
        "signed_kernel_artifact": _artifact_audit(artifact),
        "argv": remaining,
    }
    (out / "persistent_site_parameter_selection_v10_0_5_14_1.json").write_text(
        json.dumps(selection, indent=2, sort_keys=True, default=str) + "\n"
    )

    try:
        result = _base.main(remaining)
        _update_manifest(out, option, policy, artifact, completed=True)
        return result
    except BaseException:
        _update_manifest(out, option, policy, artifact, completed=False)
        raise
    finally:
        _core.load_barrier_option = saved["core_load"]
        _core.TWO_D_STATE_POLICY = saved["core_policy"]
        _core._compact_audit = saved["core_audit"]
        _legacy_registry.TWO_D_STATE_POLICY = saved["registry_policy"]
        _policy_entry.TIP_ONLY_POLICY = saved["entry_policy"]
        _v911.MovingProcessZone2DFrontEngine = saved["v911_engine"]
        _mm._j_wrapper_factory = saved["j_wrapper_factory"]
        PersistentSiteMovingProcessZoneFrontEngineV100514.clear_configuration()


if __name__ == "__main__":
    main()


__all__ = [
    "POINT_RELEASE",
    "MODEL_ID",
    "PRODUCTION_MANIFEST",
    "main",
]
