"""v10.0.5.14: PF v10.2.22 persistent-site physics on validated FEM/CZM mechanics."""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any

from . import barrier_only_response_registry_v100513 as _legacy_registry
from . import mixed_mode_first_passage_v8 as _mm
from . import mixed_mode_first_passage_v9_11 as _v911
from . import mode_i_first_passage_v10_0_5_13_barrier_only as _core
from . import mode_i_first_passage_v10_0_5_13_4_barrier_only as _policy_entry
from . import mode_i_first_passage_v10_0_5_13_5_barrier_only as _base
from .anisotropic_two_channel_drive_v100514 import augmented_j_wrapper_factory
from .persistent_site_front_engine_v100514 import (
    PersistentSiteMovingProcessZoneFrontEngineV100514,
)
from .persistent_site_registry_v100514 import (
    PARAMETER_SOURCE,
    ROWS,
    PersistentSiteRowV100514,
    select_persistent_site_row,
)
from .persistent_site_signed_mpz_v100514 import SignedShieldingKernelV100514

POINT_RELEASE = "10.0.5.14"
MODEL_ID = "FEM_CZM_full_2D_PF_v10_2_22_persistent_site_parity_v10_0_5_14"
PRODUCTION_MANIFEST = "persistent_site_production_manifest_v10_0_5_14.json"


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _fingerprint(payload: Any) -> str:
    return hashlib.sha256(
        json.dumps(
            payload, sort_keys=True, separators=(",", ":"), default=str
        ).encode()
    ).hexdigest()


def _replace_option(argv: list[str], name: str, value: str) -> None:
    while name in argv:
        index = argv.index(name)
        del argv[index : min(index + 2, len(argv))]
    argv.extend([name, str(value)])


class PersistentSiteOptionAdapterV100514:
    """Compatibility object consumed by the inherited v10.0.5.13 entry."""

    def __init__(self, candidate: PersistentSiteRowV100514):
        self.candidate = candidate
        self.option_key = candidate.option_key
        self.candidate_id = candidate.candidate_id
        self.canonical_class = "DBTT"
        self.barrier_row = candidate.barrier_row()
        self.barrier_fingerprint_sha256 = _fingerprint(self.barrier_row)
        self.source_registry_path = "PF_v10_2_22_to_FEM_CZM_physics_handoff.md"
        self.ignored_candidate_state: dict[str, Any] = {}

    def full_row(self) -> dict[str, Any]:
        c = self.candidate
        return {
            **self.barrier_row,
            "option_key": c.option_key,
            "candidate_id": c.candidate_id,
            "target_class": "DBTT",
            "selection_role": "primary",
            "parameter_source": PARAMETER_SOURCE,
            "barrier_fingerprint_sha256": self.barrier_fingerprint_sha256,
            "source_sites_per_system": c.source_sites_per_system_provenance,
            "source_sites_per_system_active": False,
            "rho_source0_m2": c.rho_source0_m2,
            "encounter_efficiency": c.encounter_efficiency,
            "retained_recovery_rate_s": 0.0,
            "source_refresh_length_um": c.source_refresh_length_um_provenance,
            "source_refresh_length_active": False,
            "c_blunt": c.c_blunt,
            "taylor_corr_rho_c_m2": c.taylor_corr_rho_c_m2,
            "taylor_corr_scale": c.taylor_corr_scale,
            "L_pz_um_recommended": c.L_pz_um_recommended,
            "n_bins_recommended": c.n_bins_recommended,
            "rho_forest_floor_m2": c.rho_forest_floor_m2,
            "peierls_stress_fraction": c.peierls_stress_fraction,
            "taylor_stress_fraction": c.taylor_stress_fraction,
            "mobile_shield_fraction": c.mobile_shield_fraction,
            "source_recovery_rate_s": 0.0,
            "reference_source_area_um2": c.reference_source_area_um2,
            "reference_front_width_um": c.reference_front_width_um,
            "source_zone_length_um": c.source_zone_length_um,
            "legacy_source_sites_active": 0,
            "legacy_source_refresh_active": 0,
            "explicit_recovery_active": 0,
            "persistent_source_inventory_active": False,
            "source_depletion_active": False,
            "source_refresh_active": False,
            "front_width_grid_independent": True,
        }

    def legacy_row(self, manifest_path: str | None = None) -> dict[str, Any]:
        row = self.full_row()
        if manifest_path is not None:
            row["parameter_manifest"] = str(manifest_path)
        row["parameter_fingerprint_sha256"] = _fingerprint(row)
        return row

    def write_barrier_csv(self, path: str | Path) -> Path:
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        row = self.full_row()
        with target.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(row))
            writer.writeheader()
            writer.writerow(row)
        return target

    def audit_payload(self) -> dict[str, Any]:
        c = self.candidate
        return {
            "point_release": POINT_RELEASE,
            "parameter_source": PARAMETER_SOURCE,
            "option_key": c.option_key,
            "candidate_id": c.candidate_id,
            "role": c.role,
            "candidate_fingerprint_sha256": c.fingerprint(),
            "barrier_fingerprint_sha256": self.barrier_fingerprint_sha256,
            "candidate_state_fields_applied": True,
            "barrier_fields_transferred": list(self.barrier_row),
            "candidate_state_fields_ignored": {},
            "rho_source0_semantics": "persistent_areal_nucleation_site_density_m2",
            "finite_source_inventory": False,
            "source_depletion": False,
            "source_refresh": False,
            "explicit_recovery": False,
            "signed_mobile_retained_state": True,
            "front_width_grid_independent": True,
            "front_width_minimum_semantics": "max(explicit_physical_minimum,b)",
            "ahead_of_tip_dx_used_as_front_width_floor": False,
            "shared_contract": self.full_row(),
        }


def persistent_site_policy(candidate: PersistentSiteRowV100514) -> dict[str, Any]:
    dx_um = candidate.L_pz_um_recommended / candidate.n_bins_recommended
    return {
        "policy_id": "PF_v10_2_22_persistent_sites_signed_MPZ_v100514",
        "bulk_plasticity_mode": "tip_only",
        "mpz_length_um": candidate.L_pz_um_recommended,
        "mpz_n_bins": candidate.n_bins_recommended,
        "source_sites_per_system": candidate.source_sites_per_system_provenance,
        "source_recovery_rate_s": 0.0,
        "source_refresh_length_um": candidate.source_refresh_length_um_provenance,
        "source_bin_count": max(
            1, int(math.ceil(candidate.source_zone_length_um / dx_um))
        ),
        "shielding_orientation_factors": (1.0, 1.0),
        "mobile_shield_fraction": 0.0,
        "shielding_core_m": 2.5e-10,
        "retained_recovery_nu0_s": 0.0,
        "retained_recovery_barrier_eV": 0.0,
        "retained_recovery_activation_volume_b3": 0.0,
        "mobile_recovery_rate_s": 0.0,
        "pair_annihilation_rate_per_count_s": 0.0,
        "blunting_length_um": 0.5,
        "blunting_slip_fraction": 1.0,
        "peierls_stress_fraction": candidate.peierls_stress_fraction,
        "taylor_stress_fraction": candidate.taylor_stress_fraction,
        "taylor_corr_rho_c_m2": candidate.taylor_corr_rho_c_m2,
        "taylor_renewal_time_s": 1.0,
        "taylor_m_exponent": 1.0,
        "taylor_m_scale": candidate.taylor_corr_scale,
        "taylor_m_cap": float("inf"),
        "encounter_efficiency": candidate.encounter_efficiency,
        "forest_density_floor_m2": candidate.rho_forest_floor_m2,
        "candidate_source_inventory_applied": False,
        "candidate_persistent_source_density_applied": True,
        "candidate_source_refresh_applied": False,
        "candidate_encounter_applied": True,
        "candidate_shielding_blunting_applied": True,
        "candidate_initial_state_applied": False,
        "state_configuration_source": "PF_v10.2.22_handoff_exact_contract",
        "state_evolution_source": "persistent_signed_moving_crack_tip_MPZ",
        "continuum_bulk_role": "elastic_fem_only",
        "uniform_bulk_mobile_retained_state_active": False,
        "trial_commit_state_active": True,
    }


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
    parser.add_argument(
        "--persistent-site-option", required=True, choices=tuple(ROWS)
    )
    parser.add_argument("--signed-shielding-kernel", type=Path, required=True)
    return parser


def _compact_audit_factory(
    option: PersistentSiteOptionAdapterV100514, policy: dict[str, Any]
):
    def compact(row):
        return {
            "model": MODEL_ID,
            "candidate_id": option.candidate_id,
            "target_class": option.canonical_class,
            "parameter_source": PARAMETER_SOURCE,
            "candidate_state_fields_applied": True,
            "persistent_site_option": option.audit_payload(),
            "two_d_state_policy": policy,
            "bulk_state_evolves_in_fem": False,
            "moving_crack_tip_mpz_active": True,
            "shielding_derived_from_signed_retained_state": True,
            "finite_source_inventory": False,
            "source_depletion": False,
            "source_refresh": False,
        }

    return compact


def _update_manifest(
    out: Path,
    option: PersistentSiteOptionAdapterV100514,
    policy: dict[str, Any],
    kernel: SignedShieldingKernelV100514,
    completed: bool,
) -> None:
    inherited = out / _base.PRODUCTION_MANIFEST
    inherited_payload: dict[str, Any] = {}
    if inherited.is_file():
        try:
            inherited_payload = json.loads(inherited.read_text())
        except Exception:
            inherited_payload = {}
    payload = {
        **inherited_payload,
        "schema": "persistent_site_production_manifest_v10_0_5_14",
        "model": MODEL_ID,
        "point_release": POINT_RELEASE,
        "status": (
            "complete" if completed else inherited_payload.get("status", "failed")
        ),
        "completed_utc": _utc_now(),
        "run_completed_without_exception": bool(completed),
        "persistent_site_option": option.audit_payload(),
        "two_d_state_policy": policy,
        "candidate_state_fields_applied": True,
        "signed_shielding_kernel": {
            "source_path": kernel.source_path,
            "active_shape": list(
                kernel.active_kernel_Pa_sqrt_m_per_signed_line.shape
            ),
            "wake_shape": list(
                kernel.wake_kernel_Pa_sqrt_m_per_signed_line.shape
            ),
            "activation_to_line_content_by_system": (
                kernel.activation_to_line_content_by_system.tolist()
            ),
            "metadata": kernel.metadata,
        },
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
        },
    }
    (out / PRODUCTION_MANIFEST).write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n"
    )


def main(argv: list[str] | None = None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    wrapper, remaining = _parser().parse_known_args(user_args)
    candidate = select_persistent_site_row(wrapper.persistent_site_option)
    option = PersistentSiteOptionAdapterV100514(candidate)
    kernel = SignedShieldingKernelV100514.from_json(
        wrapper.signed_shielding_kernel
    )
    kernel.validate(candidate.n_slip_channels, candidate.n_bins_recommended)
    policy = persistent_site_policy(candidate)

    out_value = None
    if "--out" in remaining:
        index = remaining.index("--out")
        if index + 1 < len(remaining):
            out_value = remaining[index + 1]
    if out_value is None:
        raise SystemExit("v10.0.5.14 requires --out")
    out = Path(out_value).expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)

    _replace_option(remaining, "--barrier-option", candidate.option_key)
    _replace_option(remaining, "--bulk-plasticity-mode", "tip_only")
    _replace_option(
        remaining, "--mpz-length-um", str(candidate.L_pz_um_recommended)
    )
    _replace_option(
        remaining, "--mpz-n-bins", str(candidate.n_bins_recommended)
    )
    _replace_option(remaining, "--c-blunt", str(candidate.c_blunt))

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
        candidate, kernel
    )
    _core.load_barrier_option = lambda value, path=None: option
    _core.TWO_D_STATE_POLICY = policy
    _core._compact_audit = _compact_audit_factory(option, policy)
    _legacy_registry.TWO_D_STATE_POLICY = policy
    _policy_entry.TIP_ONLY_POLICY = policy
    _v911.MovingProcessZone2DFrontEngine = (
        PersistentSiteMovingProcessZoneFrontEngineV100514
    )
    _mm._j_wrapper_factory = augmented_j_wrapper_factory(
        saved["j_wrapper_factory"]
    )

    (out / "persistent_site_parameter_selection_v10_0_5_14.json").write_text(
        json.dumps(
            {
                "schema": MODEL_ID,
                "created_utc": _utc_now(),
                "option": option.audit_payload(),
                "policy": policy,
                "kernel_source": kernel.source_path,
                "argv": remaining,
            },
            indent=2,
            sort_keys=True,
            default=str,
        )
        + "\n"
    )

    try:
        result = _base.main(remaining)
        _update_manifest(out, option, policy, kernel, completed=True)
        return result
    except BaseException:
        _update_manifest(out, option, policy, kernel, completed=False)
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
    "PersistentSiteOptionAdapterV100514",
    "persistent_site_policy",
    "main",
]
