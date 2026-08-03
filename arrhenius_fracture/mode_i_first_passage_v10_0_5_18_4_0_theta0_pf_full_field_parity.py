"""v10.0.5.18.4.0 theta=0 control with PF v10.4.1 bulk parity.

The moving-tip persistent-site MPZ, stochastic cleavage renewal, event length,
and v3.9 atomic CZM geometry are unchanged. This overlay replaces only the
surrounding continuum-bulk policy and Peierls--Taylor kinetic module with the
exact v10.4.1 detailed-balance formulation used by the PF reference case.
"""
from __future__ import annotations

from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import sys
from typing import Any, Mapping

from . import emission_derived_plasticity as _installed_bulk_pt
from . import mode_i_first_passage_v10_0_5_13_barrier_only as _core
from . import mode_i_first_passage_v10_0_5_13_1_barrier_only as _v131
from . import mode_i_first_passage_v10_0_5_13_2_barrier_only as _v132
from . import mode_i_first_passage_v10_0_5_13_4_barrier_only as _v134
from . import mode_i_first_passage_v10_0_5_14_persistent_site as _v14
from . import mode_i_first_passage_v10_0_5_18_3_9_atomic_path_corridor as _base
from . import mode_i_first_passage_v10_0_5_18_4_0_theta0_pf_parity as _control
from .bulk_pt_detailed_balance_v1041_exact import (
    EmissionDerivedPeierlsTaylorModel as ExactBulkPTModel,
    KB_EV_PER_K,
    MODEL_ID as EXACT_BULK_MODEL,
    config_from_dislocation_config as exact_bulk_config,
)

POINT_RELEASE = "10.0.5.18.4.0"
MODEL_ID = "FEM_CZM_theta0_PF_v10_4_1_full_field_bulk_parity_v10_0_5_18_4_0"
AUDIT_FILE = "theta0_pf_full_field_parity_contract_v10_0_5_18_4_0.json"
SEMANTIC_BULK_MODE = "full_field"
SOLVER_BULK_MODE = "bulk_same_pt_km"


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _exact_bulk_parameters(cfg: Any, row: Mapping[str, Any]) -> None:
    emit0 = max(float(row["emit_G00_eV"]), 1.0e-30)
    emit_gT = float(row["emit_gT_eV_per_K"])
    peierls_gT = -float(row["peierls_activation_entropy_kB"]) * KB_EV_PER_K
    taylor_gT = -float(row["taylor_activation_entropy_kB"]) * KB_EV_PER_K

    values = {
        "enable_plasticity": True,
        "use_emission_derived_pt": True,
        "bulk_kinetics_model": "emission_derived_peierls_taylor_multihit",
        "bulk_kinetics_model_detail": EXACT_BULK_MODEL,
        "plastic_update_mode": "explicit_rate",
        "thermo_consistency_mode": "time_cone",
        "bulk_mult_frac": 1.0,
        "tip_source_rho_per_emit": 0.0,
        "rho_transport_c": 0.0,
        "exhaustion_enabled": False,
        "freeze_rho": False,
        "use_static_recovery": False,
        "mobile_rho_floor": float(row["rho_forest_floor_m2"]),
        "pt_emit_G00_eV": emit0,
        "pt_emit_gT_eV_per_K": emit_gT,
        "pt_emit_sigc0_Pa": float(row["emit_sigc0_GPa"]) * 1.0e9,
        "pt_emit_sT_Pa_per_K": float(row["emit_sT_GPa_per_K"]) * 1.0e9,
        "pt_emit_exp_a": float(row["emit_exp_a"]),
        "pt_emit_exp_n": float(row["emit_exp_n"]),
        "pt_emit_floor_frac": float(row["emit_floor_frac"]),
        "pt_emit_floor_min_eV": 1.0e-4,
        "pt_emit_floor_max_frac": 0.95,
        "pt_emit_Tref_K": float(row["Tref_K"]),
        "pt_peierls_energy_ratio": float(row["peierls_H0_eV"]) / emit0,
        "pt_peierls_entropy_ratio": (
            peierls_gT / emit_gT if abs(emit_gT) > 1.0e-30 else 0.0
        ),
        "pt_peierls_stress_ratio": float(row["peierls_stress_fraction"]),
        "pt_peierls_nu0_s": float(row["peierls_nu0_s"]),
        "pt_peierls_G00_eV": float(row["peierls_H0_eV"]),
        "pt_peierls_gT_eV_per_K": peierls_gT,
        "pt_peierls_exp_a": float(row["peierls_exp_a"]),
        "pt_peierls_exp_n": float(row["peierls_exp_n"]),
        "pt_peierls_Tref_K": float(row["Tref_K"]),
        "pt_taylor_energy_ratio": float(row["taylor_H0_eV"]) / emit0,
        "pt_taylor_entropy_ratio": (
            taylor_gT / emit_gT if abs(emit_gT) > 1.0e-30 else 0.0
        ),
        "pt_taylor_stress_ratio": float(row["taylor_stress_fraction"]),
        "pt_taylor_nu0_s": float(row["taylor_nu0_s"]),
        "pt_taylor_G00_eV": float(row["taylor_H0_eV"]),
        "pt_taylor_gT_eV_per_K": taylor_gT,
        "pt_taylor_exp_a": float(row["taylor_exp_a"]),
        "pt_taylor_exp_n": float(row["taylor_exp_n"]),
        "pt_taylor_Tref_K": float(row["Tref_K"]),
        "pt_taylor_corr_rho_c": float(row["taylor_corr_rho_c_m2"]),
        "pt_taylor_renewal_time_s": 1.0e-9,
        "pt_taylor_m_exponent": 1.0,
        "pt_taylor_m_scale": float(row["taylor_corr_scale"]),
        "pt_taylor_m_cap": float("inf"),
        "pt_mobile_fraction": 0.01,
        "pt_mobile_saturation_density_m2": 1.0e14,
        "pt_mobile_density_floor_m2": float(row["rho_forest_floor_m2"]),
        "pt_jump_fraction": 1.0,
        "pt_jump_length_min_m": 2.5e-10,
        "pt_taylor_phi_max": 20.0,
    }
    for name, value in values.items():
        setattr(cfg, name, value)


def _exact_full_field_defaults(
    args: Any,
    row: Mapping[str, Any],
    bulk_mode: str,
) -> None:
    if str(bulk_mode) not in {SEMANTIC_BULK_MODE, SOLVER_BULK_MODE}:
        raise RuntimeError(
            "PF full-field parity requires the inherited full-bulk solver mode"
        )
    _exact_bulk_parameters(args, row)
    values = {
        "front_state_model": "moving_pz",
        "pz_store_to_rho_scale": 0.0,
        "tip_source_rho_per_emit": 0.0,
        "bulk_mult_frac": 1.0,
        "exhaustion": False,
        "bulk_plasticity_mode_v911": SEMANTIC_BULK_MODE,
        "bulk_plasticity_solver_token": SOLVER_BULK_MODE,
        "rho0": float(row["rho_forest_floor_m2"]),
    }
    for name, value in values.items():
        setattr(args, name, value)


def _full_field_policy(candidate) -> dict[str, Any]:
    policy = dict(_full_field_policy._original(candidate))
    policy.update(
        {
            "policy_id": "PF_v10_4_1_full_field_bulk_plus_persistent_tip_MPZ_v10051840",
            "bulk_plasticity_mode": SOLVER_BULK_MODE,
            "bulk_plasticity_semantic_mode": SEMANTIC_BULK_MODE,
            "bulk_plasticity_solver_token": SOLVER_BULK_MODE,
            "continuum_bulk_role": "exact_selected_row_Peierls_Taylor_detailed_balance",
            "bulk_state_evolves_in_fem": True,
            "uniform_bulk_mobile_retained_state_active": True,
            "tip_and_bulk_populations_distinct": True,
            "direct_tip_to_bulk_density_transfer": False,
            "bulk_mult_frac": 1.0,
            "bulk_density_transport": False,
            "bulk_static_recovery": False,
            "bulk_detailed_balance": True,
            "bulk_model": EXACT_BULK_MODEL,
        }
    )
    return policy


def _replace_option_full_field(argv: list[str], name: str, value: str) -> None:
    forced = SOLVER_BULK_MODE if name == "--bulk-plasticity-mode" else value
    _replace_option_full_field._original(argv, name, forced)


def _normalize_solver_args(argv: list[str]) -> list[str]:
    normalized = list(argv)
    name = "--bulk-plasticity-mode"
    prefix = name + "="
    for index, token in enumerate(normalized):
        if token == name and index + 1 < len(normalized):
            normalized[index + 1] = SOLVER_BULK_MODE
            return normalized
        if token.startswith(prefix):
            normalized[index] = prefix + SOLVER_BULK_MODE
            return normalized
    normalized.extend([name, SOLVER_BULK_MODE])
    return normalized


def _bypass_tip_only_policy_wrappers(argv: list[str] | None = None):
    """Retain v13.2 startup and v13.5 mesh repairs without tip-only vetoes."""
    return _v132.main(list(sys.argv[1:] if argv is None else argv))


def _full_fields() -> dict[str, Any]:
    fields = dict(_full_fields._original())
    fields.update(
        {
            "bulk_plasticity_mode_required": SEMANTIC_BULK_MODE,
            "bulk_plasticity_mode": SEMANTIC_BULK_MODE,
            "bulk_plasticity_solver_token": SOLVER_BULK_MODE,
            "bulk_state_evolves_in_fem": True,
            "bulk_model": EXACT_BULK_MODEL,
            "bulk_net_slip": "Gamma_forward_minus_Gamma_reverse",
            "zero_stress_net_plastic_rate_exactly_zero": True,
            "tip_and_bulk_populations_distinct": True,
            "direct_tip_to_bulk_density_transfer": False,
            "bulk_density_transport": False,
            "bulk_static_recovery": False,
            "PF_v10_4_1_bulk_parity_active": True,
        }
    )
    return fields


def _contract(
    argv: list[str], kernel: Path, completed: bool, error: BaseException | None
) -> dict[str, Any]:
    return {
        "schema": MODEL_ID,
        "point_release": POINT_RELEASE,
        "created_or_updated_utc": _utc_now(),
        "run_completed_without_exception": bool(completed),
        "runtime_error_type": None if error is None else type(error).__name__,
        "runtime_error": None if error is None else str(error),
        "reference_case": _control.REFERENCE_CASE,
        "parameter_option": _control.REFERENCE_OPTION,
        "temperature_K": 1000,
        "theta_deg": 0.0,
        "hazard_seed": _control.REFERENCE_HAZARD_SEED,
        "signed_kernel_family": str(kernel),
        "signed_kernel_family_sha256": _control._sha256(kernel),
        "matched_numerical_controls": {
            "nx": 36,
            "ny": 72,
            "dU_m": 2.0e-7,
            "dt_s": 8.4,
            "n_stagger": 2,
            "tip_h_fine_m": 1.0e-6,
            "tip_ratio": 1.20,
            "da_phys_m": 5.0e-6,
            "adaptive_event_target": 0.15,
        },
        "bulk_parity": {
            "PF_bulk_plasticity_mode": SEMANTIC_BULK_MODE,
            "FEM_bulk_plasticity_mode": SEMANTIC_BULK_MODE,
            "inherited_solver_bulk_mode_token": SOLVER_BULK_MODE,
            "model": EXACT_BULK_MODEL,
            "exact_selected_registry_row": True,
            "initial_bulk_density_from_rho_forest_floor": True,
            "Peierls_Taylor_parameters_from_exact_row": True,
            "detailed_balance_forward_minus_reverse": True,
            "zero_stress_net_rate_exactly_zero": True,
            "bulk_mult_frac": 1.0,
            "tip_source_rho_per_emit": 0.0,
            "bulk_density_transport": False,
            "tip_and_bulk_populations_distinct": True,
        },
        "remaining_model_form_difference": {
            "PF_crack_backend": "sharp_wake",
            "FEM_crack_backend": "adaptive_czm_atomic_path_corridor",
            "adaptive_remeshing_active": True,
        },
        "tip_physics_changed": False,
        "cleavage_hazard_changed": False,
        "event_length_changed": False,
        "argv": list(argv),
    }


def _augment_manifest(out: Path, contract: dict[str, Any]) -> None:
    path = out / _base.PRODUCTION_MANIFEST
    if not path.is_file():
        return
    try:
        payload = json.loads(path.read_text())
    except Exception:
        return
    payload["theta0_pf_full_field_parity"] = contract
    payload["point_release_overlay"] = POINT_RELEASE
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n")


def main(argv: list[str] | None = None):
    args = list(sys.argv[1:] if argv is None else argv)
    out_raw = _control._option_value(args, "--out")
    if out_raw is None:
        raise SystemExit(f"{MODEL_ID} requires --out")
    out = Path(out_raw).expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)

    _control._require_text(args, "--parameter-option", _control.REFERENCE_OPTION)
    _control._require_text(args, "--temperatures", "1000")
    _control._require_text(args, "--bulk-plasticity-mode", SEMANTIC_BULK_MODE)
    _control._require_text(args, "--j-decomposition", "cluster")
    _control._require_int(args, "--nx", 36)
    _control._require_int(args, "--ny", 72)
    _control._require_int(args, "--n-stagger", 2)
    _control._require_int(args, "--max-fronts", 1)
    _control._require_float(args, "--crystal-theta-deg", 0.0)
    _control._require_float(args, "--dU", 2.0e-7)
    _control._require_float(args, "--dt", 8.4)
    _control._require_float(args, "--tip-h-fine", 1.0e-6)
    _control._require_float(args, "--tip-ratio", 1.20)
    _control._require_float(args, "--da-phys", 5.0e-6)
    _control._require_float(args, "--adaptive-event-target", 0.15)

    seed_raw = os.environ.get("CLEAVAGE_HAZARD_SEED", "").strip()
    if seed_raw != str(_control.REFERENCE_HAZARD_SEED):
        raise SystemExit(
            f"requires CLEAVAGE_HAZARD_SEED={_control.REFERENCE_HAZARD_SEED}; "
            f"observed {seed_raw!r}"
        )
    kernel_raw = _control._option_value(args, "--signed-kernel-family")
    if kernel_raw is None:
        raise SystemExit(f"{MODEL_ID} requires --signed-kernel-family")
    kernel = Path(kernel_raw).expanduser().resolve()
    if not kernel.is_file():
        raise SystemExit(f"PF reference kernel is missing: {kernel}")
    observed_sha = _control._sha256(kernel)
    if observed_sha != _control.REFERENCE_KERNEL_SHA256:
        raise SystemExit(
            "PF reference kernel SHA-256 mismatch: "
            f"expected {_control.REFERENCE_KERNEL_SHA256}, observed {observed_sha}"
        )

    solver_args = _normalize_solver_args(args)
    saved = {
        "v14_policy": _v14.persistent_site_policy,
        "v14_replace": _v14._replace_option,
        "v134_main": _v134.main,
        "core_bulk_defaults": _core._set_bulk_barrier_defaults,
        "v131_apply": _v131._apply_barrier_pt_config,
        "bulk_model": _installed_bulk_pt.EmissionDerivedPeierlsTaylorModel,
        "bulk_config": _installed_bulk_pt.config_from_dislocation_config,
        "base_fields": _base._fields,
    }
    _full_field_policy._original = saved["v14_policy"]
    _replace_option_full_field._original = saved["v14_replace"]
    _full_fields._original = saved["base_fields"]

    _v14.persistent_site_policy = _full_field_policy
    _v14._replace_option = _replace_option_full_field
    _v134.main = _bypass_tip_only_policy_wrappers
    _core._set_bulk_barrier_defaults = _exact_full_field_defaults
    _v131._apply_barrier_pt_config = _exact_bulk_parameters
    _installed_bulk_pt.EmissionDerivedPeierlsTaylorModel = ExactBulkPTModel
    _installed_bulk_pt.config_from_dislocation_config = exact_bulk_config
    _base._fields = _full_fields

    completed = False
    error: BaseException | None = None
    contract = _contract(args, kernel, completed=False, error=None)
    (out / AUDIT_FILE).write_text(
        json.dumps(contract, indent=2, sort_keys=True, default=str) + "\n"
    )
    try:
        result = _base.main(solver_args)
        completed = True
        return result
    except BaseException as exc:
        error = exc
        raise
    finally:
        _v14.persistent_site_policy = saved["v14_policy"]
        _v14._replace_option = saved["v14_replace"]
        _v134.main = saved["v134_main"]
        _core._set_bulk_barrier_defaults = saved["core_bulk_defaults"]
        _v131._apply_barrier_pt_config = saved["v131_apply"]
        _installed_bulk_pt.EmissionDerivedPeierlsTaylorModel = saved["bulk_model"]
        _installed_bulk_pt.config_from_dislocation_config = saved["bulk_config"]
        _base._fields = saved["base_fields"]
        contract = _contract(args, kernel, completed=completed, error=error)
        (out / AUDIT_FILE).write_text(
            json.dumps(contract, indent=2, sort_keys=True, default=str) + "\n"
        )
        _augment_manifest(out, contract)


if __name__ == "__main__":
    main()


__all__ = [
    "AUDIT_FILE",
    "MODEL_ID",
    "POINT_RELEASE",
    "SEMANTIC_BULK_MODE",
    "SOLVER_BULK_MODE",
    "main",
]
