"""v10.0.5.18.4.0: theta=0 Peak-0242980 PF/FEM parity control.

This entry does not change the v10.0.5.18.3.9 atomic path-corridor geometry or
front-local constitutive law.  It locks every numerical and tip-state control
that can currently be made identical to the completed PF v10.4.1 reference
case ``T1000K_th0_seed8666`` and records the remaining model-form mismatch:
the PF reference has full-field bulk Peierls--Taylor plasticity, while this
FEM/CZM stack remains the qualified tip-only elastic-bulk formulation.
"""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import sys
from typing import Any

from . import mode_i_first_passage_v10_0_5_18_3_9_atomic_path_corridor as _base

POINT_RELEASE = "10.0.5.18.4.0"
MODEL_ID = "FEM_CZM_theta0_PF_v10_4_1_numerical_parity_control_v10_0_5_18_4_0"
AUDIT_FILE = "theta0_pf_parity_contract_v10_0_5_18_4_0.json"
REFERENCE_OPTION = "v913_paper_peak01_0242980_persistent_sites"
REFERENCE_CASE = "T1000K_th0_seed8666"
REFERENCE_HAZARD_SEED = 8666
REFERENCE_KERNEL_SHA256 = "a85b57ad9eee8331ef34ea222760ce7d4064f48ddef668f1e28a666d14946f3a"


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _option_value(argv: list[str], name: str) -> str | None:
    prefix = name + "="
    for index, token in enumerate(argv):
        if token.startswith(prefix):
            return token[len(prefix) :]
        if token == name and index + 1 < len(argv):
            return argv[index + 1]
    return None


def _require_text(argv: list[str], name: str, expected: str) -> None:
    observed = _option_value(argv, name)
    if observed != expected:
        raise SystemExit(
            f"{MODEL_ID} requires {name}={expected!r}; observed {observed!r}"
        )


def _require_int(argv: list[str], name: str, expected: int) -> None:
    observed = _option_value(argv, name)
    try:
        value = int(observed) if observed is not None else None
    except ValueError as exc:
        raise SystemExit(f"{name} must be an integer") from exc
    if value != int(expected):
        raise SystemExit(
            f"{MODEL_ID} requires {name}={expected}; observed {observed!r}"
        )


def _require_float(
    argv: list[str], name: str, expected: float, *, abs_tol: float = 1.0e-15
) -> None:
    observed = _option_value(argv, name)
    try:
        value = float(observed) if observed is not None else None
    except ValueError as exc:
        raise SystemExit(f"{name} must be numeric") from exc
    if value is None or not math.isclose(
        value, float(expected), rel_tol=1.0e-12, abs_tol=abs_tol
    ):
        raise SystemExit(
            f"{MODEL_ID} requires {name}={expected:.16g}; observed {observed!r}"
        )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _contract(argv: list[str], kernel: Path, completed: bool, error: BaseException | None):
    return {
        "schema": MODEL_ID,
        "point_release": POINT_RELEASE,
        "created_or_updated_utc": _utc_now(),
        "run_completed_without_exception": bool(completed),
        "runtime_error_type": None if error is None else type(error).__name__,
        "runtime_error": None if error is None else str(error),
        "reference": {
            "PF_model_entry": "arrhenius_fracture.sharp_front_v10_4_bulk_peierls_taylor_audited",
            "PF_release": "v10.4.1-bulk-detailed-balance",
            "case": REFERENCE_CASE,
            "parameter_option": REFERENCE_OPTION,
            "temperature_K": 1000,
            "theta_deg": 0.0,
            "hazard_seed": REFERENCE_HAZARD_SEED,
            "target_projected_extension_um": 1000.0,
        },
        "matched_controls": {
            "parameter_option": REFERENCE_OPTION,
            "temperature_K": 1000,
            "theta_deg": 0.0,
            "hazard_seed": REFERENCE_HAZARD_SEED,
            "nx": 36,
            "ny": 72,
            "dU_m": 2.0e-7,
            "dt_s": 8.4,
            "nominal_opening_rate_m_per_s": 2.0e-7 / 8.4,
            "n_stagger": 2,
            "tip_h_fine_m": 1.0e-6,
            "tip_ratio": 1.20,
            "da_phys_m": 5.0e-6,
            "adaptive_event_target": 0.15,
            "adaptive_min_frac": 1.0e-8,
            "adaptive_grow": 4.0,
            "crystal_C11_Pa": 523.0e9,
            "crystal_C12_Pa": 203.0e9,
            "crystal_C44_Pa": 160.0e9,
            "j_decomposition": "cluster",
            "maximum_fronts": 1,
            "wake_shielding": False,
            "mobile_shield_fraction": 0.0,
            "signed_kernel_family": str(kernel),
            "signed_kernel_family_sha256": _sha256(kernel),
        },
        "intentional_model_form_differences": {
            "PF_crack_backend": "sharp_wake",
            "FEM_crack_backend": "adaptive_czm_atomic_path_corridor",
            "PF_bulk_plasticity_mode": "full_field",
            "FEM_bulk_plasticity_mode": "tip_only",
            "PF_bulk_net_slip": "detailed_balance_forward_minus_reverse",
            "FEM_bulk_role": "elastic_continuum_with_moving_tip_MPZ_only",
            "bulk_plasticity_parity_complete": False,
            "theta0_run_is_remeshing_and_CZM_control_not_full_constitutive_parity": True,
        },
        "physics_guards": {
            "Peak_0242980_exact_row_required": True,
            "exact_PF_kernel_required": True,
            "exact_seed_required": True,
            "exact_stochastic_event_endpoint_preserved": True,
            "exact_stochastic_event_length_preserved": True,
            "triangle_quality_floor_relaxed": False,
            "area_ratio_floor_relaxed": False,
            "hazard_law_changed": False,
            "tip_MPZ_law_changed": False,
        },
        "argv": list(argv),
    }


def _augment_base_manifest(out: Path, contract: dict[str, Any]) -> None:
    path = out / _base.PRODUCTION_MANIFEST
    if not path.is_file():
        return
    try:
        payload = json.loads(path.read_text())
    except Exception:
        return
    payload["theta0_pf_parity_control"] = contract
    payload["point_release_overlay"] = POINT_RELEASE
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n")


def main(argv: list[str] | None = None):
    args = list(sys.argv[1:] if argv is None else argv)
    out_raw = _option_value(args, "--out")
    if out_raw is None:
        raise SystemExit(f"{MODEL_ID} requires --out")
    out = Path(out_raw).expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)

    _require_text(args, "--parameter-option", REFERENCE_OPTION)
    _require_text(args, "--temperatures", "1000")
    _require_text(args, "--bulk-plasticity-mode", "tip_only")
    _require_text(args, "--j-decomposition", "cluster")
    _require_int(args, "--nx", 36)
    _require_int(args, "--ny", 72)
    _require_int(args, "--n-stagger", 2)
    _require_int(args, "--max-fronts", 1)
    _require_float(args, "--crystal-theta-deg", 0.0)
    _require_float(args, "--dU", 2.0e-7)
    _require_float(args, "--dt", 8.4)
    _require_float(args, "--tip-h-fine", 1.0e-6)
    _require_float(args, "--tip-ratio", 1.20)
    _require_float(args, "--da-phys", 5.0e-6)
    _require_float(args, "--adaptive-event-target", 0.15)
    _require_float(args, "--adaptive-min-frac", 1.0e-8)
    _require_float(args, "--adaptive-grow", 4.0)

    seed_raw = os.environ.get("CLEAVAGE_HAZARD_SEED", "").strip()
    if seed_raw != str(REFERENCE_HAZARD_SEED):
        raise SystemExit(
            f"{MODEL_ID} requires CLEAVAGE_HAZARD_SEED={REFERENCE_HAZARD_SEED}; "
            f"observed {seed_raw!r}"
        )
    if os.environ.get("CLEAVAGE_HAZARD_MODE", "").strip().lower() != "exponential":
        raise SystemExit(f"{MODEL_ID} requires CLEAVAGE_HAZARD_MODE=exponential")
    if (
        os.environ.get("CLEAVAGE_EVENT_LENGTH_MODE", "").strip().lower()
        != "threshold_scaled"
    ):
        raise SystemExit(
            f"{MODEL_ID} requires CLEAVAGE_EVENT_LENGTH_MODE=threshold_scaled"
        )

    kernel_raw = _option_value(args, "--signed-kernel-family")
    if kernel_raw is None:
        raise SystemExit(f"{MODEL_ID} requires --signed-kernel-family")
    kernel = Path(kernel_raw).expanduser().resolve()
    if not kernel.is_file():
        raise SystemExit(f"PF reference kernel is missing: {kernel}")
    observed_sha = _sha256(kernel)
    if observed_sha != REFERENCE_KERNEL_SHA256:
        raise SystemExit(
            "PF reference kernel SHA-256 mismatch: "
            f"expected {REFERENCE_KERNEL_SHA256}, observed {observed_sha}"
        )

    error: BaseException | None = None
    completed = False
    contract = _contract(args, kernel, completed=False, error=None)
    (out / AUDIT_FILE).write_text(
        json.dumps(contract, indent=2, sort_keys=True, default=str) + "\n"
    )
    try:
        result = _base.main(args)
        completed = True
        return result
    except BaseException as exc:
        error = exc
        raise
    finally:
        contract = _contract(args, kernel, completed=completed, error=error)
        (out / AUDIT_FILE).write_text(
            json.dumps(contract, indent=2, sort_keys=True, default=str) + "\n"
        )
        _augment_base_manifest(out, contract)


if __name__ == "__main__":
    main()


__all__ = [
    "AUDIT_FILE",
    "MODEL_ID",
    "POINT_RELEASE",
    "REFERENCE_CASE",
    "REFERENCE_HAZARD_SEED",
    "REFERENCE_KERNEL_SHA256",
    "REFERENCE_OPTION",
    "main",
]
