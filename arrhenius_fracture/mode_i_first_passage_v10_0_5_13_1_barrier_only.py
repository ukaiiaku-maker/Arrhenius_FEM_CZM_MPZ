"""v10.0.5.13.1 barrier-only point release.

The v10.0.5.13 entry correctly blocked candidate-specific 1-D state fields, but
its common-state adapter still reassigned source, refresh, recovery, encounter,
correlation, shielding, and blunting values to nominal defaults.  This point
release removes those assignments.  It patches only the barrier installation
hooks and leaves every non-barrier field as constructed by the existing full
2-D solver and explicit CLI.
"""
from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path
import sys
from typing import Any, Mapping

from . import mode_i_first_passage_v10_0_5_13_barrier_only as _base
from .barrier_only_response_registry_v100513 import TWO_D_STATE_POLICY
from .mode_i_first_passage_v10_0 import _option_value
from .moving_process_zone import build_mpz_config_from_namespace

POINT_RELEASE = "10.0.5.13.1"
MODEL_ID = "FEM_CZM_full_2D_barrier_only_preserved_state_v10_0_5_13_1"
PRODUCTION_MANIFEST = _base.PRODUCTION_MANIFEST


def _apply_barrier_pt_config(cfg: Any, row: Mapping[str, Any]) -> None:
    """Install only emission/Peierls/Taylor barrier-surface parameters."""
    emit0 = max(float(row["emit_G00_eV"]), 1.0e-30)
    values = {
        "use_emission_derived_pt": True,
        "bulk_kinetics_model": "emission_derived_peierls_taylor_multihit",
        "bulk_kinetics_model_detail": "v1005131_barriers_only_preserved_2d_state",
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
        "pt_peierls_nu0_s": float(row["peierls_nu0_s"]),
        "pt_taylor_energy_ratio": float(row["taylor_H0_eV"]) / emit0,
        "pt_taylor_entropy_ratio": float(row["taylor_activation_entropy_kB"]),
        "pt_taylor_activation_entropy_kB": float(
            row["taylor_activation_entropy_kB"]
        ),
        "pt_taylor_exp_a": float(row["taylor_exp_a"]),
        "pt_taylor_exp_n": float(row["taylor_exp_n"]),
        "pt_taylor_nu0_s": float(row["taylor_nu0_s"]),
    }
    for name, value in values.items():
        setattr(cfg, name, value)


def _build_barrier_only_mpz_config(args: Any, row: Mapping[str, Any]):
    """Build the existing MPZ state, changing only resolution and barriers."""
    length_m = float(TWO_D_STATE_POLICY["mpz_length_um"]) * 1.0e-6
    cfg = build_mpz_config_from_namespace(args, default_length_m=length_m)
    cfg.length_m = length_m
    cfg.n_bins = int(TWO_D_STATE_POLICY["mpz_n_bins"])
    _apply_barrier_pt_config(cfg, row)
    return cfg


def _update_manifest(out: Path, completed: bool) -> None:
    path = out / PRODUCTION_MANIFEST
    if not path.is_file():
        return
    try:
        payload = json.loads(path.read_text())
    except Exception:
        return
    payload.update(
        {
            "schema": "barrier_only_production_manifest_v10_0_5_13_1",
            "model": MODEL_ID,
            "point_release": POINT_RELEASE,
            "nonbarrier_state_assignment_repair": {
                "active": True,
                "candidate_state_fields_applied": False,
                "common_state_values_reassigned": False,
                "state_configuration_source": (
                    "existing_full_2d_solver_and_explicit_cli"
                ),
                "barrier_hooks_changed_only": True,
                "physics_removed_or_simplified": False,
                "recorded_utc": datetime.now(timezone.utc).isoformat(),
            },
        }
    )
    if completed:
        payload["status"] = "complete"
        payload["run_completed_without_exception"] = True
    path.write_text(json.dumps(payload, indent=2, default=str))


def main(argv: list[str] | None = None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    out_value = _option_value(user_args, "--out")
    if out_value is None:
        raise SystemExit("v10.0.5.13.1 requires --out")
    out = Path(out_value).resolve()

    saved_apply = _base._apply_barrier_pt_config
    saved_build = _base._build_barrier_only_mpz_config
    _base._apply_barrier_pt_config = _apply_barrier_pt_config
    _base._build_barrier_only_mpz_config = _build_barrier_only_mpz_config
    try:
        result = _base.main(user_args)
        _update_manifest(out, completed=True)
        return result
    except BaseException:
        _update_manifest(out, completed=False)
        raise
    finally:
        _base._apply_barrier_pt_config = saved_apply
        _base._build_barrier_only_mpz_config = saved_build


if __name__ == "__main__":
    main()


__all__ = [
    "POINT_RELEASE",
    "MODEL_ID",
    "PRODUCTION_MANIFEST",
    "_apply_barrier_pt_config",
    "_build_barrier_only_mpz_config",
    "main",
]
