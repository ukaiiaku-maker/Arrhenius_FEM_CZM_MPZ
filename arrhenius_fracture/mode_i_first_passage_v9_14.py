"""v9.14 Mode-I entry point: event-localized remesh and same-load equilibrium."""
from __future__ import annotations

import csv
import json
import os
from pathlib import Path
import sys
from typing import Any

from . import crack_backend as _crack_backend
from . import fem as _fem
from . import mixed_mode_first_passage_v9_11 as _mm_v911
from . import mode_i_first_passage_v9_13 as _base_v913
from .event_equilibrium_v914 import (
    ACTIVE_CONTEXT,
    install_mechanics_recorder,
    restore_mechanics_recorder,
)
from .event_remesh_czm_v914 import build_event_remesh_backend

MODEL_ID = "FEM_CZM_Mode_I_MPZ_v9_14_event_remesh_same_load_equilibrium"


def _option_value(argv: list[str], name: str) -> str | None:
    for i, token in enumerate(argv):
        if token == name and i + 1 < len(argv):
            return argv[i + 1]
        if token.startswith(name + "="):
            return token.split("=", 1)[1]
    return None


def _recording_j_factory(original_factory):
    def factory(original_compute, context, mpz_args):
        base = original_factory(original_compute, context, mpz_args)

        def wrapped(
            mesh,
            u,
            sigma_gp,
            psi_e_gp,
            d,
            crack_tip,
            crack_direction,
            mat,
            ell,
            cfg=None,
            crack_segments=None,
            exclude_radius=0.0,
        ):
            J, KJ, info = base(
                mesh,
                u,
                sigma_gp,
                psi_e_gp,
                d,
                crack_tip,
                crack_direction,
                mat,
                ell,
                cfg=cfg,
                crack_segments=crack_segments,
                exclude_radius=exclude_radius,
            )
            ACTIVE_CONTEXT.record_j_call(
                base,
                mat,
                ell,
                cfg,
                exclude_radius,
                J,
                KJ,
            )
            return J, KJ, info

        return wrapped

    return factory


def _event_backend_factory(original_factory):
    def build(args, geom):
        enabled = os.environ.get("ARRHENIUS_EVENT_REMESH_V914", "1") != "0"
        if enabled:
            return build_event_remesh_backend(args, geom)
        return original_factory(args, geom)

    return build


def _write_equilibrium_audit(out: Path) -> None:
    records = list(ACTIVE_CONTEXT.records)
    finite_j = [
        r for r in records
        if str(r.get("J_after_event_status", "")) == "ok"
    ]
    payload: dict[str, Any] = {
        "schema": "event_equilibrium_v914",
        "model": MODEL_ID,
        "same_time_same_load_protocol": True,
        "physical_time_advanced_during_equilibrium": False,
        "hazard_action_advanced_during_equilibrium": False,
        "n_post_event_equilibria": len(records),
        "records": records,
        "all_same_time": bool(
            records and all(float(r.get("physical_time_increment_s", 1.0)) == 0.0 for r in records)
        ),
        "all_zero_hazard_increment": bool(
            records and all(float(r.get("hazard_action_increment", 1.0)) == 0.0 for r in records)
        ),
        "all_J_recomputed": bool(records and len(finite_j) == len(records)),
        "max_relative_boundary_displacement_drift": max(
            (float(r.get("max_relative_boundary_displacement_drift", float("nan"))) for r in records),
            default=float("nan"),
        ),
        "max_relative_rho_area_integral_error": max(
            (float(r.get("relative_rho_area_integral_error", float("nan"))) for r in records),
            default=float("nan"),
        ),
        "max_relative_ep_area_integral_error": max(
            (float(r.get("max_relative_ep_area_integral_error", float("nan"))) for r in records),
            default=float("nan"),
        ),
        "max_relative_total_mesh_area_error": max(
            (float(r.get("relative_total_mesh_area_error", float("nan"))) for r in records),
            default=float("nan"),
        ),
    }
    out.mkdir(parents=True, exist_ok=True)
    (out / "event_equilibrium_audit_v914.json").write_text(
        json.dumps(payload, indent=2, default=str)
    )
    if records:
        flat = []
        for row in records:
            flat.append({
                key: (json.dumps(value) if isinstance(value, (list, dict)) else value)
                for key, value in row.items()
            })
        with (out / "event_equilibrium_audit_v914.csv").open("w", newline="") as fp:
            writer = csv.DictWriter(fp, fieldnames=sorted({k for row in flat for k in row}))
            writer.writeheader()
            writer.writerows(flat)


def main(argv: list[str] | None = None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    bulk_mode = (_option_value(user_args, "--bulk-plasticity-mode") or "tip_only").strip()
    if bulk_mode != "tip_only":
        raise SystemExit(
            "v9.14 event-remesh validation is initially restricted to tip_only; "
            "bulk_same_pt_km remains available on v9.13 until its full remesh state "
            "transfer is validated independently"
        )
    backend = (_option_value(user_args, "--crack-backend") or "adaptive_czm").strip()
    if backend not in {"adaptive_czm", "event_remesh_czm"}:
        raise SystemExit(
            "v9.14 event-remesh entry point requires --crack-backend adaptive_czm "
            "(the accepted parent CLI token) or event_remesh_czm"
        )

    ACTIVE_CONTEXT.clear()
    original_assemble = install_mechanics_recorder(_fem)
    original_backend_factory = _crack_backend.build_crack_backend
    original_j_factory = _mm_v911._j_profile_wrapper_factory
    _crack_backend.build_crack_backend = _event_backend_factory(original_backend_factory)
    _mm_v911._j_profile_wrapper_factory = _recording_j_factory(original_j_factory)
    try:
        results = _base_v913.main(user_args)
    finally:
        _mm_v911._j_profile_wrapper_factory = original_j_factory
        _crack_backend.build_crack_backend = original_backend_factory
        restore_mechanics_recorder(_fem, original_assemble)

    out_value = _option_value(user_args, "--out")
    if out_value is not None:
        out = Path(out_value)
        _write_equilibrium_audit(out)
        summary = out / "anisotropic_calibrated_tip_first_passage_summary.json"
        if summary.exists():
            try:
                payload = json.loads(summary.read_text())
                payload.update({
                    "model_v914": MODEL_ID,
                    "effective_crack_backend": "event_remesh_czm",
                    "one_physical_event_per_hazard_renewal": True,
                    "conservative_parent_map_transfer": True,
                    "post_event_same_time_same_load_equilibrium": True,
                    "n_post_event_equilibria": len(ACTIVE_CONTEXT.records),
                })
                summary.write_text(json.dumps(payload, indent=2, default=str))
            except Exception:
                pass
    print("MODE_I_MPZ_V9_14_EVENT_REMESH complete")
    return results


if __name__ == "__main__":
    main()
