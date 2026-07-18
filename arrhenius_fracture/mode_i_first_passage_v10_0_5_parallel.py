"""v10.0.5 Mode-I entry point for parallel crack-opening/emission coupling.

The certified v10.0.3 progressive lifecycle is reused unchanged. This wrapper
replaces only the live engine factory, J-observer, and cohesive stepper used
inside that lifecycle.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
import sys
from typing import Any

from . import j_integral as _j_integral
from . import kinetic_progressive_2d_v1002 as _v1002
from . import mode_i_first_passage_v10_0_3_progressive as _v1003
from .kinetic_campaign_czm_v1005 import engine_factory_v1005
from .tensor_resolved_coupling_v1005 import (
    TensorResolvedDriveConfig,
    TensorResolvedKineticCohesiveStepper,
    make_tensor_resolved_J_wrapper,
    reset_tensor_drive_runtime,
    tensor_drive_runtime_payload,
)

MODEL_ID = "FEM_CZM_Mode_I_parallel_opening_tensor_emission_v10_0_5"


def _option_value(args: list[str], name: str) -> str | None:
    for i, value in enumerate(args):
        if value == name and i + 1 < len(args):
            return args[i + 1]
        prefix = name + "="
        if value.startswith(prefix):
            return value[len(prefix):]
    return None


def _write_v1005_results(
    out: Path,
    results: Any,
    audit_name: str,
) -> Path:
    rows = []
    for row in results or []:
        item = dict(row)
        item.update({
            "model": MODEL_ID,
            "point_release": "10.0.5",
            "front_state_model": "kinetic_campaign_czm",
            "front_state_model_detail": (
                "pf_v10_1_7_1_parallel_opening_tensor_emission_reset_safe_v1005"
            ),
            "tensor_resolved_parallel_coupling": True,
            "emission_drive_applied_once_inside_hazard": True,
            "directional_multiplier_applied_after_hazard": False,
            "fit_derived_shielding_cap_active": False,
            "response_classification_gate_active": False,
            "parallel_coupling_audit": audit_name,
        })
        rows.append(item)
    path = out / "mode_i_v10_0_5_results.json"
    path.write_text(json.dumps(rows, indent=2, default=str))
    return path


def main(argv: list[str] | None = None):
    args = list(sys.argv[1:] if argv is None else argv)
    theta = float(_option_value(args, "--crystal-theta-deg") or 45.0)
    probe_radius = float(
        os.environ.get("ARRHENIUS_TENSOR_DRIVE_PROBE_RADIUS_M", "1e-5")
    )
    sector = float(
        os.environ.get("ARRHENIUS_TENSOR_DRIVE_SECTOR_HALF_ANGLE_DEG", "25")
    )
    reset_tensor_drive_runtime(
        TensorResolvedDriveConfig(
            crystal_theta_deg=theta,
            probe_radius_m=probe_radius,
            sector_half_angle_deg=sector,
            damage_cutoff=float(
                os.environ.get("ARRHENIUS_TENSOR_DRIVE_DAMAGE_CUTOFF", "0.85")
            ),
            min_elements=int(
                os.environ.get("ARRHENIUS_TENSOR_DRIVE_MIN_ELEMENTS", "3")
            ),
            schmid_reference=0.5,
        )
    )

    original_J = _j_integral.compute_J_integral
    original_stepper = _v1002.KineticCohesiveStepper
    original_factory = _v1003.engine_factory_v1003
    _j_integral.compute_J_integral = make_tensor_resolved_J_wrapper(original_J)
    _v1002.KineticCohesiveStepper = TensorResolvedKineticCohesiveStepper
    _v1003.engine_factory_v1003 = engine_factory_v1005

    try:
        results = _v1003.main(args)
    finally:
        _j_integral.compute_J_integral = original_J
        _v1002.KineticCohesiveStepper = original_stepper
        _v1003.engine_factory_v1003 = original_factory

    out_value = _option_value(args, "--out")
    if out_value is None:
        raise RuntimeError("v10.0.5 requires --out")
    out = Path(out_value)
    out.mkdir(parents=True, exist_ok=True)
    audit = tensor_drive_runtime_payload()
    if not audit.get("implementation_certified", False):
        raise RuntimeError(
            "v10.0.5 tensor-resolved coupling audit failed: "
            + json.dumps(audit, default=str)
        )
    if int(audit.get("nonzero_emission_drive_capture_count", 0)) < 1:
        raise RuntimeError(
            "v10.0.5 captured no nonzero slip-system emission drive"
        )
    audit_path = out / "parallel_opening_emission_v10_0_5_audit.json"
    audit_path.write_text(json.dumps(audit, indent=2, default=str))
    results_path = _write_v1005_results(out, results, audit_path.name)

    print("V10.0.5 PARALLEL OPENING/EMISSION COUPLING COMPLETE")
    print(json.dumps({
        "tensor_drive_capture_count": audit["capture_count"],
        "nonzero_emission_drive_capture_count": (
            audit["nonzero_emission_drive_capture_count"]
        ),
        "max_emission_drive_factor": audit["max_emission_drive_factor"],
        "directional_multiplier_applied_after_hazard": False,
        "fit_derived_shielding_cap_active": False,
        "response_classification_gate_active": False,
        "audit": str(audit_path),
        "results": str(results_path),
    }, indent=2))
    return results


if __name__ == "__main__":
    main()
