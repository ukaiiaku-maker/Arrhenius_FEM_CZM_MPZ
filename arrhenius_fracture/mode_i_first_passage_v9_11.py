"""Direct Mode-I full 2-D FEM/CZM validation for the selected v9.11 MPZ law.

The ceramic, weakT, and DBTT parameterizations were calibrated in the moving
reference-frame constitutive campaign.  A pure Mode-I validation must therefore
not require a second mixed-mode boundary-response calibration.  This entry point
uses remote opening only and preserves the calibrated scalar K drive by forcing
the additional mixed-mode directional multipliers to unity.

Crystal anisotropy remains active in the FEM equilibrium, domain integral,
process-zone stress profile, and bulk plasticity.  Finite-radius FEM stresses are
retained for spatial-profile sampling and diagnostics, but they do not rescale
the calibrated Mode-I cleavage or emission drive.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path
import sys
from typing import Any

from . import mixed_mode_first_passage_v8 as v8
from . import mixed_mode_first_passage_v9_11 as v911

MODEL_ID = "FEM_CZM_Mode_I_MPZ_v9_11_direct_material_calibration"

_RESERVED_OPTIONS = {
    "--mixity-loading-angle-deg",
    "--mixity-open-coeff",
    "--mixity-shear-coeff",
    "--target-traction-phase-deg",
    "--traction-shear-sign",
    "--reference-cleavage-shape",
    "--reference-slip-shape",
    "--shear-emission-weight",
    "--directional-factor-max",
}


def unit_mode_i_directional_factors(*args: Any, **kwargs: Any) -> dict[str, Any]:
    """Return no additional directional rescaling for the calibrated Mode-I law."""
    return {
        "cleavage_factor": 1.0,
        "emission_factor": 1.0,
        "cleavage_factor_raw": 1.0,
        "emission_factor_raw": 1.0,
        "shear_excess_shape": 0.0,
        "directional_factor_cap_active": False,
        "mode_I_direct_material_calibration": True,
    }


def validate_direct_mode_args(argv: list[str]) -> None:
    """Prevent accidental reintroduction of mixed-mode calibration controls."""
    for token in argv:
        name = token.split("=", 1)[0]
        if name in _RESERVED_OPTIONS:
            raise SystemExit(
                f"{name} is fixed by the direct Mode-I v9.11 entry point; "
                "use mixed_mode_first_passage_v9_11 for a mixed-mode campaign"
            )


def _option_value(argv: list[str], name: str) -> str | None:
    for i, token in enumerate(argv):
        if token == name and i + 1 < len(argv):
            return argv[i + 1]
        if token.startswith(name + "="):
            return token.split("=", 1)[1]
    return None


def _annotate_outputs(out: Path, results: list[dict[str, Any]]) -> None:
    flags = {
        "model": MODEL_ID,
        "mode_I_direct_material_calibration": True,
        "remote_loading_open_coeff": 1.0,
        "remote_loading_shear_coeff": 0.0,
        "directional_drive_factors_forced_unity": True,
        "mixed_mode_boundary_calibration_required": False,
        "crystal_anisotropy_role": (
            "FEM equilibrium, J integral, process-zone spatial profile, and bulk plasticity"
        ),
    }
    for payload in results:
        payload.update(flags)

    summary = out / "anisotropic_calibrated_tip_first_passage_summary.json"
    if summary.exists() and results:
        summary.write_text(json.dumps(results[-1], indent=2, default=str))

    summary_csv = out / "anisotropic_calibrated_tip_first_passage_summary.csv"
    if summary_csv.exists() and results:
        with summary_csv.open("w", newline="") as fp:
            writer = csv.DictWriter(fp, fieldnames=list(results[-1]))
            writer.writeheader()
            writer.writerow(results[-1])

    audit = out / "mpz_v9_11_integration_audit.json"
    if audit.exists():
        data = json.loads(audit.read_text())
        data["mode_I_validation"] = flags
        audit.write_text(json.dumps(data, indent=2, default=str))


def main(argv: list[str] | None = None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    validate_direct_mode_args(user_args)

    fixed = [
        "--mixity-open-coeff", "1",
        "--mixity-shear-coeff", "0",
        "--target-traction-phase-deg", "0",
        "--traction-shear-sign", "1",
        "--reference-cleavage-shape", "1",
        "--reference-slip-shape", "0",
        "--shear-emission-weight", "0",
        "--directional-factor-max", "1",
    ]

    original = v8.directional_drive_factors
    v8.directional_drive_factors = unit_mode_i_directional_factors
    try:
        results = v911.main(fixed + user_args)
    finally:
        v8.directional_drive_factors = original

    out_value = _option_value(user_args, "--out")
    if out_value is not None:
        _annotate_outputs(Path(out_value), results)
    print("MODE_I_MPZ_V9_11_DIRECT complete")
    return results


if __name__ == "__main__":
    main()
