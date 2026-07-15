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
import math
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


def _derived_mpz_length_args(argv: list[str]) -> tuple[list[str], float]:
    """Keep every legacy/reporting length synchronized with the v9.11 MPZ.

    ``sharp_front.run_2d`` prints its resolution audit before the v9.11 engine is
    built.  Therefore its legacy ``L_pz`` and ``mpz_length_m`` namespace values
    must be derived from the authoritative ``--mpz-length-um`` option before the
    solver starts.  A conflicting explicit value is rejected rather than
    silently producing a misleading console audit or a split length scale.
    """
    mpz_um = float(_option_value(argv, "--mpz-length-um") or 100.0)
    if not math.isfinite(mpz_um) or mpz_um <= 0.0:
        raise SystemExit("--mpz-length-um must be finite and positive")
    mpz_m = mpz_um * 1.0e-6
    derived: list[str] = []
    for name in ("--L-pz", "--mpz-length-m"):
        raw = _option_value(argv, name)
        if raw is None:
            derived.extend([name, f"{mpz_m:.16g}"])
            continue
        value = float(raw)
        if not math.isclose(value, mpz_m, rel_tol=1.0e-12, abs_tol=1.0e-15):
            raise SystemExit(
                f"{name}={value:.16g} conflicts with --mpz-length-um={mpz_um:.16g}; "
                "v9.11 uses one authoritative moving-process-zone length"
            )
    return derived, mpz_m


def _float_or_none(value: Any) -> float | None:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if math.isfinite(out) else None


def _clock_summary(out: Path, T_K: float) -> dict[str, Any]:
    """Read the accepted solver rows and export actual clock residuals.

    The generic v8 summary previously requested ``B_final`` from a compact base
    summary that did not contain the field, producing JSON ``null`` even though
    the accepted step CSV and console contained the finite residual clock.  The
    step CSV is the authoritative accepted-state record.
    """
    path = out / f"steps_{int(round(float(T_K))):04d}K.csv"
    if not path.exists():
        return {
            "B_final": None,
            "B_first_fire_residual": None,
            "B_summary_source": "steps_csv_missing",
        }
    try:
        rows = list(csv.DictReader(path.open()))
    except Exception:
        rows = []
    if not rows:
        return {
            "B_final": None,
            "B_first_fire_residual": None,
            "B_summary_source": "steps_csv_empty",
        }

    final = rows[-1]
    first_fire = next(
        (row for row in rows if (_float_or_none(row.get("n_fire")) or 0.0) > 0.0),
        None,
    )
    return {
        "B_final": _float_or_none(final.get("B")),
        "B_first_fire_residual": (
            _float_or_none(first_fire.get("B")) if first_fire is not None else None
        ),
        "B_final_step": _float_or_none(final.get("step")),
        "B_first_fire_step": (
            _float_or_none(first_fire.get("step")) if first_fire is not None else None
        ),
        "B_summary_source": path.name,
    }


def _annotate_outputs(
    out: Path,
    results: list[dict[str, Any]],
    mpz_length_m: float,
) -> None:
    flags = {
        "model": MODEL_ID,
        "mode_I_direct_material_calibration": True,
        "remote_loading_open_coeff": 1.0,
        "remote_loading_shear_coeff": 0.0,
        "directional_drive_factors_forced_unity": True,
        "mixed_mode_boundary_calibration_required": False,
        "active_mpz_length_m": float(mpz_length_m),
        "active_mpz_length_um": float(mpz_length_m * 1.0e6),
        "sharp_front_L_pz_synchronized_to_active_mpz": True,
        "crystal_anisotropy_role": (
            "FEM equilibrium, J integral, process-zone spatial profile, and bulk plasticity"
        ),
    }
    for payload in results:
        payload.update(flags)
        payload.update(_clock_summary(out, float(payload.get("T_K", 0.0))))

    summary = out / "anisotropic_calibrated_tip_first_passage_summary.json"
    if results:
        summary.write_text(json.dumps(results[-1], indent=2, default=str))

    summary_csv = out / "anisotropic_calibrated_tip_first_passage_summary.csv"
    if results:
        with summary_csv.open("w", newline="") as fp:
            writer = csv.DictWriter(fp, fieldnames=list(results[-1]))
            writer.writeheader()
            writer.writerow(results[-1])

    # Preserve every temperature when one invocation runs a multi-temperature
    # R-curve campaign; the historical single-summary filenames remain as the
    # last-temperature compatibility view.
    if results:
        (out / "mode_i_v9_11_temperature_summaries.json").write_text(
            json.dumps(results, indent=2, default=str)
        )
        columns = sorted({key for row in results for key in row})
        with (out / "mode_i_v9_11_temperature_summaries.csv").open("w", newline="") as fp:
            writer = csv.DictWriter(fp, fieldnames=columns)
            writer.writeheader()
            writer.writerows(results)

    audit = out / "mpz_v9_11_integration_audit.json"
    if audit.exists():
        data = json.loads(audit.read_text())
        data["mode_I_validation"] = flags
        audit.write_text(json.dumps(data, indent=2, default=str))


def main(argv: list[str] | None = None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    validate_direct_mode_args(user_args)
    derived_length_args, mpz_length_m = _derived_mpz_length_args(user_args)

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
        results = v911.main(fixed + derived_length_args + user_args)
    finally:
        v8.directional_drive_factors = original

    out_value = _option_value(user_args, "--out")
    if out_value is not None:
        _annotate_outputs(Path(out_value), results, mpz_length_m)
    print("MODE_I_MPZ_V9_11_DIRECT complete")
    return results


if __name__ == "__main__":
    main()
