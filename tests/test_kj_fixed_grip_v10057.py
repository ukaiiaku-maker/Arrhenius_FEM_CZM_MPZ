import csv
import json
import math
from pathlib import Path
import subprocess
import sys

import pytest

from arrhenius_fracture.kj_audit_v10056 import SpecimenGeometryV10056
from arrhenius_fracture.kj_audit_v10057 import (
    DEFAULT_FIXED_GRIP_REFERENCE,
    build_kj_audit_row,
    fixed_grip_reference_K_Pa_sqrt_m,
)


def test_fixed_grip_reference_is_primary_and_uniform_tension_is_retained():
    geometry = SpecimenGeometryV10056()
    sigma = 100.0e6
    K = fixed_grip_reference_K_Pa_sqrt_m(sigma, geometry)
    assert K == pytest.approx(
        sigma * math.sqrt(math.pi * geometry.initial_crack_m) * 1.2003
    )
    row = build_kj_audit_row(
        Ftop_N_per_thickness=sigma * geometry.width_m,
        KJ_Pa_sqrt_m=K,
        outer_radius_m=100.0e-6,
        geometry=geometry,
    )
    assert row["K_reference_boundary_condition"] == "symmetric_fixed_grip_displacement"
    assert row["KJ_over_K_LEFM_gross"] == pytest.approx(1.0)
    assert row["KJ_over_K_fixed_grip_reference"] == pytest.approx(1.0)
    assert row["uniform_tension_edge_geometry_factor_Y"] > row["fixed_grip_geometry_factor_Y"]


def test_fixed_grip_reference_fails_closed_for_other_geometry():
    geometry = SpecimenGeometryV10056(width_m=3.0e-3)
    with pytest.raises(ValueError, match="geometry specific"):
        fixed_grip_reference_K_Pa_sqrt_m(100.0e6, geometry)


def test_generator_recovers_converged_factor(tmp_path):
    rows = tmp_path / "fixed_grip_rows.csv"
    out = tmp_path / "reference.json"
    a = 0.5e-3
    factors = [1.18, 1.195, 1.2000, 1.2004, 1.2003]
    with rows.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["mesh_label", "sigma_gross_Pa", "KJ_Pa_sqrt_m"],
        )
        writer.writeheader()
        for index, factor in enumerate(factors):
            sigma = 100.0e6
            writer.writerow(
                {
                    "mesh_label": f"m{index}",
                    "sigma_gross_Pa": sigma,
                    "KJ_Pa_sqrt_m": sigma * math.sqrt(math.pi * a) * factor,
                }
            )
    script = Path(__file__).resolve().parents[1] / "scripts" / "generate_fixed_grip_reference_v10057.py"
    completed = subprocess.run(
        [sys.executable, str(script), "--rows", str(rows), "--out", str(out)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr + completed.stdout
    payload = json.loads(out.read_text())
    assert payload["convergence_passed"]
    assert payload["geometry_factor_Y"] == pytest.approx(1.2003, rel=5.0e-4)
    assert payload["tail_max_relative_spread"] < 0.02
    assert DEFAULT_FIXED_GRIP_REFERENCE.geometry_factor_Y == pytest.approx(1.2003)
