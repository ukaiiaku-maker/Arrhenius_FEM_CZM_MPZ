#!/usr/bin/env python3
"""Generate candidate-independent PF/FEM source-drive factor maps.

Only deterministic elastic solves and read-only tensor probes are used.  The
stored factors are the normalized values consumed by the production source
engines (factor times opening stress), not ``tau/opening`` Schmid ratios.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_outputs" / "oneD_v2_predictive_model"
MAPS = ROOT / "analysis_outputs" / "oneD_v2_mechanics_maps_and_baselines"
sys.path.insert(0, str(ROOT))

from reduced_fracture_v2.production_oracles import (
    FEMCZMExactElasticFieldOracle,
    FEMCZMExactSourceProbeOperator,
    PFExactElasticFieldOracle,
    PFExactSourceProbeOperator,
)

RADII_UM = (1.0, 2.0, 4.0, 8.0, 12.0, 25.0, 50.0, 100.0)
EXTENSIONS_UM = (0.0, 25.0, 100.0, 300.0, 500.0, 700.0, 1000.0)


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _straight_events(extension_m: float) -> list[dict]:
    events: list[dict] = []
    start = np.asarray((5.0e-4, 0.0), dtype=float)
    remaining = max(float(extension_m), 0.0)
    while remaining > 1.0e-15:
        advance = min(5.0e-6, remaining)
        end = start + np.asarray((advance, 0.0))
        events.append({"p0_m": start.tolist(), "p1_m": end.tolist(), "front_id": 0})
        start = end
        remaining -= advance
    return events


def generate(backend: str) -> Path:
    OUT.mkdir(parents=True, exist_ok=True)
    if backend == "pf":
        oracle = PFExactElasticFieldOracle()
        probe = PFExactSourceProbeOperator()
        mechanics = MAPS / "oneD_v2_pf_mechanics_map.csv"
        label = "PF"
    else:
        oracle = FEMCZMExactElasticFieldOracle()
        probe = FEMCZMExactSourceProbeOperator()
        mechanics = MAPS / "oneD_v2_fem_native_mechanics_map.csv"
        label = "FEMCZM"
    rows: list[dict] = []
    for extension_um in EXTENSIONS_UM:
        extension_m = extension_um * 1.0e-6
        geometry = (
            {"events": _straight_events(extension_m), "tip_xy_m": (5.0e-4 + extension_m, 0.0)}
            if backend == "pf"
            else {"events": [], "extension_m": extension_m}
        )
        snapshot = oracle.solve(geometry)
        for radius_um in RADII_UM:
            radius_m = radius_um * 1.0e-6
            state = {
                "tip_radius_m": radius_m,
                "front_width_m": 10.0e-6,
                "source_area_m2": radius_m * 10.0e-6,
                "source_geometry_fingerprint": f"qualification:r={radius_um:g}um",
                "crystal_theta_deg": 0.0,
            }
            raw = probe.evaluate_raw(snapshot, state, snapshot.reference_opening_m)
            opening_stress = float(np.asarray(raw["opening_tensor_Pa"], dtype=float)[1, 1])
            tau = np.asarray(raw["tau_signed_Pa"], dtype=float)
            factors = np.asarray(raw["drive_factors"], dtype=float)
            rows.append({
                "backend": label,
                "extension_um": extension_um,
                "tip_radius_um": radius_um,
                "reference_opening_m": snapshot.reference_opening_m,
                "drive_factor_system_0": factors[0],
                "drive_factor_system_1": factors[1],
                "tau_signed_system_0_Pa": tau[0],
                "tau_signed_system_1_Pa": tau[1],
                "opening_stress_Pa": opening_stress,
                "tau_over_opening_system_0": tau[0] / max(abs(opening_stress), 1.0e-300),
                "tau_over_opening_system_1": tau[1] / max(abs(opening_stress), 1.0e-300),
                "reliable": bool(raw["reliable"]),
                "field_snapshot_hash": snapshot.snapshot_hash,
                "topology_fingerprint": snapshot.topology_wake_fingerprint,
                "source_commit": snapshot.source_commit,
                "source_repository": snapshot.source_repository,
            })
    path = OUT / f"oneD_v2_{backend}_source_drive_map.csv"
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=list(rows[0]),
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    manifest = {
        "schema": "oneD_v2_source_drive_map_v1",
        "backend": label,
        "status": "QUALIFIED_DETERMINISTIC_SOURCE_PROBE_MAP",
        "factor_semantics": "production_normalized_drive_factor_times_opening_stress",
        "tau_over_opening_is_diagnostic_only": True,
        "extensions_um": list(EXTENSIONS_UM),
        "tip_radii_um": list(RADII_UM),
        "elastic_solve_count": oracle.solve_count,
        "probe_query_count": probe.query_count,
        "mechanics_map": str(mechanics),
        "mechanics_map_sha256": _sha(mechanics),
        "csv": str(path),
        "csv_sha256": _sha(path),
        "no_stochastic_evolution": True,
    }
    (OUT / f"oneD_v2_{backend}_source_drive_map_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    print(path)
    return path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("backend", choices=("pf", "fem"))
    args = parser.parse_args()
    generate(args.backend)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
