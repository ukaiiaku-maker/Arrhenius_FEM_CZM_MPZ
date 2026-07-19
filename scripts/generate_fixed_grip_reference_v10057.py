#!/usr/bin/env python3
"""Generate a geometry-specific fixed-grip K-reference artifact.

Input CSV rows are elastic fixed-crack FEM calibrations ordered from coarse to
fine. Required columns are ``sigma_gross_Pa`` and ``KJ_Pa_sqrt_m``. Optional
mesh/contour columns are preserved in the convergence record.  The geometry
factor is Y=K/(sigma*sqrt(pi*a)).
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

from arrhenius_fracture.kj_audit_v10057 import REFERENCE_SCHEMA


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rows", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--width-m", type=float, default=2.0e-3)
    parser.add_argument("--height-m", type=float, default=4.0e-3)
    parser.add_argument("--initial-crack-m", type=float, default=0.5e-3)
    parser.add_argument("--tail-count", type=int, default=3)
    parser.add_argument("--relative-spread-tolerance", type=float, default=0.02)
    args = parser.parse_args()
    if args.out.exists():
        raise SystemExit(f"refusing to overwrite {args.out}")
    with args.rows.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if len(rows) < max(int(args.tail_count), 3):
        raise SystemExit("fixed-grip convergence requires at least three rows")
    derived = []
    for index, row in enumerate(rows):
        try:
            sigma = float(row["sigma_gross_Pa"])
            K = float(row["KJ_Pa_sqrt_m"])
        except Exception as exc:
            raise SystemExit(f"invalid fixed-grip row {index}: {row}") from exc
        if not (math.isfinite(sigma) and sigma > 0.0 and math.isfinite(K) and K > 0.0):
            raise SystemExit(f"non-positive or non-finite fixed-grip row {index}")
        Y = K / (sigma * math.sqrt(math.pi * float(args.initial_crack_m)))
        derived.append({**row, "geometry_factor_Y": Y})
    tail = derived[-max(int(args.tail_count), 3):]
    values = np.asarray([float(row["geometry_factor_Y"]) for row in tail])
    reference = float(np.median(values))
    spread = float(np.max(np.abs(values / reference - 1.0)))
    passed = spread <= float(args.relative_spread_tolerance)
    if not passed:
        raise SystemExit(
            "fixed-grip geometry factor did not converge: "
            f"tail spread={spread:.6g} > {args.relative_spread_tolerance:.6g}"
        )
    payload = {
        "schema": REFERENCE_SCHEMA,
        "width_m": float(args.width_m),
        "height_m": float(args.height_m),
        "initial_crack_m": float(args.initial_crack_m),
        "geometry_factor_Y": reference,
        "boundary_condition": "symmetric_fixed_grip_displacement",
        "convergence_passed": True,
        "tail_count": len(tail),
        "tail_max_relative_spread": spread,
        "relative_spread_tolerance": float(args.relative_spread_tolerance),
        "input_rows": str(args.rows.resolve()),
        "convergence_rows": derived,
        "provenance": "generated from elastic fixed-crack FEM convergence rows",
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2))
    print(json.dumps({"out": str(args.out), "geometry_factor_Y": reference, "spread": spread}, indent=2))


if __name__ == "__main__":
    main()
