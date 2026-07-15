#!/usr/bin/env python3
"""Verify that v9.12 is installed on top of the full 2-D FEM/CZM solver."""
from __future__ import annotations

import argparse
import ast
import json
from pathlib import Path

REQUIRED_FULL_FEM = [
    "arrhenius_fracture/fem.py",
    "arrhenius_fracture/j_integral.py",
    "arrhenius_fracture/crack_backend.py",
    "arrhenius_fracture/sharp_front.py",
    "arrhenius_fracture/mixed_mode_first_passage_v9_11.py",
    "arrhenius_fracture/mode_i_first_passage_v9_11.py",
]
REQUIRED_V912 = [
    "arrhenius_fracture/material_rcurve_audit_v912.py",
    "run_mpz_v9_12_tip_only_material_rcurve.py",
    "run_mpz_v9_12_tip_only_material_rcurve_700K.sh",
    "tests/test_material_rcurve_audit_v912.py",
    "tests/test_mpz_v9_12_runner.py",
]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("repo", type=Path)
    ns = ap.parse_args()
    repo = ns.repo.resolve()
    required = REQUIRED_FULL_FEM + REQUIRED_V912
    missing = [p for p in required if not (repo / p).exists()]
    syntax_errors = []
    for rel in required:
        path = repo / rel
        if path.suffix != ".py" or not path.exists():
            continue
        try:
            ast.parse(path.read_text(), filename=str(path))
        except SyntaxError as exc:
            syntax_errors.append(f"{rel}:{exc.lineno}:{exc.msg}")
    sharp = (repo / "arrhenius_fracture/sharp_front.py").read_text() if not missing else ""
    contracts = {
        "full_2d_fem_entry_present": "existing elastic-plastic FEM supplies K via the J-integral" in sharp,
        "field_snapshot_renderer_present": "def _render_field_snapshots" in sharp,
        "field_snapshot_has_stress": "sigma1 FEM (MPa)" in sharp,
        "field_snapshot_has_density": "log10 rho" in sharp,
        "field_snapshot_overlays_crack_path": "front_paths" in sharp,
    }
    result = {
        "repo": str(repo),
        "missing": missing,
        "syntax_errors": syntax_errors,
        "contracts": contracts,
        "passed": not missing and not syntax_errors and all(contracts.values()),
    }
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result["passed"] else 1)


if __name__ == "__main__":
    main()
