#!/usr/bin/env python3
"""Verify the v9.13 deterministic transfer gate on the full FEM/CZM stack."""
from __future__ import annotations

import argparse
import ast
import json
from pathlib import Path

REQUIRED = [
    "arrhenius_fracture/fem.py",
    "arrhenius_fracture/j_integral.py",
    "arrhenius_fracture/crack_backend.py",
    "arrhenius_fracture/sharp_front.py",
    "arrhenius_fracture/mode_i_first_passage_v9_13.py",
    "arrhenius_fracture/field_snapshots_v913.py",
    "arrhenius_fracture/material_rcurve_audit_v913.py",
    "run_mpz_v9_13_mode_i_rcurve.py",
    "run_mpz_v9_13_deterministic_material_transfer.py",
    "run_mpz_v9_13_deterministic_material_transfer_700K.sh",
]


def main() -> None:
    ap = argparse.ArgumentParser(); ap.add_argument("repo", type=Path); ns = ap.parse_args()
    repo = ns.repo.resolve()
    missing = [p for p in REQUIRED if not (repo / p).exists()]
    syntax_errors = []
    for rel in REQUIRED:
        path = repo / rel
        if path.suffix == ".py" and path.exists():
            try:
                ast.parse(path.read_text(), filename=str(path))
            except SyntaxError as exc:
                syntax_errors.append(f"{rel}:{exc.lineno}:{exc.msg}")
    runner = (repo / "run_mpz_v9_13_deterministic_material_transfer.py").read_text() if not missing else ""
    renderer = (repo / "arrhenius_fracture/field_snapshots_v913.py").read_text() if not missing else ""
    audit = (repo / "arrhenius_fracture/material_rcurve_audit_v913.py").read_text() if not missing else ""
    contracts = {
        "deterministic_default": 'default="deterministic"' in runner,
        "expected_emission_default": "default=False" in runner,
        "correct_root_temperature_summary": 'run_root / "rcurve_temperature_summary.csv"' in runner,
        "right_censor_gate_hardened": "completion_gate_passed" in audit and "failed_solver_cases" in audit,
        "pair_coverage_not_vacuous": "pairwise_coverage_gate_passed" in audit,
        "full_field_renderer": "field_snapshots_" in renderer,
        "tip_zoom_renderer": "field_snapshots_tip_zoom_" in renderer,
        "emitted_total_annotation": "Nemit=" in renderer,
        "resolution_aware_density": "display_coarse_grain_width_m" in renderer,
    }
    result = {
        "repo": str(repo), "missing": missing, "syntax_errors": syntax_errors,
        "contracts": contracts,
        "passed": not missing and not syntax_errors and all(contracts.values()),
    }
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result["passed"] else 1)


if __name__ == "__main__":
    main()
