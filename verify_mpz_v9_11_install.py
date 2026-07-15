#!/usr/bin/env python3
"""Static and import-level verification for the v9.11 integration."""
from __future__ import annotations

import argparse
import ast
import json
from pathlib import Path

REQUIRED_CORE = [
    "arrhenius_fracture/moving_process_zone_v9102.py",
    "arrhenius_fracture/mpz_front_engine.py",
    "arrhenius_fracture/mixed_mode_first_passage_v8.py",
    "arrhenius_fracture/plasticity.py",
    "arrhenius_fracture/sharp_front.py",
]
REQUIRED_V911 = [
    "arrhenius_fracture/mpz_parameterization_v911.py",
    "arrhenius_fracture/bulk_plasticity_v9102.py",
    "arrhenius_fracture/process_zone_2d_v911.py",
    "arrhenius_fracture/moving_process_zone_v911.py",
    "arrhenius_fracture/mpz_front_engine_v911.py",
    "arrhenius_fracture/mixed_mode_first_passage_v9_11.py",
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("repo", type=Path)
    args = ap.parse_args()
    repo = args.repo.resolve()
    missing = [p for p in REQUIRED_CORE + REQUIRED_V911 if not (repo / p).exists()]
    syntax_errors = []
    for rel in REQUIRED_V911:
        path = repo / rel
        if not path.exists():
            continue
        try:
            ast.parse(path.read_text(), filename=str(path))
        except SyntaxError as exc:
            syntax_errors.append(f"{rel}:{exc.lineno}:{exc.msg}")
    result = {
        "repo": str(repo),
        "missing": missing,
        "syntax_errors": syntax_errors,
        "passed": not missing and not syntax_errors,
    }
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result["passed"] else 1)


if __name__ == "__main__":
    main()
