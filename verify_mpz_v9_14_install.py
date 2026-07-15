#!/usr/bin/env python3
from __future__ import annotations

import ast
import json
from pathlib import Path
import sys

root = Path(sys.argv[1] if len(sys.argv) > 1 else ".").resolve()
required = [
    "arrhenius_fracture/mode_i_first_passage_v9_14.py",
    "arrhenius_fracture/remesh_audit_v914.py",
    "run_mpz_v9_14_mode_i_rcurve.py",
    "run_mpz_v9_14_event_driven_remesh.py",
    "run_mpz_v9_14_event_driven_remesh_700K.sh",
    "tests/test_event_driven_remesh_v914.py",
]
missing = [p for p in required if not (root / p).exists()]
syntax_errors = []
for rel in required:
    path = root / rel
    if path.suffix == ".py" and path.exists():
        try:
            ast.parse(path.read_text(), filename=str(path))
        except SyntaxError as exc:
            syntax_errors.append(f"{rel}:{exc.lineno}:{exc.msg}")
mode = (root / required[0]).read_text() if not missing else ""
audit = (root / required[1]).read_text() if not missing else ""
contracts = {
    "adaptive_czm_forced": '"--crack-backend", "adaptive_czm"' in mode,
    "directional_transfer_path_forced": '"--crystal-aniso"' in mode,
    # The inherited parser default is branch-off.  v9.14 must neither enable
    # --crystal-branch nor inject the nonexistent --no-crystal-branch option.
    "branching_disabled_by_default": (
        '_inject_once(user_args, "--crystal-branch")' not in mode
        and '_inject_once(user_args, "--no-crystal-branch")' not in mode
    ),
    "exact_forward_plane": 'np.array([1.0, 0.0]' in mode,
    "adaptive_events_forced": '"--adaptive-events"' in mode,
    "one_event_audit": "one_topology_event_per_accepted_solve" in audit,
    "same_load_equilibrium_audit": "same_load_post_event_reequilibration_observed" in audit,
}
result = {
    "repo": str(root), "missing": missing, "syntax_errors": syntax_errors,
    "contracts": contracts,
    "passed": not missing and not syntax_errors and all(contracts.values()),
}
print(json.dumps(result, indent=2))
raise SystemExit(0 if result["passed"] else 1)
