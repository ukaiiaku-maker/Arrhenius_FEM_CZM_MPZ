#!/usr/bin/env python3
"""Verify the v9.14 event-remeshed full FEM/CZM implementation."""
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
    "arrhenius_fracture/event_remesh_czm_v914.py",
    "arrhenius_fracture/event_equilibrium_v914.py",
    "arrhenius_fracture/event_remesh_audit_v914.py",
    "arrhenius_fracture/mode_i_first_passage_v9_14.py",
    "run_mpz_v9_14_mode_i_rcurve.py",
    "run_mpz_v9_14_event_remesh_gate.py",
    "run_mpz_v9_14_event_remesh_700K.sh",
]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("repo", type=Path)
    ns = ap.parse_args()
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

    remesh = (repo / "arrhenius_fracture/event_remesh_czm_v914.py").read_text() if not missing else ""
    eq = (repo / "arrhenius_fracture/event_equilibrium_v914.py").read_text() if not missing else ""
    mode = (repo / "arrhenius_fracture/mode_i_first_passage_v9_14.py").read_text() if not missing else ""
    audit = (repo / "arrhenius_fracture/event_remesh_audit_v914.py").read_text() if not missing else ""
    runner = (repo / "run_mpz_v9_14_event_remesh_gate.py").read_text() if not missing else ""
    contracts = {
        "full_fem_czm_stack_present": all((repo / p).exists() for p in (
            "arrhenius_fracture/fem.py", "arrhenius_fracture/j_integral.py",
            "arrhenius_fracture/crack_backend.py", "arrhenius_fracture/sharp_front.py",
        )),
        "absolute_hazard_action_localization": (
            "absolute_integrated_hazard_action" in mode
            and "_BaseEngine.predict_clock_increment_drives" in mode
            and "dB / remaining" not in mode
        ),
        "absolute_action_tolerance_default": (
            'default=0.01' in runner
            and "maximum accepted absolute integrated-hazard increment dB" in runner
        ),
        "one_physical_event_backend": "one_physical_cohesive_event" in remesh,
        "refinement_only_parent_map": "cumulative = cumulative[parent_map]" in remesh,
        "cohesive_edges_not_bisected": "if edge_key in cohesive_edges" in remesh,
        "transactional_rollback": "_transaction_rollback" in remesh,
        "fail_fast_after_rollback": (
            "fail_fast_on_event_error: bool = True" in remesh
            and "event remesh transaction rolled back" in remesh
        ),
        "structural_pre_event_mesh_match": (
            "_mesh_state_compatibility" in eq
            and "structurally_identical_rebuilt_mesh" in eq
        ),
        "same_time_equilibrium": '"physical_time_increment_s": 0.0' in eq,
        "zero_hazard_during_equilibrium": '"hazard_action_increment": 0.0' in eq,
        "same_load_boundary_reuse": "_boundary_values(pre_displacement, pre_boundary)" in eq,
        "post_event_J_recomputed": "J_after_event_equilibrium" in eq and "strict_equilibrate" in mode,
        "post_event_MPZ_profile_recomputed": (
            "mpz_profile_recomputed_after_event" in eq
            and "all_MPZ_profiles_recomputed" in audit
        ),
        "conservative_rho_and_ep_audit": "relative_rho_area_integral_error" in eq and "max_relative_ep_area_integral_error" in eq,
        "numerical_and_material_gates_separate": "numerical_event_remesh_gate_passed" in audit and "material_transfer_gate_passed_v914" in audit,
        "tip_only_initial_gate": '"bulk_plasticity_mode": "tip_only"' in runner,
        "short_50um_default": 'default=50.0' in runner,
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
