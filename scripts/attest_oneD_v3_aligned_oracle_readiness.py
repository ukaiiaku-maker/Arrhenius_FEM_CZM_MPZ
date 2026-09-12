#!/usr/bin/env python3
"""Run one bounded 2-D source-certification case before oracle generation."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from scripts.audit_oneD_v3_material_transfer import (  # noqa: E402
    UNIFIED_V5_SHA,
    _verify_v5,
    _write_json,
    _write_csv,
)


def _probe(v5_repo: Path) -> dict[str, object]:
    code = r'''
from dataclasses import replace
import json
import numpy as np

from arrhenius_fracture.finalization_v3_schema import SCIENTIFIC_ACCEPTANCE_TOLERANCES
from arrhenius_fracture.unified_fracture_material_v5 import identity_record, material_bundle
from arrhenius_fracture.voiding_production_v5 import (
    deterministic_trajectory,
    equilibrate_fixed_load_with_production_fem,
    ligament_transaction,
    refine_downstream_source,
)

bundle = material_bundle("DBTT")
preconnection, _ = deterministic_trajectory(bundle=bundle, stop_before_ligament=True)
loaded = equilibrate_fixed_load_with_production_fem(
    replace(preconnection, displacement=preconnection.displacement * 2.0)
)
connected, ligament = ligament_transaction(loaded)
qualified, audit = refine_downstream_source(
    connected,
    max_refinement_levels=3,
    refinement_region="complete_cavity_ring",
    quality_improvement="constrained_v1",
)
attempts = []
for item in audit.get("attempts", []):
    previous = item["proof"]["previous_metrics"]
    current = item["proof"]["current_metrics"]
    before = np.asarray(previous["tensor_Pa"], dtype=float)
    after = np.asarray(current["tensor_Pa"], dtype=float)
    attempts.append({
        "refinement_level": int(item["level"]),
        "qualified": bool(item["qualified"]),
        "fixed_arc_tensor_relative_change": float(
            np.linalg.norm(after - before) / max(np.linalg.norm(after), 1.0e-300)
        ),
        "normalized_cavity_traction": float(current["normalized_traction"]),
        "eta_n_max": float(current["eta_n_max"]),
        "eta_t_max": float(current["eta_t_max"]),
        "minimum_mesh_quality": float(current["minimum_quality"]),
        "tensor_Pa": current["tensor_Pa"],
        "recovery_operator": current["recovery_operator"],
        "recovery_record": current["recovery_record"],
        "state_binding": item["proof"]["current_binding"],
    })
payload = {
    "material_class": "DBTT",
    "fracture_material_row_id": bundle.fracture_material_row_id,
    "material_identity": dict(identity_record(bundle, connected.material)),
    "geometry": {
        "cavity_center_m": list(connected.void_state.cavities[0].center_m),
        "cavity_radius_m": float(connected.void_state.cavities[0].radius_m),
        "connection_entry_m": list(connected.void_state.cavities[0].connection_entry_m),
        "connection_exit_m": list(connected.void_state.cavities[0].connection_exit_m),
        "opening_scale_from_4e-7_m": 2.0,
    },
    "ligament_energy_gate": {
        "accepted": bool(ligament.accepted),
        "energy_release_J_per_m": float(ligament.energy_release_J_per_m),
        "hazard_dissipation_J_per_m": float(ligament.hazard_dissipation_J_per_m),
        "energy_margin_J_per_m": float(ligament.energy_margin_J_per_m),
    },
    "source_qualification_status": audit["status"],
    "returned_original_state_on_noncertification": qualified is connected,
    "attempts": attempts,
    "frozen_limits": {
        "fixed_arc_tensor_relative_change_max": SCIENTIFIC_ACCEPTANCE_TOLERANCES["tensor_probe_relative"],
        "normalized_cavity_traction_max": SCIENTIFIC_ACCEPTANCE_TOLERANCES["cavity_traction_normalized"],
        "eta_n_max": 0.03,
        "eta_t_max": 0.025,
        "minimum_mesh_quality": SCIENTIFIC_ACCEPTANCE_TOLERANCES["mesh_minimum_quality"],
    },
}
print(json.dumps(payload, sort_keys=True))
'''
    completed = subprocess.run(
        (sys.executable, "-c", code), cwd=v5_repo,
        env={"PYTHONPATH": str(v5_repo), "MPLCONFIGDIR": "/private/tmp/matplotlib-v5-oracle"},
        text=True, capture_output=True,
    )
    if completed.returncode:
        raise RuntimeError("bounded V5 readiness worker failed:\n" + completed.stderr)
    return json.loads(completed.stdout)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--v5-repo", type=Path, required=True)
    parser.add_argument(
        "--out", type=Path,
        default=ROOT / "analysis_outputs/oneD_v3_aligned_temperature_transfer",
    )
    args = parser.parse_args()
    v5_repo = args.v5_repo.resolve()
    head = _verify_v5(v5_repo)
    result = _probe(v5_repo)
    if result["source_qualification_status"] not in (
        "SOURCE_TENSOR_QUALIFIED", "SOURCE_TENSOR_UNQUALIFIED"
    ):
        raise RuntimeError("readiness probe returned an unregistered scientific status")
    final_attempt = result["attempts"][-1]
    payload = {
        "schema": "oneD.v3.aligned-oracle-readiness/2",
        "source_repository": "ukaiiaku-maker/PF-fracture-fatigue",
        "source_head": head,
        "source_head_matches_pinned_unified_commit": head == UNIFIED_V5_SHA,
        "bounded_readiness_cases_run": 1,
        "requested_refinement_levels": 3,
        "accepted_oracle_states": 0,
        "oracle_matrix_status": "BLOCKED_DOWNSTREAM_SOURCE_TENSOR_UNQUALIFIED",
        "source_recovery_operator": "CAVITY_FIXED_ARC_PATCH_RECOVERY_V2",
        "historical_source_recovery": {
            "CAVITY_SOURCE_RECOVERY_V1_INCIDENT_CST_MAX_PRINCIPAL": "FAIL_NONCONVERGENT"
        },
        "expected_scientific_noncertification_converted_to_unavailable": True,
        "unexpected_programming_exceptions_caught": False,
        "missing_fields_inferred_or_synthesized": False,
        "downstream_child_created": False,
        "failure_classification": [
            "TANGENTIAL_STRESS_CONVERGENCE",
            "NORMAL_DIRECTION_RESOLUTION",
        ],
        "readiness_gate_taxonomy": {
            "decision": "B_STILL_MANDATORY_V2_GATES",
            "eta_n": "MANDATORY_V2_GATE",
            "eta_t": "MANDATORY_V2_GATE",
        },
        "final_governing_predicates": {
            "fixed_arc_tensor_convergence": (
                final_attempt["fixed_arc_tensor_relative_change"]
                <= result["frozen_limits"]["fixed_arc_tensor_relative_change_max"]
            ),
            "cavity_traction": (
                final_attempt["normalized_cavity_traction"]
                <= result["frozen_limits"]["normalized_cavity_traction_max"]
            ),
            "normal_direction_resolution": (
                final_attempt["eta_n_max"] <= result["frozen_limits"]["eta_n_max"]
            ),
            "tangential_direction_resolution": (
                final_attempt["eta_t_max"] <= result["frozen_limits"]["eta_t_max"]
            ),
        },
        "probe": result,
    }
    _write_json(args.out / "aligned_oracle_readiness.json", payload)

    decision_path = args.out / "decision.json"
    decision = json.loads(decision_path.read_text())
    decision["execution_counts"].update({
        "new_2d_mechanics_solves": 1,
        "bounded_2d_oracle_readiness_cases": 1,
        "new_2d_oracle_states_accepted": 0,
        "source_refinement_levels_in_readiness_case": 3,
    })
    decision["prospective_execution_contract"]["oracle_matrix"]["status"] = payload["oracle_matrix_status"]
    decision["pilot_terminal"].update({
        "ONE_D_V3_MECHANICS_MAP_FIT": "BLOCKED_2D_DOWNSTREAM_SOURCE_UNQUALIFIED",
        "ONE_D_V3_FIXED_STATE_2D_TRANSFER": "BLOCKED_INCOMPLETE_2D_ORACLE",
        "ONE_D_V3_LOW_T_TRAJECTORY_TRANSFER": "NOT_RUN_ORACLE_GATE_BLOCKED",
        "ONE_D_V3_TRANSITION_T_TRAJECTORY_TRANSFER": "NOT_RUN_ORACLE_GATE_BLOCKED",
        "ONE_D_V3_HIGH_T_TRAJECTORY_TRANSFER": "NOT_RUN_ORACLE_GATE_BLOCKED",
        "ONE_D_V3_COMPLETE_ALIGNED_MONOTONIC_TRANSFER": "BLOCKED_2D_ORACLE_INCOMPLETE",
    })
    decision["temperature_transfer"] = {
        key: "NOT_RUN_ORACLE_GATE_BLOCKED"
        for key in decision["temperature_transfer"]
    }
    decision["blocking_reason"] = {
        "code": "DOWNSTREAM_SOURCE_TENSOR_UNQUALIFIED",
        "artifact": "aligned_oracle_readiness.json",
        "governing_predicates": payload["final_governing_predicates"],
        "failure_classification": payload["failure_classification"],
        "required_action": "retain the blocked oracle gate; any further source-recovery design requires a separate 2-D commit",
    }
    decision["next_bounded_step"] = decision["blocking_reason"]["required_action"]
    _write_json(decision_path, decision)
    ledger_path = args.out / "paired_case_ledger.csv"
    with ledger_path.open(newline="") as stream:
        ledger = list(csv.DictReader(stream))
    for row in ledger:
        for key in ("no_void_1d_status", "no_void_2d_status", "void_1d_status", "void_2d_status"):
            row[key] = "NOT_RUN_ORACLE_GATE_BLOCKED"
        row["delta_observable_status"] = "UNAVAILABLE_INCOMPLETE_2D_ORACLE"
    _write_csv(ledger_path, ledger)
    print(payload["oracle_matrix_status"])
    print("ORACLE_STATES_ACCEPTED=0/18")
    print("PAIRED_CASES_RUN=0/12")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
