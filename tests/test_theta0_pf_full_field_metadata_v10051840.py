from __future__ import annotations

import json
from pathlib import Path

from arrhenius_fracture.theta0_pf_full_field_metadata_v10051840 import (
    AUDIT_FILE,
    BULK_MODEL,
    normalize_theta0_full_field_outputs,
)


def test_metadata_normalization_preserves_old_values_in_audit(tmp_path: Path):
    (tmp_path / "summary.json").write_text(
        json.dumps(
            [
                {
                    "crystal_theta_deg": 45.0,
                    "target_class": "DBTT",
                    "bulk_state_evolves_in_fem": False,
                    "taylor_renewal_time_s": 1.0,
                }
            ]
        )
    )
    manifest_path = (
        tmp_path / "persistent_site_production_manifest_v10_0_5_18_3_9.json"
    )
    manifest_path.write_text(
        json.dumps(
            {
                "crystal_theta_deg": 45.0,
                "bulk_state_evolves_in_fem": False,
                "physics_contract": {
                    "bulk_state_evolves_in_fem": False,
                },
                "two_d_state_policy": {
                    "taylor_renewal_time_s": 1.0,
                },
            }
        )
    )
    (tmp_path / "persistent_site_parameter_selection_test.json").write_text(
        json.dumps({"policy": {"crystal_theta_deg": 45.0}})
    )
    (tmp_path / "run_args.json").write_text(
        json.dumps(
            {
                "crystal_theta_deg": 45.0,
                "bulk_plasticity_mode": "tip_only",
            }
        )
    )

    audit = normalize_theta0_full_field_outputs(tmp_path)

    summary = json.loads((tmp_path / "summary.json").read_text())[0]
    assert summary["crystal_theta_deg"] == 0.0
    assert summary["material_class"] == "peak"
    assert summary["target_class"] == "peak"
    assert summary["bulk_state_evolves_in_fem"] is True
    assert summary["taylor_renewal_time_s"] == 1.0e-9
    assert summary["bulk_model"] == BULK_MODEL

    manifest = json.loads(manifest_path.read_text())
    assert manifest["crystal_theta_deg"] == 0.0
    assert manifest["bulk_state_evolves_in_fem"] is True
    assert manifest["physics_contract"]["bulk_state_evolves_in_fem"] is True
    assert manifest["two_d_state_policy"]["taylor_renewal_time_s"] == 1.0e-9

    assert audit["metadata_only"] is True
    assert audit["physics_or_mesh_modified"] is False
    assert audit["FEM_state_modified"] is False
    assert audit["cohesive_state_modified"] is False
    previous_values = {
        (row["location"], row["field"]): row["previous"]
        for row in audit["changes"]
    }
    assert previous_values[("summary.json[0]", "crystal_theta_deg")] == 45.0
    assert previous_values[("summary.json[0]", "bulk_state_evolves_in_fem")] is False
    assert (tmp_path / AUDIT_FILE).is_file()
