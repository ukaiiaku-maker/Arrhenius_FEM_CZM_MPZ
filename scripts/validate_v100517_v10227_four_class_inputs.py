#!/usr/bin/env python3
"""Fail-closed preflight for the exact PF v10.2.27 four-class registry."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from arrhenius_fracture.audited_pf_parameter_bridge_v100517_v10227_canonical import (
    EXPECTED_ACTIVE_FINGERPRINT_SHA256,
    EXPECTED_OPTIONS,
    load_audited_parameter_option,
)


def validate(pf_root: Path, options: list[str], out: Path) -> dict:
    if options != list(EXPECTED_OPTIONS):
        raise ValueError(
            "four-class option order mismatch; expected exactly "
            + json.dumps(list(EXPECTED_OPTIONS))
        )
    records = []
    candidate_ids = set()
    row_hashes = set()
    classes = set()
    registry_hashes = set()
    for option in options:
        candidate, audit = load_audited_parameter_option(pf_root, option)
        if candidate.candidate_id in candidate_ids:
            raise ValueError(f"duplicate candidate id {candidate.candidate_id}")
        if audit["selected_row_sha256"] in row_hashes:
            raise ValueError(f"duplicate selected registry row for {option}")
        candidate_ids.add(candidate.candidate_id)
        row_hashes.add(audit["selected_row_sha256"])
        classes.add(audit["material_class"])
        registry_hashes.add(audit["registry_sha256"])
        records.append(audit)
    if classes != {"peak", "DBTT", "weakT", "ceramic"}:
        raise ValueError(f"unexpected four-class labels: {sorted(classes)}")
    if len(registry_hashes) != 1:
        raise ValueError(f"inconsistent registry hashes across options: {sorted(registry_hashes)}")
    observed_registry_sha = next(iter(registry_hashes))
    payload = {
        "schema": "v10.0.5.17_audited_PF_v10_2_27_four_class_preflight",
        "PF_repo_root": str(pf_root.resolve()),
        "n_options": len(records),
        "options_in_production_order": list(options),
        "classes": [record["material_class"] for record in records],
        "all_candidate_ids_unique": True,
        "all_selected_rows_unique": True,
        "registry_sha256_observed": observed_registry_sha,
        "registry_byte_hash_used_as_physics_identity": False,
        "active_parameter_fingerprint_sha256": EXPECTED_ACTIVE_FINGERPRINT_SHA256,
        "mechanics_changed": False,
        "source_closure_changed": False,
        "stochastic_cleavage_law_changed": False,
        "persistent_sites": True,
        "finite_source_inventory": False,
        "source_refresh": False,
        "explicit_recovery": False,
        "records": records,
    }
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2, sort_keys=True, default=str) + "\n")
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pf-repo-root", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("options", nargs="+")
    args = parser.parse_args()
    payload = validate(
        args.pf_repo_root.expanduser().resolve(),
        args.options,
        args.out.expanduser().resolve(),
    )
    print(
        json.dumps(
            {
                "n_options": payload["n_options"],
                "classes": payload["classes"],
                "registry_sha256_observed": payload["registry_sha256_observed"],
                "active_parameter_fingerprint_sha256": payload[
                    "active_parameter_fingerprint_sha256"
                ],
                "out": str(args.out),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
