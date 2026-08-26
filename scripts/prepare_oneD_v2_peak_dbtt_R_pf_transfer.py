#!/usr/bin/env python3
"""Prepare the immutable four-row control/finalist registry for bounded PF transfer."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_outputs" / "oneD_v2_peak_dbtt_rcurve_search"

ROWS = (
    ("Peak", "Peak_control", "v913_zeroD_sobol_0242980", "peak"),
    ("Peak", "Peak_R", "oneD_v2_peak_R_41f8789bcbc1f097", "peak"),
    ("DBTT", "DBTT_control", "v913_zeroD_sobol_0202500", "dbtt"),
    ("DBTT", "DBTT_R", "oneD_v2_dbtt_R_9f5160f509e713e2", "dbtt"),
)


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    frames = {
        slug: pd.read_csv(OUT / f"oneD_v2_{slug}_R_candidates.csv")
        for slug in ("peak", "dbtt")
    }
    records = []
    selection = []
    for material, option, candidate_id, slug in ROWS:
        source = frames[slug][frames[slug].candidate_id == candidate_id]
        if len(source) != 1:
            raise ValueError(f"missing unique PF transfer row {candidate_id}")
        row = source.iloc[0].copy()
        row["option_key"] = f"oneD_v2_{option}"
        row["material_class"] = "peak" if material == "Peak" else "DBTT"
        row["role"] = "focused R-curve control" if option.endswith("control") else "focused R-curve finalist"
        row["validation_status"] = "BOUNDED_DIRECT_PF_VALIDATION_REQUIRED"
        # Search-only identity fields are not accepted by the production PF
        # material registry parser and are not material coordinates.
        records.append(row)
        selection.append({
            "material_class": material,
            "option_key": row["option_key"],
            "candidate_id": candidate_id,
            "parameter_sha256": row["parameter_sha256"],
        })
    registry = pd.DataFrame(records)
    search_only = [
        "target_response_class", "search_campaign_id", "parent_or_anchor_id",
        "search_generation_method", "active_parameter_schema",
        "active_parameter_count", "canonical_parameter_json", "parameter_sha256",
    ]
    registry_path = OUT / "oneD_v2_peak_dbtt_R_pf_transfer_registry.csv"
    registry.drop(columns=search_only).to_csv(registry_path, index=False)
    selection_path = OUT / "oneD_v2_peak_dbtt_R_pf_transfer_selection.json"
    selection_path.write_text(json.dumps({
        "schema": "oneD_v2_peak_dbtt_R_pf_transfer_selection_v1",
        "analysis_only": True,
        "canonical_PF_registry_modified": False,
        "primary_candidates": selection,
        "canonical_option_order": [item["option_key"] for item in selection],
        "new_case_policy": "FINALIST_ROWS_ONLY_CONTROLS_REUSE_AUTHORITATIVE_RUNS",
        "temperatures_K": {
            "Peak": [600.0, 900.0, 1200.0],
            "DBTT": [600.0, 1100.0, 1200.0],
        },
        "seeds": {"Peak": 8666, "DBTT": 1008666},
        "target_extension_um": 100.0,
        "maximum_concurrent_heavy_workers": 2,
        "new_FEMCZM_runs": 0,
    }, indent=2, sort_keys=True) + "\n")
    manifest = {
        "schema": "oneD_v2_peak_dbtt_R_pf_transfer_inputs_v1",
        "registry_sha256": _sha(registry_path),
        "selection_sha256": _sha(selection_path),
        "row_count": len(registry),
        "source_candidate_files": {
            slug: _sha(OUT / f"oneD_v2_{slug}_R_candidates.csv")
            for slug in ("peak", "dbtt")
        },
        "canonical_PF_registry_modified": False,
    }
    (OUT / "oneD_v2_peak_dbtt_R_pf_transfer_inputs_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    print(f"PF_TRANSFER_INPUTS_COMPLETE rows={len(registry)} registry={registry_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
