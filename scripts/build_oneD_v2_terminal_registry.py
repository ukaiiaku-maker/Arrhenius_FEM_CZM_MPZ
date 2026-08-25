#!/usr/bin/env python3
"""Materialize the final shared registry and PF-transfer registry."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_outputs" / "oneD_v2_terminal_predictive_program"
SOURCE = Path(
    "/Volumes/Data/Data/Nanopillar_calculation/"
    "PF-fracture-fatigue_v10_2_21_persistent_sites_top1/"
    "arrhenius_fracture/data/materials/v10_2_27_v913_four_class_paper_registry.csv"
)
FINALISTS = {
    "Peak": "v913_zeroD_sobol_0242980",
    "DBTT": "v913_zeroD_sobol_0202500",
    "weak-T": "oneD_v2_focused_weak_T_0016",
    "ceramic-like": "oneD_v2_focused_ceramic_like_0018",
}


def main() -> int:
    source = pd.read_csv(SOURCE)
    pareto = pd.read_csv(OUT / "oneD_v2_pareto_candidates.csv").set_index("candidate_id")
    templates = {
        "Peak": source[source.candidate_id == "v913_zeroD_sobol_0242980"].iloc[0].copy(),
        "DBTT": source[source.candidate_id == "v913_zeroD_sobol_0202500"].iloc[0].copy(),
        "weak-T": source[source.candidate_id == "v913_zeroD_sobol_0129902"].iloc[0].copy(),
        "ceramic-like": source[source.candidate_id == "v913_zeroD_sobol_0077080"].iloc[0].copy(),
    }
    records = []
    transfer = []
    for material, candidate_id in FINALISTS.items():
        template = templates[material]
        if candidate_id in pareto.index:
            candidate = pareto.loc[candidate_id]
            for field in candidate.index.intersection(template.index):
                if field in {"candidate_id", "search_class", "search_source", "search_anchor"}:
                    continue
                if pd.notna(candidate[field]) and pd.api.types.is_number(candidate[field]):
                    template[field] = candidate[field]
            template["candidate_id"] = candidate_id
            template["validation_status"] = (
                "V2 terminal two-provider reduced finalist; bounded direct PF transfer required."
            )
        record = {
            "material_class": material,
            "candidate_id": candidate_id,
            "source_or_search_provenance": (
                "HISTORICAL_V913_RETAINED" if candidate_id.startswith("v913_")
                else str(pareto.loc[candidate_id, "search_source"])
            ),
            "parameter_decision": (
                "RETAINED_QUALIFIED_NATIVE_BASELINE" if material in {"Peak", "DBTT"}
                else "NEW_SHARED_TWO_PROVIDER_FINALIST"
            ),
            "reduced_selection_status": (
                "FULL_CLASS_CONTRACT_PASS" if material == "weak-T"
                else "CREDIBLE_PROVIDER_SENSITIVE" if material == "ceramic-like"
                else "NATIVE_BASELINE_TOPOLOGY_AND_ONSET_PASS"
            ),
            "same_material_row_both_providers": True,
        }
        for field in template.index:
            if field in {"candidate_id", "material_class"}:
                continue
            if isinstance(template[field], (int, float)):
                record[field] = template[field]
        records.append(record)
        transfer.append(template)
    pd.DataFrame(records).to_csv(OUT / "oneD_v2_new_four_class_registry.csv", index=False)
    transfer_frame = pd.DataFrame(transfer)[source.columns]
    transfer_frame.to_csv(OUT / "oneD_v2_pf_transfer_registry.csv", index=False)
    selection = {
        "schema": "oneD_v2_terminal_pf_transfer_selection_v1",
        "canonical_option_order": transfer_frame.option_key.tolist(),
        "primary_candidates": [
            {"option_key": row.option_key, "candidate_id": row.candidate_id,
             "material_class": material}
            for material, row in zip(FINALISTS, transfer_frame.itertuples())
        ],
        "analysis_only": True,
        "canonical_PF_registry_modified": False,
    }
    (OUT / "oneD_v2_pf_transfer_selection.json").write_text(
        json.dumps(selection, indent=2, sort_keys=True) + "\n"
    )
    print("FINAL_REGISTRY_COMPLETE")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
