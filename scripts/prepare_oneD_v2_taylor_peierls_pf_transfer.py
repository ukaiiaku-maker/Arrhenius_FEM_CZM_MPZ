#!/usr/bin/env python3
"""Build the isolated four-option PF transfer registry for kinetic finalists."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_outputs" / "oneD_v2_taylor_peierls_rcurve_search"
TEMPLATE = (
    ROOT / "analysis_outputs" / "oneD_v2_peak_dbtt_rcurve_search"
    / "oneD_v2_peak_dbtt_R_pf_transfer_registry.csv"
)
FINALISTS = {
    "Peak": ("v913_zeroD_sobol_0242980", "oneD_v2_peak_TP_6962e84b2eb78fbb"),
    "DBTT": ("v913_zeroD_sobol_0202500", "oneD_v2_dbtt_TP_6ca03f05fbae34e9"),
}
OPTION_KEYS = {
    "Peak": ("oneD_v2_Peak_control", "oneD_v2_Peak_TP"),
    "DBTT": ("oneD_v2_DBTT_control", "oneD_v2_DBTT_TP"),
}
TEMPERATURES = {"Peak": (600.0, 1000.0, 1200.0), "DBTT": (600.0, 1100.0, 1200.0)}
SEEDS = {"Peak": 8666, "DBTT": 1008666}


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    template = pd.read_csv(TEMPLATE)
    candidates = pd.read_csv(OUT / "oneD_v2_taylor_peierls_stage1_selection.csv")
    records = []
    selection = []
    for material in ("Peak", "DBTT"):
        template_control = template[
            template.candidate_id.eq(FINALISTS[material][0])
        ].iloc[0].copy()
        for index, (candidate_id, option_key) in enumerate(
            zip(FINALISTS[material], OPTION_KEYS[material])
        ):
            source = candidates[candidates.candidate_id.eq(candidate_id)]
            if len(source) != 1:
                raise RuntimeError(f"missing unique PF transfer material {candidate_id}")
            source = source.iloc[0]
            row = template_control.copy()
            row["option_key"] = option_key
            row["candidate_id"] = candidate_id
            row["material_class"] = material.lower() if material == "Peak" else material
            row["role"] = "Taylor/Peierls control" if index == 0 else "Taylor/Peierls kinetic finalist"
            row["mechanism_summary"] = (
                "Fixed Peak/DBTT cleavage and emission barriers; only source-owned "
                "Taylor/Peierls kinetics differ."
            )
            row["validation_status"] = "BOUNDED_DIRECT_PF_VALIDATION_REQUIRED"
            for field in candidates.columns:
                if field in row.index and pd.notna(source[field]):
                    row[field] = source[field]
            records.append(row)
            selection.append({
                "candidate_id": candidate_id, "material_class": material,
                "option_key": option_key,
                "full_material_sha256": str(source.full_material_sha256),
                "cleavage_barrier_sha256": str(source.cleavage_barrier_sha256),
                "emission_barrier_sha256": str(source.emission_barrier_sha256),
                "taylor_peierls_subvector_sha256": str(source.taylor_peierls_subvector_sha256),
            })
    registry = pd.DataFrame(records, columns=template.columns)
    registry_path = OUT / "oneD_v2_taylor_peierls_pf_transfer_registry.csv"
    registry.to_csv(registry_path, index=False)
    payload = {
        "schema": "oneD_v2_taylor_peierls_pf_transfer_selection_v1",
        "analysis_only": True, "canonical_PF_registry_modified": False,
        "canonical_option_order": [item for material in ("Peak", "DBTT") for item in OPTION_KEYS[material]],
        "primary_candidates": selection, "temperatures_K": TEMPERATURES,
        "seeds": SEEDS, "target_extension_um": 100.0,
        "maximum_concurrent_heavy_workers": 2, "new_FEMCZM_runs": 0,
        "control_policy": "SAME_CANONICAL_SEED_NEW_PAIRED_CONTROL_AND_VARIANT_RUNS",
        "producer_commit": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
    }
    selection_path = OUT / "oneD_v2_taylor_peierls_pf_transfer_selection.json"
    selection_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    manifest = {
        "schema": "oneD_v2_taylor_peierls_pf_transfer_inputs_v1",
        "registry_sha256": _sha(registry_path), "selection_sha256": _sha(selection_path),
        "option_count": int(len(registry)), "source_template_sha256": _sha(TEMPLATE),
    }
    (OUT / "oneD_v2_taylor_peierls_pf_transfer_inputs_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    print(registry[["option_key", "candidate_id", "material_class", "role"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
