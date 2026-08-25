#!/usr/bin/env python3
"""Run the final shared four-class rows through both reduced providers."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_outputs" / "oneD_v2_terminal_predictive_program"
sys.path.insert(0, str(ROOT))

from scripts.run_oneD_v2_predictive_campaign import TEMPERATURES, inputs, run_case, summary

CLASS_ORDER = ("Peak", "DBTT", "weak-T", "ceramic-like")


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    registry_path = OUT / "oneD_v2_pf_transfer_registry.csv"
    registry = pd.read_csv(registry_path)
    physics, _, providers = inputs()
    counts = {}
    for target_um in (100.0, 1000.0):
        records = []
        for material, (_, row) in zip(CLASS_ORDER, registry.iterrows()):
            for backend, (mechanics, drive) in providers.items():
                for temperature in TEMPERATURES:
                    records.append(summary(run_case(
                        row, material, temperature, backend, mechanics, drive,
                        physics, target_um,
                    )))
        frame = pd.DataFrame(records)
        name = (
            "oneD_v2_final_four_class_results.csv" if target_um == 100.0
            else "oneD_v2_final_four_class_1000um_results.csv"
        )
        frame.to_csv(OUT / name, index=False)
        counts[str(int(target_um))] = {
            "case_count": len(frame),
            "target_complete": int(frame.status.eq("TARGET_RIGHT_CENSORED").sum()),
            "status_counts": frame.status.value_counts().to_dict(),
            "sha256": sha(OUT / name),
        }
        print(f"FINAL_MATRIX target_um={target_um:g} complete={counts[str(int(target_um))]['target_complete']}/{len(frame)}")
    (OUT / "oneD_v2_final_matrix_manifest.json").write_text(json.dumps({
        "schema": "oneD_v2_terminal_final_matrix_v1",
        "same_material_row_both_providers": True,
        "temperatures_K": list(TEMPERATURES),
        "targets_um": counts,
        "new_2D_FEMCZM_runs": 0,
    }, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
