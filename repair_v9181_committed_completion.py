#!/usr/bin/env python3
"""Repair v9.18.1 campaign metadata without rerunning FEM cases."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from arrhenius_fracture.material_rcurve_audit_v913 import audit_campaign
from arrhenius_fracture.mpz_parameterization_v911 import normalize_class_name
from run_mpz_v9_18_2_persistent_plastic_wake import promote_committed_completion


def values(text: str, cast=str):
    return [cast(x) for x in str(text).replace(",", " ").split() if x]


def read_summary(case_dir: Path) -> dict:
    for name in (
        "v9_18_case_summary.json",
        "v9_17_case_summary.json",
        "v9_16_case_summary.json",
        "v9_13_case_summary.json",
    ):
        path = case_dir / name
        if path.exists():
            try:
                value = json.loads(path.read_text())
            except Exception:
                continue
            if isinstance(value, dict):
                return value
    return {}


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--campaign-root", type=Path, required=True)
    p.add_argument("--T-K", type=float, required=True)
    p.add_argument("--target-extension-um", type=float, required=True)
    p.add_argument("--seeds", default="1")
    p.add_argument("--classes", default="ceramic weakT DBTT")
    p.add_argument("--bulk-mode", default="tip_only")
    p.add_argument("--crystal-theta-deg", type=float, default=45.0)
    args = p.parse_args()

    root = args.campaign_root.resolve()
    rows = []
    for seed in values(args.seeds, int):
        for cls_raw in values(args.classes, str):
            cls = normalize_class_name(cls_raw)
            case_dir = (
                root / f"seed_{seed}" / args.bulk_mode / cls
                / f"T{int(round(args.T_K))}_th{args.crystal_theta_deg:g}"
            )
            row = read_summary(case_dir)
            if not row:
                print(f"MISSING {case_dir}")
                continue
            fixed = promote_committed_completion(row, args.target_extension_um)
            rows.append(fixed)
            print(
                f"{cls:7s} status={fixed.get('status')} "
                f"committed={fixed.get('analysis_committed_extension_um')} um "
                f"promoted={fixed.get('committed_completion_promoted_v9182', False)}"
            )
        audit = audit_campaign(
            root,
            seed,
            args.T_K,
            classes=[normalize_class_name(x) for x in values(args.classes, str)],
            bulk_mode=args.bulk_mode,
        )
        print(
            f"AUDIT seed={seed}: transfer_gate={audit['material_transfer_gate_passed']} "
            f"interpretation={audit['interpretation']} "
            f"failed={audit['failed_solver_cases']} incomplete={audit['incomplete_cases']}"
        )

    if rows:
        pd.DataFrame(rows).to_csv(root / "v9_18_2_campaign_summary.csv", index=False)
        (root / "v9_18_2_campaign_summary.json").write_text(
            json.dumps(rows, indent=2, default=str)
        )


if __name__ == "__main__":
    main()
