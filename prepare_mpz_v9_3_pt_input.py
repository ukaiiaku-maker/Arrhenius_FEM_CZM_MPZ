#!/usr/bin/env python3
"""Join v9.2 refined predictions with complete material rows for v9.3.

The refined analytical shortlist intentionally stores the six sampled barrier
coordinates plus refined first-passage predictions. Fixed EXP-floor shape
parameters are stored in ``mpz_analytic_shortlist_material_rows.csv``. The
Peierls--Taylor search needs both. This utility joins them by analytical
candidate ID and fails early with a clear schema report when a required field
is still absent.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


REQUIRED_EMISSION_COLUMNS = (
    "emit_G00_eV",
    "emit_gT_eV_per_K",
    "emit_sigc0_GPa",
    "emit_exp_a",
    "emit_exp_n",
    "emit_floor_frac",
)


def resolve_material_rows(atlas_shortlist: Path, material_rows: Path | None) -> Path:
    if material_rows is not None:
        return material_rows
    return atlas_shortlist.parent / "mpz_analytic_shortlist_material_rows.csv"


def prepare_input(
    atlas_shortlist: Path,
    material_rows: Path,
    output: Path,
) -> pd.DataFrame:
    if not atlas_shortlist.exists():
        raise FileNotFoundError(f"atlas shortlist not found: {atlas_shortlist}")
    if not material_rows.exists():
        raise FileNotFoundError(
            "complete material rows not found: "
            f"{material_rows}. Expected the v9.2 atlas output "
            "mpz_analytic_shortlist_material_rows.csv."
        )

    atlas = pd.read_csv(atlas_shortlist)
    materials = pd.read_csv(material_rows)
    if "candidate_id" not in atlas.columns:
        raise ValueError(
            f"{atlas_shortlist} lacks candidate_id; columns={list(atlas.columns)}"
        )

    rename = {}
    if "analytic_candidate_id" in materials.columns:
        rename["analytic_candidate_id"] = "candidate_id"
    if "analytic_region" in materials.columns and "region" not in materials.columns:
        rename["analytic_region"] = "region"
    if "analytic_shape_family" in materials.columns and "shape_family" not in materials.columns:
        rename["analytic_shape_family"] = "shape_family"
    materials = materials.rename(columns=rename)
    if "candidate_id" not in materials.columns:
        raise ValueError(
            f"{material_rows} lacks analytic_candidate_id/candidate_id; "
            f"columns={list(materials.columns)}"
        )

    # Candidate IDs include family and loading-rate tags and are unique in the
    # atlas. Keep one complete material row per ID and preserve refined atlas
    # predictions and sampled coordinates as the authoritative values.
    materials = materials.drop_duplicates("candidate_id", keep="last")
    joined = atlas.merge(
        materials,
        on="candidate_id",
        how="left",
        suffixes=("", "__material"),
        validate="many_to_one",
        indicator=True,
    )
    unmatched = joined.loc[joined["_merge"] != "both", "candidate_id"].astype(str)
    if not unmatched.empty:
        examples = ", ".join(unmatched.head(8))
        raise ValueError(
            f"{len(unmatched)} atlas candidates lack complete material rows; "
            f"examples: {examples}"
        )
    joined = joined.drop(columns=["_merge"])

    material_suffix = "__material"
    for col in list(joined.columns):
        if not col.endswith(material_suffix):
            continue
        base = col[: -len(material_suffix)]
        if base in joined.columns:
            joined[base] = joined[base].where(joined[base].notna(), joined[col])
            joined = joined.drop(columns=[col])
        else:
            joined = joined.rename(columns={col: base})

    missing = [c for c in REQUIRED_EMISSION_COLUMNS if c not in joined.columns]
    null = [c for c in REQUIRED_EMISSION_COLUMNS if c in joined.columns and joined[c].isna().any()]
    if missing or null:
        raise ValueError(
            "incomplete emission EXP-floor schema after joining atlas and material rows; "
            f"missing={missing}, contains_null={null}, columns={list(joined.columns)}"
        )

    output.parent.mkdir(parents=True, exist_ok=True)
    joined.to_csv(output, index=False)
    return joined


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--atlas-shortlist", type=Path, required=True)
    ap.add_argument("--material-rows", type=Path, default=None)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    material_rows = resolve_material_rows(args.atlas_shortlist, args.material_rows)
    joined = prepare_input(args.atlas_shortlist, material_rows, args.out)
    print(
        "Prepared v9.3 PT input: "
        f"rows={len(joined)} columns={len(joined.columns)} "
        f"atlas={args.atlas_shortlist} materials={material_rows} out={args.out}",
        flush=True,
    )


if __name__ == "__main__":
    main()
