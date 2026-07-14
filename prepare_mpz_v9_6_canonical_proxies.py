#!/usr/bin/env python3
"""Attach exact historical barrier references to evaluated atlas K(T) curves.

The exact four-class reference table contains barrier and reduced-state
parameters but not the complete refined first-passage trajectory.  This utility
matches each reference to the nearest fully evaluated atlas row in normalized
barrier-coordinate space, copies the requested K(T) columns, and records the
proxy candidate ID and distance.  The matching is explicit and auditable; it is
not treated as an exact identity.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


MATCH_COLUMNS = [
    "emit_G00_eV", "emit_gT_eV_per_K", "emit_sigc0_GPa",
    "emit_sT_GPa_per_K", "emit_exp_a", "emit_exp_n", "emit_floor_frac",
    "cleave_G00_eV", "cleave_gT_eV_per_K", "cleave_sigc0_GPa",
    "cleave_sT_GPa_per_K", "cleave_exp_a", "cleave_exp_n",
    "cleave_floor_frac", "cleave_S_hs_kB",
]


def floats(text: str) -> list[float]:
    return [float(x) for x in str(text).replace(",", " ").split() if x]


def prepare(
    atlas_path: Path,
    canonical_path: Path,
    temperatures: list[float],
    output: Path,
) -> pd.DataFrame:
    atlas = pd.read_csv(atlas_path)
    refs = pd.read_csv(canonical_path).copy()
    if "candidate_id" not in atlas:
        atlas["candidate_id"] = [f"atlas_{i:06d}" for i in range(len(atlas))]

    kc_columns: list[str] = []
    for T in temperatures:
        tag = str(int(round(T)))
        name = next(
            (c for c in (f"refined_Kc_T{tag}", f"Kc_T{tag}") if c in atlas),
            None,
        )
        if name is None:
            raise ValueError(f"atlas has no first-passage K column for {T:g} K")
        kc_columns.append(name)

    pool = atlas.dropna(subset=kc_columns).copy()
    common = [c for c in MATCH_COLUMNS if c in pool and c in refs]
    if pool.empty or not common:
        raise ValueError("no complete atlas rows or common barrier coordinates")

    X = pool[common].to_numpy(float)
    scale = np.nanstd(X, axis=0)
    scale = np.where(np.isfinite(scale) & (scale > 1.0e-12), scale, 1.0)

    refs["candidate_id"] = "canonical_" + refs.target_class.astype(str)
    refs["region"] = refs.target_class.map({
        "ceramic": "ceramic_reference",
        "peak": "peak_reference",
        "weakT": "weakT_reference",
        "DBTT": "DBTT_reference",
    })
    refs["candidate_source"] = "prior_first_passage_reference"

    for idx, ref in refs.iterrows():
        y = ref[common].to_numpy(float)
        distance = np.sqrt(np.nanmean(((X - y) / scale) ** 2, axis=1))
        nearest_index = int(np.nanargmin(distance))
        nearest = pool.iloc[nearest_index]
        for name in kc_columns:
            refs.loc[idx, name] = float(nearest[name])
        refs.loc[idx, "canonical_kc_proxy_candidate_id"] = str(
            nearest.candidate_id
        )
        refs.loc[idx, "canonical_kc_proxy_distance"] = float(
            distance[nearest_index]
        )

    output.parent.mkdir(parents=True, exist_ok=True)
    refs.to_csv(output, index=False)
    return refs


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--atlas", type=Path, required=True)
    ap.add_argument("--canonical", type=Path, required=True)
    ap.add_argument("--temperatures", default="300 700 900 1200")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    result = prepare(
        args.atlas,
        args.canonical,
        floats(args.temperatures),
        args.out,
    )
    print(
        result[[
            "target_class", "candidate_id",
            "canonical_kc_proxy_candidate_id",
            "canonical_kc_proxy_distance",
        ]].to_string(index=False),
        flush=True,
    )
    print(f"Prepared canonical proxy table: {args.out}", flush=True)


if __name__ == "__main__":
    main()
