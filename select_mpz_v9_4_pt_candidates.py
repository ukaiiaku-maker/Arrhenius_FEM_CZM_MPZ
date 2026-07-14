#!/usr/bin/env python3
"""Select strict and cross-class PT closures from a completed v9.4 search.

The raw search reports both topology-only and strict-strength rows.  Developed
MPZ calculations must start from the strict subset.  Because the Peierls--
Taylor closure represents bulk/process-zone transport rather than the intrinsic
cleavage/emission class itself, this utility also identifies transport closures
that pass for all requested intrinsic regions.  A common closure is preferred
for the first developed-state comparison so class differences remain traceable
to the intrinsic EXP-floor surfaces rather than unrelated bulk plasticity laws.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


PT_COLUMNS = [
    "pt_peierls_energy_ratio",
    "pt_taylor_energy_ratio",
    "pt_entropy_multiplier",
    "pt_peierls_entropy_ratio",
    "pt_taylor_entropy_ratio",
    "pt_taylor_corr_rho_c",
    "pt_taylor_renewal_time_s",
    "pt_taylor_m_exponent",
    "pt_taylor_m_scale",
    "pt_taylor_m_cap",
    "pt_mobile_saturation_density_m2",
    "pt_mobile_fraction",
]


def _as_bool(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series
    return series.astype(str).str.lower().isin({"1", "true", "yes"})


def select(search_root: Path, out: Path) -> dict[str, object]:
    accepted_path = search_root / "peierls_taylor_search_accepted.csv"
    if not accepted_path.exists():
        raise FileNotFoundError(f"missing completed search table: {accepted_path}")
    df = pd.read_csv(accepted_path)
    missing = [c for c in PT_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"accepted table lacks PT columns: {missing}")

    mask = _as_bool(df["accepted"]) & _as_bool(df["strict_strength_window"])
    if "barrier_valid" in df:
        mask &= _as_bool(df["barrier_valid"])
    if "detailed_balance_valid" in df:
        mask &= _as_bool(df["detailed_balance_valid"])
    strict = df.loc[mask].copy()
    if strict.empty:
        raise RuntimeError("no strict detailed-balance PT rows were found")

    regions = sorted(strict["region"].astype(str).unique())
    closure = (
        strict.groupby(PT_COLUMNS, dropna=False)
        .agg(
            n_rows=("candidate_id", "size"),
            n_regions=("region", "nunique"),
            n_candidates=("candidate_id", "nunique"),
            mean_pt_score=("pt_score", "mean"),
            max_pt_score=("pt_score", "max"),
            min_sigma_ref_700K_1e14_GPa=(
                "sigma_ref_700K_1e14_GPa", "min"
            ),
            max_sigma_ref_700K_1e14_GPa=(
                "sigma_ref_700K_1e14_GPa", "max"
            ),
            max_sigma_over_grid_GPa=("sigma_max_GPa", "max"),
            min_raw_scaled_G0_eV=("min_raw_scaled_G0_eV", "min"),
            max_zero_stress_rate_s=("zero_stress_rate_max_s", "max"),
        )
        .reset_index()
        .sort_values(
            ["n_regions", "n_candidates", "mean_pt_score", "max_sigma_over_grid_GPa"],
            ascending=[False, False, True, True],
        )
        .reset_index(drop=True)
    )
    closure.insert(0, "closure_rank", np.arange(1, len(closure) + 1))
    closure["covers_all_regions"] = closure["n_regions"] == len(regions)

    all_region = closure[closure["covers_all_regions"]].copy()
    if all_region.empty:
        raise RuntimeError(
            f"no single PT closure passed all strict regions: {regions}"
        )
    recommended = all_region.iloc[[0]].copy()
    key = recommended[PT_COLUMNS]
    recommended_rows = strict.merge(key, on=PT_COLUMNS, how="inner")
    recommended_rows = (
        recommended_rows.sort_values(["region", "pt_score"])
        .groupby("region", as_index=False, sort=True)
        .head(1)
        .reset_index(drop=True)
    )
    recommended_rows["status"] = (
        "PT_STRICT_COMMON_CLOSURE_REQUIRES_DEVELOPED_MPZ_VALIDATION"
    )

    out.mkdir(parents=True, exist_ok=True)
    strict.to_csv(out / "pt_v9_4_strict_rows.csv", index=False)
    closure.to_csv(out / "pt_v9_4_common_closure_ranking.csv", index=False)
    recommended.to_csv(out / "pt_v9_4_recommended_common_closure.csv", index=False)
    recommended_rows.to_csv(
        out / "pt_v9_4_recommended_intrinsic_rows.csv", index=False
    )

    report = {
        "search_root": str(search_root),
        "regions": regions,
        "n_accepted_input": int(len(df)),
        "n_strict_rows": int(len(strict)),
        "n_unique_common_closures": int(len(all_region)),
        "recommended_closure_rank": int(recommended.iloc[0]["closure_rank"]),
        "recommended_n_candidates": int(recommended.iloc[0]["n_candidates"]),
        "recommended_mean_pt_score": float(
            recommended.iloc[0]["mean_pt_score"]
        ),
        "recommended_max_sigma_over_grid_GPa": float(
            recommended.iloc[0]["max_sigma_over_grid_GPa"]
        ),
        "recommended_min_raw_scaled_G0_eV": float(
            recommended.iloc[0]["min_raw_scaled_G0_eV"]
        ),
        "recommended_candidate_ids": recommended_rows[
            ["region", "candidate_id"]
        ].to_dict(orient="records"),
    }
    (out / "pt_v9_4_selection_summary.json").write_text(
        json.dumps(report, indent=2)
    )
    return report


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--search-root",
        type=Path,
        default=Path("runs/mpz_v9_4_peierls_taylor_search_v1"),
    )
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args()
    out = args.out or args.search_root / "strict_common_selection"
    report = select(args.search_root, out)
    print(json.dumps(report, indent=2), flush=True)


if __name__ == "__main__":
    main()
