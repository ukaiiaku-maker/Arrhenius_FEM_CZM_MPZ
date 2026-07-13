#!/usr/bin/env python3
"""Convergence audit for fitted v9.1 three-class moving-PZ parameters."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

from fit_mpz_three_classes import parse_float_list, simulate
from mpz_run_utils import check_parameter_status


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--parameters", required=True)
    ap.add_argument("--classes", default="ceramic weakT DBTT")
    ap.add_argument("--temperatures", default="300 500 700 900 1100")
    ap.add_argument("--dK-values", default="0.20 0.10 0.05")
    ap.add_argument("--bin-counts", default="50 100 200")
    ap.add_argument("--da-values-um", default="10 5")
    ap.add_argument("--target-extension-um", type=float, default=1000.0)
    ap.add_argument("--Kdot", type=float, default=0.005)
    ap.add_argument("--Kmax", type=float, default=65.0)
    ap.add_argument("--out", default="runs/mpz_v9_1_three_class_convergence")
    ap.add_argument("--strict", action="store_true")
    ap.add_argument("--K-init-range-max", type=float, default=1.0)
    ap.add_argument("--K-plateau-range-max", type=float, default=1.5)
    ap.add_argument("--delta-KR-range-max", type=float, default=1.5)
    a = ap.parse_args()

    params = pd.read_csv(a.parameters).set_index("target_class", drop=False)
    check_parameter_status(params, a.parameters, require_fitted=False)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    dK_values = parse_float_list(a.dK_values)
    bins = [int(x) for x in parse_float_list(a.bin_counts)]
    da_values = parse_float_list(a.da_values_um)
    temps = parse_float_list(a.temperatures)
    classes = a.classes.replace(",", " ").split()

    variants: list[dict[str, float | int | str]] = []
    # One-factor-at-a-time around the central variant.
    dK_ref = dK_values[len(dK_values) // 2]
    bins_ref = bins[len(bins) // 2]
    da_ref = da_values[-1]
    for dK in dK_values:
        variants.append({"variant": f"dK_{dK:g}", "dK": dK,
                         "bins": bins_ref, "da_um": da_ref})
    for nb in bins:
        if nb != bins_ref:
            variants.append({"variant": f"bins_{nb}", "dK": dK_ref,
                             "bins": nb, "da_um": da_ref})
    for da in da_values:
        if da != da_ref:
            variants.append({"variant": f"da_{da:g}um", "dK": dK_ref,
                             "bins": bins_ref, "da_um": da})

    rows: list[dict] = []
    events: list[dict] = []
    for klass in classes:
        base = params.loc[klass].copy()
        for T in temps:
            for variant in variants:
                row = base.copy()
                row["mpz_n_bins"] = int(variant["bins"])
                da_um = float(variant["da_um"])
                n_adv = int(round(a.target_extension_um / da_um)) + 1
                opt = SimpleNamespace(
                    dK=float(variant["dK"]), Kdot=a.Kdot,
                    n_advances=n_adv, Kmax=a.Kmax, da_um=da_um,
                    early_window_um=(20.0, 220.0),
                    plateau_window_um=(700.0, 1000.0),
                    target_dB_substep=0.25,
                    target_emission_hazard_substep=1.0,
                    source_active_fraction_min=1.0e-4,
                    min_substep_fraction=1.0e-8,
                    max_substeps=2_000_000,
                )
                q = simulate(row, T, opt)
                metrics = {k: v for k, v in q.items() if k != "events"}
                rows.append({"target_class": klass, "variant": variant["variant"],
                             "dK": variant["dK"], "mpz_n_bins": variant["bins"],
                             "da_um": da_um, **metrics})
                for ev in q["events"]:
                    events.append({"target_class": klass, "T_K": T,
                                   "variant": variant["variant"], **ev})
                print(klass, T, variant["variant"], q["K_init"], q["K_plateau"])

    df = pd.DataFrame(rows)
    edf = pd.DataFrame(events)
    df.to_csv(out / "convergence_metrics.csv", index=False)
    edf.to_csv(out / "convergence_event_Rcurves.csv", index=False)

    summary = []
    failed = False
    for (klass, T), g in df.groupby(["target_class", "T_K"]):
        rec = {"target_class": klass, "T_K": T, "n_variants": len(g)}
        for col, limit in [
            ("K_init", a.K_init_range_max),
            ("K_plateau", a.K_plateau_range_max),
            ("delta_KR", a.delta_KR_range_max),
        ]:
            finite = g[col].to_numpy(dtype=float)
            finite = finite[np.isfinite(finite)]
            spread = float(np.ptp(finite)) if finite.size else float("inf")
            rec[f"{col}_range"] = spread
            rec[f"{col}_pass"] = bool(spread <= limit)
            failed = failed or spread > limit
        summary.append(rec)
    sdf = pd.DataFrame(summary)
    sdf.to_csv(out / "convergence_summary.csv", index=False)
    (out / "run_config.json").write_text(json.dumps(vars(a), indent=2))
    print(sdf.to_string(index=False))
    if a.strict and failed:
        raise SystemExit("convergence thresholds failed")


if __name__ == "__main__":
    main()
