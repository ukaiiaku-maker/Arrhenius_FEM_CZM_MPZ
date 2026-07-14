#!/usr/bin/env python3
"""Audit virgin and seeded MPZ branches after the v9.5 spatial correction.

The audit does not optimize parameters.  It asks whether a physically seeded
near-tip forest-density profile decays, persists, or grows during crack advance,
and whether a virgin state approaches the same branch.  This directly tests the
original DBTT premise: a developed shielding state may exist even when virgin
first passage cannot create it instantaneously.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pandas as pd

from arrhenius_fracture import sharp_front as sf
from fit_mpz_three_classes import simulate
from search_mpz_v9_4_developed_state import (
    REGION_TO_CLASS,
    SHARED_COLUMNS,
    apply_shared,
    parse_temperatures,
)


def parse_float_list(text: str) -> list[float]:
    return [float(x) for x in str(text).replace(",", " ").split() if x]


def run_seeded(
    row: pd.Series,
    T_K: float,
    opt: SimpleNamespace,
    seed_rho_m2: float,
    seed_decay_length_m: float,
    seed_available_site_fraction: float,
) -> tuple[dict[str, Any], Any]:
    original_build = sf.build_engine
    holder: dict[str, Any] = {}

    def seeded_build(args, material):
        eng = original_build(args, material)
        holder["engine"] = eng
        floor = float(eng.mpz_config.pt_forest_density_floor_m2)
        if seed_rho_m2 > floor * (1.0 + 1.0e-12):
            if not hasattr(eng.mpz_state, "initialize_forest_profile"):
                raise RuntimeError(
                    "active MPZ state lacks v9.5 initialize_forest_profile"
                )
            eng.mpz_state.initialize_forest_profile(
                seed_rho_m2,
                decay_length_m=seed_decay_length_m,
                available_site_fraction=seed_available_site_fraction,
            )
            eng._sync_compat()
        return eng

    sf.build_engine = seeded_build
    try:
        result = simulate(row, T_K, opt)
    finally:
        sf.build_engine = original_build
    return result, holder.get("engine")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--selected-rows", type=Path,
        default=Path(
            "runs/mpz_v9_4_peierls_taylor_search_v1/strict_common_selection/"
            "pt_v9_4_recommended_intrinsic_rows.csv"
        ),
    )
    ap.add_argument(
        "--state-search", type=Path,
        default=Path(
            "runs/mpz_v9_4_developed_state_search_v1/"
            "developed_state_search_all.csv"
        ),
    )
    ap.add_argument("--candidate-indices", default="99 68 100")
    ap.add_argument(
        "--temperatures",
        default="ceramic:300,1200;weakT:300,1200;DBTT:300,900,1200",
    )
    ap.add_argument(
        "--seed-rho-values-m2",
        default="5e12 1e13 3e13 1e14 3e14",
    )
    ap.add_argument("--seed-decay-length-um", type=float, default=5.0)
    ap.add_argument(
        "--seed-available-site-fraction", type=float, default=0.5
    )
    ap.add_argument("--target-extension-um", type=float, default=200.0)
    ap.add_argument("--da-um", type=float, default=5.0)
    ap.add_argument("--dK", type=float, default=0.25)
    ap.add_argument("--Kdot", type=float, default=0.005)
    ap.add_argument("--Kmax", type=float, default=65.0)
    ap.add_argument(
        "--out", type=Path,
        default=Path("runs/mpz_v9_5_state_continuation_audit_v1"),
    )
    a = ap.parse_args()

    selected = pd.read_csv(a.selected_rows)
    selected["target_class"] = selected["region"].map(REGION_TO_CLASS)
    base_rows = {
        klass: g.sort_values("pt_score").iloc[0].copy()
        for klass, g in selected.groupby("target_class")
    }
    states = pd.read_csv(a.state_search)
    requested = [int(x) for x in parse_float_list(a.candidate_indices)]
    state_rows = states[states.candidate_index.astype(int).isin(requested)]
    missing = sorted(set(requested) - set(state_rows.candidate_index.astype(int)))
    if missing:
        raise SystemExit(f"state candidates not found: {missing}")

    temperature_map = parse_temperatures(a.temperatures)
    seed_values = parse_float_list(a.seed_rho_values_m2)
    opt = SimpleNamespace(
        dK=float(a.dK),
        Kdot=float(a.Kdot),
        n_advances=int(round(a.target_extension_um / a.da_um)) + 1,
        Kmax=float(a.Kmax),
        da_um=float(a.da_um),
        early_window_um=(20.0, min(100.0, 0.5 * a.target_extension_um)),
        plateau_window_um=(0.70 * a.target_extension_um, a.target_extension_um),
        target_dB_substep=0.25,
        target_emission_hazard_substep=1.0,
        source_active_fraction_min=1.0e-4,
        min_substep_fraction=1.0e-8,
        max_substeps=2_000_000,
        objective_mode="rcurve",
    )

    out = a.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    events: list[dict[str, Any]] = []
    total = len(state_rows) * sum(len(v) for v in temperature_map.values()) * len(seed_values)
    count = 0

    for _, state_row in state_rows.sort_values("candidate_index").iterrows():
        shared = {c: float(state_row[c]) for c in SHARED_COLUMNS}
        for klass in ("ceramic", "weakT", "DBTT"):
            material_row = apply_shared(base_rows[klass], shared)
            for T_K in temperature_map[klass]:
                for seed_rho in seed_values:
                    count += 1
                    result, eng = run_seeded(
                        material_row,
                        float(T_K),
                        opt,
                        float(seed_rho),
                        a.seed_decay_length_um * 1.0e-6,
                        a.seed_available_site_fraction,
                    )
                    if eng is None:
                        raise RuntimeError("seeded engine was not captured")
                    final_rho = eng.mpz_state.local_forest_density_m2()
                    final_Ksh = eng.mpz_state.shielding_K(
                        eng.G, eng.nu, eng.b
                    ) / 1.0e6
                    floor = float(eng.mpz_config.pt_forest_density_floor_m2)
                    initial_excess = max(float(seed_rho) - floor, 0.0)
                    final_excess = max(float(np.max(final_rho)) - floor, 0.0)
                    persistence = (
                        final_excess / initial_excess
                        if initial_excess > 0.0 else float("nan")
                    )
                    rec = {
                        "candidate_index": int(state_row.candidate_index),
                        "target_class": klass,
                        "T_K": float(T_K),
                        "seed_rho_tip_m2": float(seed_rho),
                        "seed_branch": "virgin" if seed_rho <= floor else "developed",
                        **shared,
                        **{k: v for k, v in result.items() if k != "events"},
                        "final_rho_forest_max_m2": float(np.max(final_rho)),
                        "final_rho_forest_median_m2": float(np.median(final_rho)),
                        "final_retained_count": float(eng.mpz_state.retained_count),
                        "final_mobile_count": float(eng.mpz_state.mobile_count),
                        "final_K_shield_MPa_sqrt_m": float(final_Ksh),
                        "seed_excess_density_persistence": float(persistence),
                    }
                    rows.append(rec)
                    for event in result.get("events", []):
                        events.append({
                            "candidate_index": int(state_row.candidate_index),
                            "target_class": klass,
                            "T_K": float(T_K),
                            "seed_rho_tip_m2": float(seed_rho),
                            **event,
                        })
                    print(
                        f"evaluated {count}/{total} candidate={int(state_row.candidate_index)} "
                        f"class={klass} T={T_K:g} seed={seed_rho:.3g} "
                        f"Kp={result['K_plateau']:.4g} "
                        f"rho_final={np.max(final_rho):.3g} "
                        f"Ksh_final={final_Ksh:.4g}",
                        flush=True,
                    )

    df = pd.DataFrame(rows)
    edf = pd.DataFrame(events)
    df.to_csv(out / "state_continuation_metrics.csv", index=False)
    edf.to_csv(out / "state_continuation_events.csv", index=False)

    summary = []
    for keys, g in df.groupby(["candidate_index", "target_class", "T_K"]):
        ordered = g.sort_values("seed_rho_tip_m2")
        virgin = ordered.iloc[0]
        developed = ordered.iloc[-1]
        summary.append({
            "candidate_index": int(keys[0]),
            "target_class": keys[1],
            "T_K": float(keys[2]),
            "virgin_K_plateau": float(virgin.K_plateau),
            "developed_K_plateau": float(developed.K_plateau),
            "branch_gap_K_MPa_sqrt_m": float(
                developed.K_plateau - virgin.K_plateau
            ),
            "virgin_final_rho_max_m2": float(
                virgin.final_rho_forest_max_m2
            ),
            "developed_final_rho_max_m2": float(
                developed.final_rho_forest_max_m2
            ),
            "developed_seed_persistence": float(
                developed.seed_excess_density_persistence
            ),
            "developed_final_K_shield_MPa_sqrt_m": float(
                developed.final_K_shield_MPa_sqrt_m
            ),
        })
    sdf = pd.DataFrame(summary)
    sdf.to_csv(out / "state_continuation_summary.csv", index=False)

    config = vars(a).copy()
    for key, value in list(config.items()):
        if isinstance(value, Path):
            config[key] = str(value)
    config.update({
        "temperature_map": temperature_map,
        "seed_values": seed_values,
        "status": "V9_5_REDUCED_BRANCH_AUDIT_NOT_2D_VALIDATED",
    })
    (out / "state_continuation_config.json").write_text(
        json.dumps(config, indent=2)
    )
    print(sdf.to_string(index=False), flush=True)
    print(f"Outputs: {out}", flush=True)


if __name__ == "__main__":
    main()
