#!/usr/bin/env python3
"""Broad analytical DBTT-capacity map using the uncapped v9.6 PT model.

Unlike the v9.4 screen, this stage does not preselect five intrinsic rows, does
not require one common PT closure before developed-state behavior is assessed,
and does not rank candidates by an artificially capped flow-stress curve.

The full refined first-passage candidate table is combined with the exact prior
four-class EXP-floor references.  For each intrinsic row and PT sample, a
transparent reduced state balance estimates whether emission can build and
retain a dislocation population on the scale that previously generated a DBTT:

    N_ref = 1488.7579558862624
    chi_ref = 0.7899746000194677

These are benchmark coordinates from the prior first-passage DBTT model, not
production caps.  Production saturation must emerge from finite sources,
escape, recovery, transport, and back stress in the moving process zone.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import qmc

from arrhenius_fracture.emission_derived_plasticity import (
    CorrelatedTaylorConfig,
    EmissionDerivedPeierlsTaylorConfig,
    ExpFloorSurface,
    MechanismScale,
)
from arrhenius_fracture.emission_derived_plasticity_v96 import (
    EmissionDerivedPeierlsTaylorModel,
)


def floats(text: str) -> list[float]:
    return [float(x) for x in str(text).replace(",", " ").split() if x]


def log_scale(u: np.ndarray, lo: float, hi: float) -> np.ndarray:
    return 10.0 ** (
        math.log10(lo) + u * (math.log10(hi) - math.log10(lo))
    )


def sample_pt(n: int, seed: int) -> pd.DataFrame:
    sampler = qmc.Sobol(7, scramble=True, seed=seed)
    if n > 0 and (n & (n - 1)) == 0:
        u = sampler.random_base2(int(math.log2(n)))
    else:
        u = sampler.random(n)
    rows = pd.DataFrame({
        "pt_peierls_energy_ratio": log_scale(u[:, 0], 0.002, 0.02),
        "pt_taylor_energy_ratio": log_scale(u[:, 1], 0.005, 0.08),
        "pt_entropy_multiplier": log_scale(u[:, 2], 0.25, 8.0),
        "pt_correlation_rho_c_m2": log_scale(u[:, 3], 1.0e10, 1.0e16),
        "pt_correlation_scale": log_scale(u[:, 4], 0.2, 5.0),
        "pt_correlation_exponent": 0.7 + 0.6 * u[:, 5],
        "pt_mobile_fraction": log_scale(u[:, 6], 1.0e-4, 0.1),
    })
    anchors = pd.DataFrame([
        {
            "pt_peierls_energy_ratio": 0.005,
            "pt_taylor_energy_ratio": 0.02,
            "pt_entropy_multiplier": e,
            "pt_correlation_rho_c_m2": 1.0e14,
            "pt_correlation_scale": 1.0,
            "pt_correlation_exponent": 1.0,
            "pt_mobile_fraction": 0.01,
        }
        for e in (0.5, 1.0, 2.0, 4.0)
    ])
    return pd.concat([anchors, rows], ignore_index=True)


def surface(row: pd.Series) -> ExpFloorSurface:
    return ExpFloorSurface(
        G00_eV=float(row.emit_G00_eV),
        gT_eV_per_K=float(row.emit_gT_eV_per_K),
        sigc0_Pa=float(row.emit_sigc0_GPa) * 1.0e9,
        sT_Pa_per_K=float(row.get("emit_sT_GPa_per_K", 0.0)) * 1.0e9,
        Tref_K=float(row.get("emit_Tref_K", 481.33)),
        a=float(row.emit_exp_a),
        n=float(row.emit_exp_n),
        floor_fraction=float(row.emit_floor_frac),
        floor_min_eV=float(row.get("emit_floor_min_eV", 1.0e-4)),
        floor_max_fraction=float(row.get("emit_floor_max_frac", 0.95)),
    )


def model(row: pd.Series, p: pd.Series) -> EmissionDerivedPeierlsTaylorModel:
    ep = float(p.pt_peierls_energy_ratio)
    et = float(p.pt_taylor_energy_ratio)
    em = float(p.pt_entropy_multiplier)
    return EmissionDerivedPeierlsTaylorModel(
        EmissionDerivedPeierlsTaylorConfig(
            parent=surface(row),
            peierls=MechanismScale(ep, ep * em, 1.0, 1.0e12),
            taylor=MechanismScale(et, et * em, 1.0, 1.0e11),
            correlated_taylor=CorrelatedTaylorConfig(
                rho_c_m2=float(p.pt_correlation_rho_c_m2),
                renewal_time_s=1.0,
                m_exponent=float(p.pt_correlation_exponent),
                m_scale=float(p.pt_correlation_scale),
                m_cap=float("inf"),
            ),
            mobile_fraction_low_density=float(p.pt_mobile_fraction),
            mobile_saturation_density_m2=float("inf"),
            mobile_density_floor_m2=0.0,
            jump_fraction_of_forest_spacing=1.0,
            jump_length_min_m=0.0,
            taylor_phi_max=float("inf"),
            rate_cap_s=float("inf"),
        )
    )


def kc_value(row: pd.Series, T: float) -> float:
    tag = f"{int(round(T))}"
    for name in (f"refined_Kc_T{tag}", f"Kc_T{tag}"):
        if name in row and pd.notna(row[name]):
            return float(row[name])
    return float("nan")


def target_curve(path: Path, temperatures: list[float]) -> dict[float, float]:
    if not path.exists():
        return {300.0: 15.0, 700.0: 18.0, 900.0: 43.0, 1200.0: 50.0}
    df = pd.read_csv(path)
    class_col = "target_class" if "target_class" in df else "class"
    g = df[df[class_col].astype(str).str.upper() == "DBTT"].copy()
    if g.empty:
        return {300.0: 15.0, 700.0: 18.0, 900.0: 43.0, 1200.0: 50.0}
    value_col = next(
        (c for c in (
            "target_K_plateau", "target_Kc_MPa_sqrtm",
            "target_K_MPa_sqrt_m", "target_K_init"
        ) if c in g),
        None,
    )
    if value_col is None or "T_K" not in g:
        raise ValueError(f"cannot identify DBTT target columns in {path}")
    source_T = g.T_K.to_numpy(float)
    source_K = g[value_col].to_numpy(float)
    return {
        float(T): float(np.interp(T, source_T, source_K))
        for T in temperatures
    }


def normalize_canonical(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df = df.copy()
    df["candidate_id"] = "canonical_" + df.target_class.astype(str)
    df["region"] = df.target_class.map({
        "ceramic": "ceramic_reference",
        "peak": "peak_reference",
        "weakT": "weakT_reference",
        "DBTT": "DBTT_reference",
    })
    df["candidate_source"] = "prior_first_passage_reference"
    return df


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--atlas", type=Path,
        default=Path(
            "runs/mpz_v9_2_analytic_first_passage_atlas/"
            "analytic_first_passage_atlas_shortlist_refined.csv"
        ),
    )
    ap.add_argument(
        "--prepared-atlas", type=Path,
        default=Path(
            "runs/mpz_v9_4_peierls_taylor_search_v1/"
            "pt_search_input_joined.csv"
        ),
    )
    ap.add_argument(
        "--canonical", type=Path,
        default=Path("mpz_v9_6_canonical_first_passage_references.csv"),
    )
    ap.add_argument(
        "--targets", type=Path,
        default=Path("mpz_three_class_design_targets.csv"),
    )
    ap.add_argument("--temperatures", default="300 700 900 1200")
    ap.add_argument("--pt-samples", type=int, default=64)
    ap.add_argument("--seed", type=int, default=96061)
    ap.add_argument("--max-intrinsic", type=int, default=0)
    ap.add_argument("--r-pz-m", type=float, default=1.0e-6)
    ap.add_argument("--b-m", type=float, default=2.74e-10)
    ap.add_argument("--nu0-emission", type=float, default=1.0e12)
    ap.add_argument("--recovery-rate-s", type=float, default=1.0e-5)
    ap.add_argument("--backstress-unit-MPa-sqrt-m-per-sqrt-N", type=float, default=1.4)
    ap.add_argument("--top-count", type=int, default=200)
    ap.add_argument(
        "--out", type=Path,
        default=Path("runs/mpz_v9_6_broad_dbtt_map_v1"),
    )
    a = ap.parse_args()

    source = a.prepared_atlas if a.prepared_atlas.exists() else a.atlas
    if not source.exists():
        raise SystemExit(f"analytical atlas input not found: {source}")
    intrinsic = pd.read_csv(source).copy()
    intrinsic["candidate_source"] = "refined_analytic_atlas"
    if "candidate_id" not in intrinsic:
        intrinsic["candidate_id"] = [f"atlas_{i:06d}" for i in range(len(intrinsic))]
    canonical = normalize_canonical(a.canonical)
    candidates = pd.concat([intrinsic, canonical], ignore_index=True, sort=False)
    if a.max_intrinsic > 0:
        candidates = pd.concat([
            candidates[candidates.candidate_source == "prior_first_passage_reference"],
            candidates[candidates.candidate_source != "prior_first_passage_reference"].head(a.max_intrinsic),
        ], ignore_index=True)

    temperatures = floats(a.temperatures)
    targets = target_curve(a.targets, temperatures)
    dbtt_ref = canonical[canonical.target_class == "DBTT"].iloc[0]
    N_ref = float(dbtt_ref.N_sat)
    chi_ref = float(dbtt_ref.chi_shield)
    rho_ref = N_ref / (math.pi * a.r_pz_m ** 2)
    pt = sample_pt(a.pt_samples, a.seed)
    out = a.out.resolve()
    out.mkdir(parents=True, exist_ok=True)

    summary_rows = []
    temp_rows = []
    total = len(candidates) * len(pt)
    count = 0
    for _, row in candidates.iterrows():
        for _, p in pt.iterrows():
            count += 1
            m = model(row, p)
            predicted = []
            valid = True
            for T in temperatures:
                K0 = kc_value(row, T)
                if not np.isfinite(K0):
                    valid = False
                    break
                sigma_tip = K0 * 1.0e6 / math.sqrt(2.0 * math.pi * a.r_pz_m)
                G_emit = m.barrier_eV("emission", sigma_tip, T)
                emit_rate = float(m._arrhenius_rate(G_emit, T, a.nu0_emission))
                rates = m.rates(sigma_tip, rho_ref, T, a.b_m)
                escape_rate = float(np.asarray(rates["series_rate_s"]))
                denom = emit_rate + escape_rate + max(a.recovery_rate_s, 0.0)
                retained_fraction = emit_rate / denom if denom > 0.0 else 0.0
                N_developed = N_ref * retained_fraction
                K_shield = (
                    chi_ref
                    * a.backstress_unit_MPa_sqrt_m_per_sqrt_N
                    * math.sqrt(max(N_developed, 0.0))
                )
                K_developed = K0 + K_shield
                predicted.append(K_developed)
                temp_rows.append({
                    "candidate_id": row.candidate_id,
                    "candidate_source": row.candidate_source,
                    "region": row.get("region", ""),
                    **p.to_dict(),
                    "T_K": T,
                    "K_intrinsic": K0,
                    "K_target_DBTT": targets[T],
                    "K_developed_proxy": K_developed,
                    "K_shield_proxy": K_shield,
                    "emission_rate_s": emit_rate,
                    "escape_rate_s": escape_rate,
                    "retained_fraction": retained_fraction,
                    "N_developed_proxy": N_developed,
                    "rho_reference_m2": rho_ref,
                    "taylor_m_eff": float(np.asarray(rates["taylor_m_eff"])),
                    "taylor_amplification": float(np.asarray(rates["taylor_amplification"])),
                    "caps_active": False,
                })
            if valid:
                pred = np.asarray(predicted)
                target = np.asarray([targets[T] for T in temperatures])
                rmse = float(np.sqrt(np.mean((pred - target) ** 2)))
                rise = float(pred[-1] - pred[0])
                intrinsic_rise = float(
                    kc_value(row, temperatures[-1]) - kc_value(row, temperatures[0])
                )
                score = rmse + 0.5 * max(20.0 - rise, 0.0) + 0.25 * max(intrinsic_rise, 0.0)
                summary_rows.append({
                    "candidate_id": row.candidate_id,
                    "candidate_source": row.candidate_source,
                    "region": row.get("region", ""),
                    **p.to_dict(),
                    "dbtt_proxy_score": score,
                    "dbtt_proxy_rmse": rmse,
                    "developed_rise_MPa_sqrt_m": rise,
                    "intrinsic_rise_MPa_sqrt_m": intrinsic_rise,
                    "K_developed_300": pred[0],
                    "K_developed_1200": pred[-1],
                    "N_reference": N_ref,
                    "chi_reference": chi_ref,
                    "rho_reference_m2": rho_ref,
                    "caps_active": False,
                    "status": "ANALYTICAL_DBTT_CAPACITY_PROXY_NOT_MPZ_CALIBRATION",
                })
            if count % 1000 == 0 or count == total:
                print(f"evaluated {count}/{total}", flush=True)

    summary = pd.DataFrame(summary_rows).sort_values("dbtt_proxy_score")
    detail = pd.DataFrame(temp_rows)
    summary.to_csv(out / "broad_dbtt_map_all.csv.gz", index=False, compression="gzip")
    detail.to_csv(out / "broad_dbtt_map_temperature_detail.csv.gz", index=False, compression="gzip")
    summary.head(a.top_count).to_csv(out / "broad_dbtt_map_shortlist.csv", index=False)
    canonical_hits = summary[
        summary.candidate_source == "prior_first_passage_reference"
    ]
    canonical_hits.to_csv(out / "canonical_reference_results.csv", index=False)

    report = {
        "n_intrinsic_candidates": int(len(candidates)),
        "n_pt_closures": int(len(pt)),
        "n_evaluated": int(len(summary)),
        "N_reference": N_ref,
        "chi_reference": chi_ref,
        "rho_reference_m2": rho_ref,
        "benchmark_not_cap": True,
        "common_closure_required": False,
        "preselection_top_per_region": False,
        "production_caps_or_saturations": [],
        "best": summary.head(1).to_dict(orient="records"),
    }
    (out / "broad_dbtt_map_summary.json").write_text(json.dumps(report, indent=2))
    print(json.dumps({k: v for k, v in report.items() if k != "best"}, indent=2), flush=True)
    print(f"Outputs: {out}", flush=True)


if __name__ == "__main__":
    main()
