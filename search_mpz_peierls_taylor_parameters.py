#!/usr/bin/env python3
"""Screen emission-derived Peierls--Taylor parameters after the intrinsic atlas.

The input rows come from the v9.2 analytical cleavage/emission atlas. This
stage does not run a crack-growth simulation. It asks whether the emission
surface can generate a physically admissible bulk/process-zone transport law
when Peierls and Taylor barriers are scaled from that surface and Taylor
completion is correlated at high forest density.

Candidates are rejected if the fixed-rate flow stress is unresolved or
decreases materially with increasing forest density. Total forest density is
never capped constitutively.
"""
from __future__ import annotations

import argparse
from dataclasses import replace
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import qmc

from arrhenius_fracture.emission_derived_plasticity import (
    CorrelatedTaylorConfig,
    EmissionDerivedPeierlsTaylorConfig,
    EmissionDerivedPeierlsTaylorModel,
    ExpFloorSurface,
    MechanismScale,
)


def floats(text: str) -> list[float]:
    return [float(x) for x in str(text).replace(",", " ").split() if x]


def strings(text: str) -> list[str]:
    return [x for x in str(text).replace(",", " ").split() if x]


def choose_intrinsic_rows(
    df: pd.DataFrame,
    regions: list[str],
    intrinsic_rate: float,
    top_per_region: int,
    min_virgin_K: float,
) -> pd.DataFrame:
    out = df.copy()
    if "Kdot_MPa_sqrt_m_per_s" in out:
        vals = out["Kdot_MPa_sqrt_m_per_s"].to_numpy(float)
        out = out[np.isclose(vals, intrinsic_rate, rtol=1e-8, atol=1e-12)]
    if "region" in out:
        out = out[out.region.astype(str).isin(regions)]

    kc_cols = [c for c in out.columns if c.startswith("refined_Kc_T")]
    if not kc_cols:
        kc_cols = [c for c in out.columns if c.startswith("Kc_T")]
    if kc_cols:
        out = out.copy()
        out["min_virgin_K_MPa_sqrt_m"] = out[kc_cols].min(
            axis=1, skipna=True
        )
        out = out[out.min_virgin_K_MPa_sqrt_m >= min_virgin_K]

    score_map = {
        "ceramic_intrinsic": "ceramic_score",
        "weakT_intrinsic": "weakT_score",
        "DBTT_precursor": "DBTT_precursor_score",
    }
    chunks = []
    for region in regions:
        g = out[out.region.astype(str) == region].copy()
        score = score_map.get(region)
        if score in g:
            g = g.sort_values(score)
        elif "shortlist_score" in g:
            g = g.sort_values("shortlist_score")
        chunks.append(g.head(top_per_region))
    ans = pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame()
    if ans.empty:
        raise SystemExit("No intrinsic candidates passed the requested filters")
    return ans


def log_scale(u: np.ndarray, lo: float, hi: float) -> np.ndarray:
    return 10.0 ** (
        math.log10(lo) + u * (math.log10(hi) - math.log10(lo))
    )


def sample_transport_parameters(n: int, seed: int) -> pd.DataFrame:
    """Sobol sample of the post-emission kinetic closure."""
    sampler = qmc.Sobol(10, scramble=True, seed=seed)
    if n > 0 and (n & (n - 1)) == 0:
        u = sampler.random_base2(int(math.log2(n)))
    else:
        u = sampler.random(n)
    sampled = pd.DataFrame({
        "pt_peierls_energy_ratio": log_scale(u[:, 0], 0.0025, 0.010),
        "pt_taylor_energy_ratio": log_scale(u[:, 1], 0.010, 0.040),
        "pt_entropy_multiplier": log_scale(u[:, 2], 0.25, 128.0),
        "pt_taylor_corr_rho_c": log_scale(u[:, 3], 1.0e9, 1.0e14),
        "pt_taylor_renewal_time_s": log_scale(u[:, 4], 1.0e-18, 1.0e-8),
        "pt_taylor_m_exponent": 0.5 + 1.5 * u[:, 5],
        "pt_taylor_m_scale": log_scale(u[:, 6], 0.10, 10.0),
        "pt_taylor_m_cap": 6.0 + 42.0 * u[:, 7],
        "pt_mobile_saturation_density_m2": log_scale(
            u[:, 8], 1.0e12, 1.0e16
        ),
        "pt_mobile_fraction": log_scale(u[:, 9], 1.0e-4, 5.0e-2),
    })

    anchors = []
    for entropy_mult in (1.0, 10.0, 50.0):
        for m_cap in (15.0, 22.0, 30.0):
            anchors.append({
                "pt_peierls_energy_ratio": 0.005,
                "pt_taylor_energy_ratio": 0.020,
                "pt_entropy_multiplier": entropy_mult,
                "pt_taylor_corr_rho_c": 1.0e11,
                "pt_taylor_renewal_time_s": 1.0e-10,
                "pt_taylor_m_exponent": 1.0,
                "pt_taylor_m_scale": 1.0,
                "pt_taylor_m_cap": m_cap,
                "pt_mobile_saturation_density_m2": 1.0e14,
                "pt_mobile_fraction": 0.01,
            })
    return pd.concat([pd.DataFrame(anchors), sampled], ignore_index=True)


def surface_from_row(row: pd.Series) -> ExpFloorSurface:
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


def model_from_rows(
    parent: ExpFloorSurface, p: pd.Series
) -> EmissionDerivedPeierlsTaylorModel:
    pe = float(p.pt_peierls_energy_ratio)
    te = float(p.pt_taylor_energy_ratio)
    em = float(p.pt_entropy_multiplier)
    cfg = EmissionDerivedPeierlsTaylorConfig(
        parent=parent,
        peierls=MechanismScale(pe, pe * em, 1.0, 1.0e12),
        taylor=MechanismScale(te, te * em, 1.0, 1.0e11),
        correlated_taylor=CorrelatedTaylorConfig(
            rho_c_m2=float(p.pt_taylor_corr_rho_c),
            renewal_time_s=float(p.pt_taylor_renewal_time_s),
            m_exponent=float(p.pt_taylor_m_exponent),
            m_scale=float(p.pt_taylor_m_scale),
            m_cap=float(p.pt_taylor_m_cap),
        ),
        mobile_fraction_low_density=float(p.pt_mobile_fraction),
        mobile_saturation_density_m2=float(
            p.pt_mobile_saturation_density_m2
        ),
        mobile_density_floor_m2=1.0e6,
        jump_fraction_of_forest_spacing=1.0,
        jump_length_min_m=2.5e-10,
        rate_cap_s=float("inf"),
    )
    return EmissionDerivedPeierlsTaylorModel(cfg)


def evaluate_one(
    model: EmissionDerivedPeierlsTaylorModel,
    rho: np.ndarray,
    temperatures: list[float],
    strain_rates: list[float],
    b: float,
    min_reference_stress_GPa: float,
    max_reference_stress_GPa: float,
    max_stress_GPa: float,
    slope_tol_GPa_decade: float,
    drop_tol_fraction: float,
    zero_stress_threshold_GPa: float,
) -> dict[str, float | bool | str]:
    curves = []
    slopes = []
    drops = []
    resolved = True
    for T in temperatures:
        for edot in strain_rates:
            sig = model.flow_stress(rho, T, edot, b) / 1.0e9
            curves.append(sig)
            if not np.all(np.isfinite(sig)):
                resolved = False
                continue
            dlog = np.diff(np.log10(rho))
            ds = np.diff(sig)
            slopes.append(float(np.min(ds / dlog)))
            drops.append(float(np.max(
                np.maximum(-ds, 0.0) / np.maximum(sig[:-1], 1.0e-6)
            )))
    arr = np.asarray(curves, dtype=float)
    finite = np.isfinite(arr)
    smin = float(np.nanmin(arr)) if np.any(finite) else float("nan")
    smax = float(np.nanmax(arr)) if np.any(finite) else float("nan")
    zero_fraction = (
        float(np.mean(arr[finite] <= zero_stress_threshold_GPa))
        if np.any(finite) else 1.0
    )
    min_slope = min(slopes) if slopes else float("-inf")
    max_drop = max(drops) if drops else float("inf")
    monotonic = (
        min_slope >= -abs(slope_tol_GPa_decade)
        and max_drop <= drop_tol_fraction
    )

    iref = int(np.argmin(np.abs(np.log10(rho) - 14.0)))
    Tref = min(temperatures, key=lambda x: abs(x - 700.0))
    eref = min(strain_rates, key=lambda x: abs(math.log10(x) + 5.0))
    k = (
        temperatures.index(Tref) * len(strain_rates)
        + strain_rates.index(eref)
    )
    sigma_ref = float(arr[k, iref]) if k < len(arr) else float("nan")

    mvals = model.cfg.correlated_taylor.hit_order(rho)
    mcap = float(model.cfg.correlated_taylor.m_cap)
    cap_active_fraction = (
        float(np.mean(mvals >= mcap * (1.0 - 1.0e-9)))
        if np.isfinite(mcap) and mcap > 1.0 else 0.0
    )

    admissible = (
        resolved
        and monotonic
        and np.isfinite(smax)
        and smax <= max_stress_GPa
    )
    strength_window = (
        admissible
        and np.isfinite(sigma_ref)
        and min_reference_stress_GPa
        <= sigma_ref
        <= max_reference_stress_GPa
    )
    if strength_window:
        status = "strict_strength_window"
    elif admissible:
        status = "monotonic_topology_only"
    elif not resolved:
        status = "unresolved_at_sigma_limit"
    elif not monotonic:
        status = "high_density_downturn"
    else:
        status = "stress_above_limit"

    return {
        "resolved": bool(resolved),
        "monotonic": bool(monotonic),
        "accepted": bool(admissible),
        "strict_strength_window": bool(strength_window),
        "pt_screen_status": status,
        "sigma_min_GPa": smin,
        "sigma_max_GPa": smax,
        "sigma_ref_700K_1e14_GPa": sigma_ref,
        "zero_stress_fraction": zero_fraction,
        "min_slope_GPa_per_decade": min_slope,
        "max_step_drop_fraction": max_drop,
        "m_eff_at_rho_min": float(mvals[0]),
        "m_eff_at_rho_max": float(mvals[-1]),
        "m_cap_active_fraction": cap_active_fraction,
    }


def independent_peak_diagnostic(
    model: EmissionDerivedPeierlsTaylorModel,
    rho: np.ndarray,
    T: float,
    edot: float,
    b: float,
) -> dict[str, float]:
    independent = EmissionDerivedPeierlsTaylorModel(
        replace(
            model.cfg,
            correlated_taylor=replace(
                model.cfg.correlated_taylor,
                rho_c_m2=1.0e300,
                m_exponent=1.0,
                m_scale=0.0,
                m_cap=1.0,
            ),
        )
    )
    sig = independent.flow_stress(rho, T, edot, b) / 1.0e9
    if not np.any(np.isfinite(sig)):
        return {
            "independent_peak_rho_m2": float("nan"),
            "independent_peak_stress_GPa": float("nan"),
            "independent_min_slope_GPa_per_decade": float("nan"),
        }
    imax = int(np.nanargmax(sig))
    slope = np.diff(sig) / np.diff(np.log10(rho))
    return {
        "independent_peak_rho_m2": float(rho[imax]),
        "independent_peak_stress_GPa": float(sig[imax]),
        "independent_min_slope_GPa_per_decade": float(np.nanmin(slope)),
    }


def make_plots(
    short: pd.DataFrame,
    rho: np.ndarray,
    temperatures: list[float],
    edot: float,
    out: Path,
) -> None:
    if short.empty:
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    for _, row in short.iterrows():
        parent = surface_from_row(row)
        model = model_from_rows(parent, row)
        fig, ax = plt.subplots(figsize=(7.5, 5.4))
        for T in temperatures:
            sig = model.flow_stress(rho, T, edot, 2.74e-10) / 1.0e9
            ax.plot(rho, sig, label=f"{T:g} K")
        ax.set_xscale("log")
        ax.set_xlabel(r"Forest density $\rho_f$ [m$^{-2}$]")
        ax.set_ylabel("Flow stress [GPa]")
        ax.set_title(
            f"{row.get('region','')} {row.get('candidate_id','')} "
            f"PT rank {int(row.pt_rank)}"
        )
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=8)
        fig.tight_layout()
        tag = str(row.get("candidate_id", "candidate")).replace("/", "_")
        fig.savefig(
            out / f"pt_strength_{tag}_r{int(row.pt_rank):02d}.png",
            dpi=180,
        )
        plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--atlas-shortlist", type=Path,
        default=Path(
            "runs/mpz_v9_2_analytic_first_passage_atlas/"
            "analytic_first_passage_atlas_shortlist_refined.csv"
        ),
    )
    ap.add_argument(
        "--regions",
        default="ceramic_intrinsic weakT_intrinsic DBTT_precursor",
    )
    ap.add_argument("--intrinsic-rate", type=float, default=0.005)
    ap.add_argument("--intrinsic-top-per-region", type=int, default=5)
    ap.add_argument("--min-virgin-K", type=float, default=3.0)
    ap.add_argument("--transport-samples", type=int, default=256)
    ap.add_argument("--seed", type=int, default=93017)
    ap.add_argument("--temperatures", default="300 700 900 1200")
    ap.add_argument("--strain-rates", default="1e-5 1e-3")
    ap.add_argument("--rho-min", type=float, default=5.0e12)
    ap.add_argument("--rho-max", type=float, default=1.0e18)
    ap.add_argument("--rho-points", type=int, default=65)
    ap.add_argument("--min-reference-stress-GPa", type=float, default=0.05)
    ap.add_argument("--max-reference-stress-GPa", type=float, default=20.0)
    ap.add_argument("--max-stress-GPa", type=float, default=80.0)
    ap.add_argument("--zero-stress-threshold-GPa", type=float, default=1.0e-4)
    ap.add_argument("--slope-tol-GPa-decade", type=float, default=0.02)
    ap.add_argument("--drop-tol-fraction", type=float, default=0.02)
    ap.add_argument("--top-per-intrinsic", type=int, default=3)
    ap.add_argument(
        "--out", type=Path,
        default=Path("runs/mpz_v9_3_peierls_taylor_search"),
    )
    a = ap.parse_args()

    source = a.atlas_shortlist
    if not source.exists() and source.name.endswith("_refined.csv"):
        source = source.with_name(source.name.replace("_refined", ""))
    if not source.exists():
        raise SystemExit(f"Atlas shortlist not found: {a.atlas_shortlist}")
    df = pd.read_csv(source)
    intrinsic = choose_intrinsic_rows(
        df,
        strings(a.regions),
        a.intrinsic_rate,
        a.intrinsic_top_per_region,
        a.min_virgin_K,
    )
    samples = sample_transport_parameters(a.transport_samples, a.seed)
    rho = np.logspace(
        math.log10(a.rho_min), math.log10(a.rho_max), a.rho_points
    )
    temperatures = floats(a.temperatures)
    strain_rates = floats(a.strain_rates)
    out = a.out.resolve()
    out.mkdir(parents=True, exist_ok=True)

    rows = []
    total = len(intrinsic) * len(samples)
    count = 0
    for _, base in intrinsic.iterrows():
        parent = surface_from_row(base)
        for _, p in samples.iterrows():
            count += 1
            model = model_from_rows(parent, p)
            metrics = evaluate_one(
                model,
                rho,
                temperatures,
                strain_rates,
                2.74e-10,
                a.min_reference_stress_GPa,
                a.max_reference_stress_GPa,
                a.max_stress_GPa,
                a.slope_tol_GPa_decade,
                a.drop_tol_fraction,
                a.zero_stress_threshold_GPa,
            )
            prior = (
                abs(math.log(float(p.pt_peierls_energy_ratio) / 0.005))
                + abs(math.log(float(p.pt_taylor_energy_ratio) / 0.02))
                + 0.10 * abs(math.log(float(p.pt_entropy_multiplier)))
                + 0.10 * abs(math.log(float(p.pt_mobile_fraction) / 0.01))
            )
            sigma_ref = max(
                float(metrics["sigma_ref_700K_1e14_GPa"]), 1.0e-9
            )
            stress_penalty = abs(math.log(sigma_ref / 2.0))
            topology_penalty = (
                max(-float(metrics["min_slope_GPa_per_decade"]), 0.0)
                * 20.0
                + 2.0 * float(metrics["max_step_drop_fraction"])
                + 0.25 * float(metrics["zero_stress_fraction"])
            )
            rejection = 0.0 if metrics["accepted"] else 1000.0
            strict_bonus = (
                -1.0 if metrics["strict_strength_window"] else 0.0
            )
            score = (
                rejection + prior + stress_penalty
                + topology_penalty + strict_bonus
            )
            rec = {
                **base.to_dict(),
                **p.to_dict(),
                **metrics,
                "pt_score": score,
            }
            rec["pt_peierls_entropy_ratio"] = (
                float(p.pt_peierls_energy_ratio)
                * float(p.pt_entropy_multiplier)
            )
            rec["pt_taylor_entropy_ratio"] = (
                float(p.pt_taylor_energy_ratio)
                * float(p.pt_entropy_multiplier)
            )
            rows.append(rec)
            if count % 250 == 0 or count == total:
                print(f"evaluated {count}/{total}", flush=True)

    all_df = pd.DataFrame(rows)
    accepted = all_df[all_df.accepted.astype(bool)].copy().sort_values(
        ["strict_strength_window", "pt_score"],
        ascending=[False, True],
    )
    ranked = []
    for _, g in accepted.groupby("candidate_id", sort=False):
        h = g.sort_values("pt_score").head(a.top_per_intrinsic).copy()
        h["pt_rank"] = np.arange(1, len(h) + 1)
        ranked.append(h)
    short = (
        pd.concat(ranked, ignore_index=True)
        if ranked else pd.DataFrame(columns=accepted.columns)
    )

    if not short.empty:
        diags = []
        for _, row in short.iterrows():
            model = model_from_rows(surface_from_row(row), row)
            diags.append(
                independent_peak_diagnostic(
                    model, rho, 700.0, 1.0e-5, 2.74e-10
                )
            )
        short = pd.concat(
            [short.reset_index(drop=True), pd.DataFrame(diags)], axis=1
        )

    all_df.to_csv(
        out / "peierls_taylor_search_all.csv.gz",
        index=False,
        compression="gzip",
    )
    accepted.to_csv(
        out / "peierls_taylor_search_accepted.csv", index=False
    )
    short.to_csv(
        out / "peierls_taylor_search_shortlist.csv", index=False
    )
    if not short.empty:
        material = short.copy()
        material["status"] = (
            "PT_MONOTONIC_SCREEN_PASSED_NOT_TRANSIENTLY_VALIDATED"
        )
        material.to_csv(
            out / "mpz_v9_3_material_shortlist.csv", index=False
        )
        make_plots(
            short.groupby("region", sort=False).head(2),
            rho,
            temperatures,
            1.0e-5,
            out,
        )

    counts = pd.DataFrame({
        "n_intrinsic_candidates": [len(intrinsic)],
        "n_transport_samples_per_candidate": [len(samples)],
        "n_total": [len(all_df)],
        "n_accepted": [len(accepted)],
        "acceptance_fraction": [len(accepted) / max(len(all_df), 1)],
    })
    counts.to_csv(
        out / "peierls_taylor_search_summary.csv", index=False
    )
    config = vars(a).copy()
    for k, v in list(config.items()):
        if isinstance(v, Path):
            config[k] = str(v)
    config.update({
        "source_resolved": str(source),
        "rho_grid": rho.tolist(),
        "temperatures_resolved": temperatures,
        "strain_rates_resolved": strain_rates,
        "production_ratios": {
            "peierls_over_emission": 0.005,
            "taylor_over_emission": 0.02,
        },
        "total_density_cap_used": False,
        "taylor_hit_order_cap_interpretation": (
            "finite obstacle count in one correlation domain; every candidate "
            "is rejected if a high-density flow-stress downturn reappears"
        ),
    })
    (out / "peierls_taylor_search_config.json").write_text(
        json.dumps(config, indent=2)
    )
    print(counts.to_string(index=False), flush=True)
    print(f"Outputs: {out}", flush=True)


if __name__ == "__main__":
    main()
