#!/usr/bin/env python3
"""Audit the uncapped v9.6 Arrhenius Peierls--Taylor closure.

The audit plots flow stress together with every density-dependent ingredient so
that shoulders cannot be hidden inside the inversion.  It accepts either the
canonical first-passage references or any joined analytical-atlas table.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

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


def make_model(row: pd.Series, a: argparse.Namespace) -> EmissionDerivedPeierlsTaylorModel:
    cfg = EmissionDerivedPeierlsTaylorConfig(
        parent=surface(row),
        peierls=MechanismScale(
            a.peierls_energy_ratio,
            a.peierls_entropy_ratio,
            a.peierls_stress_ratio,
            a.nu0_peierls,
        ),
        taylor=MechanismScale(
            a.taylor_energy_ratio,
            a.taylor_entropy_ratio,
            a.taylor_stress_ratio,
            a.nu0_taylor,
        ),
        correlated_taylor=CorrelatedTaylorConfig(
            rho_c_m2=a.correlation_rho_c_m2,
            renewal_time_s=1.0,
            m_exponent=a.correlation_exponent,
            m_scale=a.correlation_scale,
            m_cap=float("inf"),
        ),
        mobile_fraction_low_density=a.mobile_fraction,
        mobile_saturation_density_m2=float("inf"),
        mobile_density_floor_m2=0.0,
        jump_fraction_of_forest_spacing=a.jump_fraction,
        jump_length_min_m=0.0,
        taylor_phi_max=float("inf"),
        rate_cap_s=float("inf"),
    )
    return EmissionDerivedPeierlsTaylorModel(cfg)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--materials", type=Path,
        default=Path("mpz_v9_6_canonical_first_passage_references.csv"),
    )
    ap.add_argument("--classes", default="ceramic weakT DBTT")
    ap.add_argument("--temperatures", default="300 700 900 1200")
    ap.add_argument("--strain-rates", default="1e-5 1e-3")
    ap.add_argument("--rho-min", type=float, default=1.0e10)
    ap.add_argument("--rho-max", type=float, default=1.0e18)
    ap.add_argument("--rho-points", type=int, default=129)
    ap.add_argument("--b-m", type=float, default=2.74e-10)
    ap.add_argument("--peierls-energy-ratio", type=float, default=0.005)
    ap.add_argument("--peierls-entropy-ratio", type=float, default=0.005)
    ap.add_argument("--peierls-stress-ratio", type=float, default=1.0)
    ap.add_argument("--taylor-energy-ratio", type=float, default=0.02)
    ap.add_argument("--taylor-entropy-ratio", type=float, default=0.02)
    ap.add_argument("--taylor-stress-ratio", type=float, default=1.0)
    ap.add_argument("--nu0-peierls", type=float, default=1.0e12)
    ap.add_argument("--nu0-taylor", type=float, default=1.0e11)
    ap.add_argument("--correlation-rho-c-m2", type=float, default=1.0e14)
    ap.add_argument("--correlation-scale", type=float, default=1.0)
    ap.add_argument("--correlation-exponent", type=float, default=1.0)
    ap.add_argument("--mobile-fraction", type=float, default=0.01)
    ap.add_argument("--jump-fraction", type=float, default=1.0)
    ap.add_argument("--sigma-max-GPa", type=float, default=200.0)
    ap.add_argument(
        "--out", type=Path,
        default=Path("runs/mpz_v9_6_uncapped_pt_audit"),
    )
    a = ap.parse_args()

    df = pd.read_csv(a.materials)
    class_names = str(a.classes).replace(",", " ").split()
    class_col = "target_class" if "target_class" in df else "region"
    selected = df[df[class_col].astype(str).isin(class_names)].copy()
    if selected.empty:
        raise SystemExit(f"no requested classes in {a.materials}")

    rho = np.logspace(math.log10(a.rho_min), math.log10(a.rho_max), a.rho_points)
    temps = floats(a.temperatures)
    rates = floats(a.strain_rates)
    out = a.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, float | str]] = []

    for _, material in selected.iterrows():
        model = make_model(material, a)
        for T in temps:
            for target_rate in rates:
                stress = model.flow_stress(
                    rho, T, target_rate, a.b_m,
                    sigma_max_Pa=a.sigma_max_GPa * 1.0e9,
                )
                solved = model.rates(stress, rho, T, a.b_m)
                for i, density in enumerate(rho):
                    rows.append({
                        "target_class": str(material[class_col]),
                        "candidate_id": str(material.get("candidate_id", material[class_col])),
                        "T_K": T,
                        "target_rate_s": target_rate,
                        "rho_forest_m2": density,
                        "flow_stress_GPa": float(stress[i] / 1.0e9),
                        "rho_mobile_m2": float(solved["rho_mobile_m2"][i]),
                        "forest_spacing_m": float(solved["forest_spacing_m"][i]),
                        "jump_length_m": float(solved["jump_length_m"][i]),
                        "taylor_amplification": float(solved["taylor_amplification"][i]),
                        "taylor_m_eff": float(solved["taylor_m_eff"][i]),
                        "peierls_rate_s": float(solved["peierls_rate_s"][i]),
                        "taylor_single_hit_rate_s": float(solved["taylor_single_hit_rate_s"][i]),
                        "taylor_completion_rate_s": float(solved["taylor_completion_rate_s"][i]),
                        "series_rate_s": float(solved["series_rate_s"][i]),
                        "equivalent_plastic_rate_s": float(solved["equivalent_plastic_rate_s"][i]),
                        "constitutive_caps_active": False,
                    })

    audit = pd.DataFrame(rows)
    audit.to_csv(out / "uncapped_pt_audit.csv", index=False)

    summaries = []
    for keys, g in audit.groupby(["target_class", "T_K", "target_rate_s"]):
        g = g.sort_values("rho_forest_m2")
        slope = np.diff(g.flow_stress_GPa) / np.diff(np.log10(g.rho_forest_m2))
        summaries.append({
            "target_class": keys[0],
            "T_K": keys[1],
            "target_rate_s": keys[2],
            "stress_min_GPa": float(g.flow_stress_GPa.min()),
            "stress_max_GPa": float(g.flow_stress_GPa.max()),
            "min_slope_GPa_per_decade": float(np.nanmin(slope)),
            "max_slope_GPa_per_decade": float(np.nanmax(slope)),
            "m_start": float(g.taylor_m_eff.iloc[0]),
            "m_end": float(g.taylor_m_eff.iloc[-1]),
            "phi_start": float(g.taylor_amplification.iloc[0]),
            "phi_end": float(g.taylor_amplification.iloc[-1]),
            "caps_active": False,
        })
    pd.DataFrame(summaries).to_csv(out / "uncapped_pt_summary.csv", index=False)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        for klass, gclass in audit.groupby("target_class"):
            for target_rate, grate in gclass.groupby("target_rate_s"):
                fig, axes = plt.subplots(3, 2, figsize=(11, 11))
                for T, g in grate.groupby("T_K"):
                    g = g.sort_values("rho_forest_m2")
                    label = f"{T:g} K"
                    axes[0, 0].plot(g.rho_forest_m2, g.flow_stress_GPa, label=label)
                    axes[0, 1].plot(g.rho_forest_m2, g.taylor_m_eff, label=label)
                    axes[1, 0].plot(g.rho_forest_m2, g.taylor_amplification, label=label)
                    axes[1, 1].plot(g.rho_forest_m2, g.rho_mobile_m2, label=label)
                    axes[2, 0].plot(g.rho_forest_m2, g.jump_length_m, label=label)
                    axes[2, 1].plot(g.rho_forest_m2, g.series_rate_s, label=label)
                labels = [
                    "Flow stress [GPa]", "Taylor hit order", "Taylor amplification",
                    "Mobile density [m^-2]", "Jump length [m]", "PT series rate [s^-1]",
                ]
                for ax, ylabel in zip(axes.flat, labels):
                    ax.set_xscale("log")
                    ax.set_yscale("log")
                    ax.set_xlabel("Forest density [m^-2]")
                    ax.set_ylabel(ylabel)
                    ax.grid(True, alpha=0.25)
                axes[0, 0].legend(fontsize=8)
                fig.suptitle(f"{klass}: uncapped PT audit, target rate={target_rate:g} s^-1")
                fig.tight_layout()
                fig.savefig(out / f"uncapped_pt_{klass}_{target_rate:g}.png", dpi=180)
                plt.close(fig)
    except Exception as exc:
        print(f"plotting skipped: {type(exc).__name__}: {exc}", flush=True)

    config = vars(a).copy()
    for key, value in list(config.items()):
        if isinstance(value, Path):
            config[key] = str(value)
    config["production_caps_or_saturations"] = []
    config["taylor_completion"] = "mean gamma waiting time lambda_T/m"
    (out / "uncapped_pt_config.json").write_text(json.dumps(config, indent=2))
    print(f"Wrote {out}", flush=True)


if __name__ == "__main__":
    main()
