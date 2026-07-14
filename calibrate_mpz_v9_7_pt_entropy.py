#!/usr/bin/env python3
"""Calibrate PT magnitude first, then independent activation entropy.

The reference-temperature strength cannot be repaired by changing activation
entropy because the entropy term vanishes at Tref.  This script therefore uses
two transparent stages:

1. Select Peierls/Taylor energy-ratio pairs that produce a finite GPa-scale
   flow stress at the chosen magnitude temperature.
2. Hold each selected energy pair fixed and sweep independent Peierls and
   Taylor activation entropies to map thermal softening, near-athermal response,
   and thermal hardening.

No common closure across material classes is imposed.  No density, stress,
hit-order, mobile-density, jump-length, or rate cap is introduced.
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
)
from arrhenius_fracture.emission_derived_plasticity_v97 import (
    EmissionDerivedPeierlsTaylorModel,
    IndependentEntropyMechanismScale,
    KB_EV_PER_K,
)


def floats(text: str) -> list[float]:
    return [float(x) for x in str(text).replace(",", " ").split() if x]


def parent_surface(row: pd.Series) -> ExpFloorSurface:
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


def make_model(
    row: pd.Series,
    peierls_energy_ratio: float,
    taylor_energy_ratio: float,
    peierls_entropy_kB: float,
    taylor_entropy_kB: float,
    args: argparse.Namespace,
) -> EmissionDerivedPeierlsTaylorModel:
    cfg = EmissionDerivedPeierlsTaylorConfig(
        parent=parent_surface(row),
        peierls=IndependentEntropyMechanismScale(
            peierls_energy_ratio,
            peierls_entropy_kB,
            args.peierls_stress_ratio,
            args.nu0_peierls,
        ),
        taylor=IndependentEntropyMechanismScale(
            taylor_energy_ratio,
            taylor_entropy_kB,
            args.taylor_stress_ratio,
            args.nu0_taylor,
        ),
        correlated_taylor=CorrelatedTaylorConfig(
            rho_c_m2=args.correlation_rho_c_m2,
            renewal_time_s=1.0,
            m_exponent=args.correlation_exponent,
            m_scale=args.correlation_scale,
            m_cap=float("inf"),
        ),
        mobile_fraction_low_density=args.mobile_fraction,
        mobile_saturation_density_m2=float("inf"),
        mobile_density_floor_m2=0.0,
        jump_fraction_of_forest_spacing=args.jump_fraction,
        jump_length_min_m=0.0,
        taylor_phi_max=float("inf"),
        rate_cap_s=float("inf"),
    )
    return EmissionDerivedPeierlsTaylorModel(cfg)


def entropy_samples(n: int, seed: int, lo: float, hi: float) -> pd.DataFrame:
    sampler = qmc.Sobol(2, scramble=True, seed=seed)
    if n > 0 and (n & (n - 1)) == 0:
        u = sampler.random_base2(int(math.log2(n)))
    else:
        u = sampler.random(n)
    sampled = pd.DataFrame(
        {
            "peierls_activation_entropy_kB": lo + (hi - lo) * u[:, 0],
            "taylor_activation_entropy_kB": lo + (hi - lo) * u[:, 1],
        }
    )
    anchors = pd.DataFrame(
        [
            {
                "peierls_activation_entropy_kB": p,
                "taylor_activation_entropy_kB": t,
            }
            for p, t in (
                (0.0, 0.0),
                (-10.0, -10.0),
                (-20.0, -20.0),
                (-30.0, -30.0),
                (-20.0, -10.0),
                (-10.0, -20.0),
                (10.0, 10.0),
            )
        ]
    )
    return pd.concat([anchors, sampled], ignore_index=True)


def solve_stresses(
    model: EmissionDerivedPeierlsTaylorModel,
    T_K: float,
    rates: list[float],
    rho: np.ndarray,
    args: argparse.Namespace,
) -> dict[float, np.ndarray]:
    return {
        rate: model.flow_stress(
            rho,
            T_K,
            rate,
            args.b_m,
            sigma_max_Pa=args.sigma_max_GPa * 1.0e9,
            iterations=args.iterations,
        )
        for rate in rates
    }


def classify_thermal_ratio(ratio: float) -> str:
    if not np.isfinite(ratio):
        return "unresolved"
    if ratio < 0.25:
        return "strong_softening"
    if ratio < 0.70:
        return "moderate_softening"
    if ratio <= 1.30:
        return "near_athermal"
    return "thermal_hardening"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--materials",
        type=Path,
        default=Path("mpz_v9_6_canonical_first_passage_references.csv"),
    )
    ap.add_argument("--classes", default="ceramic weakT DBTT")
    ap.add_argument("--temperatures", default="300 700 900 1200")
    ap.add_argument("--strain-rates", default="1e-5 1e-3")
    ap.add_argument("--rho-values-m2", default="5e12 1e14 1e16")
    ap.add_argument("--reference-rho-m2", type=float, default=1.0e14)
    ap.add_argument("--magnitude-temperature-K", type=float, default=700.0)
    ap.add_argument("--target-reference-stress-GPa", type=float, default=2.0)
    ap.add_argument("--reference-stress-min-GPa", type=float, default=0.1)
    ap.add_argument("--reference-stress-max-GPa", type=float, default=10.0)
    ap.add_argument("--global-stress-max-GPa", type=float, default=40.0)
    ap.add_argument("--energy-ratio-min", type=float, default=0.02)
    ap.add_argument("--energy-ratio-max", type=float, default=2.0)
    ap.add_argument("--energy-ratio-points", type=int, default=17)
    ap.add_argument("--magnitude-top-per-class", type=int, default=16)
    ap.add_argument("--entropy-min-kB", type=float, default=-60.0)
    ap.add_argument("--entropy-max-kB", type=float, default=30.0)
    ap.add_argument("--entropy-samples", type=int, default=256)
    ap.add_argument("--seed", type=int, default=97017)
    ap.add_argument("--peierls-stress-ratio", type=float, default=1.0)
    ap.add_argument("--taylor-stress-ratio", type=float, default=1.0)
    ap.add_argument("--nu0-peierls", type=float, default=1.0e12)
    ap.add_argument("--nu0-taylor", type=float, default=1.0e11)
    ap.add_argument("--correlation-rho-c-m2", type=float, default=1.0e14)
    ap.add_argument("--correlation-scale", type=float, default=1.0)
    ap.add_argument("--correlation-exponent", type=float, default=1.0)
    ap.add_argument("--mobile-fraction", type=float, default=0.01)
    ap.add_argument("--jump-fraction", type=float, default=1.0)
    ap.add_argument("--b-m", type=float, default=2.74e-10)
    ap.add_argument("--sigma-max-GPa", type=float, default=200.0)
    ap.add_argument("--iterations", type=int, default=64)
    ap.add_argument("--shortlist-count", type=int, default=200)
    ap.add_argument(
        "--out",
        type=Path,
        default=Path("runs/mpz_v9_7_pt_entropy_calibration_v1"),
    )
    a = ap.parse_args()

    materials = pd.read_csv(a.materials)
    class_col = "target_class" if "target_class" in materials else "region"
    requested = str(a.classes).replace(",", " ").split()
    materials = materials[materials[class_col].astype(str).isin(requested)].copy()
    if materials.empty:
        raise SystemExit(f"no requested classes found in {a.materials}")

    temperatures = floats(a.temperatures)
    strain_rates = sorted(floats(a.strain_rates))
    rho_values = np.asarray(sorted(floats(a.rho_values_m2)), dtype=float)
    ref_index = int(np.argmin(np.abs(np.log(rho_values / a.reference_rho_m2))))
    ref_rho = float(rho_values[ref_index])
    ratios = np.logspace(
        math.log10(a.energy_ratio_min),
        math.log10(a.energy_ratio_max),
        a.energy_ratio_points,
    )
    out = a.out.resolve()
    out.mkdir(parents=True, exist_ok=True)

    # Stage A: strength magnitude at fixed zero activation entropy.
    magnitude_rows: list[dict[str, float | str | bool]] = []
    for _, material in materials.iterrows():
        klass = str(material[class_col])
        for p_ratio in ratios:
            for t_ratio in ratios:
                model = make_model(material, p_ratio, t_ratio, 0.0, 0.0, a)
                solved = solve_stresses(
                    model,
                    a.magnitude_temperature_K,
                    strain_rates,
                    rho_values,
                    a,
                )
                low = solved[strain_rates[0]] / 1.0e9
                high = solved[strain_rates[-1]] / 1.0e9
                ref_low = float(low[ref_index])
                ref_high = float(high[ref_index])
                finite = bool(np.all(np.isfinite(low)) and np.all(np.isfinite(high)))
                rate_sensitive = bool(
                    finite and np.all(high >= low * (1.0 - 1.0e-8))
                )
                in_window = bool(
                    finite
                    and a.reference_stress_min_GPa <= ref_low
                    <= a.reference_stress_max_GPa
                )
                score = (
                    abs(
                        math.log10(
                            max(ref_low, 1.0e-30)
                            / a.target_reference_stress_GPa
                        )
                    )
                    if finite and ref_low > 0.0
                    else 1.0e6
                )
                if not rate_sensitive:
                    score += 100.0
                if not in_window:
                    score += 10.0
                magnitude_rows.append(
                    {
                        "target_class": klass,
                        "candidate_id": str(material.get("candidate_id", klass)),
                        "peierls_energy_ratio": p_ratio,
                        "taylor_energy_ratio": t_ratio,
                        "magnitude_temperature_K": a.magnitude_temperature_K,
                        "reference_rho_m2": ref_rho,
                        "reference_stress_low_rate_GPa": ref_low,
                        "reference_stress_high_rate_GPa": ref_high,
                        "stress_min_GPa": (
                            float(np.nanmin(low))
                            if np.any(np.isfinite(low))
                            else np.nan
                        ),
                        "stress_max_GPa": (
                            float(np.nanmax(high))
                            if np.any(np.isfinite(high))
                            else np.nan
                        ),
                        "finite": finite,
                        "rate_sensitive": rate_sensitive,
                        "in_reference_window": in_window,
                        "magnitude_score": score,
                    }
                )
    magnitude = pd.DataFrame(magnitude_rows)
    magnitude.to_csv(out / "pt_magnitude_grid.csv", index=False)
    selected_magnitude = (
        magnitude[
            magnitude.finite
            & magnitude.rate_sensitive
            & magnitude.in_reference_window
        ]
        .sort_values(["target_class", "magnitude_score"])
        .groupby("target_class", as_index=False)
        .head(a.magnitude_top_per_class)
        .reset_index(drop=True)
    )
    selected_magnitude.to_csv(out / "pt_magnitude_selected.csv", index=False)
    missing_classes = sorted(
        set(materials[class_col].astype(str))
        - set(selected_magnitude.target_class.astype(str))
    )
    if missing_classes:
        raise SystemExit(
            "no finite GPa-scale magnitude pair found for class(es): "
            + ", ".join(missing_classes)
            + "; inspect pt_magnitude_grid.csv or widen the energy-ratio range"
        )

    # Stage B: independent entropy sweep around accepted magnitude families.
    entropy = entropy_samples(
        a.entropy_samples, a.seed, a.entropy_min_kB, a.entropy_max_kB
    )
    thermal_rows: list[dict[str, float | str | bool]] = []
    material_by_class = {
        str(row[class_col]): row for _, row in materials.iterrows()
    }
    total = len(selected_magnitude) * len(entropy)
    count = 0
    for _, mag in selected_magnitude.iterrows():
        material = material_by_class[str(mag.target_class)]
        for _, ent in entropy.iterrows():
            count += 1
            pS = float(ent.peierls_activation_entropy_kB)
            tS = float(ent.taylor_activation_entropy_kB)
            model = make_model(
                material,
                float(mag.peierls_energy_ratio),
                float(mag.taylor_energy_ratio),
                pS,
                tS,
                a,
            )
            values: dict[tuple[float, float], np.ndarray] = {}
            finite = True
            for T in temperatures:
                solved = solve_stresses(model, T, strain_rates, rho_values, a)
                for rate, stress in solved.items():
                    values[(T, rate)] = stress / 1.0e9
                    finite = finite and bool(np.all(np.isfinite(stress)))
            low_rate = strain_rates[0]
            ref_by_T = {
                T: float(values[(T, low_rate)][ref_index]) for T in temperatures
            }
            reference_min = float(np.nanmin(list(ref_by_T.values())))
            reference_max = float(np.nanmax(list(ref_by_T.values())))
            reference_window_all_T = bool(
                finite
                and reference_min >= a.reference_stress_min_GPa
                and reference_max <= a.reference_stress_max_GPa
            )
            stress_all = np.concatenate(list(values.values()))
            raw_barriers = [
                model.raw_zero_stress_barrier_eV(mech, T)
                for mech in ("peierls", "taylor")
                for T in temperatures
            ]
            min_raw_barrier = float(np.min(raw_barriers))
            rate_sensitive = (
                all(
                    np.all(
                        values[(T, strain_rates[-1])]
                        >= values[(T, strain_rates[0])] * (1.0 - 1.0e-8)
                    )
                    for T in temperatures
                )
                if finite
                else False
            )
            stress_min = (
                float(np.nanmin(stress_all))
                if np.any(np.isfinite(stress_all))
                else np.nan
            )
            stress_max = (
                float(np.nanmax(stress_all))
                if np.any(np.isfinite(stress_all))
                else np.nan
            )
            plausible = bool(
                finite
                and rate_sensitive
                and min_raw_barrier > 0.0
                and reference_window_all_T
                and stress_max <= a.global_stress_max_GPa
            )
            T_lo, T_hi = min(temperatures), max(temperatures)
            ratio = ref_by_T[T_hi] / max(ref_by_T[T_lo], 1.0e-30)
            ref_mid = ref_by_T[
                min(
                    temperatures,
                    key=lambda x: abs(x - a.magnitude_temperature_K),
                )
            ]
            score = abs(
                math.log10(
                    max(ref_mid, 1.0e-30) / a.target_reference_stress_GPa
                )
            )
            if not plausible:
                score += 100.0
            thermal_rows.append(
                {
                    "target_class": str(mag.target_class),
                    "candidate_id": str(mag.candidate_id),
                    "peierls_energy_ratio": float(mag.peierls_energy_ratio),
                    "taylor_energy_ratio": float(mag.taylor_energy_ratio),
                    "peierls_activation_entropy_kB": pS,
                    "taylor_activation_entropy_kB": tS,
                    "peierls_gT_eV_per_K": -pS * KB_EV_PER_K,
                    "taylor_gT_eV_per_K": -tS * KB_EV_PER_K,
                    "reference_rho_m2": ref_rho,
                    **{
                        f"stress_ref_{int(T)}K_GPa": ref_by_T[T]
                        for T in temperatures
                    },
                    "reference_stress_min_over_T_GPa": reference_min,
                    "reference_stress_max_over_T_GPa": reference_max,
                    "reference_window_all_T": reference_window_all_T,
                    "stress_min_GPa": stress_min,
                    "stress_max_GPa": stress_max,
                    "min_raw_zero_stress_barrier_eV": min_raw_barrier,
                    "rate_sensitive": rate_sensitive,
                    "finite": finite,
                    "plausible": plausible,
                    "thermal_ratio_high_over_low": ratio,
                    "thermal_regime": classify_thermal_ratio(ratio),
                    "thermal_score": score,
                }
            )
            if count % 250 == 0 or count == total:
                print(f"evaluated {count}/{total}", flush=True)

    thermal = pd.DataFrame(thermal_rows)
    thermal.to_csv(
        out / "pt_entropy_map_all.csv.gz",
        index=False,
        compression="gzip",
    )
    plausible = thermal[thermal.plausible.astype(bool)].copy()
    plausible.to_csv(out / "pt_entropy_map_plausible.csv", index=False)
    shortlist = (
        plausible.sort_values(
            ["target_class", "thermal_regime", "thermal_score"]
        )
        .groupby(["target_class", "thermal_regime"], as_index=False)
        .head(max(1, a.shortlist_count // max(1, len(requested) * 4)))
        .sort_values(["thermal_score", "target_class"])
        .head(a.shortlist_count)
        .reset_index(drop=True)
    )
    shortlist.to_csv(out / "pt_entropy_shortlist.csv", index=False)

    summary = {
        "n_magnitude_grid": int(len(magnitude)),
        "n_magnitude_selected": int(len(selected_magnitude)),
        "n_thermal_evaluations": int(len(thermal)),
        "n_plausible": int(len(plausible)),
        "plausible_by_class": plausible.target_class.value_counts().to_dict(),
        "plausible_by_regime": plausible.thermal_regime.value_counts().to_dict(),
        "reference_temperature_note": (
            "Entropy is independent of emission. Energy ratios set the reference "
            "strength; activation entropy sets the temperature slope."
        ),
        "output": str(out),
    }
    (out / "pt_entropy_calibration_summary.json").write_text(
        json.dumps(summary, indent=2)
    )
    config = vars(a).copy()
    for key, value in list(config.items()):
        if isinstance(value, Path):
            config[key] = str(value)
    config["constitutive_caps_or_saturations"] = []
    config["production_activation_status"] = (
        "CALIBRATION_ONLY_NOT_YET_ACTIVATED_IN_MPZ_OR_BULK_FEM"
    )
    (out / "pt_entropy_calibration_config.json").write_text(
        json.dumps(config, indent=2)
    )
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
