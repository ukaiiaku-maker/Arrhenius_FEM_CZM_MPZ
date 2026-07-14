#!/usr/bin/env python3
"""Analytical virgin-tip first-passage atlas for the moving-PZ model.

The atlas evaluates the exact EXP-floor cleavage and emission surfaces under a
monotonic K ramp without advancing the MPZ state. It predicts the cleavage
first-passage condition B_c=1 and the emission exposure accumulated before that
condition. The result identifies intrinsic ceramic-like, weak-T, and
DBTT-precursor regions before transient process-zone simulations are run.

This is intentionally not a steady-state or R-curve solver. A DBTT precursor
means that emission exposure changes strongly from low to high temperature; a
subsequent moving-process-zone calculation must establish whether the emitted
state is transported, trapped, retained, and mechanically effective.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.special import gammainc
from scipy.stats import qmc

EV_TO_J = 1.602176634e-19
KB = 1.380649e-23

PARAM_NAMES = (
    "cleave_G00_eV",
    "cleave_gT_eV_per_K",
    "cleave_sigc0_GPa",
    "emit_G00_eV",
    "emit_gT_eV_per_K",
    "emit_sigc0_GPa",
)

DEFAULT_BOUNDS = {
    "cleave_G00_eV": (0.6, 4.5),
    "cleave_gT_eV_per_K": (-0.0015, 0.0025),
    "cleave_sigc0_GPa": (1.5, 8.0),
    "emit_G00_eV": (0.6, 3.5),
    "emit_gT_eV_per_K": (-0.0015, 0.0030),
    "emit_sigc0_GPa": (0.5, 8.0),
}


def parse_float_list(text: str) -> list[float]:
    return [float(x) for x in str(text).replace(",", " ").split() if x]


def parse_str_list(text: str) -> list[str]:
    return [x for x in str(text).replace(",", " ").split() if x]


def temperature_tag(T: float) -> str:
    return f"{float(T):g}".replace("-", "m").replace(".", "p")


def rate_tag(rate: float) -> str:
    return (
        f"{float(rate):.8g}"
        .replace("-", "m")
        .replace(".", "p")
        .replace("+", "p")
    )


def load_shape_rows(path: Path, families: Iterable[str]) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "target_class" not in df.columns:
        raise SystemExit(f"{path} does not contain target_class")
    wanted = list(families)
    missing = [x for x in wanted if x not in set(df.target_class.astype(str))]
    if missing:
        raise SystemExit(f"shape families missing from {path}: {missing}")
    return (
        df[df.target_class.astype(str).isin(wanted)]
        .drop_duplicates("target_class", keep="last")
        .set_index("target_class", drop=False)
    )


def sobol_parameters(n: int, seed: int) -> pd.DataFrame:
    if n <= 0:
        raise ValueError("samples per family must be positive")
    sampler = qmc.Sobol(d=len(PARAM_NAMES), scramble=True, seed=seed)
    m = int(math.floor(math.log2(n)))
    if 2**m == n:
        u = sampler.random_base2(m)
    else:
        u = sampler.random(n)
    lo = np.array([DEFAULT_BOUNDS[p][0] for p in PARAM_NAMES], dtype=float)
    hi = np.array([DEFAULT_BOUNDS[p][1] for p in PARAM_NAMES], dtype=float)
    x = qmc.scale(u, lo, hi)
    return pd.DataFrame(x, columns=PARAM_NAMES)


def exp_floor_barrier_eV(
    sigma_Pa: float | np.ndarray,
    T: np.ndarray,
    G00: np.ndarray,
    gT: np.ndarray,
    sigc0_GPa: np.ndarray,
    sT_GPa_per_K: float,
    a: float,
    n: float,
    floor_frac: float,
    Tref_K: float,
    floor_min_eV: float,
    floor_max_frac: float,
) -> np.ndarray:
    """Vectorized EXP-floor barrier matching FractureBarrier._exp_floor."""
    G0 = np.maximum(G00 + gT * (T - Tref_K), 1.0e-9)
    sigc = np.maximum(
        (sigc0_GPa + sT_GPa_per_K * (T - Tref_K)) * 1.0e9,
        1.0,
    )
    raw_floor = np.maximum(floor_min_eV, floor_frac * G0)
    floor = np.minimum(floor_max_frac * G0, raw_floor)
    sigma_arr = np.maximum(np.asarray(sigma_Pa, dtype=float), 0.0)
    x = sigma_arr / sigc
    return floor + (G0 - floor) * np.exp(-float(a) * np.power(x, float(n)))


def cleavage_effective_rate(
    G_eV: np.ndarray,
    T: np.ndarray,
    nu0_s: float,
    m_hits: float,
    tau_s: float,
) -> np.ndarray:
    raw = float(nu0_s) * np.exp(
        np.clip(-G_eV * EV_TO_J / (KB * T), -700.0, 0.0)
    )
    if m_hits <= 1.0 + 1.0e-12:
        return raw
    return gammainc(
        float(m_hits),
        np.minimum(raw * float(tau_s), 1.0e12),
    ) / float(tau_s)


def emission_rate(G_eV: np.ndarray, T: np.ndarray, nu0_s: float) -> np.ndarray:
    return float(nu0_s) * np.exp(
        np.clip(-G_eV * EV_TO_J / (KB * T), -700.0, 0.0)
    )


def evaluate_candidates(
    params: pd.DataFrame,
    shape: pd.Series,
    temperatures: list[float],
    Kdot: float,
    dK: float,
    Kmax: float,
    nu0_c: float,
    nu0_e: float,
    m_hits: float,
    tau_c: float,
    Tref_K: float,
    floor_min_eV: float,
    floor_max_frac: float,
    progress_every: int = 200,
) -> dict[str, np.ndarray]:
    """Integrate cleavage and emission hazards analytically on a K grid."""
    ncase = len(params)
    T = np.asarray(temperatures, dtype=float)[None, :]
    nt = T.shape[1]

    def col(name: str) -> np.ndarray:
        return params[name].to_numpy(dtype=float)[:, None]

    Gc00 = col("cleave_G00_eV")
    gcT = col("cleave_gT_eV_per_K")
    sc0 = col("cleave_sigc0_GPa")
    Ge00 = col("emit_G00_eV")
    geT = col("emit_gT_eV_per_K")
    se0 = col("emit_sigc0_GPa")

    c_sT = float(shape.get("cleave_sT_GPa_per_K", 0.0))
    c_a = float(shape["cleave_exp_a"])
    c_n = float(shape["cleave_exp_n"])
    c_ff = float(shape["cleave_floor_frac"])
    e_sT = float(shape.get("emit_sT_GPa_per_K", 0.0))
    e_a = float(shape["emit_exp_a"])
    e_n = float(shape["emit_exp_n"])
    e_ff = float(shape["emit_floor_frac"])
    r0 = max(float(shape["r_pz_m"]), 1.0e-30)

    B = np.zeros((ncase, nt), dtype=float)
    H = np.zeros_like(B)
    Kc = np.full_like(B, np.nan)
    H_at_Kc = np.full_like(B, np.nan)
    Gc_at_Kc = np.full_like(B, np.nan)
    Ge_at_Kc = np.full_like(B, np.nan)
    sigma_at_Kc_GPa = np.full_like(B, np.nan)

    nsystems = max(int(shape.get("mpz_n_systems", 1)), 1)
    sites_per_system = max(
        float(shape.get("mpz_source_sites_per_system", 1.0)), 0.0
    )
    ntotal = max(nsystems * sites_per_system, 1.0e-300)
    target_expected_one = -math.log(
        max(1.0 - min(1.0 / ntotal, 1.0 - 1.0e-15), 1.0e-15)
    )
    hazard_targets = {
        "K_emit_one_expected": target_expected_one,
        "K_emit_1pct": -math.log(0.99),
        "K_emit_10pct": -math.log(0.90),
        "K_emit_50pct": -math.log(0.50),
        "K_emit_90pct": -math.log(0.10),
    }
    K_emit = {name: np.full_like(B, np.nan) for name in hazard_targets}

    Gc0 = exp_floor_barrier_eV(
        0.0, T, Gc00, gcT, sc0, c_sT, c_a, c_n, c_ff,
        Tref_K, floor_min_eV, floor_max_frac,
    )
    Ge0 = exp_floor_barrier_eV(
        0.0, T, Ge00, geT, se0, e_sT, e_a, e_n, e_ff,
        Tref_K, floor_min_eV, floor_max_frac,
    )
    lc_prev = cleavage_effective_rate(Gc0, T, nu0_c, m_hits, tau_c)
    le_prev = emission_rate(Ge0, T, nu0_e)

    nsteps = int(math.ceil(Kmax / dK))
    Kprev = 0.0
    H_full = np.zeros_like(B)
    for istep in range(1, nsteps + 1):
        K = min(istep * dK, Kmax)
        step = K - Kprev
        sigma = K * 1.0e6 / math.sqrt(2.0 * math.pi * r0)
        Gc = exp_floor_barrier_eV(
            sigma, T, Gc00, gcT, sc0, c_sT, c_a, c_n, c_ff,
            Tref_K, floor_min_eV, floor_max_frac,
        )
        Ge = exp_floor_barrier_eV(
            sigma, T, Ge00, geT, se0, e_sT, e_a, e_n, e_ff,
            Tref_K, floor_min_eV, floor_max_frac,
        )
        lc = cleavage_effective_rate(Gc, T, nu0_c, m_hits, tau_c)
        le = emission_rate(Ge, T, nu0_e)
        dt = step / Kdot
        dB = 0.5 * (lc_prev + lc) * dt
        dH = 0.5 * (le_prev + le) * dt

        Bnew = B + dB
        Hnew = H + dH
        active = ~np.isfinite(Kc)
        crossed = active & (Bnew >= 1.0)
        if np.any(crossed):
            denom = np.maximum(Bnew - B, 1.0e-300)
            frac = np.clip((1.0 - B) / denom, 0.0, 1.0)
            Kcross = Kprev + frac * step
            Hcross = H + frac * dH
            Kc[crossed] = Kcross[crossed]
            H_at_Kc[crossed] = Hcross[crossed]
            sigma_cross = Kcross * 1.0e6 / math.sqrt(2.0 * math.pi * r0)
            Gccross = exp_floor_barrier_eV(
                sigma_cross, T, Gc00, gcT, sc0, c_sT, c_a, c_n, c_ff,
                Tref_K, floor_min_eV, floor_max_frac,
            )
            Gecross = exp_floor_barrier_eV(
                sigma_cross, T, Ge00, geT, se0, e_sT, e_a, e_n, e_ff,
                Tref_K, floor_min_eV, floor_max_frac,
            )
            Gc_at_Kc[crossed] = Gccross[crossed]
            Ge_at_Kc[crossed] = Gecross[crossed]
            sigma_at_Kc_GPa[crossed] = sigma_cross[crossed] / 1.0e9

        B[active] = Bnew[active]
        H[active] = Hnew[active]

        H_full += dH
        for name, Htarget in hazard_targets.items():
            unresolved = ~np.isfinite(K_emit[name])
            if np.any(unresolved):
                prev_full = H_full - dH
                hc = unresolved & (H_full >= Htarget)
                if np.any(hc):
                    denomH = np.maximum(dH, 1.0e-300)
                    fracH = np.clip((Htarget - prev_full) / denomH, 0.0, 1.0)
                    K_emit[name][hc] = (Kprev + fracH * step)[hc]

        lc_prev = lc
        le_prev = le
        Kprev = K
        if progress_every > 0 and (
            istep % progress_every == 0 or istep == nsteps
        ):
            resolved = int(np.isfinite(Kc).sum())
            print(
                f"  K={K:7.3f}/{Kmax:g} MPa√m: "
                f"resolved {resolved}/{ncase * nt}",
                flush=True,
            )
        if np.all(np.isfinite(Kc)) and all(
            np.all(np.isfinite(v)) for v in K_emit.values()
        ):
            break

    source_fraction = 1.0 - np.exp(-np.clip(H_at_Kc, 0.0, 700.0))
    expected_emitted = ntotal * source_fraction
    return {
        "Kc": Kc,
        "H_emit_at_Kc": H_at_Kc,
        "source_fraction_at_Kc": source_fraction,
        "expected_emitted_at_Kc": expected_emitted,
        "Gc_at_Kc_eV": Gc_at_Kc,
        "Ge_at_Kc_eV": Ge_at_Kc,
        "sigma_at_Kc_GPa": sigma_at_Kc_GPa,
        **K_emit,
    }


def add_temperature_columns(
    df: pd.DataFrame,
    result: dict[str, np.ndarray],
    temperatures: list[float],
) -> pd.DataFrame:
    cols: dict[str, np.ndarray] = {}
    for j, T in enumerate(temperatures):
        tag = temperature_tag(T)
        for key, arr in result.items():
            cols[f"{key}_T{tag}"] = arr[:, j]
    return pd.concat([df.reset_index(drop=True), pd.DataFrame(cols)], axis=1)


def finite_ptp(a: np.ndarray, axis: int = 1) -> np.ndarray:
    return np.nanmax(a, axis=axis) - np.nanmin(a, axis=axis)


def summarize_regions(
    df: pd.DataFrame,
    result: dict[str, np.ndarray],
    temperatures: list[float],
    shape: pd.Series,
) -> pd.DataFrame:
    Kc = result["Kc"]
    frac = result["source_fraction_at_Kc"]
    emitted = result["expected_emitted_at_Kc"]
    T = np.asarray(temperatures, dtype=float)
    low = T <= 500.0
    high = T >= 900.0

    Klow = np.nanmedian(Kc[:, low], axis=1)
    Khigh = np.nanmedian(Kc[:, high], axis=1)
    flow = np.nanmedian(frac[:, low], axis=1)
    fhigh = np.nanmedian(frac[:, high], axis=1)
    Elow = np.nanmedian(emitted[:, low], axis=1)
    Ehigh = np.nanmedian(emitted[:, high], axis=1)
    Krange = finite_ptp(Kc)
    dK = Khigh - Klow
    finite = np.all(np.isfinite(Kc), axis=1)
    dKsteps = np.diff(Kc, axis=1)
    dec_fraction = np.nanmean(dKsteps <= 0.25, axis=1)

    cross_T = np.full(len(df), np.nan)
    for i in range(len(df)):
        hit = np.where(frac[i] >= 0.10)[0]
        if hit.size:
            cross_T[i] = T[hit[0]]

    eps = 1.0e-12
    crossover_decades = np.log10((Ehigh + eps) / (Elow + eps))
    total_sites = max(
        float(shape.get("mpz_n_systems", 1))
        * float(shape.get("mpz_source_sites_per_system", 1.0)),
        eps,
    )

    ceramic = finite & (dK <= -3.0) & (fhigh <= 0.01) & (dec_fraction >= 0.8)
    weakT = finite & (Krange <= 3.0) & (fhigh >= 1.0e-4) & (fhigh <= 0.85)
    dbtt = (
        finite
        & (flow <= 0.01)
        & (fhigh >= 0.10)
        & (crossover_decades >= 2.0)
        & np.isfinite(cross_T)
        & (cross_T >= 600.0)
        & (cross_T <= 1000.0)
    )
    saturated = finite & ((flow >= 0.50) | (fhigh >= 0.95))

    region = np.full(len(df), "mixed_unclassified", dtype=object)
    region[ceramic] = "ceramic_intrinsic"
    region[weakT] = "weakT_intrinsic"
    region[dbtt] = "DBTT_precursor"
    region[saturated] = "emission_saturated"

    T300_idx = int(np.argmin(np.abs(T - 300.0)))
    T1200_idx = int(np.argmin(np.abs(T - 1200.0)))
    K300 = Kc[:, T300_idx]
    K1200 = Kc[:, T1200_idx]
    ceramic_score = (
        np.abs(K300 - 18.0) / 2.0
        + np.abs(K1200 - 6.5) / 2.0
        + 20.0 * fhigh
        + np.maximum(dK + 3.0, 0.0)
    )
    weakT_score = (
        np.abs(np.nanmean(Kc, axis=1) - 15.0) / 1.5
        + Krange / 2.5
        + np.maximum(1.0e-4 - fhigh, 0.0) * 1.0e4
        + np.maximum(fhigh - 0.70, 0.0) * 10.0
    )
    dbtt_score = (
        np.abs(Klow - 15.0) / 2.0
        + np.maximum(flow - 0.01, 0.0) * 100.0
        + np.maximum(0.10 - fhigh, 0.0) * 20.0
        + np.abs(np.nan_to_num(cross_T, nan=1400.0) - 800.0) / 200.0
        + np.maximum(2.0 - crossover_decades, 0.0)
        + np.maximum(fhigh - 0.90, 0.0) * 10.0
    )

    out = df.copy()
    out["K_low_median_MPa_sqrt_m"] = Klow
    out["K_high_median_MPa_sqrt_m"] = Khigh
    out["delta_K_high_minus_low_MPa_sqrt_m"] = dK
    out["K_range_MPa_sqrt_m"] = Krange
    out["source_fraction_lowT"] = flow
    out["source_fraction_highT"] = fhigh
    out["expected_emitted_lowT"] = Elow
    out["expected_emitted_highT"] = Ehigh
    out["emission_crossover_decades"] = crossover_decades
    out["emission_cross_10pct_T_K"] = cross_T
    out["temperature_decrease_fraction"] = dec_fraction
    out["all_temperatures_resolved"] = finite
    out["region"] = region
    out["ceramic_score"] = ceramic_score
    out["weakT_score"] = weakT_score
    out["DBTT_precursor_score"] = dbtt_score
    out["source_sites_total"] = total_sites
    return out


def shortlist(regions: pd.DataFrame, top_n: int) -> pd.DataFrame:
    rows = []
    specs = (
        ("ceramic_intrinsic", "ceramic_score"),
        ("weakT_intrinsic", "weakT_score"),
        ("DBTT_precursor", "DBTT_precursor_score"),
    )
    for region, score in specs:
        g = regions[regions.region == region].sort_values(score).head(top_n).copy()
        if not g.empty:
            g["shortlist_rank"] = np.arange(1, len(g) + 1)
            g["shortlist_score_name"] = score
            g["shortlist_score"] = g[score]
            rows.append(g)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=list(regions.columns))


def material_rows_from_shortlist(short: pd.DataFrame, shape_rows: pd.DataFrame) -> pd.DataFrame:
    if short.empty:
        return pd.DataFrame()
    class_map = {
        "ceramic_intrinsic": "ceramic",
        "weakT_intrinsic": "weakT",
        "DBTT_precursor": "DBTT",
    }
    rows = []
    for _, rec in short.iterrows():
        region = str(rec.get("region", ""))
        if region not in class_map:
            continue
        family = str(rec["shape_family"])
        if family not in shape_rows.index:
            continue
        row = shape_rows.loc[family].copy()
        for name in PARAM_NAMES:
            row[name] = float(rec[name])
        row["target_class"] = class_map[region]
        row["analytic_region"] = region
        row["analytic_shape_family"] = family
        row["analytic_candidate_id"] = rec["candidate_id"]
        row["analytic_shortlist_rank"] = int(rec["shortlist_rank"])
        row["analytic_shortlist_score"] = float(rec["shortlist_score"])
        row["Kdot_MPa_sqrt_m_per_s"] = float(rec["Kdot_MPa_sqrt_m_per_s"])
        row["status"] = "ANALYTIC_ATLAS_CANDIDATE_NOT_TRANSIENTLY_VALIDATED"
        rows.append(dict(row))
    return pd.DataFrame(rows)


def refine_shortlist(
    short: pd.DataFrame,
    shape_rows: pd.DataFrame,
    temperatures: list[float],
    args: argparse.Namespace,
) -> pd.DataFrame:
    if short.empty or args.refine_dK <= 0:
        return pd.DataFrame()
    frames = []
    for (Kdot, family), g in short.groupby(
        ["Kdot_MPa_sqrt_m_per_s", "shape_family"], sort=False
    ):
        if family not in shape_rows.index:
            continue
        p = g[list(PARAM_NAMES)].copy().reset_index(drop=True)
        res = evaluate_candidates(
            p, shape_rows.loc[family], temperatures, float(Kdot),
            args.refine_dK, args.Kmax, args.nu0_cleave, args.nu0_emit,
            args.multihit_m, args.multihit_tau, args.Tref_K,
            args.floor_min_eV, args.floor_max_frac, progress_every=0,
        )
        base = g.reset_index(drop=True).copy()
        cols: dict[str, np.ndarray] = {}
        for j, T in enumerate(temperatures):
            tag = temperature_tag(T)
            for key, arr in res.items():
                cols[f"refined_{key}_T{tag}"] = arr[:, j]
        base = pd.concat([base, pd.DataFrame(cols)], axis=1)
        base["refine_dK_MPa_sqrt_m"] = args.refine_dK
        frames.append(base)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def evaluate_anchor_table(
    path: Path,
    shape_rows: pd.DataFrame,
    temperatures: list[float],
    Kdot: float,
    args: argparse.Namespace,
) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    src = pd.read_csv(path)
    if "target_class" not in src.columns:
        return pd.DataFrame()
    rows = []
    for klass, g in src.groupby(src.target_class.astype(str), sort=False):
        if klass not in shape_rows.index:
            continue
        row = g.iloc[-1]
        if any(p not in row.index for p in PARAM_NAMES):
            continue
        p = pd.DataFrame([{name: float(row[name]) for name in PARAM_NAMES}])
        res = evaluate_candidates(
            p, shape_rows.loc[klass], temperatures, Kdot, args.dK, args.Kmax,
            args.nu0_cleave, args.nu0_emit, args.multihit_m, args.multihit_tau,
            args.Tref_K, args.floor_min_eV, args.floor_max_frac, progress_every=0,
        )
        rec = {
            "anchor_file": str(path),
            "target_class": klass,
            **{name: float(row[name]) for name in PARAM_NAMES},
        }
        for j, T in enumerate(temperatures):
            tag = temperature_tag(T)
            for key in (
                "Kc", "H_emit_at_Kc", "source_fraction_at_Kc",
                "expected_emitted_at_Kc",
            ):
                rec[f"{key}_T{tag}"] = float(res[key][0, j])
        rows.append(rec)
    return pd.DataFrame(rows)


def make_plots(
    regions: pd.DataFrame,
    short: pd.DataFrame,
    temperatures: list[float],
    out: Path,
) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return

    fig, ax = plt.subplots(figsize=(8.2, 6.2))
    groups = [
        ("mixed_unclassified", 0.15, 8),
        ("ceramic_intrinsic", 0.8, 16),
        ("weakT_intrinsic", 0.8, 16),
        ("DBTT_precursor", 0.9, 20),
        ("emission_saturated", 0.25, 10),
    ]
    for name, alpha, size in groups:
        g = regions[regions.region == name]
        if g.empty:
            continue
        ax.scatter(
            g["delta_K_high_minus_low_MPa_sqrt_m"],
            g["emission_crossover_decades"],
            s=size, alpha=alpha, label=name,
        )
    ax.axhline(2.0, lw=1.0, ls="--")
    ax.axvline(0.0, lw=1.0, ls=":")
    ax.set_xlabel(r"$K_{high}-K_{low}$ [MPa$\sqrt{m}$]")
    ax.set_ylabel("Emission crossover [decades]")
    ax.set_title("Analytical intrinsic first-passage atlas")
    ax.legend(fontsize=8, loc="best")
    ax.grid(True, alpha=0.2)
    fig.tight_layout()
    fig.savefig(out / "analytic_atlas_region_map.png", dpi=200)
    plt.close(fig)

    if short.empty:
        return
    for ykey, ylabel, filename in (
        ("Kc", r"$K_c$ [MPa$\sqrt{m}$]", "analytic_atlas_shortlist_Kc_T.png"),
        (
            "source_fraction_at_Kc",
            "Source fraction emitted before cleavage",
            "analytic_atlas_shortlist_emission_T.png",
        ),
    ):
        fig, ax = plt.subplots(figsize=(8.6, 6.2))
        for _, r in short.iterrows():
            y = [r.get(f"{ykey}_T{temperature_tag(T)}", np.nan) for T in temperatures]
            ax.plot(
                temperatures, y, marker="o", ms=3,
                label=f"{r['region']} {r['candidate_id']}",
            )
        ax.set_xlabel("Temperature [K]")
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=7, ncol=2)
        fig.tight_layout()
        fig.savefig(out / filename, dpi=200)
        plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--initial", type=Path, default=Path("mpz_three_class_initial_guesses.csv"))
    ap.add_argument("--shape-families", default="ceramic weakT DBTT")
    ap.add_argument("--temperatures", default="300 400 500 600 700 800 900 1000 1100 1200")
    ap.add_argument("--Kdot-values", default="0.005")
    ap.add_argument("--samples-per-family", type=int, default=16384)
    ap.add_argument("--seed", type=int, default=92031)
    ap.add_argument("--dK", type=float, default=0.05)
    ap.add_argument("--Kmax", type=float, default=80.0)
    ap.add_argument("--nu0-cleave", type=float, default=1.0e12)
    ap.add_argument("--nu0-emit", type=float, default=1.0e11)
    ap.add_argument("--multihit-m", type=float, default=3.0)
    ap.add_argument("--multihit-tau", type=float, default=1.0e-6)
    ap.add_argument("--Tref-K", type=float, default=481.33)
    ap.add_argument("--floor-min-eV", type=float, default=1.0e-4)
    ap.add_argument("--floor-max-frac", type=float, default=0.95)
    ap.add_argument("--top-per-region", type=int, default=20)
    ap.add_argument(
        "--refine-dK", type=float, default=0.01,
        help="Re-evaluate shortlisted rows at this finer K increment; <=0 disables.",
    )
    ap.add_argument("--anchor-tables", default="")
    ap.add_argument(
        "--out", type=Path,
        default=Path("runs/mpz_v9_2_analytic_first_passage_atlas"),
    )
    ap.add_argument("--progress-every", type=int, default=200)
    a = ap.parse_args()

    if a.dK <= 0 or a.Kmax <= 0 or a.samples_per_family <= 0:
        raise SystemExit("dK, Kmax, and samples-per-family must be positive")
    temperatures = parse_float_list(a.temperatures)
    rates = parse_float_list(a.Kdot_values)
    families = parse_str_list(a.shape_families)
    if not temperatures or not rates or not families:
        raise SystemExit("temperatures, Kdot values, and shape families may not be empty")

    out = a.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    shapes = load_shape_rows(a.initial, families)

    all_candidates = []
    all_regions = []
    all_short = []
    anchors = []
    run_index = 0

    for irate, Kdot in enumerate(rates):
        if Kdot <= 0:
            raise SystemExit("Kdot values must be positive")
        for ifam, family in enumerate(families):
            run_index += 1
            print(
                f"[{run_index}/{len(rates) * len(families)}] "
                f"family={family} Kdot={Kdot:g} samples={a.samples_per_family}",
                flush=True,
            )
            shape = shapes.loc[family]
            params = sobol_parameters(
                a.samples_per_family, a.seed + 1009 * irate + 7919 * ifam
            )
            params.insert(0, "candidate_index", np.arange(len(params), dtype=int))
            params.insert(
                0, "candidate_id",
                [f"{family}_r{rate_tag(Kdot)}_{i:07d}" for i in range(len(params))],
            )
            params.insert(0, "shape_family", family)
            params["Kdot_MPa_sqrt_m_per_s"] = Kdot
            result = evaluate_candidates(
                params, shape, temperatures, Kdot, a.dK, a.Kmax,
                a.nu0_cleave, a.nu0_emit, a.multihit_m, a.multihit_tau,
                a.Tref_K, a.floor_min_eV, a.floor_max_frac,
                progress_every=a.progress_every,
            )
            wide = add_temperature_columns(params, result, temperatures)
            regions = summarize_regions(params, result, temperatures, shape)
            short = shortlist(regions, a.top_per_region)
            if not short.empty:
                temp_cols = [c for c in wide.columns if c not in params.columns]
                short = short.merge(
                    wide[["candidate_id", *temp_cols]],
                    on="candidate_id", how="left", validate="one_to_one",
                )
            all_candidates.append(wide)
            all_regions.append(regions)
            all_short.append(short)

        for ptext in parse_str_list(a.anchor_tables):
            p = Path(ptext)
            adf = evaluate_anchor_table(p, shapes, temperatures, Kdot, a)
            if not adf.empty:
                adf["Kdot_MPa_sqrt_m_per_s"] = Kdot
                anchors.append(adf)

    candidates_df = pd.concat(all_candidates, ignore_index=True)
    regions_df = pd.concat(all_regions, ignore_index=True)
    short_df = pd.concat(all_short, ignore_index=True) if all_short else pd.DataFrame()
    anchor_df = pd.concat(anchors, ignore_index=True) if anchors else pd.DataFrame()
    refined_df = refine_shortlist(short_df, shapes, temperatures, a)
    material_df = material_rows_from_shortlist(short_df, shapes)

    candidates_df.to_csv(
        out / "analytic_first_passage_atlas_candidates.csv.gz",
        index=False, compression="gzip",
    )
    regions_df.to_csv(out / "analytic_first_passage_atlas_regions.csv", index=False)
    short_df.to_csv(out / "analytic_first_passage_atlas_shortlist.csv", index=False)
    if not refined_df.empty:
        refined_df.to_csv(out / "analytic_first_passage_atlas_shortlist_refined.csv", index=False)
    if not material_df.empty:
        material_df.to_csv(out / "mpz_analytic_shortlist_material_rows.csv", index=False)
    if not anchor_df.empty:
        anchor_df.to_csv(out / "analytic_first_passage_anchor_predictions.csv", index=False)

    counts = (
        regions_df.groupby(["Kdot_MPa_sqrt_m_per_s", "shape_family", "region"])
        .size().rename("n_candidates").reset_index()
    )
    counts.to_csv(out / "analytic_first_passage_region_counts.csv", index=False)
    make_plots(regions_df, short_df, temperatures, out)

    config = vars(a).copy()
    for k, v in list(config.items()):
        if isinstance(v, Path):
            config[k] = str(v)
    config.update({
        "temperatures_resolved": temperatures,
        "Kdot_values_resolved": rates,
        "shape_families_resolved": families,
        "parameter_bounds": {k: list(v) for k, v in DEFAULT_BOUNDS.items()},
        "interpretation": {
            "ceramic_intrinsic": "cleavage first passage with negligible emission exposure",
            "weakT_intrinsic": "nearly flat virgin first passage with finite but nonsaturated emission exposure",
            "DBTT_precursor": "low-T cleavage dominance and high-T emission exposure; requires transient MPZ validation",
            "emission_saturated": "source inventory mostly exhausted before cleavage; reject or audit separately",
        },
    })
    (out / "analytic_first_passage_atlas_config.json").write_text(
        json.dumps(config, indent=2)
    )

    print("\nAnalytical first-passage atlas complete", flush=True)
    print(counts.to_string(index=False), flush=True)
    print(f"Outputs: {out}", flush=True)


if __name__ == "__main__":
    main()
