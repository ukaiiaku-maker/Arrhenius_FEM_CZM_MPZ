#!/usr/bin/env python3
"""Assemble the first-pass V1 Figure 1 atlas (Panels A, C, F).

This script is deliberately a data-preserving assembly layer. It does not
re-implement fracture, fatigue, or strength physics. It reads solver outputs,
extracts explicit metrics, writes panel-specific CSVs, and produces provisional
publication-style PNG/PDF plots that can later be re-styled from the saved data.

Panels in the first pass:
  A - four canonical Kc(T) response classes
  C - endurance-like vs non-endurance-like fatigue response
  F - DBTT amplitude index vs endurance index across six curated V1 cases

Metric definitions:
  I_DBTT_amp = [max_T Kc(T) - min_T Kc(T)] / Kc(T_min)
  T_inflection = T at max |dKc/dT|
  DeltaKth(T) = rate-defined threshold for da/dN = criterion
  I_endurance = DeltaKth(T_ref) / Kc(T_ref)
  T_pers_50 = highest T with DeltaKth(T) >= 0.5 DeltaKth(T_ref)

By default the cross-case metric scenario is S_emit=-40 kB, S_cleave=0 kB,
which is shared by the six-case refined core study.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Dict, Iterable, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DEFAULT_CASE_ORDER = [
    "FCC_like_case29",
    "shifted_ductile_case64",
    "steep_cleavage_case35",
    "slow_threshold_case101",
    "higher_barrier_case171",
    "plastic_shielded_case64_M1",
]

DEFAULT_DISPLAY = {
    "FCC_like_case29": "FCC-like / non-endurance",
    "shifted_ductile_case64": "Shifted ductile",
    "steep_cleavage_case35": "Steep cleavage",
    "slow_threshold_case101": "Slow threshold",
    "higher_barrier_case171": "Higher barrier",
    "plastic_shielded_case64_M1": "Plastic-shielded / endurance-like",
}

DEFAULT_PLOT_LABEL = {
    "FCC_like_case29": "FCC-like",
    "shifted_ductile_case64": "Shifted ductile",
    "steep_cleavage_case35": "Steep cleavage",
    "slow_threshold_case101": "Slow threshold",
    "higher_barrier_case171": "Higher barrier",
    "plastic_shielded_case64_M1": "Plastic-shielded",
}

LABEL_OFFSET_POINTS = {
    "FCC_like_case29": (5, -2),
    "shifted_ductile_case64": (5, -13),
    "steep_cleavage_case35": (5, 4),
    "slow_threshold_case101": (5, 4),
    "higher_barrier_case171": (5, 4),
    "plastic_shielded_case64_M1": (5, 9),
}

A_ORDER = ["ceramic", "peak", "weakT", "dbtt"]
A_DISPLAY = {
    "ceramic": "Ceramic-like",
    "peak": "Peak",
    "weakT": "Weak-T",
    "dbtt": "DBTT-like",
}


def finite_or_nan(x) -> float:
    try:
        v = float(x)
    except Exception:
        return float("nan")
    return v if math.isfinite(v) else float("nan")


def interp_at(x: np.ndarray, y: np.ndarray, x0: float) -> float:
    mask = np.isfinite(x) & np.isfinite(y)
    x = np.asarray(x[mask], float)
    y = np.asarray(y[mask], float)
    if x.size == 0:
        return float("nan")
    order = np.argsort(x)
    x, y = x[order], y[order]
    if x0 < x[0] or x0 > x[-1]:
        return float("nan")
    return float(np.interp(float(x0), x, y))


def read_manifest(path: Optional[Path]) -> pd.DataFrame:
    if path is None:
        return pd.DataFrame()
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path)


def display_map(manifest: pd.DataFrame) -> Dict[str, str]:
    out = dict(DEFAULT_DISPLAY)
    if not manifest.empty and {"case_id", "display_label"}.issubset(manifest.columns):
        for _, r in manifest.dropna(subset=["case_id", "display_label"]).iterrows():
            out[str(r["case_id"])] = str(r["display_label"])
    return out


def load_panel_a(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    aliases = {
        "regime_key": ["regime_key", "regime_id", "key"],
        "regime": ["regime", "display", "class_label"],
        "T_K": ["T_K", "T", "temperature_K"],
        "Kc_pred_MPa_sqrtm": ["Kc_pred_MPa_sqrtm", "Kc_model_MPa_sqrt_m", "Kc_MPa_sqrtm", "Kc_first_MPa_sqrtm"],
    }
    ren = {}
    for target, names in aliases.items():
        for n in names:
            if n in df.columns:
                ren[n] = target
                break
    df = df.rename(columns=ren)
    need = ["regime_key", "T_K", "Kc_pred_MPa_sqrtm"]
    miss = [c for c in need if c not in df.columns]
    if miss:
        raise ValueError(f"Panel A source missing columns {miss}; got {list(df.columns)}")
    if "regime" not in df.columns:
        df["regime"] = df["regime_key"].map(A_DISPLAY).fillna(df["regime_key"])
    out = df[df["regime_key"].isin(A_ORDER)].copy()
    out = out[[c for c in ["regime_key", "regime", "T_K", "Kc_pred_MPa_sqrtm",
                          "cleave_H0_eV", "chi_shield", "N_sat"] if c in out.columns]]
    return out.sort_values(["regime_key", "T_K"]).reset_index(drop=True)


def scenario_filter(df: pd.DataFrame, Se: float, Sc: float) -> pd.DataFrame:
    out = df.copy()
    if "S_emit_kB" in out.columns:
        out = out[np.isclose(pd.to_numeric(out["S_emit_kB"], errors="coerce"), Se)]
    if "S_cleave_kB" in out.columns:
        out = out[np.isclose(pd.to_numeric(out["S_cleave_kB"], errors="coerce"), Sc)]
    return out


def extract_panel_c(fatigue: pd.DataFrame, cases: list[str], T: float, Se: float, Sc: float,
                    disp: Dict[str, str]) -> pd.DataFrame:
    df = scenario_filter(fatigue, Se, Sc)
    df = df[df["case_label"].isin(cases)].copy()
    df = df[np.isclose(pd.to_numeric(df["T_K"], errors="coerce"), T)]
    df["display_label"] = df["case_label"].map(disp).fillna(df["case_label"])
    keep = [
        "case_label", "display_label", "S_emit_kB", "S_cleave_kB", "T_K",
        "Kmax_MPa_sqrtm", "DeltaK_MPa_sqrtm", "cycles_total",
        "cycles_to_first_fire", "n_adv", "a_adv_m", "da_dN_m_per_cycle",
        "da_dN_upper_bound_m_per_cycle", "status", "direct_lt_1_cycle",
        "mu_emit_initial_per_cycle", "mu_cleave_initial_per_cycle",
        "log10_mu_emit_over_cleave_initial",
    ]
    keep = [c for c in keep if c in df.columns]
    return df[keep].sort_values(["case_label", "DeltaK_MPa_sqrtm"]).reset_index(drop=True)


def threshold_from_rate_points(g: pd.DataFrame, criterion: float) -> Dict[str, float | str]:
    g = g.sort_values("DeltaK_MPa_sqrtm").copy()
    dk = pd.to_numeric(g["DeltaK_MPa_sqrtm"], errors="coerce").to_numpy(float)
    measured = pd.to_numeric(g.get("da_dN_m_per_cycle", np.nan), errors="coerce").to_numpy(float)
    upper = pd.to_numeric(g.get("da_dN_upper_bound_m_per_cycle", np.nan), errors="coerce").to_numpy(float)

    low_dk: list[float] = []
    high_dk: list[float] = []
    for x, r, ub in zip(dk, measured, upper):
        if not math.isfinite(x):
            continue
        if math.isfinite(r):
            (high_dk if r >= criterion else low_dk).append(float(x))
        elif math.isfinite(ub) and ub < criterion:
            low_dk.append(float(x))

    lo = max(low_dk) if low_dk else float("nan")
    hi_candidates = [x for x in high_dk if not math.isfinite(lo) or x >= lo]
    hi = min(hi_candidates) if hi_candidates else (min(high_dk) if high_dk else float("nan"))

    if math.isfinite(lo) and math.isfinite(hi) and hi >= lo:
        est = 0.5 * (lo + hi)
        cls = "bracketed_midpoint"
    elif math.isfinite(hi):
        est = hi
        cls = "first_measured_above_criterion"
    else:
        est = float("nan")
        cls = "unresolved_above_scan"
    return {
        "DeltaK_threshold_lower_MPa_sqrtm": lo,
        "DeltaK_threshold_upper_MPa_sqrtm": hi,
        "DeltaK_threshold_estimate_MPa_sqrtm": est,
        "threshold_class": cls,
    }


def thresholds_from_fatigue(fatigue: pd.DataFrame, criterion: float, Se: float, Sc: float) -> pd.DataFrame:
    df = scenario_filter(fatigue, Se, Sc)
    rows = []
    for (case, T), g in df.groupby(["case_label", "T_K"], sort=False):
        rec = threshold_from_rate_points(g, criterion)
        rec.update({
            "case_label": case,
            "S_emit_kB": Se,
            "S_cleave_kB": Sc,
            "T_K": float(T),
            "da_dN_criterion_m_per_cycle": criterion,
            "threshold_source": "derived_from_fatigue_points",
        })
        rows.append(rec)
    return pd.DataFrame(rows)


def load_thresholds(path: Optional[Path], fatigue: pd.DataFrame, criterion: float, Se: float, Sc: float) -> pd.DataFrame:
    if path is None or not path.exists():
        return thresholds_from_fatigue(fatigue, criterion, Se, Sc)
    thr = pd.read_csv(path)
    thr = scenario_filter(thr, Se, Sc)
    ccol = "da_dN_criterion_m_per_cycle"
    if ccol in thr.columns:
        thr = thr[np.isclose(pd.to_numeric(thr[ccol], errors="coerce"), criterion, rtol=1e-8, atol=0.0)]
    if "DeltaK_threshold_estimate_MPa_sqrtm" not in thr.columns:
        return thresholds_from_fatigue(fatigue, criterion, Se, Sc)
    thr = thr.copy()
    thr["threshold_source"] = "rate_defined_threshold_table"
    return thr


def monotonic_metrics(monotonic: pd.DataFrame, Se: float, Sc: float,
                      T_ref: float, T_hi_ref: float) -> pd.DataFrame:
    df = scenario_filter(monotonic, Se, Sc)
    rows = []
    for case, g in df.groupby("case_label", sort=False):
        g = g.sort_values("T_K")
        T = pd.to_numeric(g["T_K"], errors="coerce").to_numpy(float)
        K = pd.to_numeric(g["Kc_first_MPa_sqrtm"], errors="coerce").to_numpy(float)
        mask = np.isfinite(T) & np.isfinite(K)
        T, K = T[mask], K[mask]
        if len(T) < 2:
            continue
        order = np.argsort(T)
        T, K = T[order], K[order]
        dKdT = np.gradient(K, T)
        i = int(np.nanargmax(np.abs(dKdT)))
        K_Tmin = float(K[0])
        K_Tmax = float(K[-1])
        Kc_ref = interp_at(T, K, T_ref)
        Kc_hi = interp_at(T, K, T_hi_ref)
        amp = (float(np.nanmax(K)) - float(np.nanmin(K))) / max(abs(K_Tmin), 1e-12)
        slope_index = float(np.nanmax(np.abs(dKdT))) * 100.0 / max(abs(K_Tmin), 1e-12)

        # Transition midpoint temperature: crossing of the half-amplitude level
        # between the observed low and high Kc asymptotes. This is robust to
        # either increasing or decreasing Kc(T).
        K_hi_curve = float(np.nanmax(K))
        K_lo_curve = float(np.nanmin(K))
        K_half = K_lo_curve + 0.5 * (K_hi_curve - K_lo_curve)
        T_half = float("nan")
        for j in range(len(T) - 1):
            y0, y1 = K[j] - K_half, K[j + 1] - K_half
            if y0 == 0:
                T_half = float(T[j]); break
            if y0 * y1 <= 0 and K[j + 1] != K[j]:
                f = (K_half - K[j]) / (K[j + 1] - K[j])
                T_half = float(T[j] + f * (T[j + 1] - T[j]))
                break

        rows.append({
            "case_label": case,
            "T_min_K": float(T[0]),
            "T_max_K": float(T[-1]),
            "Kc_Tmin_MPa_sqrtm": K_Tmin,
            "Kc_Tmax_MPa_sqrtm": K_Tmax,
            "Kc_ref_MPa_sqrtm": Kc_ref,
            "Kc_hi_ref_MPa_sqrtm": Kc_hi,
            "Kc_hi_minus_ref_MPa_sqrtm": Kc_hi - Kc_ref if math.isfinite(Kc_hi) and math.isfinite(Kc_ref) else float("nan"),
            "I_DBTT_amplitude": amp,
            "I_DBTT_max_abs_slope_per_100K": slope_index,
            "T_transition_half_K": T_half,
            "Kc_transition_half_MPa_sqrtm": K_half,
            "T_inflection_K": float(T[i]),
            "dKc_dT_at_inflection_MPa_sqrtm_per_K": float(dKdT[i]),
        })
    return pd.DataFrame(rows)


def endurance_metrics(thr: pd.DataFrame, T_ref: float, disp: Dict[str, str]) -> pd.DataFrame:
    rows = []
    for case, g in thr.groupby("case_label", sort=False):
        g = g.sort_values("T_K")
        T = pd.to_numeric(g["T_K"], errors="coerce").to_numpy(float)
        D = pd.to_numeric(g["DeltaK_threshold_estimate_MPa_sqrtm"], errors="coerce").to_numpy(float)
        Dref = interp_at(T, D, T_ref)
        Tpers = float("nan")
        if math.isfinite(Dref) and Dref > 0:
            good = np.isfinite(T) & np.isfinite(D) & (D >= 0.5 * Dref)
            if np.any(good):
                Tpers = float(np.nanmax(T[good]))
        rows.append({
            "case_label": case,
            "display_label": disp.get(case, case),
            "DeltaKth_ref_MPa_sqrtm": Dref,
            "T_pers_50_K": Tpers,
            "DeltaKth_reference_temperature_K": T_ref,
            "persistence_fraction": 0.5,
        })
    return pd.DataFrame(rows)


def cross_case_metrics(monotonic: pd.DataFrame, thr: pd.DataFrame, cases: list[str],
                       Se: float, Sc: float, T_ref: float, T_hi_ref: float,
                       disp: Dict[str, str]) -> pd.DataFrame:
    m = monotonic_metrics(monotonic, Se, Sc, T_ref, T_hi_ref)
    e = endurance_metrics(thr, T_ref, disp)
    out = m.merge(e, on="case_label", how="outer")
    out = out[out["case_label"].isin(cases)].copy()
    out["display_label"] = out["case_label"].map(disp).fillna(out.get("display_label", out["case_label"]))
    out["I_endurance_DKth_over_Kc_ref"] = out["DeltaKth_ref_MPa_sqrtm"] / out["Kc_ref_MPa_sqrtm"]
    out["S_emit_kB"] = Se
    out["S_cleave_kB"] = Sc
    out["metric_temperature_ref_K"] = T_ref
    order_map = {c: i for i, c in enumerate(cases)}
    out["_order"] = out["case_label"].map(order_map)
    out = out.sort_values("_order").drop(columns="_order")
    return out.reset_index(drop=True)


def save_fig(fig, base: Path) -> None:
    fig.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def plot_panel_a(df: pd.DataFrame, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.6, 4.8))
    for key in A_ORDER:
        g = df[df["regime_key"] == key].sort_values("T_K")
        if g.empty:
            continue
        ax.plot(g["T_K"], g["Kc_pred_MPa_sqrtm"], lw=2.2, label=A_DISPLAY.get(key, key))
    ax.set_xlabel("Temperature, T (K)")
    ax.set_ylabel(r"Fracture toughness, $K_c$ (MPa $\sqrt{m}$)")
    ax.set_title("A  Toughness classes")
    ax.legend(frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(direction="out")
    fig.tight_layout()
    save_fig(fig, out / "panel_A_toughness_classes")


def plot_panel_c(df: pd.DataFrame, out: Path, case_order: list[str]) -> None:
    fig, ax = plt.subplots(figsize=(6.6, 4.8))
    for case in case_order:
        g = df[df["case_label"] == case].sort_values("DeltaK_MPa_sqrtm")
        if g.empty:
            continue
        label = str(g.iloc[0]["display_label"])
        m = np.isfinite(pd.to_numeric(g.get("da_dN_m_per_cycle", np.nan), errors="coerce"))
        gm = g[m]
        if not gm.empty:
            ax.plot(gm["DeltaK_MPa_sqrtm"], gm["da_dN_m_per_cycle"], marker="o", lw=2.0, label=label)
        gc = g[~m].copy()
        if not gc.empty and "da_dN_upper_bound_m_per_cycle" in gc.columns:
            ub = pd.to_numeric(gc["da_dN_upper_bound_m_per_cycle"], errors="coerce")
            ok = np.isfinite(ub)
            ax.scatter(gc.loc[ok, "DeltaK_MPa_sqrtm"], ub[ok], marker="v", facecolors="none", s=55)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"Cyclic driving range, $\Delta K$ (MPa $\sqrt{m}$)")
    ax.set_ylabel(r"Crack-growth rate, $da/dN$ (m cycle$^{-1}$)")
    ax.set_title("C  Fatigue classes")
    ax.legend(frameon=False, fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(direction="out")
    fig.tight_layout()
    save_fig(fig, out / "panel_C_fatigue_classes")


def plot_panel_f(metrics: pd.DataFrame, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.6, 4.8))
    x = pd.to_numeric(metrics["T_transition_half_K"], errors="coerce")
    y = pd.to_numeric(metrics["DeltaKth_ref_MPa_sqrtm"], errors="coerce")
    ok = np.isfinite(x) & np.isfinite(y)
    ax.scatter(x[ok], y[ok], s=70)
    for _, r in metrics.loc[ok].iterrows():
        case = str(r["case_label"])
        dx, dy = LABEL_OFFSET_POINTS.get(case, (5, 5))
        ax.annotate(DEFAULT_PLOT_LABEL.get(case, str(r["display_label"])),
                    (float(r["T_transition_half_K"]), float(r["DeltaKth_ref_MPa_sqrtm"])),
                    xytext=(dx, dy), textcoords="offset points", fontsize=8)
    ax.set_xlabel(r"DBTT transition temperature, $T_{50}^{K_c}$ (K)")
    ax.set_ylabel(r"Rate-defined endurance threshold, $\Delta K_{th}$ (MPa $\sqrt{m}$)")
    ax.set_title("F  DBTT-endurance association")
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(direction="out")
    fig.tight_layout()
    save_fig(fig, out / "panel_F_DBTT_endurance_association")


def combined_figure(a: pd.DataFrame, c: pd.DataFrame, f: pd.DataFrame, out: Path,
                    c_cases: list[str]) -> None:
    fig, axs = plt.subplots(1, 3, figsize=(15.5, 4.7))
    ax = axs[0]
    for key in A_ORDER:
        g = a[a["regime_key"] == key].sort_values("T_K")
        if not g.empty:
            ax.plot(g["T_K"], g["Kc_pred_MPa_sqrtm"], lw=2.0, label=A_DISPLAY.get(key, key))
    ax.set_xlabel("T (K)")
    ax.set_ylabel(r"$K_c$ (MPa $\sqrt{m}$)")
    ax.set_title("A  Toughness classes")
    ax.legend(frameon=False, fontsize=8)

    ax = axs[1]
    for case in c_cases:
        g = c[c["case_label"] == case].sort_values("DeltaK_MPa_sqrtm")
        if g.empty:
            continue
        label = str(g.iloc[0]["display_label"])
        m = np.isfinite(pd.to_numeric(g.get("da_dN_m_per_cycle", np.nan), errors="coerce"))
        gm = g[m]
        if not gm.empty:
            ax.plot(gm["DeltaK_MPa_sqrtm"], gm["da_dN_m_per_cycle"], marker="o", lw=1.8, label=label)
        gc = g[~m]
        if not gc.empty and "da_dN_upper_bound_m_per_cycle" in gc.columns:
            ub = pd.to_numeric(gc["da_dN_upper_bound_m_per_cycle"], errors="coerce")
            ok2 = np.isfinite(ub)
            ax.scatter(gc.loc[ok2, "DeltaK_MPa_sqrtm"], ub[ok2], marker="v", facecolors="none", s=38)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"$\Delta K$ (MPa $\sqrt{m}$)")
    ax.set_ylabel(r"$da/dN$ (m cycle$^{-1}$)")
    ax.set_title("C  Fatigue classes")
    ax.legend(frameon=False, fontsize=7)

    ax = axs[2]
    x = pd.to_numeric(f["T_transition_half_K"], errors="coerce")
    y = pd.to_numeric(f["DeltaKth_ref_MPa_sqrtm"], errors="coerce")
    ok = np.isfinite(x) & np.isfinite(y)
    ax.scatter(x[ok], y[ok], s=55)
    for _, r in f.loc[ok].iterrows():
        case = str(r["case_label"])
        dx, dy = LABEL_OFFSET_POINTS.get(case, (4, 4))
        ax.annotate(DEFAULT_PLOT_LABEL.get(case, str(r["display_label"])),
                    (float(r["T_transition_half_K"]), float(r["DeltaKth_ref_MPa_sqrtm"])),
                    xytext=(dx, dy), textcoords="offset points", fontsize=6.7)
    ax.set_xlabel(r"$T_{50}^{K_c}$ (K)")
    ax.set_ylabel(r"$\Delta K_{th}$ (MPa $\sqrt{m}$)")
    ax.set_title("F  DBTT-endurance association")

    for ax in axs:
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(direction="out")
    fig.tight_layout(w_pad=2.2)
    save_fig(fig, out / "figure1_first_pass_ACF")


def write_readme(out: Path, args, metrics: pd.DataFrame) -> None:
    text = "# Figure 1 V1 atlas - first pass\n\n"
    text += "This directory contains a data-preserving first-pass assembly of Panels A, C, and F.\n\n"
    text += "## Metric definitions\n\n"
    text += "- primary Panel F DBTT metric: `T_transition_half_K`, the half-amplitude Kc transition temperature\n"
    text += "- secondary DBTT amplitude: `I_DBTT_amplitude = (max(Kc)-min(Kc))/Kc(Tmin)`\n"
    text += "- `T_inflection_K`: temperature at maximum absolute finite-difference slope `|dKc/dT|`\n"
    text += f"- primary fatigue threshold: `da/dN = {args.primary_rate:g} m/cycle`\n"
    text += f"- primary Panel F endurance metric: rate-defined `DeltaKth({args.metric_T_ref:g} K)`\n"
    text += f"- secondary normalized endurance index: `I_endurance = DeltaKth({args.metric_T_ref:g} K)/Kc({args.metric_T_ref:g} K)`\n"
    text += f"- `T_pers_50`: highest T retaining at least 50% of DeltaKth({args.metric_T_ref:g} K)\n\n"
    text += "## Scenario used for cross-case metrics\n\n"
    text += f"`S_emit = {args.S_emit:g} kB`, `S_cleave = {args.S_cleave:g} kB`.\n\n"
    text += "## Source-family note\n\n"
    text += ("Panel A uses the four canonical monotonic fracture-regime definitions from the forward Kc(T) workflow. "
             "Panels C and F use the six curated V1 fatigue/material-response cases from the refined two-barrier workflow. "
             "These are retained as distinct source families instead of being relabeled as the same parameter sets.\n\n")
    text += "## Replot-ready data\n\n"
    text += "- `panel_A_Kc_temperature.csv`\n- `panel_C_fatigue_contrast.csv`\n- `panel_F_DBTT_endurance_metrics.csv`\n"
    text += "- `thresholds_used_for_metrics.csv`\n- `case_metrics_summary.csv`\n- `atlas_settings.json`\n\n"
    text += "## Provisional plots\n\n"
    text += "Each panel and the combined A-C-F figure are exported as both PNG and vector PDF.\n\n"
    text += "## Current cross-case metric values\n\n"
    if not metrics.empty:
        text += metrics[[c for c in ["display_label", "T_transition_half_K", "I_DBTT_amplitude", "DeltaKth_ref_MPa_sqrtm",
                                      "I_endurance_DKth_over_Kc_ref", "T_inflection_K", "T_pers_50_K"]
                         if c in metrics.columns]].to_markdown(index=False)
        text += "\n"
    (out / "README_FIGURE1_V1_ATLAS_FIRST_PASS.md").write_text(text, encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--panel-a-csv", type=Path, required=True)
    ap.add_argument("--fatigue-csv", type=Path, required=True)
    ap.add_argument("--monotonic-csv", type=Path, required=True)
    ap.add_argument("--threshold-csv", type=Path, default=None)
    ap.add_argument("--manifest", type=Path, default=None)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--S-emit", dest="S_emit", type=float, default=-40.0)
    ap.add_argument("--S-cleave", dest="S_cleave", type=float, default=0.0)
    ap.add_argument("--panel-c-T", type=float, default=300.0)
    ap.add_argument("--metric-T-ref", type=float, default=300.0)
    ap.add_argument("--metric-T-hi-ref", type=float, default=900.0)
    ap.add_argument("--primary-rate", type=float, default=1e-10)
    ap.add_argument("--panel-c-cases", nargs="+", default=["FCC_like_case29", "plastic_shielded_case64_M1"])
    ap.add_argument("--metric-cases", nargs="+", default=DEFAULT_CASE_ORDER)
    args = ap.parse_args()

    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    manifest = read_manifest(args.manifest)
    disp = display_map(manifest)

    panel_a = load_panel_a(args.panel_a_csv)
    fatigue = pd.read_csv(args.fatigue_csv)
    monotonic = pd.read_csv(args.monotonic_csv)
    thresholds = load_thresholds(args.threshold_csv, fatigue, args.primary_rate, args.S_emit, args.S_cleave)

    panel_c = extract_panel_c(fatigue, args.panel_c_cases, args.panel_c_T, args.S_emit, args.S_cleave, disp)
    metrics = cross_case_metrics(monotonic, thresholds, args.metric_cases,
                                 args.S_emit, args.S_cleave, args.metric_T_ref,
                                 args.metric_T_hi_ref, disp)

    panel_a.to_csv(out / "panel_A_Kc_temperature.csv", index=False)
    panel_c.to_csv(out / "panel_C_fatigue_contrast.csv", index=False)
    metrics.to_csv(out / "panel_F_DBTT_endurance_metrics.csv", index=False)
    thresholds.to_csv(out / "thresholds_used_for_metrics.csv", index=False)
    metrics.to_csv(out / "case_metrics_summary.csv", index=False)

    settings = {
        "panel_A_source": str(args.panel_a_csv),
        "fatigue_source": str(args.fatigue_csv),
        "monotonic_source": str(args.monotonic_csv),
        "threshold_source": str(args.threshold_csv) if args.threshold_csv else None,
        "S_emit_kB": args.S_emit,
        "S_cleave_kB": args.S_cleave,
        "panel_C_T_K": args.panel_c_T,
        "metric_T_ref_K": args.metric_T_ref,
        "metric_T_hi_ref_K": args.metric_T_hi_ref,
        "primary_da_dN_criterion_m_per_cycle": args.primary_rate,
        "panel_C_cases": args.panel_c_cases,
        "metric_cases": args.metric_cases,
        "metric_definitions": {
            "T_transition_half_K": "temperature where Kc crosses half of its observed amplitude",
            "I_DBTT_amplitude": "(max(Kc)-min(Kc))/Kc(Tmin)",
            "T_inflection_K": "temperature of max abs(dKc/dT)",
            "DeltaKth_ref_MPa_sqrtm": "rate-defined DeltaK threshold at Tref",
            "I_endurance_DKth_over_Kc_ref": "DeltaKth(Tref)/Kc(Tref)",
            "T_pers_50_K": "highest T with DeltaKth(T)>=0.5*DeltaKth(Tref)",
        },
    }
    (out / "atlas_settings.json").write_text(json.dumps(settings, indent=2), encoding="utf-8")

    plot_panel_a(panel_a, out)
    plot_panel_c(panel_c, out, args.panel_c_cases)
    plot_panel_f(metrics, out)
    combined_figure(panel_a, panel_c, metrics, out, args.panel_c_cases)
    write_readme(out, args, metrics)

    print("\n=== Figure 1 first-pass cross-case metrics ===")
    cols = ["case_label", "T_transition_half_K", "I_DBTT_amplitude", "DeltaKth_ref_MPa_sqrtm",
            "I_endurance_DKth_over_Kc_ref", "T_inflection_K", "T_pers_50_K"]
    print(metrics[cols].to_string(index=False))
    print(f"\nWrote Figure 1 atlas to: {out}")


if __name__ == "__main__":
    main()
