#!/usr/bin/env python3
"""Rate-defined fatigue-threshold analysis for the two-barrier DBTT/fatigue map.

This script re-analyzes fatigue da/dN data using an operational fatigue threshold
criterion rather than an event/no-event definition.  For each material case,
emission entropy, cleavage entropy, and temperature, it estimates DeltaK_th from
log10(da/dN) vs DeltaK and then merges that result with the monotonic Kc(T)
calculation.

The primary paper-facing outputs are:
  * rate_defined_thresholds.csv
  * rate_defined_DBTT_fatigue_link.csv
  * correlation_summary.csv
  * fatigue_excess_residuals.csv
  * Kc_vs_rate_threshold_trajectories.png
  * threshold_vs_temperature_panels.png
  * fatigue_excess_residuals.png

The script never modifies the simulation CSV files.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

GROUP = ["case_label", "S_emit_kB", "S_cleave_kB", "T_K"]


def finite(v) -> bool:
    try:
        return math.isfinite(float(v))
    except Exception:
        return False


def load_map_dir(map_dir: Path) -> Tuple[pd.DataFrame, pd.DataFrame]:
    fatigue_path = map_dir / "fatigue_paris_points.csv"
    mono_path = map_dir / "monotonic_DBTT_points.csv"
    if not fatigue_path.exists():
        raise FileNotFoundError(f"missing {fatigue_path}")
    if not mono_path.exists():
        raise FileNotFoundError(f"missing {mono_path}")
    return pd.read_csv(fatigue_path), pd.read_csv(mono_path)


def effective_rate_rows(g: pd.DataFrame) -> pd.DataFrame:
    """Build a conservative rate table suitable for threshold bracketing.

    Measured points use da/dN.  Horizon-censored no-growth points use their
    reported upper bound.  Block-limited unresolved points are retained but
    excluded from interpolation unless no better information exists.
    """
    rows: List[Dict] = []
    for _, r in g.iterrows():
        rate = float("nan")
        source = "unresolved"
        if finite(r.get("da_dN_m_per_cycle", np.nan)) and float(r["da_dN_m_per_cycle"]) > 0:
            rate = float(r["da_dN_m_per_cycle"])
            source = "measured"
        elif finite(r.get("da_dN_upper_bound_m_per_cycle", np.nan)) and float(r["da_dN_upper_bound_m_per_cycle"]) > 0:
            # Treat upper bounds as inequality information.  The estimator below
            # may use them only as the lower side of a crossing bracket.
            rate = float(r["da_dN_upper_bound_m_per_cycle"])
            source = "upper_bound"
        rows.append({
            "DeltaK": float(r["DeltaK_MPa_sqrtm"]),
            "rate": rate,
            "source": source,
            "status": str(r.get("status", "")),
            "direct_lt_1_cycle": bool(r.get("direct_lt_1_cycle", False)),
        })
    out = pd.DataFrame(rows).sort_values("DeltaK").drop_duplicates("DeltaK", keep="last")
    return out


def interpolate_crossing(g: pd.DataFrame, criterion: float) -> Dict:
    tab = effective_rate_rows(g)
    valid = tab[np.isfinite(tab["rate"]) & (tab["rate"] > 0)].copy()
    result = {
        "threshold_class": "unresolved",
        "DeltaK_threshold_lower_MPa_sqrtm": np.nan,
        "DeltaK_threshold_upper_MPa_sqrtm": np.nan,
        "DeltaK_threshold_estimate_MPa_sqrtm": np.nan,
        "lower_rate_m_per_cycle": np.nan,
        "upper_rate_m_per_cycle": np.nan,
        "lower_source": "",
        "upper_source": "",
        "n_rate_points": int(len(valid)),
    }
    if valid.empty:
        return result

    below = valid[valid["rate"] < criterion]
    above = valid[valid["rate"] >= criterion]

    # Find adjacent bracketing pairs in DeltaK order.
    pair = None
    vals = valid.sort_values("DeltaK").to_dict("records")
    for a, b in zip(vals[:-1], vals[1:]):
        if a["rate"] < criterion <= b["rate"]:
            pair = (a, b)
            break

    if pair is not None:
        lo, hi = pair
        x0, x1 = float(lo["DeltaK"]), float(hi["DeltaK"])
        y0, y1 = math.log10(float(lo["rate"])), math.log10(float(hi["rate"]))
        yc = math.log10(criterion)
        if abs(y1 - y0) < 1e-14:
            x = 0.5 * (x0 + x1)
        else:
            x = x0 + (x1 - x0) * (yc - y0) / (y1 - y0)
        x = min(max(x, x0), x1)
        result.update({
            "threshold_class": "bracketed",
            "DeltaK_threshold_lower_MPa_sqrtm": x0,
            "DeltaK_threshold_upper_MPa_sqrtm": x1,
            "DeltaK_threshold_estimate_MPa_sqrtm": x,
            "lower_rate_m_per_cycle": float(lo["rate"]),
            "upper_rate_m_per_cycle": float(hi["rate"]),
            "lower_source": str(lo["source"]),
            "upper_source": str(hi["source"]),
        })
        return result

    min_dk = float(valid["DeltaK"].min())
    max_dk = float(valid["DeltaK"].max())
    low_row = valid.loc[valid["DeltaK"].idxmin()]
    high_row = valid.loc[valid["DeltaK"].idxmax()]
    if float(low_row["rate"]) >= criterion:
        result.update({
            "threshold_class": "below_grid",
            "DeltaK_threshold_upper_MPa_sqrtm": min_dk,
            "upper_rate_m_per_cycle": float(low_row["rate"]),
            "upper_source": str(low_row["source"]),
        })
    elif float(high_row["rate"]) < criterion:
        result.update({
            "threshold_class": "above_grid",
            "DeltaK_threshold_lower_MPa_sqrtm": max_dk,
            "lower_rate_m_per_cycle": float(high_row["rate"]),
            "lower_source": str(high_row["source"]),
        })
    return result


def build_thresholds(fatigue: pd.DataFrame, criteria: Iterable[float]) -> pd.DataFrame:
    rows: List[Dict] = []
    for keys, g in fatigue.groupby(GROUP, sort=False):
        metadata = {}
        for col in ["response_regime", "chi_shield", "N_sat", "scenario"]:
            if col in g.columns:
                metadata[col] = g.iloc[0][col]
        for crit in criteria:
            rec = dict(zip(GROUP, keys))
            rec.update(metadata)
            rec["da_dN_criterion_m_per_cycle"] = float(crit)
            rec.update(interpolate_crossing(g, float(crit)))
            rows.append(rec)
    return pd.DataFrame(rows)


def correlation_rows(link: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for crit, g0 in link.groupby("da_dN_criterion_m_per_cycle"):
        g = g0[(g0["threshold_class"] == "bracketed") & np.isfinite(g0["Kc_first_MPa_sqrtm"]) & np.isfinite(g0["DeltaK_threshold_estimate_MPa_sqrtm"])].copy()
        if len(g) >= 2:
            pearson = g["Kc_first_MPa_sqrtm"].corr(g["DeltaK_threshold_estimate_MPa_sqrtm"], method="pearson")
            spearman = g["Kc_first_MPa_sqrtm"].corr(g["DeltaK_threshold_estimate_MPa_sqrtm"], method="spearman")
        else:
            pearson = spearman = np.nan
        rows.append({"criterion": crit, "scope": "pooled", "T_K": np.nan, "n": len(g), "pearson_r": pearson, "spearman_rho": spearman})
        for T, gt in g.groupby("T_K"):
            if len(gt) >= 2:
                p = gt["Kc_first_MPa_sqrtm"].corr(gt["DeltaK_threshold_estimate_MPa_sqrtm"], method="pearson")
                s = gt["Kc_first_MPa_sqrtm"].corr(gt["DeltaK_threshold_estimate_MPa_sqrtm"], method="spearman")
            else:
                p = s = np.nan
            rows.append({"criterion": crit, "scope": "fixed_temperature", "T_K": T, "n": len(gt), "pearson_r": p, "spearman_rho": s})
    return pd.DataFrame(rows)


def add_common_trend_residuals(link: pd.DataFrame, primary_criterion: float, exclude_cases: List[str]) -> Tuple[pd.DataFrame, Dict]:
    g = link[np.isclose(link["da_dN_criterion_m_per_cycle"], primary_criterion)].copy()
    g = g[(g["threshold_class"] == "bracketed") & np.isfinite(g["Kc_first_MPa_sqrtm"]) & np.isfinite(g["DeltaK_threshold_estimate_MPa_sqrtm"])].copy()
    train = g[~g["case_label"].isin(exclude_cases)].copy()
    meta = {"primary_criterion": primary_criterion, "excluded_from_common_trend": exclude_cases}
    if len(train) < 2:
        g["DeltaK_common_trend_MPa_sqrtm"] = np.nan
        g["fatigue_excess_residual_MPa_sqrtm"] = np.nan
        meta.update({"slope": np.nan, "intercept": np.nan, "n_fit": len(train)})
        return g, meta
    x = train["Kc_first_MPa_sqrtm"].to_numpy(float)
    y = train["DeltaK_threshold_estimate_MPa_sqrtm"].to_numpy(float)
    slope, intercept = np.polyfit(x, y, 1)
    pred = slope * g["Kc_first_MPa_sqrtm"].to_numpy(float) + intercept
    g["DeltaK_common_trend_MPa_sqrtm"] = pred
    g["fatigue_excess_residual_MPa_sqrtm"] = g["DeltaK_threshold_estimate_MPa_sqrtm"].to_numpy(float) - pred
    meta.update({"slope": float(slope), "intercept": float(intercept), "n_fit": int(len(train))})
    return g, meta


def plot_trajectories(link: pd.DataFrame, primary: float, out: Path):
    g = link[np.isclose(link["da_dN_criterion_m_per_cycle"], primary)].copy()
    g = g[(g["threshold_class"] == "bracketed") & np.isfinite(g["Kc_first_MPa_sqrtm"]) & np.isfinite(g["DeltaK_threshold_estimate_MPa_sqrtm"])].copy()
    if g.empty:
        return
    cases = list(dict.fromkeys(g["case_label"].tolist()))
    markers = ["o", "s", "^", "D", "v", "P", "X", "*"]
    Tmin, Tmax = float(g["T_K"].min()), float(g["T_K"].max())
    norm = plt.Normalize(Tmin, Tmax)
    cmap = plt.get_cmap("viridis")
    fig, ax = plt.subplots(figsize=(8.4, 6.5))
    for i, case in enumerate(cases):
        gc = g[g["case_label"] == case].sort_values(["S_emit_kB", "S_cleave_kB", "T_K"])
        # Connect each thermodynamic scenario independently so lines have a
        # clear meaning: increasing temperature along one material/scenario.
        for _, gs in gc.groupby(["S_emit_kB", "S_cleave_kB"], sort=False):
            gs = gs.sort_values("T_K")
            ax.plot(gs["Kc_first_MPa_sqrtm"], gs["DeltaK_threshold_estimate_MPa_sqrtm"], linewidth=1.0, alpha=0.28)
        ax.scatter(gc["Kc_first_MPa_sqrtm"], gc["DeltaK_threshold_estimate_MPa_sqrtm"],
                   c=gc["T_K"], cmap=cmap, norm=norm, marker=markers[i % len(markers)],
                   s=62, edgecolors="black", linewidths=0.35, label=case)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label("Temperature (K)")
    ax.set_xlabel(r"Monotonic $K_c$ (MPa $\sqrt{m}$)")
    ax.set_ylabel(r"Rate-defined $\Delta K_{th}$ (MPa $\sqrt{m}$)")
    ax.set_title(rf"Two-barrier toughness–fatigue trajectories, $da/dN={primary:.0e}$ m/cycle")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=7.5, loc="best")
    fig.tight_layout()
    fig.savefig(out / "Kc_vs_rate_threshold_trajectories.png", dpi=260)
    plt.close(fig)


def plot_temperature_panels(link: pd.DataFrame, primary: float, out: Path):
    g = link[np.isclose(link["da_dN_criterion_m_per_cycle"], primary)].copy()
    g = g[g["threshold_class"] == "bracketed"]
    cases = list(dict.fromkeys(g["case_label"].tolist()))
    if not cases:
        return
    ncols = 3
    nrows = int(math.ceil(len(cases) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(15, 4.2 * nrows), squeeze=False)
    for ax, case in zip(axes.flat, cases):
        gc = g[g["case_label"] == case]
        for (Se, Sc), gs in gc.groupby(["S_emit_kB", "S_cleave_kB"], sort=False):
            gs = gs.sort_values("T_K")
            ax.plot(gs["T_K"], gs["DeltaK_threshold_estimate_MPa_sqrtm"], marker="o", label=f"Se={Se:g}, Sc={Sc:g}")
        ax.set_title(case)
        ax.set_xlabel("Temperature (K)")
        ax.set_ylabel(r"$\Delta K_{th}$ (MPa $\sqrt{m}$)")
        ax.grid(True, alpha=0.25)
        ax.legend(fontsize=6.5)
    for ax in axes.flat[len(cases):]:
        ax.axis("off")
    fig.suptitle(rf"Rate-defined fatigue thresholds, $da/dN={primary:.0e}$ m/cycle", y=0.995)
    fig.tight_layout()
    fig.savefig(out / "threshold_vs_temperature_panels.png", dpi=240)
    plt.close(fig)


def plot_residuals(resid: pd.DataFrame, out: Path):
    if resid.empty:
        return
    fig, ax = plt.subplots(figsize=(9.0, 5.8))
    cases = list(dict.fromkeys(resid["case_label"].tolist()))
    xmap = {c: i for i, c in enumerate(cases)}
    rng = np.random.default_rng(7)
    for case, g in resid.groupby("case_label", sort=False):
        x = np.full(len(g), xmap[case], dtype=float) + rng.normal(0, 0.045, len(g))
        sc = ax.scatter(x, g["fatigue_excess_residual_MPa_sqrtm"], c=g["T_K"], cmap="viridis", s=44, alpha=0.8)
    ax.axhline(0.0, linewidth=1.2)
    ax.set_xticks(range(len(cases)))
    ax.set_xticklabels(cases, rotation=28, ha="right")
    ax.set_ylabel(r"Fatigue excess $\delta_{fatigue}$ (MPa $\sqrt{m}$)")
    ax.set_title("Departure from the cleavage-dominated toughness–fatigue trend")
    ax.grid(True, axis="y", alpha=0.25)
    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label("Temperature (K)")
    fig.tight_layout()
    fig.savefig(out / "fatigue_excess_residuals.png", dpi=250)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--map-dir", default="", help="Directory containing fatigue_paris_points.csv and monotonic_DBTT_points.csv")
    ap.add_argument("--fatigue-csv", default="", help="Optional explicit fatigue points CSV; overrides --map-dir for fatigue input")
    ap.add_argument("--monotonic-csv", default="", help="Optional explicit monotonic points CSV; overrides --map-dir for monotonic input")
    ap.add_argument("--out", required=True)
    ap.add_argument("--criteria", nargs="+", type=float, default=[1e-12, 1e-10])
    ap.add_argument("--primary-criterion", type=float, default=1e-10)
    ap.add_argument("--exclude-common-trend-cases", nargs="*", default=["plastic_shielded_case64_M1"],
                    help="Cases excluded when fitting the cleavage-dominated common trend.")
    args = ap.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    if args.fatigue_csv and args.monotonic_csv:
        fatigue = pd.read_csv(args.fatigue_csv)
        mono = pd.read_csv(args.monotonic_csv)
    elif args.map_dir:
        fatigue, mono = load_map_dir(Path(args.map_dir))
    else:
        raise SystemExit("provide --map-dir or both --fatigue-csv and --monotonic-csv")

    thr = build_thresholds(fatigue, args.criteria)
    link = thr.merge(mono, on=GROUP, how="left", suffixes=("", "_mono"))
    corr = correlation_rows(link)
    resid, fit_meta = add_common_trend_residuals(link, args.primary_criterion, args.exclude_common_trend_cases)

    thr.to_csv(out / "rate_defined_thresholds.csv", index=False)
    link.to_csv(out / "rate_defined_DBTT_fatigue_link.csv", index=False)
    corr.to_csv(out / "correlation_summary.csv", index=False)
    resid.to_csv(out / "fatigue_excess_residuals.csv", index=False)
    with open(out / "common_trend_fit.json", "w") as f:
        json.dump(fit_meta, f, indent=2)

    plot_trajectories(link, args.primary_criterion, out)
    plot_temperature_panels(link, args.primary_criterion, out)
    plot_residuals(resid, out)

    print(f"wrote {out / 'rate_defined_thresholds.csv'}")
    print(f"wrote {out / 'rate_defined_DBTT_fatigue_link.csv'}")
    print(f"wrote {out / 'correlation_summary.csv'}")
    print(f"wrote {out / 'fatigue_excess_residuals.csv'}")
    print(corr.to_string(index=False))


if __name__ == "__main__":
    main()
