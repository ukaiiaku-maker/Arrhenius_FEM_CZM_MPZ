#!/usr/bin/env python3
"""Plot seeded R-curves, extension-aligned means, and activation-volume/Gumbel diagnostics."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import re

import numpy as np
import pandas as pd

CLASSES_DEFAULT = ["ceramic", "peak", "weakT", "DBTT"]
EULER_GAMMA = 0.5772156649015329


def parse_list(text: str, cast=str):
    return [cast(x) for x in text.replace(",", " ").split() if x]


def replicate_metadata(rep_dir: Path) -> tuple[int | None, int | None]:
    m = re.search(r"replicate_(\d+)_seed(-?\d+)", rep_dir.name)
    if not m:
        return None, None
    return int(m.group(1)), int(m.group(2))


def load_summary(case_dir: Path) -> dict:
    sf = case_dir / "summary.json"
    if not sf.exists():
        return {}
    try:
        d = json.loads(sf.read_text())
        return d[0] if isinstance(d, list) and d else d
    except Exception:
        return {}


def load_curve(case_dir: Path, target_ext_um: float) -> tuple[pd.DataFrame, dict]:
    f = case_dir / "R_curve_event_sampled.csv"
    if not f.exists():
        return pd.DataFrame(), {}
    c = pd.read_csv(f)
    required = {"crack_extension_um", "KJ_MPa_sqrt_m"}
    if not required.issubset(c.columns):
        return pd.DataFrame(), {}
    c = c.replace([np.inf, -np.inf], np.nan).dropna(
        subset=["crack_extension_um", "KJ_MPa_sqrt_m"]
    )
    c = c[(c.crack_extension_um >= 0.0) & (c.crack_extension_um <= target_ext_um + 5.0)]
    if c.empty:
        return c, {}
    c = c.sort_values("crack_extension_um").drop_duplicates(
        "crack_extension_um", keep="last"
    )

    summary = load_summary(case_dir)
    try:
        Kinit = float(summary.get("Kc_first_MPa_sqrt_m", np.nan))
    except Exception:
        Kinit = np.nan
    if np.isfinite(Kinit):
        if c.crack_extension_um.min() > 1e-9:
            c = pd.concat([
                pd.DataFrame({"crack_extension_um": [0.0], "KJ_MPa_sqrt_m": [Kinit]}),
                c,
            ], ignore_index=True)
        else:
            c.loc[c.crack_extension_um.idxmin(), "KJ_MPa_sqrt_m"] = Kinit
    return c.sort_values("crack_extension_um").reset_index(drop=True), summary


def interpolate_column(curve: pd.DataFrame, grid: np.ndarray, column: str) -> np.ndarray:
    out = np.full_like(grid, np.nan, dtype=float)
    if column not in curve.columns:
        return out
    x = pd.to_numeric(curve.crack_extension_um, errors="coerce").to_numpy(float)
    y = pd.to_numeric(curve[column], errors="coerce").to_numpy(float)
    valid = np.isfinite(x) & np.isfinite(y)
    x, y = x[valid], y[valid]
    if len(x) < 2:
        return out
    order = np.argsort(x)
    x, y = x[order], y[order]
    ux, idx = np.unique(x, return_index=True)
    y = y[idx]
    if len(ux) < 2:
        return out
    inside = (grid >= ux.min()) & (grid <= ux.max())
    out[inside] = np.interp(grid[inside], ux, y)
    return out


def column_stats(matrix: np.ndarray) -> dict[str, np.ndarray]:
    n = np.sum(np.isfinite(matrix), axis=0)
    mean = np.full(matrix.shape[1], np.nan)
    median = np.full(matrix.shape[1], np.nan)
    p10 = np.full(matrix.shape[1], np.nan)
    p90 = np.full(matrix.shape[1], np.nan)
    std = np.full(matrix.shape[1], np.nan)
    sem = np.full(matrix.shape[1], np.nan)
    for j in range(matrix.shape[1]):
        vals = matrix[:, j]
        vals = vals[np.isfinite(vals)]
        if len(vals):
            mean[j] = np.mean(vals)
            median[j] = np.median(vals)
            p10[j], p90[j] = np.percentile(vals, [10, 90])
        if len(vals) >= 2:
            std[j] = np.std(vals, ddof=1)
            sem[j] = std[j] / np.sqrt(len(vals))
    return {"n": n, "mean": mean, "median": median, "p10": p10,
            "p90": p90, "std": std, "sem": sem}


def interval_mean(values: np.ndarray, grid: np.ndarray, lo: float, hi: float) -> float:
    m = (grid >= lo) & (grid <= hi) & np.isfinite(values)
    return float(np.mean(values[m])) if np.any(m) else np.nan


def discover_cases(root: Path, klass: str, T: int):
    out = []
    class_dir = root / klass
    if not class_dir.exists():
        return out
    for rep_dir in sorted(class_dir.glob("replicate_*_seed*")):
        matches = sorted(rep_dir.glob(f"T{T}_th*"))
        if matches:
            rep, seed = replicate_metadata(rep_dir)
            out.append((matches[0], rep, seed))
    return out


def min_gumbel_cdf(x: np.ndarray, mu: float, beta: float) -> np.ndarray:
    z = np.clip((np.asarray(x, float) - mu) / beta, -700.0, 700.0)
    return 1.0 - np.exp(-np.exp(z))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--root", required=True)
    ap.add_argument("--classes", default=" ".join(CLASSES_DEFAULT))
    ap.add_argument("--temperature", type=int, default=500)
    ap.add_argument("--target-ext-um", type=float, default=1000.0)
    ap.add_argument("--grid-step-um", type=float, default=5.0)
    ap.add_argument("--burgers-vector-m", type=float, default=2.74e-10)
    args = ap.parse_args()

    root = Path(args.root)
    classes = parse_list(args.classes)
    grid = np.arange(0.0, args.target_ext_um + 0.5 * args.grid_step_um,
                     args.grid_step_um)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    replicate_rows = []
    gumbel_rows = []
    class_stats = {}
    class_curves = {}

    for klass in classes:
        curves = []
        K_mats = []
        V_mats = []
        betaK_mats = []
        pred_sd_mats = []
        initiation = []

        for case_dir, rep, seed in discover_cases(root, klass, args.temperature):
            curve, summary = load_curve(case_dir, args.target_ext_um)
            if curve.empty:
                continue
            K_i = interpolate_column(curve, grid, "KJ_MPa_sqrt_m")
            V_i = interpolate_column(curve, grid, "vstar_cleave_b3")
            beta_i = interpolate_column(curve, grid, "gumbel_beta_K_MPa_sqrt_m")
            pred_sd_i = interpolate_column(curve, grid, "gumbel_sd_K_MPa_sqrt_m")
            curves.append((curve, rep, seed))
            K_mats.append(K_i)
            V_mats.append(V_i)
            betaK_mats.append(beta_i)
            pred_sd_mats.append(pred_sd_i)

            try:
                Kinit = float(summary.get("Kc_first_MPa_sqrt_m", np.nan))
            except Exception:
                Kinit = np.nan
            try:
                final_ext = (float(summary.get("a_final_mm")) - 0.5) * 1000.0
            except Exception:
                final_ext = float(curve.crack_extension_um.max())

            event_curve = curve[curve.crack_extension_um > 0.0]
            first = event_curve.iloc[0] if len(event_curve) else curve.iloc[0]
            beta_init = float(first.get("gumbel_beta_K_MPa_sqrt_m", np.nan))
            v_init = float(first.get("vstar_cleave_b3", np.nan))
            initiation.append((Kinit, beta_init, v_init, rep, seed))

            late_mean = interval_mean(K_i, grid, 500.0, args.target_ext_um)
            replicate_rows.append({
                "class": klass,
                "replicate": rep,
                "solver_seed": seed,
                "T_K": args.temperature,
                "K_init_MPa_sqrt_m": Kinit,
                "first_event_vstar_b3": v_init,
                "first_event_gumbel_beta_K_MPa_sqrt_m": beta_init,
                "first_event_gumbel_sd_K_MPa_sqrt_m": (
                    math.pi / math.sqrt(6.0) * beta_init
                    if np.isfinite(beta_init) else np.nan
                ),
                "final_crack_extension_um": final_ext,
                "K_mean_200_1000um_MPa_sqrt_m": interval_mean(
                    K_i, grid, 200.0, args.target_ext_um
                ),
                "K_mean_500_1000um_MPa_sqrt_m": late_mean,
                "delta_K_late_minus_init_MPa_sqrt_m": (
                    late_mean - Kinit
                    if np.isfinite(late_mean) and np.isfinite(Kinit) else np.nan
                ),
                "n_event_points": len(event_curve),
                "case_dir": str(case_dir.relative_to(root)),
            })

        class_curves[klass] = curves
        if not K_mats:
            continue

        Kmat = np.vstack(K_mats)
        Vmat = np.vstack(V_mats)
        Bmat = np.vstack(betaK_mats)
        Pmat = np.vstack(pred_sd_mats)
        Ks, Vs, Bs, Ps = map(column_stats, [Kmat, Vmat, Bmat, Pmat])
        stats = pd.DataFrame({
            "class": klass,
            "T_K": args.temperature,
            "crack_extension_um": grid,
            "n_replicates": Ks["n"],
            "K_mean_MPa_sqrt_m": Ks["mean"],
            "K_std_MPa_sqrt_m": Ks["std"],
            "K_sem_MPa_sqrt_m": Ks["sem"],
            "K_median_MPa_sqrt_m": Ks["median"],
            "K_p10_MPa_sqrt_m": Ks["p10"],
            "K_p90_MPa_sqrt_m": Ks["p90"],
            "vstar_mean_b3": Vs["mean"],
            "vstar_std_b3": Vs["std"],
            "gumbel_beta_K_mean_MPa_sqrt_m": Bs["mean"],
            "gumbel_predicted_local_sd_K_mean_MPa_sqrt_m": Ps["mean"],
        })
        class_stats[klass] = stats
        (root / klass).mkdir(parents=True, exist_ok=True)
        stats.to_csv(root / klass / "R_curve_ensemble_and_gumbel_statistics.csv", index=False)

        # Per-class R-curves plus empirical and activation-volume-predicted scatter.
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.8, 8.0), sharex=True)
        for curve, rep, seed in curves:
            ax1.plot(curve.crack_extension_um, curve.KJ_MPa_sqrt_m,
                     linewidth=0.9, alpha=0.50,
                     label=f"replicate {rep}, seed {seed}")
        g = stats[stats.n_replicates > 0]
        ax1.plot(g.crack_extension_um, g.K_mean_MPa_sqrt_m,
                 linewidth=2.5, label="extension-aligned mean")
        band = g[g.n_replicates >= 2]
        if not band.empty:
            ax1.fill_between(
                band.crack_extension_um,
                band.K_mean_MPa_sqrt_m - band.K_std_MPa_sqrt_m,
                band.K_mean_MPa_sqrt_m + band.K_std_MPa_sqrt_m,
                alpha=0.18,
                label=r"empirical mean $\pm$ 1 SD",
            )
        ax1.set_ylabel(r"$K_J$ (MPa$\sqrt{m}$)")
        ax1.set_title(f"{klass}, {args.temperature} K: seeded first-passage realizations")
        ax1.grid(alpha=0.25)
        ax1.legend(frameon=False, fontsize=7, ncol=2)

        ax2.plot(g.crack_extension_um, g.K_std_MPa_sqrt_m,
                 linewidth=1.8, label="empirical SD across five runs")
        ax2.plot(g.crack_extension_um,
                 g.gumbel_predicted_local_sd_K_mean_MPa_sqrt_m,
                 linewidth=1.8,
                 label=r"local Gumbel SD from $V^*=-\partial\Delta G/\partial\sigma$")
        ax2.set_xlim(0.0, args.target_ext_um)
        ax2.set_xlabel(r"Projected crack extension $\Delta a_x$ ($\mu$m)")
        ax2.set_ylabel(r"Scatter scale in $K_J$ (MPa$\sqrt{m}$)")
        ax2.grid(alpha=0.25)
        ax2.legend(frameon=False, fontsize=8)
        fig.tight_layout()
        fig.savefig(root / klass / "R_curve_replicates_mean_and_gumbel_scatter.png", dpi=240)
        plt.close(fig)

        # Initiation-level local minimum-Gumbel audit.
        Kvals = np.array([x[0] for x in initiation], float)
        betas = np.array([x[1] for x in initiation], float)
        vvals = np.array([x[2] for x in initiation], float)
        goodK = Kvals[np.isfinite(Kvals)]
        goodB = betas[np.isfinite(betas) & (betas > 0)]
        beta_theory = float(np.median(goodB)) if len(goodB) else np.nan
        emp_mean = float(np.mean(goodK)) if len(goodK) else np.nan
        emp_std = float(np.std(goodK, ddof=1)) if len(goodK) >= 2 else np.nan
        mu_from_mean = (
            emp_mean + EULER_GAMMA * beta_theory
            if np.isfinite(emp_mean) and np.isfinite(beta_theory) else np.nan
        )
        gumbel_rows.append({
            "class": klass,
            "T_K": args.temperature,
            "n_replicates": len(goodK),
            "K_init_empirical_mean_MPa_sqrt_m": emp_mean,
            "K_init_empirical_std_MPa_sqrt_m": emp_std,
            "first_event_vstar_median_b3": float(np.nanmedian(vvals)) if np.isfinite(vvals).any() else np.nan,
            "gumbel_beta_K_theory_median_MPa_sqrt_m": beta_theory,
            "gumbel_std_K_theory_MPa_sqrt_m": (
                math.pi / math.sqrt(6.0) * beta_theory
                if np.isfinite(beta_theory) else np.nan
            ),
            "minimum_gumbel_location_mu_from_empirical_mean_MPa_sqrt_m": mu_from_mean,
            "interpretation": "local-linearized minimum-type Gumbel; n=5 is diagnostic, not a distribution fit",
        })

    rep_df = pd.DataFrame(replicate_rows)
    if not rep_df.empty:
        rep_df = rep_df.sort_values(["class", "replicate", "solver_seed"], na_position="last")
    rep_df.to_csv(root / "four_class_500K_replicate_summary.csv", index=False)

    gumbel_df = pd.DataFrame(gumbel_rows)
    gumbel_df.to_csv(root / "four_class_500K_activation_volume_gumbel_summary.csv", index=False)

    if class_stats:
        all_stats = pd.concat(class_stats.values(), ignore_index=True)
        all_stats.to_csv(
            root / "four_class_500K_R_curve_ensemble_and_gumbel_statistics.csv",
            index=False,
        )

        fig, axes = plt.subplots(2, 2, figsize=(11.2, 8.2), sharex=True)
        for ax, klass in zip(axes.ravel(), classes):
            for curve, rep, seed in class_curves.get(klass, []):
                ax.plot(curve.crack_extension_um, curve.KJ_MPa_sqrt_m,
                        linewidth=0.75, alpha=0.38)
            stats = class_stats.get(klass)
            if stats is not None:
                g = stats[stats.n_replicates > 0]
                ax.plot(g.crack_extension_um, g.K_mean_MPa_sqrt_m,
                        linewidth=2.4, label="mean")
                band = g[g.n_replicates >= 2]
                if not band.empty:
                    ax.fill_between(
                        band.crack_extension_um,
                        band.K_mean_MPa_sqrt_m - band.K_std_MPa_sqrt_m,
                        band.K_mean_MPa_sqrt_m + band.K_std_MPa_sqrt_m,
                        alpha=0.16,
                    )
            ax.set_title(klass)
            ax.set_xlim(0.0, args.target_ext_um)
            ax.grid(alpha=0.25)
            ax.set_ylabel(r"$K_J$ (MPa$\sqrt{m}$)")
        for ax in axes[-1, :]:
            ax.set_xlabel(r"Projected crack extension $\Delta a_x$ ($\mu$m)")
        axes[0, 0].legend(frameon=False)
        fig.suptitle(f"Seeded FEM/CZM R-curve-like responses at {args.temperature} K")
        fig.tight_layout(rect=(0, 0, 1, 0.965))
        fig.savefig(root / "four_class_500K_individual_and_mean_R_curves.png", dpi=240)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(8.0, 5.4))
        for klass in classes:
            stats = class_stats.get(klass)
            if stats is None:
                continue
            g = stats[stats.n_replicates > 0]
            ax.plot(g.crack_extension_um, g.K_mean_MPa_sqrt_m,
                    linewidth=2.0, label=klass)
        ax.set_xlim(0.0, args.target_ext_um)
        ax.set_xlabel(r"Projected crack extension $\Delta a_x$ ($\mu$m)")
        ax.set_ylabel(r"Mean $K_J$ (MPa$\sqrt{m}$)")
        ax.set_title(f"Mean R-curve-like response by class, {args.temperature} K")
        ax.grid(alpha=0.25)
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(root / "four_class_500K_mean_R_curve_comparison.png", dpi=240)
        plt.close(fig)

        fig, axes = plt.subplots(2, 2, figsize=(11.0, 8.0), sharex=True)
        for ax, klass in zip(axes.ravel(), classes):
            stats = class_stats.get(klass)
            if stats is None:
                continue
            g = stats[stats.n_replicates > 0]
            ax.plot(g.crack_extension_um, g.vstar_mean_b3, linewidth=2.0)
            ax.set_title(klass)
            ax.set_ylabel(r"Mean effective activation volume $V^*/b^3$")
            ax.grid(alpha=0.25)
        for ax in axes[-1, :]:
            ax.set_xlabel(r"Projected crack extension $\Delta a_x$ ($\mu$m)")
        fig.suptitle(
            r"Activation volume along the mean crack-growth history, $V^*=-\partial\Delta G/\partial\sigma$"
        )
        fig.tight_layout(rect=(0, 0, 1, 0.965))
        fig.savefig(root / "four_class_500K_activation_volume_vs_extension.png", dpi=240)
        plt.close(fig)

    print(f"WROTE {root / 'four_class_500K_replicate_summary.csv'}")
    print(f"WROTE {root / 'four_class_500K_activation_volume_gumbel_summary.csv'}")
    if class_stats:
        print(f"WROTE {root / 'four_class_500K_R_curve_ensemble_and_gumbel_statistics.csv'}")
        print(f"WROTE {root / 'four_class_500K_individual_and_mean_R_curves.png'}")
        print(f"WROTE {root / 'four_class_500K_mean_R_curve_comparison.png'}")
        print(f"WROTE {root / 'four_class_500K_activation_volume_vs_extension.png'}")
    else:
        print("WARNING: no completed replicate R-curve files were found")


if __name__ == "__main__":
    main()
