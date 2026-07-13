#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def _k_from_dir(name: str) -> float:
    m = re.fullmatch(r"K(.+)", name)
    if not m:
        return float("nan")
    return float(m.group(1).replace("p", "."))


def load_points(root: Path, R: float, cycles_max: float) -> pd.DataFrame:
    rows = []
    for case_dir in sorted(p for p in root.iterdir() if p.is_dir()):
        for kdir in sorted(p for p in case_dir.iterdir() if p.is_dir() and p.name.startswith("K")):
            Kmax = _k_from_dir(kdir.name)
            for tdir in sorted(p for p in kdir.iterdir() if p.is_dir() and p.name.startswith("T") and p.name.endswith("K")):
                try:
                    T = float(tdir.name[1:-1])
                except Exception:
                    continue
                hpath = tdir / "fatigue_v1_history.csv"
                spath = tdir / "summary.json"
                if not hpath.exists() or not spath.exists():
                    continue
                hist = pd.read_csv(hpath)
                with spath.open() as fp:
                    summ = json.load(fp)
                last = hist.iloc[-1]
                cycles_total = float(last.get("cycles_total", summ.get("cycles_total", float("nan"))))
                n_adv = int(summ.get("n_adv", last.get("n_adv", 0)))
                a_adv_m = float(summ.get("a_adv_m", last.get("a_adv_m", 0.0)))
                fired = hist[hist.get("n_fire", pd.Series(index=hist.index, data=0)).astype(float) > 0]
                cycles_first = float(fired.iloc[0]["cycles_total"]) if len(fired) else float("nan")
                measured = n_adv > 0 and cycles_total > 0 and a_adv_m > 0
                da_dN = a_adv_m / cycles_total if measured else float("nan")
                # V1 default advance quantum in the production calibration is 20 um.
                ub = 20e-6 / cycles_total if (not measured and cycles_total > 0) else float("nan")
                if measured:
                    status = "measured"
                elif math.isfinite(cycles_max) and cycles_total >= 0.999 * cycles_max:
                    status = "censored_cycle_horizon"
                else:
                    status = "censored_block_limited"
                rows.append({
                    "case_label": case_dir.name,
                    "Kmax_MPa_sqrtm": Kmax,
                    "DeltaK_MPa_sqrtm": (1.0 - R) * Kmax,
                    "T_K": T,
                    "cycles_total": cycles_total,
                    "cycles_to_first_fire": cycles_first,
                    "n_adv": n_adv,
                    "a_adv_m": a_adv_m,
                    "da_dN_m_per_cycle": da_dN,
                    "da_dN_upper_bound_m_per_cycle": ub,
                    "direct_lt_1_cycle": math.isfinite(cycles_first) and cycles_first < 1.0,
                    "status": status,
                    "B_final": float(summ.get("B", last.get("B", float("nan")))),
                    "N_em_final": float(summ.get("N_em", last.get("N_em", float("nan")))),
                    "sigma_back_Pa": float(summ.get("sigma_back_Pa", last.get("sigma_back", float("nan")))),
                    "dG_emb_eV": float(summ.get("dG_emb_eV", last.get("dG_emb_eV", float("nan")))),
                    "history_csv": str(hpath),
                })
    return pd.DataFrame(rows)


def make_case_plots(points: pd.DataFrame, outdir: Path) -> None:
    for case, g in points.groupby("case_label", sort=False):
        fig, ax = plt.subplots(figsize=(7.3, 5.2))
        for T, t in g.groupby("T_K"):
            t = t.sort_values("DeltaK_MPa_sqrtm")
            m = t[t["status"] == "measured"]
            if not m.empty:
                ax.plot(m["DeltaK_MPa_sqrtm"], m["da_dN_m_per_cycle"], marker="o", label=f"{T:g} K")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"$\Delta K$ (MPa $\sqrt{m}$)")
        ax.set_ylabel(r"V1 $da/dN$ proxy (m/cycle)")
        ax.set_title(case)
        ax.grid(True, which="both", alpha=0.3)
        ax.legend()
        fig.tight_layout()
        fig.savefig(outdir / f"{case}_da_dN_vs_DeltaK_T.png", dpi=220)
        plt.close(fig)


def summarize_thresholds(points: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (case, T), g in points.groupby(["case_label", "T_K"]):
        m = g[g["status"] == "measured"].sort_values("DeltaK_MPa_sqrtm")
        c = g[g["status"] != "measured"].sort_values("DeltaK_MPa_sqrtm")
        direct = g[g["direct_lt_1_cycle"]].sort_values("DeltaK_MPa_sqrtm")
        rows.append({
            "case_label": case,
            "T_K": T,
            "DeltaK_first_measured_MPa_sqrtm": float(m.iloc[0]["DeltaK_MPa_sqrtm"]) if len(m) else float("nan"),
            "DeltaK_highest_no_growth_MPa_sqrtm": float(c.iloc[-1]["DeltaK_MPa_sqrtm"]) if len(c) else float("nan"),
            "DeltaK_first_direct_lt_1_cycle_MPa_sqrtm": float(direct.iloc[0]["DeltaK_MPa_sqrtm"]) if len(direct) else float("nan"),
            "n_measured": len(m),
            "n_censored": len(c),
        })
    return pd.DataFrame(rows)


def plot_threshold_map(thr: pd.DataFrame, outdir: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.8, 5.4))
    for case, g in thr.groupby("case_label", sort=False):
        g = g.sort_values("T_K")
        ax.plot(g["T_K"], g["DeltaK_first_measured_MPa_sqrtm"], marker="o", label=case)
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel(r"First measured-growth $\Delta K$ (MPa $\sqrt{m}$)")
    ax.set_title("V1 apparent crack-growth threshold versus temperature")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(outdir / "V1_threshold_vs_temperature.png", dpi=220)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--R", type=float, default=0.1)
    ap.add_argument("--cycles-max", type=float, default=float("nan"))
    args = ap.parse_args()

    root = Path(args.root)
    points = load_points(root, args.R, args.cycles_max)
    points.to_csv(root / "v1_temperature_paris_points.csv", index=False)
    thr = summarize_thresholds(points) if not points.empty else pd.DataFrame()
    thr.to_csv(root / "v1_temperature_threshold_summary.csv", index=False)
    if not points.empty:
        make_case_plots(points, root)
        plot_threshold_map(thr, root)
    print(f"wrote {root / 'v1_temperature_paris_points.csv'}")
    print(f"wrote {root / 'v1_temperature_threshold_summary.csv'}")


if __name__ == "__main__":
    main()
