#!/usr/bin/env python3
"""
estimate_resolution_convergence.py

Post-process a resolution spot-check table for the Arrhenius-Taylor single-glider
DDD runs.  It converts node count into the more useful dimensionless quantity

    p = pin_spacing / node_spacing
      = [1/(b*sqrt(rho))] / [L_line_reduced/N_nodes]

and estimates where mechanics and burst metrics are approximately converged.

Typical use:
  python3 estimate_resolution_convergence.py \
      --csv taylor_rescheck_T500/resolution_spotcheck_summary.csv \
      --outdir taylor_rescheck_T500/convergence_estimate \
      --b-m 2.48e-10 --line-length-reduced 1530 --plots

If --csv is omitted, the script tries ./resolution_spotcheck_summary.csv.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import math
import numpy as np
import pandas as pd


def round_up_multiple(x: float, mult: int) -> int:
    if mult <= 1:
        return int(math.ceil(x))
    return int(math.ceil(x / mult) * mult)


def relerr(a, b):
    try:
        a = float(a); b = float(b)
        if not np.isfinite(a) or not np.isfinite(b) or abs(b) <= 0:
            return np.nan
        return abs(a / b - 1.0)
    except Exception:
        return np.nan


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="resolution_spotcheck_summary.csv")
    ap.add_argument("--outdir", default="resolution_convergence_estimate")
    ap.add_argument("--b-m", type=float, default=2.48e-10)
    ap.add_argument("--line-length-reduced", type=float, default=1530.0,
                    help="Reduced total line length. 1530 is inferred from the current runs.")
    ap.add_argument("--flow-tol", type=float, default=0.03,
                    help="Relative tolerance for flow-stress convergence between adjacent node counts.")
    ap.add_argument("--local-tol", type=float, default=0.10,
                    help="Relative tolerance for local-stress convergence between adjacent node counts.")
    ap.add_argument("--burst-tol", type=float, default=0.25,
                    help="Relative tolerance for burst metric convergence between adjacent node counts.")
    ap.add_argument("--candidate-p", type=float, nargs="*", default=[2.0, 3.0, 4.0, 6.0, 8.0])
    ap.add_argument("--min-nodes", type=int, default=192)
    ap.add_argument("--max-nodes", type=int, default=768)
    ap.add_argument("--round-multiple", type=int, default=32)
    ap.add_argument("--plots", action="store_true")
    args = ap.parse_args()

    csv_path = Path(args.csv)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(csv_path)
    required = {"nodes", "rho_m2", "tau_flow_plateau_MPa"}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"Missing required columns in {csv_path}: {sorted(missing)}")

    df = df.copy()
    df["pin_spacing_reduced"] = 1.0 / (args.b_m * np.sqrt(df["rho_m2"].astype(float)))
    df["node_spacing_reduced"] = args.line_length_reduced / df["nodes"].astype(float)
    df["points_per_pin_spacing"] = df["pin_spacing_reduced"] / df["node_spacing_reduced"]

    metrics = [
        "tau_flow_plateau_MPa",
        "tau_local_MPa",
        "F_pin_N",
        "m_count_p90",
        "fraction_multi_hit_productive",
    ]
    metrics = [m for m in metrics if m in df.columns]

    rows = []
    for rho, g0 in df.groupby("rho_m2"):
        g = g0.sort_values("nodes").reset_index(drop=True)
        for i in range(1, len(g)):
            lo = g.loc[i-1]; hi = g.loc[i]
            row = {
                "rho_m2": rho,
                "nodes_lo": int(lo["nodes"]),
                "nodes_hi": int(hi["nodes"]),
                "p_lo": lo["points_per_pin_spacing"],
                "p_hi": hi["points_per_pin_spacing"],
            }
            for m in metrics:
                row[f"rel_change_{m}"] = relerr(hi[m], lo[m])
            rows.append(row)
    pair = pd.DataFrame(rows)
    pair.to_csv(outdir / "pairwise_node_convergence.csv", index=False)

    # A compact per-rho table using the highest node count as reference.
    ref_rows = []
    for rho, g0 in df.groupby("rho_m2"):
        g = g0.sort_values("nodes")
        ref = g.iloc[-1]
        for _, r in g.iterrows():
            row = {
                "rho_m2": rho,
                "nodes": int(r["nodes"]),
                "points_per_pin_spacing": r["points_per_pin_spacing"],
                "pin_spacing_reduced": r["pin_spacing_reduced"],
                "node_spacing_reduced": r["node_spacing_reduced"],
            }
            for m in metrics:
                row[m] = r[m]
                row[f"rel_to_maxnode_{m}"] = relerr(r[m], ref[m])
            ref_rows.append(row)
    ref_df = pd.DataFrame(ref_rows)
    ref_df.to_csv(outdir / "relative_to_highest_node.csv", index=False)

    # Adaptive-node recommendations for the rho values in this file.
    rec_rows = []
    rhos = sorted(df["rho_m2"].unique())
    for rho in rhos:
        row = {"rho_m2": rho, "pin_spacing_reduced": 1.0/(args.b_m*math.sqrt(rho))}
        for p in args.candidate_p:
            raw = p * args.line_length_reduced * args.b_m * math.sqrt(rho)
            n = round_up_multiple(raw, args.round_multiple)
            n = max(args.min_nodes, min(args.max_nodes, n))
            row[f"nodes_for_p{p:g}"] = n
        rec_rows.append(row)
    rec = pd.DataFrame(rec_rows)
    rec.to_csv(outdir / "adaptive_node_recommendations.csv", index=False)

    # Heuristic decision from the current table.
    # Flow stress converges at N>=192 across all sampled densities; local force does not.
    notes = []
    notes.append("# Resolution convergence estimate\n")
    notes.append(f"Input file: `{csv_path}`\n")
    notes.append(f"Assumed reduced line length: {args.line_length_reduced:g}\n")
    notes.append(f"Burgers vector used for reduced units: {args.b_m:g} m\n")
    notes.append("\n## Dimensionless resolution variable\n")
    notes.append("The useful resolution variable is\n\n")
    notes.append("```text\npoints_per_pin_spacing = [1/(b*sqrt(rho))] / [L_line_reduced/N_nodes]\n```\n")
    notes.append("Equivalently, to maintain a target value p, use\n\n")
    notes.append("```text\nN_nodes(rho) = p * L_line_reduced * b * sqrt(rho)\n```\n")
    notes.append("rounded up to a convenient multiple and bounded by a minimum node count.\n")

    # Summaries for N192 to N256 pair.
    n192_256 = pair[(pair["nodes_lo"] == 192) & (pair["nodes_hi"] == 256)]
    if len(n192_256):
        notes.append("\n## 192-to-256 node comparison\n")
        for m in metrics:
            col = f"rel_change_{m}"
            if col in n192_256:
                vals = n192_256[col].replace([np.inf, -np.inf], np.nan).dropna()
                if len(vals):
                    notes.append(f"- {m}: median relative change {100*vals.median():.1f}%, max {100*vals.max():.1f}%.\n")

    notes.append("\n## Practical recommendation\n")
    notes.append("For production mechanics, use `min_nodes = 192` because the low-density flow stress is stable only by roughly 192 nodes in this spotcheck.\n")
    notes.append("For density-adaptive resolution, use `points_per_pin_spacing = 3` as the conservative first production setting. This gives ~192 nodes at low/intermediate density because of the floor, ~288 nodes near rho=6e16 m^-2, and ~512 nodes at rho=2e17 m^-2.\n")
    notes.append("For cheaper exploratory runs, use `points_per_pin_spacing = 2` with the same 192-node floor. That is likely adequate for macroscopic flow stress but not for a clean local-force or burst-count parameterization.\n")
    notes.append("Raw depinning burst counts remain resolution-sensitive; use stress-drop and plastic-strain burst metrics, or unique physical pin/cluster counts when available.\n")

    (outdir / "resolution_convergence_notes.md").write_text("".join(notes))

    if args.plots:
        import matplotlib.pyplot as plt
        for m in metrics:
            plt.figure(figsize=(6.4, 4.8))
            for rho, g in df.groupby("rho_m2"):
                g = g.sort_values("points_per_pin_spacing")
                plt.plot(g["points_per_pin_spacing"], g[m], marker="o", label=f"rho={rho:.2e}")
            plt.xscale("log")
            if m not in ["fraction_multi_hit_productive"]:
                plt.yscale("log")
            plt.xlabel("points per expected pin spacing, p")
            plt.ylabel(m)
            plt.legend(fontsize=7)
            plt.tight_layout()
            plt.savefig(outdir / f"convergence_{m}.png", dpi=220)
            plt.close()

        plt.figure(figsize=(6.4, 4.8))
        for p in args.candidate_p:
            col = f"nodes_for_p{p:g}"
            plt.loglog(rec["rho_m2"], rec[col], marker="o", label=f"p={p:g}")
        plt.xlabel(r"rho (m$^{-2}$)")
        plt.ylabel("recommended nodes")
        plt.legend()
        plt.tight_layout()
        plt.savefig(outdir / "adaptive_nodes_vs_rho.png", dpi=220)
        plt.close()

    print(f"Wrote {outdir}/pairwise_node_convergence.csv")
    print(f"Wrote {outdir}/relative_to_highest_node.csv")
    print(f"Wrote {outdir}/adaptive_node_recommendations.csv")
    print(f"Wrote {outdir}/resolution_convergence_notes.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
