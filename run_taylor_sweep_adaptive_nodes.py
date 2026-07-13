#!/usr/bin/env python3
"""
run_taylor_sweep_adaptive_nodes.py

Wrapper around run_taylor_sweep.py that uses a density-dependent line resolution:

    N_nodes(rho) = ceil_to_multiple(p * L_line_reduced * b * sqrt(rho))

bounded by --min-nodes and --max-nodes.  This keeps the number of line nodes per
expected forest spacing approximately constant at high density, while preserving
a minimum resolution at low density where contact detection/geometry still needs
an adequate line discretization.

This wrapper runs one density at a time so each density can use a different
--nodes value.  It also writes the plan and exact commands before executing.

Dry run example:
  python3 run_taylor_sweep_adaptive_nodes.py --outroot taylor_T500_adaptN_q3 \
    --temperature-K 500 --rho-min 2e14 --rho-max 2e17 --n-density 24 \
    --points-per-pin-spacing 3 --min-nodes 192 --max-nodes 640 --round-multiple 32 \
    --dry-run -- \
    --crossing-drive-mode force_work --capture-radius 8 --dt 2e-9 --target-strain 0.05 \
    --cross-sigc0-MPa 1800 --peierls-cross-ratio 100 --avalanche-eps-min 0.0015 \
    --avalanche-n-boot 0 --no-video

Execute by replacing --dry-run with --execute.
"""
from __future__ import annotations
import argparse
import csv
import math
import os
from pathlib import Path
import shlex
import subprocess
import sys


def round_up_multiple(x: float, mult: int) -> int:
    if mult <= 1:
        return int(math.ceil(x))
    return int(math.ceil(x / mult) * mult)


def safe_rho_label(rho: float) -> str:
    return f"rho{rho:.3e}".replace("+", "").replace(".", "p")


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Run run_taylor_sweep.py density-by-density with adaptive --nodes."
    )
    ap.add_argument("--outroot", required=True)
    ap.add_argument("--driver", default="run_taylor_sweep.py")
    ap.add_argument("--temperature-K", type=float, default=500.0)
    ap.add_argument("--rho-min", type=float, required=True)
    ap.add_argument("--rho-max", type=float, required=True)
    ap.add_argument("--n-density", type=int, required=True)
    ap.add_argument("--b-m", type=float, default=2.48e-10)
    ap.add_argument("--line-length-reduced", type=float, default=1530.0)
    ap.add_argument("--points-per-pin-spacing", type=float, default=3.0,
                    help="Target pin-spacing/node-spacing ratio p.")
    ap.add_argument("--min-nodes", type=int, default=192)
    ap.add_argument("--max-nodes", type=int, default=640)
    ap.add_argument("--round-multiple", type=int, default=32)
    ap.add_argument("--execute", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--continue-on-error", action="store_true")
    ap.add_argument("--combine-trends", action="store_true", default=True)
    args, extras = ap.parse_known_args()

    if "--" in extras:
        extras = [x for x in extras if x != "--"]

    outroot = Path(args.outroot)
    outroot.mkdir(parents=True, exist_ok=True)

    if args.n_density == 1:
        rhos = [args.rho_min]
    else:
        lo = math.log10(args.rho_min); hi = math.log10(args.rho_max)
        rhos = [10 ** (lo + i * (hi - lo) / (args.n_density - 1)) for i in range(args.n_density)]

    plan_rows = []
    commands = []
    for i, rho in enumerate(rhos):
        pin_spacing_reduced = 1.0 / (args.b_m * math.sqrt(rho))
        raw_nodes = args.points_per_pin_spacing * args.line_length_reduced * args.b_m * math.sqrt(rho)
        nodes = round_up_multiple(raw_nodes, args.round_multiple)
        nodes = max(args.min_nodes, min(args.max_nodes, nodes))
        actual_p = pin_spacing_reduced / (args.line_length_reduced / nodes)
        run_out = outroot / f"run_{i:02d}_{safe_rho_label(rho)}_N{nodes}"
        cmd = [
            sys.executable, args.driver,
            "--outroot", str(run_out),
            "--temperatures-K", f"{args.temperature_K:g}",
            "--rho-min", f"{rho:.17g}",
            "--rho-max", f"{rho:.17g}",
            "--n-density", "1",
            "--nodes", str(nodes),
        ] + extras
        plan_rows.append({
            "index": i,
            "rho_m2": rho,
            "sqrt_rho": math.sqrt(rho),
            "pin_spacing_reduced": pin_spacing_reduced,
            "target_points_per_pin_spacing": args.points_per_pin_spacing,
            "raw_nodes": raw_nodes,
            "nodes": nodes,
            "actual_points_per_pin_spacing": actual_p,
            "outroot": str(run_out),
            "command": " ".join(shlex.quote(c) for c in cmd),
        })
        commands.append(cmd)

    with open(outroot / "adaptive_node_plan.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(plan_rows[0].keys()))
        w.writeheader(); w.writerows(plan_rows)
    with open(outroot / "adaptive_node_commands.sh", "w") as f:
        f.write("#!/usr/bin/env bash\nset -euo pipefail\n\n")
        for row in plan_rows:
            f.write(row["command"] + "\n")
    os.chmod(outroot / "adaptive_node_commands.sh", 0o755)

    print(f"Wrote {outroot/'adaptive_node_plan.csv'}")
    print(f"Wrote {outroot/'adaptive_node_commands.sh'}")
    print("\nPlan:")
    for r in plan_rows:
        print(f"  {r['index']:02d}: rho={r['rho_m2']:.3e}, nodes={r['nodes']}, p_actual={r['actual_points_per_pin_spacing']:.2f}")

    if args.dry_run or not args.execute:
        print("\nDry run only. Re-run with --execute, or run the generated adaptive_node_commands.sh.")
        return 0

    failures = []
    for row, cmd in zip(plan_rows, commands):
        print(f"\n[run {row['index']:02d}] rho={row['rho_m2']:.3e}, nodes={row['nodes']}")
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            failures.append((row, e.returncode))
            print(f"[FAILED] rho={row['rho_m2']:.3e}, nodes={row['nodes']}, returncode={e.returncode}")
            if not args.continue_on_error:
                return e.returncode

    # Combine one-row taylor_trend.csv files when available.
    if args.combine_trends:
        rows = []
        for r in plan_rows:
            td = Path(r["outroot"])
            candidates = list(td.rglob("taylor_trend.csv"))
            for c in candidates:
                try:
                    with open(c, newline="") as f:
                        reader = csv.DictReader(f)
                        for row in reader:
                            row["adaptive_nodes"] = r["nodes"]
                            row["adaptive_actual_points_per_pin_spacing"] = f"{r['actual_points_per_pin_spacing']:.8g}"
                            row["source_taylor_trend"] = str(c)
                            rows.append(row)
                except Exception:
                    pass
        if rows:
            fieldnames = []
            for row in rows:
                for k in row.keys():
                    if k not in fieldnames:
                        fieldnames.append(k)
            with open(outroot / "adaptive_taylor_trend_combined.csv", "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=fieldnames)
                w.writeheader(); w.writerows(rows)
            print(f"\nWrote {outroot/'adaptive_taylor_trend_combined.csv'}")

    if failures:
        with open(outroot / "adaptive_node_failures.txt", "w") as f:
            for row, rc in failures:
                f.write(f"rho={row['rho_m2']:.17g}, nodes={row['nodes']}, returncode={rc}\n")
        print(f"\nCompleted with {len(failures)} failure(s).")
        return 1

    print("\nAll adaptive-node runs completed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
