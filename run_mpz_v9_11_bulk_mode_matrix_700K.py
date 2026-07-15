#!/usr/bin/env python3
"""Run ceramic, weakT and DBTT at 700 K in both v9.11 bulk modes."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import subprocess
import sys

import pandas as pd

from arrhenius_fracture.bulk_state_v911 import VALID_BULK_MODES
from arrhenius_fracture.mpz_parameterization_v911 import normalize_class_name

CLASSES = ("ceramic", "weakT", "DBTT")


def run_one(args, cls: str, mode: str, root: Path) -> dict:
    out = root / mode
    cmd = [
        sys.executable,
        "run_mpz_v9_11_mode_i_rcurve_3T.py",
        "--parameter-root", str(args.parameter_root),
        "--material-class", cls,
        "--bulk-plasticity-mode", mode,
        "--temperatures", str(args.T_K),
        "--outroot", str(out),
        "--target-extension-um", str(args.target_extension_um),
        "--steps", str(args.steps),
        "--nx", str(args.nx),
        "--ny", str(args.ny),
        "--tip-h-fine", str(args.tip_h_fine),
        "--tip-ratio", str(args.tip_ratio),
        "--dU", str(args.dU),
        "--dt", str(args.dt),
        "--n-stagger", str(args.n_stagger),
        "--print-every", str(args.print_every),
        "--adaptive-event-target", str(args.adaptive_event_target),
        "--da-phys-um", str(args.da_phys_um),
        "--mpz-length-um", str(args.mpz_length_um),
        "--mpz-n-bins", str(args.mpz_n_bins),
        "--crystal-theta-deg", str(args.crystal_theta_deg),
        "--save-snapshots", str(args.save_snapshots),
        "--snapshot-cols", str(args.snapshot_cols),
        "--snapshot-by-extension-um", str(args.snapshot_by_extension_um),
    ]
    cmd.append("--skip-existing" if args.skip_existing else "--no-skip-existing")
    cmd.append("--make-solver-plots" if args.make_solver_plots else "--no-make-solver-plots")

    log_dir = root / "matrix_logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log = log_dir / f"{cls}_{mode}_{int(args.T_K)}K.log"
    print(f"START {cls:7s} {mode:18s} T={args.T_K:g} K")
    with log.open("w") as fp:
        cp = subprocess.run(cmd, stdout=fp, stderr=subprocess.STDOUT, text=True)

    summary_path = out / cls / "rcurve_temperature_summary.csv"
    if summary_path.exists():
        frame = pd.read_csv(summary_path)
        row = frame.iloc[0].to_dict() if not frame.empty else {}
    else:
        row = {}
    row.update({
        "class": normalize_class_name(cls),
        "bulk_plasticity_mode": mode,
        "T_K": float(args.T_K),
        "matrix_returncode": int(cp.returncode),
        "matrix_log": str(log),
    })
    print(
        f"DONE  {cls:7s} {mode:18s} rc={cp.returncode} "
        f"status={row.get('status')} ext={row.get('final_extension_um')} "
        f"Kinit={row.get('K_init_MPa_sqrt_m')}"
    )
    return row


def plot_mode_comparisons(root: Path, T_K: float, target_extension_um: float) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        print(f"plotting unavailable: {exc}")
        return

    for cls in CLASSES:
        fig, ax = plt.subplots(figsize=(7.8, 5.2))
        n = 0
        for mode in VALID_BULK_MODES:
            case_root = root / mode / normalize_class_name(cls)
            matches = sorted(case_root.glob(f"T{int(round(T_K))}_th*"))
            if not matches:
                continue
            path = matches[0] / "R_curve_event_sampled.csv"
            if not path.exists():
                continue
            rc = pd.read_csv(path)
            if rc.empty:
                continue
            ax.plot(
                rc["crack_extension_um"],
                rc["KJ_MPa_sqrt_m"],
                linewidth=1.2,
                label=mode,
            )
            n += 1
        if n:
            ax.set_xlabel(r"Projected crack extension $\Delta a_x$ ($\mu$m)")
            ax.set_ylabel(r"$K_J$ (MPa$\sqrt{m}$)")
            ax.set_title(f"{normalize_class_name(cls)}, {T_K:g} K: bulk-mode comparison")
            ax.set_xlim(left=0.0, right=float(target_extension_um))
            ax.grid(alpha=0.25)
            ax.legend(frameon=False)
            fig.tight_layout()
            fig.savefig(root / f"{normalize_class_name(cls)}_bulk_mode_comparison_{int(T_K)}K.png", dpi=220)
        plt.close(fig)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--parameter-root", type=Path, default=Path("mpz_v9_11_parameters"))
    p.add_argument("--outroot", type=Path, default=Path("runs/mpz_v9_11_bulk_mode_matrix_700K_v1"))
    p.add_argument("--T-K", type=float, default=700.0)
    p.add_argument("--target-extension-um", type=float, default=500.0)
    p.add_argument("--steps", type=int, default=6000)
    p.add_argument("--nx", type=int, default=36)
    p.add_argument("--ny", type=int, default=72)
    p.add_argument("--tip-h-fine", type=float, default=1.0e-6)
    p.add_argument("--tip-ratio", type=float, default=1.20)
    p.add_argument("--dU", type=float, default=2.0e-7)
    p.add_argument("--dt", type=float, default=8.4)
    p.add_argument("--n-stagger", type=int, default=2)
    p.add_argument("--print-every", type=int, default=25)
    p.add_argument("--adaptive-event-target", type=float, default=0.15)
    p.add_argument("--da-phys-um", type=float, default=5.0)
    p.add_argument("--mpz-length-um", type=float, default=100.0)
    p.add_argument("--mpz-n-bins", type=int, default=200)
    p.add_argument("--crystal-theta-deg", type=float, default=45.0)
    p.add_argument("--save-snapshots", type=int, default=12)
    p.add_argument("--snapshot-cols", type=int, default=4)
    p.add_argument("--snapshot-by-extension-um", type=float, default=50.0)
    p.add_argument("--make-solver-plots", action=argparse.BooleanOptionalAction, default=False)
    p.add_argument("--skip-existing", action=argparse.BooleanOptionalAction, default=True)
    args = p.parse_args()

    root = args.outroot.resolve()
    root.mkdir(parents=True, exist_ok=True)
    rows = []
    for cls in CLASSES:
        for mode in VALID_BULK_MODES:
            rows.append(run_one(args, cls, mode, root))

    frame = pd.DataFrame(rows)
    frame.to_csv(root / "bulk_mode_matrix_700K_summary.csv", index=False)
    (root / "bulk_mode_matrix_700K_summary.json").write_text(
        json.dumps(rows, indent=2, default=str)
    )
    plot_mode_comparisons(root, args.T_K, args.target_extension_um)
    print(frame.to_string(index=False))
    print("wrote", root / "bulk_mode_matrix_700K_summary.csv")
    if any(int(row.get("matrix_returncode", 1)) != 0 for row in rows):
        raise SystemExit(
            "one or more matrix cases returned nonzero; inspect matrix_logs and "
            "the per-case run.log before interpreting the comparison"
        )


if __name__ == "__main__":
    main()
