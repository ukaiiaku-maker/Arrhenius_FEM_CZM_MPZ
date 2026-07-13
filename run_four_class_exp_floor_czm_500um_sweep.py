#!/usr/bin/env python3
"""Production temperature sweep for four tuned fully EXP-floor fracture classes.

Runs the adaptive FEM/CZM sharp-front model for ceramic, peak, weakT, and DBTT
parameterizations from 300 to 1200 K at 100 K increments. Branch formation is
disabled; max_fronts=1 is also imposed as a defensive guard. Each case runs to
500 um projected crack extension by default.

The emission design CSV stores the pre-scale values used by the V1 search. The
V1 calibration applies 0.75 to emission G00 and gT; this driver therefore passes
the effective values 0.75*G00 and 0.75*gT to the FEM/CZM emission EXP-floor
surface so the 2-D model receives the same local emission barrier.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
import shlex
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd

CLASSES_DEFAULT = ["ceramic", "peak", "weakT", "DBTT"]
TEMPS_DEFAULT = list(range(300, 1201, 100))
A0_MM = 0.5


def parse_list(text: str, cast=str):
    return [cast(x) for x in text.replace(",", " ").split() if x]


def load_parameters(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    required = {
        "target_class", "exp_G00_eV", "exp_gT_eV_per_K", "exp_sigc0_GPa",
        "exp_sT_MPa_per_K", "exp_a", "exp_n", "exp_floor_frac",
        "cleave_G00_eV", "cleave_gT_eV_per_K", "cleave_sigc0_GPa",
        "cleave_sT_MPa_per_K", "cleave_exp_a", "cleave_exp_n",
        "cleave_floor_frac", "cleave_S_hs_kB", "chi_shield", "N_sat",
    }
    missing = sorted(required.difference(df.columns))
    if missing:
        raise SystemExit(f"parameter table missing columns: {missing}")
    return df.set_index("target_class", drop=False)


def resolve_python(args) -> str:
    if args.python_bin:
        p = Path(args.python_bin).expanduser().resolve()
        if not p.exists():
            raise SystemExit(f"PYTHON_BIN does not exist: {p}")
        return str(p)
    cp = subprocess.run(
        ["conda", "run", "-n", args.conda_env, "python", "-c", "import sys; print(sys.executable)"],
        text=True, capture_output=True,
    )
    if cp.returncode != 0:
        raise SystemExit(f"could not resolve Python in Conda env {args.conda_env!r}:\n{cp.stderr}")
    lines = [x.strip() for x in cp.stdout.splitlines() if x.strip()]
    if not lines:
        raise SystemExit(f"Conda env {args.conda_env!r} returned no Python path")
    return lines[-1]


def preflight(py: str) -> None:
    cp = subprocess.run(
        [py, "-m", "arrhenius_fracture.sharp_front", "--help"],
        text=True, capture_output=True,
    )
    if cp.returncode != 0:
        raise SystemExit(f"sharp_front preflight failed:\n{cp.stderr}")
    help_text = cp.stdout + cp.stderr
    required = [
        "--emit-barrier-kind", "--emit-G00-eV", "--emit-exp-a",
        "--cleave-barrier-kind", "--target-crack-extension-um",
        "--crack-backend",
    ]
    missing = [x for x in required if x not in help_text]
    if missing:
        raise SystemExit(
            "sharp_front is missing required fully EXP-floor/adaptive-CZM options: "
            + ", ".join(missing)
            + "\nApply the four-class production sweep patch before launching the batch."
        )
    cp = subprocess.run(
        [py, "-c", "import numpy, scipy, arrhenius_fracture; "
                    "print('preflight OK', numpy.__version__, scipy.__version__)"],
        text=True, capture_output=True,
    )
    if cp.returncode != 0:
        raise SystemExit(f"Python package preflight failed:\n{cp.stderr}")
    print(cp.stdout.strip())


def completion_status(case_dir: Path, target_um: float) -> tuple[bool, float | None]:
    summary = case_dir / "summary.json"
    if not summary.exists():
        return False, None
    try:
        data = json.loads(summary.read_text())
        if isinstance(data, list):
            d = data[0] if data else {}
        else:
            d = data
        afinal = float(d["a_final_mm"])
        ext_um = (afinal - A0_MM) * 1000.0
        return ext_um >= target_um - 0.5, ext_um
    except Exception:
        return False, None


def fstr(x) -> str:
    x = float(x)
    if math.isinf(x):
        return "inf"
    return f"{x:.16g}"


def build_command(py: str, row: pd.Series, T: int, case_dir: Path, args) -> list[str]:
    # Match the V1 tuning convention: scale emission G00 and gT by 0.75.
    emit_G00_eff = 0.75 * float(row.exp_G00_eV)
    emit_gT_eff = 0.75 * float(row.exp_gT_eV_per_K)

    cmd = [
        py, "-m", "arrhenius_fracture.sharp_front",
        "--mode", "2d",
        "--nx", str(args.nx), "--ny", str(args.ny),
        "--tip-h-fine", fstr(args.tip_h_fine), "--tip-ratio", fstr(args.tip_ratio),
        "--dU", fstr(args.dU), "--dt", fstr(args.dt),
        "--steps", str(args.long_steps), "--n-stagger", str(args.n_stagger),
        "--print-every", str(args.print_every),
        "--target-crack-extension-um", fstr(args.target_ext_um),
        "--crystal-aniso", "--crystal-compete", "--crystal-material", args.crystal_material,
        "--cleave-gamma-aniso", fstr(args.cleave_gamma_aniso),
        "--multihit-m", "3", "--multihit-tau", "1e-6",
        "--emb-sat-frac", "1",
        "--adaptive-events", "--adaptive-event-target", "0.35",
        "--adaptive-min-frac", "1e-8", "--adaptive-grow", "4.0",
        # Branching is OFF: no --crystal-branch flag, plus a hard single-front guard.
        "--max-fronts", "1",
        "--da-phys", "5e-6",
        "--j-decomposition", "cluster", "--rJ-cluster", "20e-6", "--rJ-outer", "25e-6",
        "--temperatures", str(T), "--crystal-theta-deg", fstr(args.theta),
        "--crack-backend", "adaptive_czm", "--czm-max-angle-error-deg", "35",
        # Fully EXP-floor local emission hazard.
        "--emit-barrier-kind", "exp_floor",
        "--emit-G00-eV", fstr(emit_G00_eff),
        "--emit-gT-eV-per-K", fstr(emit_gT_eff),
        "--emit-sigc0-GPa", fstr(row.exp_sigc0_GPa),
        "--emit-sT-GPa-per-K", fstr(float(row.exp_sT_MPa_per_K) / 1000.0),
        "--emit-exp-a", fstr(row.exp_a), "--emit-exp-n", fstr(row.exp_n),
        "--emit-floor-frac", fstr(row.exp_floor_frac), "--emit-Tref-K", "300",
        # Fully EXP-floor crack-opening hazard.
        "--cleave-barrier-kind", "exp_floor", "--cleave-exp-T-mode", "linear",
        "--cleave-G00-eV", fstr(row.cleave_G00_eV),
        "--cleave-gT-eV-per-K", fstr(row.cleave_gT_eV_per_K),
        "--cleave-sigc0-GPa", fstr(row.cleave_sigc0_GPa),
        "--cleave-sT-GPa-per-K", fstr(float(row.cleave_sT_MPa_per_K) / 1000.0),
        "--cleave-exp-a", fstr(row.cleave_exp_a), "--cleave-exp-n", fstr(row.cleave_exp_n),
        "--cleave-floor-frac", fstr(row.cleave_floor_frac),
        "--cleave-S-hs-kB", fstr(row.cleave_S_hs_kB),
        "--cleave-sigma-S-GPa", "6", "--cleave-S-hs-power", "2",
        "--cleave-S-hs-Tref-K", "300", "--cleave-Tref-K", "300",
        "--cleave-shield-chi", fstr(row.chi_shield),
        "--n-sat", fstr(row.N_sat),
        "--out", str(case_dir),
    ]

    if args.save_snapshots > 0:
        cmd += ["--save-snapshots", str(args.save_snapshots),
                "--snapshot-cols", str(args.snapshot_cols)]
        if args.snapshot_by_ext_um > 0:
            cmd += ["--snapshot-by-crack-extension-um", fstr(args.snapshot_by_ext_um)]
    else:
        cmd += ["--save-snapshots", "0", "--no-plots"]
    return cmd


def run_case(py: str, row: pd.Series, klass: str, T: int, root: Path, args) -> dict:
    case_dir = root / klass / f"T{T}_th{args.theta:g}"
    case_dir.mkdir(parents=True, exist_ok=True)
    complete, ext = completion_status(case_dir, args.target_ext_um)
    if complete and not args.force:
        rc = process_case_r_curve(case_dir, klass, T, args.target_ext_um)
        print(f"SKIP {klass:8s} T={T:4d} K: completed extension={ext:.1f} um; R-points={len(rc)}")
        return {"class": klass, "T_K": T, "status": "skipped_complete",
                "extension_um": ext, "n_r_curve_points": len(rc)}

    cmd = build_command(py, row, T, case_dir, args)
    (case_dir / "command.txt").write_text(shlex.join(cmd) + "\n")
    param_payload = row.to_dict()
    param_payload.update({
        "target_class": klass,
        "T_K": T,
        "emission_G00_effective_eV": 0.75 * float(row.exp_G00_eV),
        "emission_gT_effective_eV_per_K": 0.75 * float(row.exp_gT_eV_per_K),
        "branching_enabled": False,
        "target_extension_um": args.target_ext_um,
    })
    (case_dir / "resolved_parameters.json").write_text(json.dumps(param_payload, indent=2, default=str))

    print(f"START {klass:8s} T={T:4d} K -> {case_dir}")
    with (case_dir / "run.log").open("w") as log:
        cp = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, text=True)
    complete, ext = completion_status(case_dir, args.target_ext_um)
    if cp.returncode == 0 and complete:
        (case_dir / ".long_growth_complete").touch()
        rc = process_case_r_curve(case_dir, klass, T, args.target_ext_um)
        print(f"DONE  {klass:8s} T={T:4d} K extension={ext:.1f} um; R-points={len(rc)}")
        status = "complete"
    elif cp.returncode == 0:
        print(f"INCOMPLETE {klass:8s} T={T:4d} K extension={ext}")
        status = "incomplete"
    else:
        print(f"FAILED {klass:8s} T={T:4d} K rc={cp.returncode}")
        status = "failed"
    if cp.returncode != 0 or not complete:
        rc = process_case_r_curve(case_dir, klass, T, args.target_ext_um)
    return {"class": klass, "T_K": T, "status": status, "extension_um": ext,
            "returncode": cp.returncode, "n_r_curve_points": len(rc),
            "case_dir": str(case_dir.relative_to(root))}



def find_steps_file(case_dir: Path, T: int) -> Path | None:
    exact = case_dir / f"steps_{int(T):04d}K.csv"
    if exact.exists():
        return exact
    matches = sorted(case_dir.glob("steps_*K.csv"))
    return matches[0] if matches else None


def extract_r_curve(case_dir: Path, T: int) -> pd.DataFrame:
    """Extract event-sampled KJ(delta-a) data from the accepted step history.

    The raw step table contains loading steps at fixed crack length as well as
    crack-growth events.  For the R-curve-like output we retain only rows where
    accepted projected crack extension increased (or n_fire>0 as a fallback).
    This gives one ordered propagation-resistance history without vertical
    loading segments at fixed delta-a.  The full raw steps CSV remains intact.
    """
    sf = find_steps_file(case_dir, T)
    if sf is None:
        return pd.DataFrame()
    try:
        st = pd.read_csv(sf)
    except Exception as exc:
        print(f"WARNING could not read {sf}: {exc}")
        return pd.DataFrame()
    required = {"KJ_Pa_sqrtm", "a_tip_m"}
    if not required.issubset(st.columns):
        print(f"WARNING {sf} missing R-curve columns {sorted(required.difference(st.columns))}")
        return pd.DataFrame()

    if "crack_extension_m" in st.columns:
        ext_m = pd.to_numeric(st["crack_extension_m"], errors="coerce").to_numpy(float)
    else:
        ext_m = pd.to_numeric(st["a_tip_m"], errors="coerce").to_numpy(float) - A0_MM * 1e-3
    ext_m = np.maximum(ext_m, 0.0)

    if "da_block_m" in st.columns:
        da_m = pd.to_numeric(st["da_block_m"], errors="coerce").fillna(0.0).to_numpy(float)
    else:
        da_m = np.r_[0.0, np.maximum(np.diff(ext_m), 0.0)]
    n_fire = (pd.to_numeric(st["n_fire"], errors="coerce").fillna(0.0).to_numpy(float)
              if "n_fire" in st.columns else np.zeros(len(st)))
    growth = (da_m > 1e-12) | (n_fire > 0.0)
    idx = np.flatnonzero(growth)
    if len(idx) == 0:
        return pd.DataFrame()

    out = pd.DataFrame({
        "growth_event_id": np.arange(1, len(idx) + 1, dtype=int),
        "step": pd.to_numeric(st.iloc[idx]["step"], errors="coerce").to_numpy() if "step" in st.columns else idx,
        "crack_extension_um": ext_m[idx] * 1e6,
        "a_tip_mm": pd.to_numeric(st.iloc[idx]["a_tip_m"], errors="coerce").to_numpy(float) * 1e3,
        "da_block_um": da_m[idx] * 1e6,
        "KJ_MPa_sqrt_m": pd.to_numeric(st.iloc[idx]["KJ_Pa_sqrtm"], errors="coerce").to_numpy(float) / 1e6,
    })
    optional = {
        "N_em": ("N_em", 1.0),
        "B": ("B", 1.0),
        "n_fire": ("n_fire", 1.0),
        "sigma_tip_GPa": ("sigma_tip_Pa", 1e-9),
        "sigma_back_GPa": ("sigma_back_Pa", 1e-9),
        "lambda_c_per_s": ("lambda_c", 1.0),
        "lambda_e_per_s": ("lambda_e", 1.0),
        "G_cleave_eff_eV": ("G_cleave_eff_eV", 1.0),
    }
    for out_name, (src_name, scale) in optional.items():
        if src_name in st.columns:
            out[out_name] = pd.to_numeric(st.iloc[idx][src_name], errors="coerce").to_numpy(float) * scale

    # Keep the authoritative event order.  Remove only exact duplicate rows
    # introduced by pathological repeated writes, not repeated extensions from
    # legitimate multi-hit event bookkeeping.
    out = out.replace([np.inf, -np.inf], np.nan)
    out = out.dropna(subset=["crack_extension_um", "KJ_MPa_sqrt_m"]).reset_index(drop=True)
    out["growth_event_id"] = np.arange(1, len(out) + 1, dtype=int)
    return out


def propagation_metrics(rc: pd.DataFrame, Kinit: float | None = None) -> dict:
    if rc.empty:
        return {
            "Kprop_200_500um_median": np.nan,
            "Kprop_200_500um_mean": np.nan,
            "Kprop_200_500um_p10": np.nan,
            "Kprop_200_500um_p90": np.nan,
            "delta_KR_median_minus_init": np.nan,
            "n_growth_events": 0,
        }
    late = rc[(rc.crack_extension_um >= 200.0) & (rc.crack_extension_um <= 500.5)]
    vals = late.KJ_MPa_sqrt_m.to_numpy(float)
    if len(vals) == 0:
        med = mean = p10 = p90 = np.nan
    else:
        med = float(np.nanmedian(vals)); mean = float(np.nanmean(vals))
        p10 = float(np.nanpercentile(vals, 10)); p90 = float(np.nanpercentile(vals, 90))
    dkr = med - float(Kinit) if Kinit is not None and np.isfinite(med) else np.nan
    return {
        "Kprop_200_500um_median": med,
        "Kprop_200_500um_mean": mean,
        "Kprop_200_500um_p10": p10,
        "Kprop_200_500um_p90": p90,
        "delta_KR_median_minus_init": dkr,
        "n_growth_events": int(len(rc)),
    }


def plot_case_r_curve(rc: pd.DataFrame, case_dir: Path, klass: str, T: int,
                      target_ext_um: float) -> None:
    if rc.empty:
        return
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7.2, 4.8))
        ax.plot(rc.crack_extension_um, rc.KJ_MPa_sqrt_m, marker="o", markersize=2.5,
                linewidth=1.2)
        ax.set_xlabel(r"Projected crack extension $\Delta a_x$ ($\mu$m)")
        ax.set_ylabel(r"$K_J$ (MPa$\sqrt{m}$)")
        ax.set_title(f"{klass}, {T} K")
        ax.set_xlim(left=0.0, right=max(float(target_ext_um), float(rc.crack_extension_um.max())))
        ax.grid(alpha=0.25)
        fig.tight_layout()
        fig.savefig(case_dir / "R_curve_K_vs_crack_extension.png", dpi=220)
        plt.close(fig)
    except Exception as exc:
        print(f"WARNING R-curve plotting failed for {klass} T={T} K: {exc}")


def process_case_r_curve(case_dir: Path, klass: str, T: int,
                         target_ext_um: float) -> pd.DataFrame:
    rc = extract_r_curve(case_dir, T)
    if rc.empty:
        return rc
    rc.insert(0, "class", klass)
    rc.insert(1, "T_K", int(T))
    rc.to_csv(case_dir / "R_curve_event_sampled.csv", index=False)
    plot_case_r_curve(rc, case_dir, klass, T, target_ext_um)
    return rc


def make_class_r_curve_overlays(root: Path, classes: list[str], temps: list[int],
                                target_ext_um: float) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        for klass in classes:
            class_dir = root / klass
            if not class_dir.exists():
                continue
            fig, ax = plt.subplots(figsize=(7.8, 5.2))
            n = 0
            for T in temps:
                matches = sorted(class_dir.glob(f"T{T}_th*"))
                if not matches:
                    continue
                f = matches[0] / "R_curve_event_sampled.csv"
                if not f.exists():
                    continue
                rc = pd.read_csv(f)
                if rc.empty:
                    continue
                ax.plot(rc.crack_extension_um, rc.KJ_MPa_sqrt_m, linewidth=1.1,
                        label=f"{T} K")
                n += 1
            if n:
                ax.set_xlabel(r"Projected crack extension $\Delta a_x$ ($\mu$m)")
                ax.set_ylabel(r"$K_J$ (MPa$\sqrt{m}$)")
                ax.set_title(f"{klass}: propagation resistance curves")
                ax.set_xlim(left=0.0, right=float(target_ext_um))
                ax.grid(alpha=0.25)
                ax.legend(frameon=False, ncol=2, fontsize=8)
                fig.tight_layout()
                fig.savefig(class_dir / "R_curves_all_temperatures.png", dpi=220)
            plt.close(fig)
    except Exception as exc:
        print(f"WARNING class R-curve overlay plotting failed: {exc}")


def make_init_vs_propagation_plot(df: pd.DataFrame, root: Path) -> None:
    if df.empty or "Kprop_200_500um_median" not in df.columns:
        return
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        classes = list(df["class"].drop_duplicates())
        fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.6), sharex=True)
        axes = axes.ravel()
        for ax, klass in zip(axes, classes):
            g = df[df["class"] == klass].sort_values("T_K")
            ax.plot(g.T_K, g.Kc_first_MPa_sqrt_m, marker="o", label="initiation")
            ax.plot(g.T_K, g.Kprop_200_500um_median, marker="s",
                    label="median, 200-500 μm")
            ax.set_title(klass)
            ax.set_ylabel(r"$K$ (MPa$\sqrt{m}$)")
            ax.grid(alpha=0.25)
        for ax in axes[-2:]:
            ax.set_xlabel("Temperature (K)")
        axes[0].legend(frameon=False)
        fig.tight_layout()
        fig.savefig(root / "four_class_init_vs_propagation_K_vs_T.png", dpi=220)
        plt.close(fig)
    except Exception as exc:
        print(f"WARNING initiation/propagation plotting failed: {exc}")


def collect_summary(root: Path, classes: list[str], temps: list[int]) -> pd.DataFrame:
    rows = []
    for klass in classes:
        for T in temps:
            dirs = list((root / klass).glob(f"T{T}_th*")) if (root / klass).exists() else []
            if not dirs:
                continue
            case_dir = dirs[0]
            sf = case_dir / "summary.json"
            if not sf.exists():
                continue
            try:
                data = json.loads(sf.read_text())
                d = data[0] if isinstance(data, list) else data
                Kinit = d.get("Kc_first_MPa_sqrt_m")
                rcf = case_dir / "R_curve_event_sampled.csv"
                if rcf.exists():
                    rc = pd.read_csv(rcf)
                else:
                    rc = process_case_r_curve(case_dir, klass, T, 500.0)
                metrics = propagation_metrics(rc, Kinit)
                row = {
                    "class": klass,
                    "T_K": T,
                    "Kc_first_MPa_sqrt_m": Kinit,
                    "a_final_mm": d.get("a_final_mm"),
                    "crack_extension_um": (float(d.get("a_final_mm")) - A0_MM) * 1000.0 if d.get("a_final_mm") is not None else None,
                    "N_em_final": d.get("N_em_final"),
                    "deflection_deg": d.get("deflection_deg"),
                    "n_fronts": d.get("n_fronts"),
                    "branched": d.get("branched"),
                    "mode": d.get("mode"),
                    "case_dir": str(case_dir.relative_to(root)),
                }
                row.update(metrics)
                rows.append(row)
            except Exception as exc:
                rows.append({"class": klass, "T_K": T, "summary_error": str(exc),
                             "case_dir": str(case_dir.relative_to(root))})
    return pd.DataFrame(rows).sort_values(["class", "T_K"]) if rows else pd.DataFrame()


def make_plot(df: pd.DataFrame, root: Path) -> None:
    if df.empty or "Kc_first_MPa_sqrt_m" not in df:
        return
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(7.4, 5.0))
        for klass, g in df.groupby("class", sort=False):
            g = g.sort_values("T_K")
            ax.plot(g.T_K, g.Kc_first_MPa_sqrt_m, marker="o", label=klass)
        ax.set_xlabel("Temperature (K)")
        ax.set_ylabel(r"First-passage $K_c$ (MPa$\sqrt{m}$)")
        ax.grid(alpha=0.25)
        ax.legend(frameon=False)
        fig.tight_layout()
        fig.savefig(root / "four_class_CZM_Kc_vs_T.png", dpi=220)
        plt.close(fig)
    except Exception as exc:
        print(f"WARNING plotting failed: {exc}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--parameters", default="four_class_exp_floor_exact_model_inputs.csv")
    ap.add_argument("--outroot", default="runs/four_class_exp_floor_CZM_no_branch_500um_theta45")
    ap.add_argument("--classes", default=" ".join(CLASSES_DEFAULT))
    ap.add_argument("--temps", default=" ".join(map(str, TEMPS_DEFAULT)))
    ap.add_argument("--theta", type=float, default=45.0)
    ap.add_argument("--target-ext-um", type=float, default=500.0)
    ap.add_argument("--long-steps", type=int, default=20000)
    ap.add_argument("--max-jobs", type=int, default=1)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--conda-env", default="arrhenius-fem-czm")
    ap.add_argument("--python-bin", default="")
    ap.add_argument("--nx", type=int, default=12)
    ap.add_argument("--ny", type=int, default=24)
    ap.add_argument("--tip-h-fine", type=float, default=5e-6)
    ap.add_argument("--tip-ratio", type=float, default=1.30)
    ap.add_argument("--dU", type=float, default=2e-7)
    ap.add_argument("--dt", type=float, default=8.4)
    ap.add_argument("--n-stagger", type=int, default=2)
    ap.add_argument("--print-every", type=int, default=500)
    ap.add_argument("--crystal-material", default="branchy")
    ap.add_argument("--cleave-gamma-aniso", type=float, default=2.0)
    ap.add_argument("--save-snapshots", type=int, default=0)
    ap.add_argument("--snapshot-cols", type=int, default=5)
    ap.add_argument("--snapshot-by-ext-um", type=float, default=50.0)
    args = ap.parse_args()

    classes = parse_list(args.classes, str)
    temps = parse_list(args.temps, int)
    params = load_parameters(Path(args.parameters))
    missing = [c for c in classes if c not in params.index]
    if missing:
        raise SystemExit(f"classes absent from parameter table: {missing}")

    root = Path(args.outroot)
    root.mkdir(parents=True, exist_ok=True)
    py = resolve_python(args)
    print(f"python: {py}")
    preflight(py)

    # Preserve exact parameter table and sweep configuration in the output root.
    root.mkdir(parents=True, exist_ok=True)
    Path(root / "four_class_exp_floor_exact_model_inputs.csv").write_bytes(Path(args.parameters).read_bytes())
    config = vars(args).copy()
    config.update({"resolved_classes": classes, "resolved_temperatures_K": temps,
                   "branching_enabled": False, "emission_energy_scale": 0.75})
    (root / "sweep_config.json").write_text(json.dumps(config, indent=2))

    tasks = [(c, T) for c in classes for T in temps]
    results = []
    if args.max_jobs <= 1:
        for c, T in tasks:
            results.append(run_case(py, params.loc[c], c, T, root, args))
    else:
        with ThreadPoolExecutor(max_workers=args.max_jobs) as ex:
            futs = {ex.submit(run_case, py, params.loc[c], c, T, root, args): (c, T) for c, T in tasks}
            for fut in as_completed(futs):
                results.append(fut.result())

    pd.DataFrame(results).sort_values(["class", "T_K"]).to_csv(root / "sweep_status.csv", index=False)
    summary = collect_summary(root, classes, temps)
    summary.to_csv(root / "four_class_temperature_summary.csv", index=False)
    make_plot(summary, root)
    make_class_r_curve_overlays(root, classes, temps, args.target_ext_um)
    make_init_vs_propagation_plot(summary, root)

    bad = [r for r in results if r["status"] in ("failed", "incomplete")]
    print(f"WROTE {root / 'four_class_temperature_summary.csv'}")
    if bad:
        print(f"WARNING: {len(bad)} cases failed or did not reach {args.target_ext_um:g} um")
        raise SystemExit(1)
    print("All requested cases completed.")


if __name__ == "__main__":
    main()
