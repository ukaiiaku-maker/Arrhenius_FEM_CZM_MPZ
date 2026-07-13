#!/usr/bin/env python3
"""Compare the V1 K-controlled sharp-front fatigue reduction against the v8 2-D
sharp-front multifront fatigue adapter over a small Kmax sweep.

The 2-D driver does not prescribe K directly.  This script uses a calibrated
linear displacement-to-local-K scale for the smoke-test geometry, ramps to that
amplitude on the first accepted step, then holds the displacement amplitude for
additional fatigue blocks using --fatigue-hold-load.  The comparison is therefore
near-constant local Kmax in both V1 and v8.
"""
from __future__ import annotations

import argparse
import csv
import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List

import pandas as pd


def run_cmd(cmd: List[str], log_path: Path) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w") as log:
        log.write("$ " + " ".join(cmd) + "\n\n")
        log.flush()
        proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        raise RuntimeError(f"command failed with code {proc.returncode}: {' '.join(cmd)}; see {log_path}")


def first_fire_cycles_from_v1(hist: pd.DataFrame) -> float | None:
    rows = hist.loc[hist["n_fire"] > 0]
    if rows.empty:
        return None
    return float(rows.iloc[0]["cycles_total"])


def first_fire_cycles_from_2d(steps: pd.DataFrame) -> float | None:
    rows = steps.loc[steps["n_fire"] > 0]
    if rows.empty:
        return None
    idx = rows.index[0]
    return float(steps.loc[:idx, "fatigue_cycles"].sum())


def first_fire_row(df: pd.DataFrame) -> pd.Series | None:
    rows = df.loc[df["n_fire"] > 0]
    if rows.empty:
        return None
    return rows.iloc[0]


def _as_float(x, default=float("nan")) -> float:
    try:
        if x is None:
            return default
        return float(x)
    except Exception:
        return default


def _read_first_K_from_steps(run_dir: Path, T: float) -> float:
    csv_path = run_dir / f"steps_{int(T):04d}K.csv"
    df = pd.read_csv(csv_path)
    if df.empty or "KJ_Pa_sqrtm" not in df.columns:
        return float("nan")
    return float(df.iloc[0].KJ_Pa_sqrtm) / 1e6


def calibrate_2d_displacement(args, out: Path, K_target: float, label: str) -> dict:
    """Calibrate dU so the first-step local J-derived K matches K_target.

    The 2-D driver prescribes displacement, not K.  A small geometry-dependent
    mismatch of only a few percent changes the Arrhenius hazard by orders of
    magnitude, so the comparison harness should solve for the displacement that
    gives the requested local K before running fatigue.  Calibration runs use
    dt=0 and one mechanical step, so they measure the J handoff without
    accumulating fracture/plasticity clocks.
    """
    dU = float(K_target) * float(args.dU_per_K_MPa)
    hist = []
    if not args.calibrate_2d_K:
        return {
            "dU": dU, "K_calibrated": float("nan"), "rel_error": float("nan"),
            "iters": 0, "history": hist, "enabled": False,
        }
    tol = float(args.K_calib_tol)
    for it in range(max(int(args.K_calib_iters), 1)):
        cdir = out / "_K_calibration" / f"K{label}_iter{it:02d}"
        cmd = [
            sys.executable, "-m", "arrhenius_fracture.sharp_front",
            "--mode", "2d",
            "--temperatures", str(args.T),
            "--steps", "1",
            "--dt", "0",
            "--nx", str(args.nx), "--ny", str(args.ny),
            "--tip-h-fine", str(args.tip_h_fine),
            "--tip-ratio", str(args.tip_ratio),
            "--n-stagger", "1",
            "--crystal-aniso",
            "--crystal-branch",
            "--j-decomposition", "cluster",
            "--sigma-cap-GPa", "0",
            "--dU", str(dU),
            "--save-snapshots", "0",
            "--no-plots",
            "--out", str(cdir),
        ]
        run_cmd(cmd, out / "_K_calibration" / f"K{label}_iter{it:02d}.log")
        K_meas = _read_first_K_from_steps(cdir, args.T)
        rel = (K_meas / K_target - 1.0) if K_target > 0 and K_meas == K_meas and K_meas > 0 else float("nan")
        hist.append({"iter": it, "dU": dU, "K_measured": K_meas, "rel_error": rel})
        if K_meas == K_meas and K_meas > 0 and abs(rel) <= tol:
            break
        if K_meas == K_meas and K_meas > 0:
            ratio = max(min(K_target / K_meas, 2.0), 0.5)
            dU *= ratio
        else:
            dU *= 1.25
    last = hist[-1] if hist else {"K_measured": float("nan"), "rel_error": float("nan"), "iter": -1}
    return {
        "dU": dU, "K_calibrated": float(last["K_measured"]),
        "rel_error": float(last["rel_error"]), "iters": int(last["iter"])+1,
        "history": hist, "enabled": True,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default="runs/v8_compare_1d_2d_Ksweep_300K")
    ap.add_argument("--Kmax-MPa-sqrt-m", type=float, nargs="+", default=[5.0, 10.0, 15.0, 20.0])
    ap.add_argument("--T", type=float, default=300.0)
    ap.add_argument("--R", type=float, default=0.1)
    ap.add_argument("--frequency-Hz", type=float, default=1000.0)
    ap.add_argument("--blocks", type=int, default=6, help="Maximum accepted blocks/steps. With --cycle-block-mode hazard_limited, use this as a safety bound; --cycles-max is the physical cycle horizon.")
    ap.add_argument("--cycles-max", type=float, default=1.0e12,
                    help="Physical fatigue-cycle horizon per K condition. Both V1 and v8 stop at this horizon unless first fire occurs earlier.")
    ap.add_argument("--calibrate-2d-K", action=argparse.BooleanOptionalAction, default=True,
                    help="Pre-calibrate the 2-D displacement amplitude so the initial local J-derived K matches the requested K. Enabled by default because small K offsets cause large Arrhenius hazard offsets.")
    ap.add_argument("--K-calib-iters", type=int, default=3, help="Maximum displacement calibration iterations per K.")
    ap.add_argument("--K-calib-tol", type=float, default=5e-3, help="Relative tolerance for initial 2-D local K calibration.")
    ap.add_argument("--cyclic-mechanics-phases", type=int, default=4,
                    help="Representative 2-D cyclic-mechanics phases per held fatigue block. Four is enough for this quick smoke comparison; use 8-24 for production convergence checks.")
    ap.add_argument("--block-cycles", type=float, default=1.0e-3,
                    help="Requested cycles per smoke block. Keep small so the full 2-D cyclic-mechanics block is fast at low K.")
    ap.add_argument("--min-block-cycles", type=float, default=1e-6)
    ap.add_argument("--max-block-cycles", type=float, default=1.0e6,
                    help="Upper bound on adaptive cycle-block size passed through to both drivers.")
    ap.add_argument("--cycle-block-mode", choices=["requested_cap", "hazard_limited"], default="requested_cap",
                    help="requested_cap keeps --block-cycles as a hard upper bound. hazard_limited lets the controller grow up to --max-block-cycles when hazards are small and shrink when any fracture/plasticity increment is too large.")
    ap.add_argument("--target-dB", type=float, default=0.2)
    ap.add_argument("--target-dN-store", type=float, default=0.25)
    ap.add_argument("--target-dN-emit", type=float, default=float("inf"))
    ap.add_argument("--target-dN-mobile", type=float, default=float("inf"))
    ap.add_argument("--target-dN-escape", type=float, default=float("inf"))
    ap.add_argument("--target-dN-peierls", type=float, default=float("inf"))
    ap.add_argument("--target-dN-taylor", type=float, default=float("inf"))
    ap.add_argument("--storage-model", default="escape_limited", choices=["escape_limited", "all_retained", "fixed_fraction"])
    ap.add_argument("--fixed-retained-fraction", type=float, default=0.1)
    ap.add_argument("--dU-per-K-MPa", type=float, default=5.65e-7,
                    help="2-D displacement amplitude per target MPa*sqrt(m) for the smoke geometry.")
    ap.add_argument("--nx", type=int, default=24)
    ap.add_argument("--ny", type=int, default=48)
    ap.add_argument("--tip-h-fine", type=float, default=1e-6)
    ap.add_argument("--tip-ratio", type=float, default=1.25)
    ap.add_argument("--save-snapshots", type=int, default=8,
                    help="Number of 2-D field snapshots rendered by the v8 sharp-front run. Use 0 for CSV-only runs.")
    ap.add_argument("--snapshot-cols", type=int, default=4,
                    help="Maximum number of columns in the rendered 2-D snapshot panel.")
    ap.add_argument("--make-2d-plots", action=argparse.BooleanOptionalAction, default=True,
                    help="Render v8 2-D diagnostic figures. Enabled by default for material-atlas runs; disable for fast CSV-only screening.")
    ap.add_argument("--min-global-forward", type=float, default=0.05,
                    help="Global +x crack-growth admissibility gate for through-ligament edge-crack geometry. Set to -1 only for geometries where back-branching is intended.")
    ap.add_argument("--da-phys", type=float, default=None,
                    help="Physical crack-advance quantum passed to v8 [m]. Supplying this makes da/dN bounds and target-da conversion explicit.")
    ap.add_argument("--crack-backend", choices=["sharp_wake", "edge_split_czm", "adaptive_czm"], default="sharp_wake",
                    help="2-D crack geometry backend passed to sharp_front.")
    ap.add_argument("--czm-penalty-normal", type=float, default=1.0e18)
    ap.add_argument("--czm-penalty-tangent", type=float, default=1.0e18)
    ap.add_argument("--czm-event-damage", type=float, default=1.0)
    ap.add_argument("--czm-max-angle-error-deg", type=float, default=35.0)
    ap.add_argument("--target-crack-extension-um", type=float, default=float("inf"),
                    help="Stop each 2-D run once leading crack extension reaches this value [um]. Use for long-growth morphology/Paris runs.")
    ap.add_argument("--snapshot-by-crack-extension-um", type=float, default=0.0,
                    help="Ask v8 to save snapshots when crack extension crosses this interval [um], in addition to block/snapshot count triggers.")
    ap.add_argument("--max-da-per-block-um", type=float, default=float("inf"),
                    help="Audit warning if one accepted 2-D fatigue block advances the leading crack by more than this [um].")
    ap.add_argument("--target-da-per-block-um", type=float, default=None,
                    help="Convert a desired expected crack extension per cycle block into a tighter --target-dB using da_phys. This improves long-growth resolution when kinetics accelerate.")
    ap.add_argument("--stop-after-first-2d-fire", action=argparse.BooleanOptionalAction, default=True,
                    help="Stop each v8 2-D comparison run after the first accepted crack advance. This makes the comparison match V1, which normally stops after its first advance. Use --no-stop-after-first-2d-fire for propagation-after-init diagnostics.")
    # Cleavage free-energy surface controls, passed through to both V1 and v8.
    ap.add_argument("--cleave-barrier-kind", choices=["classic", "exp_floor"], default=None)
    ap.add_argument("--cleave-G00-eV", type=float, default=None)
    ap.add_argument("--cleave-gT-eV-per-K", type=float, default=None)
    ap.add_argument("--cleave-sigc0-GPa", type=float, default=None)
    ap.add_argument("--cleave-sT-GPa-per-K", type=float, default=None)
    ap.add_argument("--cleave-exp-a", type=float, default=None)
    ap.add_argument("--cleave-exp-n", type=float, default=None)
    ap.add_argument("--cleave-floor-frac", type=float, default=None)
    ap.add_argument("--cleave-floor-min-eV", type=float, default=None)
    ap.add_argument("--cleave-floor-max-frac", type=float, default=None)
    ap.add_argument("--cleave-Tref-K", type=float, default=None)
    ap.add_argument("--cleave-exp-T-mode", choices=["linear", "mu_scale"], default=None)
    ap.add_argument("--cleave-mu-dlnmu-dT-per-K", type=float, default=None)
    ap.add_argument("--cleave-G0-mu-power", type=float, default=None)
    ap.add_argument("--cleave-sigc-mu-power", type=float, default=None)
    ap.add_argument("--cleave-S-hs-kB", type=float, default=None)
    ap.add_argument("--cleave-sigma-S-GPa", type=float, default=None)
    ap.add_argument("--cleave-S-hs-power", type=float, default=None)
    ap.add_argument("--cleave-S-hs-dT-per-K", type=float, default=None)
    ap.add_argument("--cleave-S-hs-Tref-K", type=float, default=None)
    # Plasticity EXP-floor mechanism controls, passed through to both V1 and v8.
    ap.add_argument("--exp-system", default=None, choices=["W[100]", "Ta[111]", "Al0.7CoCrFeNi-BCC", "Al0.7CoCrFeNi-FCC", "Cu"])
    ap.add_argument("--exp-a", type=float, default=None)
    ap.add_argument("--exp-n", type=float, default=None)
    ap.add_argument("--nu0-emit-pz", type=float, default=None)
    ap.add_argument("--nu0-peierls", type=float, default=None)
    ap.add_argument("--nu0-taylor", type=float, default=None)
    ap.add_argument("--emit-energy-scale", type=float, default=None)
    ap.add_argument("--emit-entropy-scale", type=float, default=None)
    ap.add_argument("--emit-stress-scale", type=float, default=None)
    ap.add_argument("--peierls-energy-scale", type=float, default=None)
    ap.add_argument("--peierls-entropy-scale", type=float, default=None)
    ap.add_argument("--peierls-stress-scale", type=float, default=None)
    ap.add_argument("--taylor-energy-scale", type=float, default=None)
    ap.add_argument("--taylor-entropy-scale", type=float, default=None)
    ap.add_argument("--taylor-stress-scale", type=float, default=None)
    ap.add_argument("--keep-existing", action="store_true")
    args = ap.parse_args()

    out = Path(args.out)
    if out.exists() and not args.keep_existing:
        shutil.rmtree(out)
    out.mkdir(parents=True, exist_ok=True)

    # Long-growth convenience: --target-dB limits expected cleavage renewals per
    # block. If the user supplies a desired crack-extension budget and da_phys,
    # translate it to an equivalent renewal-clock budget. This complements the
    # state audit inside sharp_front.py and keeps fast points from jumping too far
    # in one accepted cycle block.
    if args.target_da_per_block_um is not None:
        if args.da_phys is None or args.da_phys <= 0:
            raise ValueError("--target-da-per-block-um requires --da-phys so the driver can convert extension to expected renewal events")
        target_events = max(float(args.target_da_per_block_um) * 1e-6 / max(float(args.da_phys), 1e-300), 1e-12)
        if target_events < float(args.target_dB):
            print(f"tightening target-dB from {args.target_dB:g} to {target_events:g} based on target-da-per-block")
            args.target_dB = target_events

    common = [
        "--R", str(args.R),
        "--frequency-Hz", str(args.frequency_Hz),
        "--block-cycles", str(args.block_cycles),
        "--min-block-cycles", str(args.min_block_cycles),
        "--max-block-cycles", str(args.max_block_cycles),
        "--cycle-block-mode", str(args.cycle_block_mode),
        "--target-dB", str(args.target_dB),
        "--target-dN-store", str(args.target_dN_store),
        "--target-dN-emit", str(args.target_dN_emit),
        "--target-dN-mobile", str(args.target_dN_mobile),
        "--target-dN-escape", str(args.target_dN_escape),
        "--target-dN-peierls", str(args.target_dN_peierls),
        "--target-dN-taylor", str(args.target_dN_taylor),
        "--storage-model", str(args.storage_model),
        "--dN-cap", "inf",
    ]
    if args.storage_model == "fixed_fraction":
        common += ["--fixed-retained-fraction", str(args.fixed_retained_fraction)]
    cleave_passthrough = [
        ("--cleave-barrier-kind", args.cleave_barrier_kind),
        ("--cleave-G00-eV", args.cleave_G00_eV),
        ("--cleave-gT-eV-per-K", args.cleave_gT_eV_per_K),
        ("--cleave-sigc0-GPa", args.cleave_sigc0_GPa),
        ("--cleave-sT-GPa-per-K", args.cleave_sT_GPa_per_K),
        ("--cleave-exp-a", args.cleave_exp_a),
        ("--cleave-exp-n", args.cleave_exp_n),
        ("--cleave-floor-frac", args.cleave_floor_frac),
        ("--cleave-floor-min-eV", args.cleave_floor_min_eV),
        ("--cleave-floor-max-frac", args.cleave_floor_max_frac),
        ("--cleave-Tref-K", args.cleave_Tref_K),
        ("--cleave-exp-T-mode", args.cleave_exp_T_mode),
        ("--cleave-mu-dlnmu-dT-per-K", args.cleave_mu_dlnmu_dT_per_K),
        ("--cleave-G0-mu-power", args.cleave_G0_mu_power),
        ("--cleave-sigc-mu-power", args.cleave_sigc_mu_power),
        ("--cleave-S-hs-kB", args.cleave_S_hs_kB),
        ("--cleave-sigma-S-GPa", args.cleave_sigma_S_GPa),
        ("--cleave-S-hs-power", args.cleave_S_hs_power),
        ("--cleave-S-hs-dT-per-K", args.cleave_S_hs_dT_per_K),
        ("--cleave-S-hs-Tref-K", args.cleave_S_hs_Tref_K),
    ]
    for flag, val in cleave_passthrough:
        if val is not None:
            common += [flag, str(val)]

    plasticity_passthrough = [
        ("--exp-system", args.exp_system),
        ("--exp-a", args.exp_a),
        ("--exp-n", args.exp_n),
        ("--nu0-emit-pz", args.nu0_emit_pz),
        ("--nu0-peierls", args.nu0_peierls),
        ("--nu0-taylor", args.nu0_taylor),
        ("--emit-energy-scale", args.emit_energy_scale),
        ("--emit-entropy-scale", args.emit_entropy_scale),
        ("--emit-stress-scale", args.emit_stress_scale),
        ("--peierls-energy-scale", args.peierls_energy_scale),
        ("--peierls-entropy-scale", args.peierls_entropy_scale),
        ("--peierls-stress-scale", args.peierls_stress_scale),
        ("--taylor-energy-scale", args.taylor_energy_scale),
        ("--taylor-entropy-scale", args.taylor_entropy_scale),
        ("--taylor-stress-scale", args.taylor_stress_scale),
    ]
    for flag, val in plasticity_passthrough:
        if val is not None:
            common += [flag, str(val)]

    rows: List[Dict[str, object]] = []
    paired: List[Dict[str, object]] = []

    for K in args.Kmax_MPa_sqrt_m:
        label = (f"{K:g}").replace(".", "p")
        # V1 direct K-controlled reduction.
        v1out = out / f"v1_K{label}"
        cmd = [
            sys.executable, "-m", "arrhenius_fracture.fatigue_sharp_front",
            "--temperatures", str(args.T),
            "--Kmax-MPa-sqrt-m", str(K),
            "--cycles-max", str(args.cycles_max),
            "--max-blocks", str(args.blocks),
            "--no-plots",
            "--out", str(v1out),
        ] + common
        run_cmd(cmd, out / f"v1_K{label}.log")
        hist = pd.read_csv(v1out / f"T{int(args.T)}K" / "fatigue_v1_history.csv")
        r = hist.iloc[-1]
        v1_first_row = first_fire_row(hist)
        v1_first = first_fire_cycles_from_v1(hist)
        v1_K_initial = float(hist.iloc[0].Kmax_Pa_sqrt_m) / 1e6
        v1_K_first = (float(v1_first_row.Kmax_Pa_sqrt_m) / 1e6) if v1_first_row is not None else float("nan")
        v1_K_final = float(r.Kmax_Pa_sqrt_m) / 1e6
        rows.append({
            "model": "V1_1D",
            "target_Kmax_MPa_sqrtm": K,
            "actual_Kmax_MPa_sqrtm": v1_K_final,
            "K_initial_MPa_sqrtm": v1_K_initial,
            "K_first_fire_MPa_sqrtm": v1_K_first,
            "K_final_MPa_sqrtm": v1_K_final,
            "DeltaK_MPa_sqrtm": float(r.DeltaK_Pa_sqrt_m) / 1e6,
            "blocks_completed": int(len(hist)),
            "cycles_total": float(r.cycles_total),
            "cycles_to_first_fire": v1_first,
            "mu_cleave_per_cycle": float(r.mu_cleave_pred),
            "mu_emit_per_cycle": float(r.mu_emit),
            "dB_last_block": float(r.dB_block),
            "cycle_unlimited": _as_float(getattr(r, "cycle_unlimited", float("nan"))),
            "cycle_limiter": str(getattr(r, "cycle_limiter", "")),
            "B_final": float(r.B),
            "N_em_final": float(r.N_em),
            "dN_store_last_block": float(r.dN_store_block),
            "dN_mobile_last_block": float(r.dN_mobile_block),
            "dN_escape_last_block": float(r.dN_escape_block),
            "storage_fraction": float(r.storage_fraction),
            "n_adv_or_fire_total": int(hist["n_fire"].sum()),
            "a_adv_um": float(r.a_adv_m) * 1e6,
            "da_dN_m_per_cycle": (float(r.a_adv_m) / max(float(r.cycles_total), 1e-300)),
            "da_dN_upper_bound_m_per_cycle": (float(r.da_m) / max(float(r.cycles_total), 1e-300)) if int(hist["n_fire"].sum()) == 0 and hasattr(r, "da_m") else float("nan"),
            "sigma_tip_GPa": float(r.sigma_tip) / 1e9,
            "G_cleave_raw_eV": _as_float(getattr(r, "G_cleave_raw_eV", float("nan"))),
            "G_cleave_eff_eV": _as_float(getattr(r, "G_cleave_eff_eV", float("nan"))),
            "S_cleave_kB": _as_float(getattr(r, "S_cleave_kB", float("nan"))),
            "dGcleave_dsigma_eV_per_GPa": _as_float(getattr(r, "dGcleave_dsigma_eV_per_GPa", float("nan"))),
            "vstar_cleave_b3": _as_float(getattr(r, "vstar_cleave_b3", float("nan"))),
            "G_emit_eV": _as_float(getattr(r, "G_emit_eV", float("nan"))),
            "S_emit_kB": _as_float(getattr(r, "S_emit_kB", float("nan"))),
            "cyclic_plastic_work_J": "",
            "run_dir": str(v1out),
        })

        # v8 2-D held-amplitude sharp-front run.
        v8out = out / f"v8_2d_K{label}"
        k_cal = calibrate_2d_displacement(args, out, float(K), label)
        dU = float(k_cal["dU"])
        cmd = [
            sys.executable, "-m", "arrhenius_fracture.sharp_front",
            "--mode", "2d",
            "--fatigue-cycles",
            "--fatigue-hold-load",
            "--cycles-max", str(args.cycles_max),
            "--temperatures", str(args.T),
            "--steps", str(args.blocks),
            "--nx", str(args.nx), "--ny", str(args.ny),
            "--tip-h-fine", str(args.tip_h_fine),
            "--tip-ratio", str(args.tip_ratio),
            "--n-stagger", "1",
            "--crystal-aniso",
            "--crystal-branch",
            "--j-decomposition", "cluster",
            "--sigma-cap-GPa", "0",
            "--dU", str(dU),
            "--save-snapshots", str(args.save_snapshots),
            "--snapshot-cols", str(args.snapshot_cols),
            "--min-global-forward", str(args.min_global_forward),
            "--cyclic-mechanics-phases", str(args.cyclic_mechanics_phases),
            "--target-crack-extension-um", str(args.target_crack_extension_um),
            "--snapshot-by-crack-extension-um", str(args.snapshot_by_crack_extension_um),
            "--max-da-per-block-um", str(args.max_da_per_block_um),
            "--crack-backend", str(args.crack_backend),
            "--czm-penalty-normal", str(args.czm_penalty_normal),
            "--czm-penalty-tangent", str(args.czm_penalty_tangent),
            "--czm-event-damage", str(args.czm_event_damage),
            "--czm-max-angle-error-deg", str(args.czm_max_angle_error_deg),
            "--out", str(v8out),
        ]
        if args.da_phys is not None:
            cmd += ["--da-phys", str(args.da_phys)]
        if not args.make_2d_plots:
            cmd += ["--no-plots"]
        if args.stop_after_first_2d_fire:
            cmd += ["--stop-after-first-fire"]
        cmd += common
        run_cmd(cmd, out / f"v8_2d_K{label}.log")
        steps = pd.read_csv(v8out / f"steps_{int(args.T):04d}K.csv")
        s = steps.iloc[-1]
        v8_first_row = first_fire_row(steps)
        v8_first = first_fire_cycles_from_2d(steps)
        v8_K_initial = float(steps.iloc[0].KJ_Pa_sqrtm) / 1e6
        v8_K_first = (float(v8_first_row.KJ_Pa_sqrtm) / 1e6) if v8_first_row is not None else float("nan")
        v8_K_final = float(s.KJ_Pa_sqrtm) / 1e6
        rows.append({
            "model": "v8_2D",
            "target_Kmax_MPa_sqrtm": K,
            "actual_Kmax_MPa_sqrtm": v8_K_final,
            "K_initial_MPa_sqrtm": v8_K_initial,
            "K_first_fire_MPa_sqrtm": v8_K_first,
            "K_final_MPa_sqrtm": v8_K_final,
            "dU_calibrated_m": dU,
            "K_calibration_MPa_sqrtm": float(k_cal.get("K_calibrated", float("nan"))),
            "K_calibration_rel_error": float(k_cal.get("rel_error", float("nan"))),
            "K_calibration_iters": int(k_cal.get("iters", 0)),
            "DeltaK_MPa_sqrtm": (1.0 - args.R) * float(s.KJ_Pa_sqrtm) / 1e6,
            "blocks_completed": int(len(steps)),
            "cycles_total": float(steps["fatigue_cycles"].sum()),
            "cycles_to_first_fire": v8_first,
            "mu_cleave_per_cycle": float(s.mu_cleave_pred_per_cycle),
            "mu_emit_per_cycle": float(s.mu_emit_per_cycle),
            "dB_last_block": float(s.dB_block),
            "cycle_unlimited": _as_float(getattr(s, "cycle_unlimited", float("nan"))),
            "cycle_limiter_code": _as_float(getattr(s, "cycle_limiter_code", float("nan"))),
            "B_final": float(s.B),
            "N_em_final": float(s.N_em),
            "dN_store_last_block": float(s.dN_store_block),
            "dN_mobile_last_block": float(s.dN_mobile_block),
            "dN_escape_last_block": float(s.dN_escape_block),
            "storage_fraction": float(s.storage_fraction),
            "n_adv_or_fire_total": int(steps["n_fire"].sum()),
            "a_adv_um": (_as_float(getattr(s, "crack_extension_m", float("nan")), (float(s.a_tip_m) - 5e-4)) * 1e6),
            "da_dN_m_per_cycle": (_as_float(getattr(s, "crack_extension_m", float("nan")), (float(s.a_tip_m) - 5e-4)) / max(float(steps["fatigue_cycles"].sum()), 1e-300)),
            "da_dN_upper_bound_m_per_cycle": ((float(args.da_phys) if args.da_phys is not None else 5.0e-6) / max(float(steps["fatigue_cycles"].sum()), 1e-300)) if int(steps["n_fire"].sum()) == 0 else float("nan"),
            "sigma_tip_GPa": float(s.sigma_tip_Pa) / 1e9,
            "G_cleave_raw_eV": _as_float(getattr(s, "G_cleave_raw_eV", float("nan"))),
            "G_cleave_eff_eV": _as_float(getattr(s, "G_cleave_eff_eV", float("nan"))),
            "S_cleave_kB": _as_float(getattr(s, "S_cleave_kB", float("nan"))),
            "dGcleave_dsigma_eV_per_GPa": _as_float(getattr(s, "dGcleave_dsigma_eV_per_GPa", float("nan"))),
            "vstar_cleave_b3": _as_float(getattr(s, "vstar_cleave_b3", float("nan"))),
            "cyclic_plastic_work_J": float(s.cyclic_plastic_work_J),
            "run_dir": str(v8out),
        })

        v1 = rows[-2]
        v8 = rows[-1]
        paired.append({
            "target_Kmax_MPa_sqrtm": K,
            "v8_K_initial_MPa_sqrtm": v8["K_initial_MPa_sqrtm"],
            "v8_K_first_fire_MPa_sqrtm": v8["K_first_fire_MPa_sqrtm"],
            "v8_K_final_MPa_sqrtm": v8["K_final_MPa_sqrtm"],
            "v8_dU_calibrated_m": v8.get("dU_calibrated_m", float("nan")),
            "v8_K_calibration_rel_error": v8.get("K_calibration_rel_error", float("nan")),
            "K_initial_ratio_v8_over_v1": v8["K_initial_MPa_sqrtm"] / max(v1["K_initial_MPa_sqrtm"], 1e-300),
            "K_final_ratio_v8_over_v1": v8["K_final_MPa_sqrtm"] / max(v1["K_final_MPa_sqrtm"], 1e-300),
            "mu_cleave_ratio_v8_over_v1": v8["mu_cleave_per_cycle"] / max(v1["mu_cleave_per_cycle"], 1e-300),
            "mu_emit_ratio_v8_over_v1": v8["mu_emit_per_cycle"] / max(v1["mu_emit_per_cycle"], 1e-300),
            "G_cleave_eff_v8_minus_v1_eV": v8.get("G_cleave_eff_eV", float("nan")) - v1.get("G_cleave_eff_eV", float("nan")),
            "S_cleave_v8_minus_v1_kB": v8.get("S_cleave_kB", float("nan")) - v1.get("S_cleave_kB", float("nan")),
            "sigma_tip_ratio_v8_over_v1": v8["sigma_tip_GPa"] / max(v1["sigma_tip_GPa"], 1e-300),
            "cycles_to_first_fire_V1": v1["cycles_to_first_fire"],
            "cycles_to_first_fire_v8": v8["cycles_to_first_fire"],
            "final_B_V1": v1["B_final"],
            "final_B_v8": v8["B_final"],
        })

    df = pd.DataFrame(rows)
    pr = pd.DataFrame(paired)
    df.to_csv(out / "compare_summary.csv", index=False)
    pr.to_csv(out / "paired_ratios.csv", index=False)
    with (out / "comparison_settings.json").open("w") as f:
        json.dump(vars(args), f, indent=2)

    print("\n=== comparison summary ===")
    print(df[["model", "target_Kmax_MPa_sqrtm", "K_initial_MPa_sqrtm", "K_first_fire_MPa_sqrtm", "K_final_MPa_sqrtm", "blocks_completed", "cycles_total", "cycles_to_first_fire", "mu_cleave_per_cycle", "B_final", "n_adv_or_fire_total", "a_adv_um"]].to_string(index=False))
    print("\n=== paired ratios ===")
    print(pr.to_string(index=False))
    print(f"\nwrote {out / 'compare_summary.csv'}")


if __name__ == "__main__":
    main()
