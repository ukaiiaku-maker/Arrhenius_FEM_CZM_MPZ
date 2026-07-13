#!/usr/bin/env python3
"""Run seeded 500 K FEM/CZM R-curve realizations for the four EXP-floor classes."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import math
from pathlib import Path
import shlex
import subprocess
from types import SimpleNamespace

import numpy as np
import pandas as pd

import run_four_class_exp_floor_czm_500um_sweep as base

CLASSES_DEFAULT = ["ceramic", "peak", "weakT", "DBTT"]
KB = 1.380649e-23


def parse_list(text: str, cast=str):
    return [cast(x) for x in text.replace(",", " ").split() if x]


def make_base_args(args) -> SimpleNamespace:
    return SimpleNamespace(
        nx=args.nx,
        ny=args.ny,
        tip_h_fine=args.tip_h_fine,
        tip_ratio=args.tip_ratio,
        dU=args.dU,
        dt=args.dt,
        long_steps=args.long_steps,
        n_stagger=args.n_stagger,
        print_every=args.print_every,
        target_ext_um=args.target_ext_um,
        crystal_material=args.crystal_material,
        cleave_gamma_aniso=args.cleave_gamma_aniso,
        theta=args.theta,
        save_snapshots=args.save_snapshots,
        snapshot_cols=args.snapshot_cols,
        snapshot_by_ext_um=args.snapshot_by_ext_um,
    )


def find_steps_file(case_dir: Path, T: int) -> Path | None:
    exact = case_dir / f"steps_{T:04d}K.csv"
    if exact.exists():
        return exact
    files = sorted(case_dir.glob("steps_*K.csv"))
    return files[0] if files else None


def enrich_curve_with_barrier_statistics(
    case_dir: Path,
    T: int,
    rc: pd.DataFrame,
    burgers_vector_m: float,
) -> pd.DataFrame:
    """Attach activation-volume and local Gumbel-scale diagnostics to events."""
    if rc.empty or "step" not in rc.columns:
        return rc
    sf = find_steps_file(case_dir, T)
    if sf is None:
        return rc
    st = pd.read_csv(sf)
    if "step" not in st.columns:
        return rc

    diag_cols = [
        "vstar_cleave_b3",
        "dGcleave_dsigma_eV_per_GPa",
        "sigma_cleave_eff_Pa",
        "S_cleave_kB",
        "G_cleave_raw_eV",
        "G_cleave_eff_eV",
    ]
    keep = ["step", *[c for c in diag_cols if c in st.columns]]
    diag = st[keep].copy()
    diag["step"] = pd.to_numeric(diag["step"], errors="coerce")
    diag = diag.dropna(subset=["step"]).drop_duplicates("step", keep="last")

    out = rc.merge(diag, on="step", how="left", suffixes=("", "_step"))
    for c in diag_cols:
        step_name = c + "_step"
        if step_name in out.columns:
            if c in out.columns:
                out[c] = pd.to_numeric(out[c], errors="coerce").fillna(
                    pd.to_numeric(out[step_name], errors="coerce")
                )
                out = out.drop(columns=[step_name])
            else:
                out = out.rename(columns={step_name: c})

    if "sigma_cleave_eff_Pa" in out.columns:
        out["sigma_cleave_eff_GPa"] = (
            pd.to_numeric(out["sigma_cleave_eff_Pa"], errors="coerce") * 1e-9
        )

    if "vstar_cleave_b3" in out.columns:
        v_b3 = pd.to_numeric(out["vstar_cleave_b3"], errors="coerce").to_numpy(float)
        v_m3 = v_b3 * float(burgers_vector_m) ** 3
        out["activation_volume_m3"] = v_m3
        out["activation_volume_nm3"] = v_m3 * 1e27

        beta_sigma = np.full(len(out), np.nan)
        valid_v = np.isfinite(v_m3) & (v_m3 > 0.0)
        beta_sigma[valid_v] = KB * float(T) / v_m3[valid_v]
        out["gumbel_beta_sigma_GPa"] = beta_sigma * 1e-9
        out["gumbel_sd_sigma_GPa"] = (math.pi / math.sqrt(6.0)) * beta_sigma * 1e-9

        K_pa = pd.to_numeric(out["KJ_MPa_sqrt_m"], errors="coerce").to_numpy(float) * 1e6
        if "sigma_cleave_eff_Pa" in out.columns:
            sigma_eff = pd.to_numeric(
                out["sigma_cleave_eff_Pa"], errors="coerce"
            ).to_numpy(float)
        elif "sigma_tip_GPa" in out.columns:
            sigma_eff = pd.to_numeric(
                out["sigma_tip_GPa"], errors="coerce"
            ).to_numpy(float) * 1e9
        else:
            sigma_eff = np.full(len(out), np.nan)

        # Simplest local mapping: sigma_eff is proportional to K over the small
        # load interval in which one hazard event is selected. This gives
        # d sigma_eff / dK ~= sigma_eff/K.
        dsigma_dK = np.full(len(out), np.nan)
        good_map = np.isfinite(K_pa) & (K_pa > 0.0) & np.isfinite(sigma_eff) & (sigma_eff > 0.0)
        dsigma_dK[good_map] = sigma_eff[good_map] / K_pa[good_map]

        beta_K_pa = np.full(len(out), np.nan)
        good = valid_v & np.isfinite(dsigma_dK) & (dsigma_dK > 0.0)
        beta_K_pa[good] = KB * float(T) / (v_m3[good] * dsigma_dK[good])
        out["gumbel_beta_K_MPa_sqrt_m"] = beta_K_pa / 1e6
        out["gumbel_sd_K_MPa_sqrt_m"] = (
            math.pi / math.sqrt(6.0)
        ) * beta_K_pa / 1e6
        out["local_dsigma_dK_inv_sqrt_m"] = dsigma_dK

    return out


def annotate_curve(
    case_dir: Path,
    klass: str,
    T: int,
    target: float,
    replicate: int,
    solver_seed: int,
    burgers_vector_m: float,
) -> pd.DataFrame:
    rc = base.process_case_r_curve(case_dir, klass, T, target)
    if rc.empty:
        return rc
    rc = enrich_curve_with_barrier_statistics(case_dir, T, rc, burgers_vector_m)
    for column in ["replicate", "solver_seed", "seed"]:
        if column in rc.columns:
            rc = rc.drop(columns=[column])
    rc.insert(2, "replicate", int(replicate))
    rc.insert(3, "solver_seed", int(solver_seed))
    rc.to_csv(case_dir / "R_curve_event_sampled.csv", index=False)
    return rc


def run_one(
    py: str,
    wrapper: Path,
    row: pd.Series,
    klass: str,
    replicate: int,
    solver_seed: int,
    root: Path,
    args,
    base_args,
) -> dict:
    T = int(args.temperature)
    rep_dir = root / klass / f"replicate_{replicate:02d}_seed{solver_seed}"
    case_dir = rep_dir / f"T{T}_th{args.theta:g}"
    case_dir.mkdir(parents=True, exist_ok=True)

    complete, ext = base.completion_status(case_dir, args.target_ext_um)
    if complete and not args.force:
        rc = annotate_curve(
            case_dir, klass, T, args.target_ext_um, replicate, solver_seed,
            args.burgers_vector_m,
        )
        print(
            f"SKIP {klass:8s} rep={replicate} seed={solver_seed}: "
            f"extension={ext:.1f} um"
        )
        return {
            "class": klass,
            "replicate": replicate,
            "solver_seed": solver_seed,
            "T_K": T,
            "status": "skipped_complete",
            "extension_um": ext,
            "n_r_curve_points": len(rc),
            "case_dir": str(case_dir.relative_to(root)),
        }

    cmd = base.build_command(py, row, T, case_dir, base_args)
    if len(cmd) < 3 or cmd[1:3] != ["-m", "arrhenius_fracture.sharp_front"]:
        raise RuntimeError(f"Unexpected base command prefix: {cmd[:4]}")
    cmd = [py, str(wrapper), "--solver-seed", str(solver_seed), *cmd[3:]]

    (case_dir / "command.txt").write_text(shlex.join(cmd) + "\n")
    payload = row.to_dict()
    payload.update({
        "target_class": klass,
        "replicate": int(replicate),
        "solver_seed": int(solver_seed),
        "temperature_K": T,
        "theta_deg": float(args.theta),
        "target_extension_um": float(args.target_ext_um),
        "random_process": "seeded nonhomogeneous Poisson cleavage renewal",
        "latent_first_passage_action": "H_n=-ln(U_n), U_n~Uniform(0,1)",
        "threshold_values_are_model_realizations_not_parameters": True,
        "common_random_number_seed_reused_across_classes": True,
        "mesh_seed_fixed": True,
        "branching_enabled": False,
        "burgers_vector_m_for_activation_volume_analysis": float(args.burgers_vector_m),
        "emission_G00_effective_eV": 0.75 * float(row.exp_G00_eV),
        "emission_gT_effective_eV_per_K": 0.75 * float(row.exp_gT_eV_per_K),
    })
    (case_dir / "resolved_parameters.json").write_text(
        json.dumps(payload, indent=2, default=str)
    )

    print(
        f"START {klass:8s} rep={replicate} seed={solver_seed} -> {case_dir}"
    )
    with (case_dir / "run.log").open("w") as log:
        cp = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, text=True)

    complete, ext = base.completion_status(case_dir, args.target_ext_um)
    rc = annotate_curve(
        case_dir, klass, T, args.target_ext_um, replicate, solver_seed,
        args.burgers_vector_m,
    )
    if cp.returncode == 0 and complete:
        (case_dir / ".long_growth_complete").touch()
        status = "complete"
        print(
            f"DONE  {klass:8s} rep={replicate} seed={solver_seed}: "
            f"extension={ext:.1f} um"
        )
    elif cp.returncode == 0:
        status = "incomplete"
        print(
            f"INCOMPLETE {klass:8s} rep={replicate} seed={solver_seed}: "
            f"extension={ext}"
        )
    else:
        status = "failed"
        print(
            f"FAILED {klass:8s} rep={replicate} seed={solver_seed}: rc={cp.returncode}"
        )

    return {
        "class": klass,
        "replicate": replicate,
        "solver_seed": solver_seed,
        "T_K": T,
        "status": status,
        "returncode": cp.returncode,
        "extension_um": ext,
        "n_r_curve_points": len(rc),
        "case_dir": str(case_dir.relative_to(root)),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--parameters", default="four_class_exp_floor_exact_model_inputs.csv")
    ap.add_argument(
        "--outroot",
        default="runs/four_class_exp_floor_CZM_500K_5rep_1000um_theta45",
    )
    ap.add_argument("--classes", default=" ".join(CLASSES_DEFAULT))
    ap.add_argument("--temperature", type=int, default=500)
    ap.add_argument("--solver-seeds", default="1101 1102 1103 1104 1105")
    ap.add_argument("--theta", type=float, default=45.0)
    ap.add_argument("--target-ext-um", type=float, default=1000.0)
    ap.add_argument("--long-steps", type=int, default=50000)
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
    ap.add_argument("--snapshot-by-ext-um", type=float, default=100.0)
    ap.add_argument("--grid-step-um", type=float, default=5.0)
    ap.add_argument("--burgers-vector-m", type=float, default=2.74e-10)
    args = ap.parse_args()

    classes = parse_list(args.classes)
    seeds = parse_list(args.solver_seeds, int)
    if not seeds:
        raise SystemExit("At least one solver seed is required")

    params = base.load_parameters(Path(args.parameters))
    missing = [c for c in classes if c not in params.index]
    if missing:
        raise SystemExit(f"classes absent from parameter table: {missing}")

    py = base.resolve_python(args)
    print(f"python: {py}")
    base.preflight(py)

    wrapper = Path(__file__).with_name("run_seeded_sharp_front.py").resolve()
    plotter = Path(__file__).with_name(
        "plot_four_class_500K_seeded_rcurves.py"
    ).resolve()
    for p in [wrapper, plotter, Path("run_four_class_exp_floor_czm_500um_sweep.py")]:
        if not p.exists():
            raise SystemExit(f"required file not found: {p}")

    check = subprocess.run(
        [
            py,
            "-c",
            (
                "import inspect; from arrhenius_fracture.sharp_front import FrontEngine; "
                "s=inspect.getsource(FrontEngine.step); "
                "assert 'max_advances_per_step' in s, "
                "'adaptive-CZM one-event FrontEngine is required'; "
                "print('seeded first-passage preflight OK')"
            ),
        ],
        text=True,
        capture_output=True,
    )
    if check.returncode != 0:
        raise SystemExit(check.stdout + check.stderr)
    print(check.stdout.strip())

    root = Path(args.outroot)
    root.mkdir(parents=True, exist_ok=True)
    (root / "four_class_exp_floor_exact_model_inputs.csv").write_bytes(
        Path(args.parameters).read_bytes()
    )
    config = vars(args).copy()
    config.update({
        "resolved_classes": classes,
        "resolved_solver_seeds": seeds,
        "n_replicates_per_class": len(seeds),
        "branching_enabled": False,
        "physical_first_passage_scatter": True,
        "latent_action_distribution": "Exp(mean=1)",
        "common_random_numbers_across_classes": True,
        "mesh_seed_policy": "fixed; not part of physical replicate scatter",
        "activation_volume_definition": "V_eff=-d(DeltaG)/d(sigma)",
        "local_gumbel_scale_sigma": "beta_sigma=kBT/V_eff",
        "local_gumbel_scale_K": "beta_K=kBT/[V_eff*(d sigma_eff/dK)]",
    })
    (root / "replicate_campaign_config.json").write_text(
        json.dumps(config, indent=2)
    )

    base_args = make_base_args(args)
    tasks = [
        (klass, i + 1, seed)
        for klass in classes
        for i, seed in enumerate(seeds)
    ]
    results = []
    if args.max_jobs <= 1:
        for klass, rep, seed in tasks:
            results.append(run_one(
                py, wrapper, params.loc[klass], klass, rep, seed,
                root, args, base_args,
            ))
    else:
        with ThreadPoolExecutor(max_workers=args.max_jobs) as ex:
            futs = {
                ex.submit(
                    run_one, py, wrapper, params.loc[klass], klass, rep, seed,
                    root, args, base_args,
                ): (klass, rep, seed)
                for klass, rep, seed in tasks
            }
            for fut in as_completed(futs):
                results.append(fut.result())

    status = pd.DataFrame(results).sort_values(["class", "replicate"])
    status.to_csv(root / "replicate_campaign_status.csv", index=False)

    plot_cmd = [
        py,
        str(plotter),
        "--root",
        str(root),
        "--classes",
        " ".join(classes),
        "--temperature",
        str(args.temperature),
        "--target-ext-um",
        str(args.target_ext_um),
        "--grid-step-um",
        str(args.grid_step_um),
        "--burgers-vector-m",
        str(args.burgers_vector_m),
    ]
    plot_cp = subprocess.run(plot_cmd, text=True)

    bad = status[~status.status.isin(["complete", "skipped_complete"])]
    print(f"WROTE {root / 'replicate_campaign_status.csv'}")
    if len(bad):
        print(f"WARNING: {len(bad)} cases failed or were incomplete")
    if plot_cp.returncode != 0:
        raise SystemExit(plot_cp.returncode)
    if len(bad):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
