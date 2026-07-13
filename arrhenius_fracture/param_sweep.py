"""
Higher-level diagnostic parameter sweeps for the emergent Arrhenius fracture model.

This module is intentionally separate from sweep.py.  It is designed for overnight
screening of a few hundred controlled cases after the single-case diagnostics have
identified plausible physics windows.

Typical use:

python -m arrhenius_fracture.param_sweep \
  --name W100_Gc63_barrier_memory \
  --temperatures 425 450 475 500 525 550 \
  --Gc-list 63 \
  --barrier-systems "W[100]" \
  --barrier-scales 0.1 0.2 0.3 \
  --entropy-scales 0.5 1.0 1.5 \
  --stress-scales 0.8 1.0 1.25 \
  --flow-epsdot-refs 1e-5 1e-4 1e-3 \
  --memory-gains 0.25 1.0 \
  --M-max-values 3 4 5 \
  --max-cases 250
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import os
import random
import time
from copy import deepcopy
from dataclasses import asdict
from typing import Any, Dict, Iterable, List

import numpy as np

from .config import make_emergent_config
from .main import (
    run_simulation,
    _load_exp_floor_plastic_barrier,
    _available_exp_floor_systems,
    _safe_system_tag,
    _safe_float_tag,
)
from .diagnostics import history_summary


def _as_list(x):
    if x is None:
        return []
    return list(x)


def _expand_systems(json_path: str, requested: List[str]) -> List[str]:
    out: List[str] = []
    for r in requested:
        rl = str(r).lower()
        if rl in ("all", "all_including_si"):
            out.extend(_available_exp_floor_systems(json_path, exclude_si=False))
        elif rl in ("all_non_si", "all_nonsi", "all_metals"):
            out.extend(_available_exp_floor_systems(json_path, exclude_si=True))
        else:
            out.append(r)
    # preserve order, remove duplicates
    seen = set()
    uniq = []
    for s in out:
        if s not in seen:
            uniq.append(s); seen.add(s)
    return uniq


def _case_tag(case: Dict[str, Any]) -> str:
    parts = [
        _safe_system_tag(case["system"]),
        f"Gc{_safe_float_tag(case['Gc'])}",
        f"Eb{_safe_float_tag(case['barrier_scale'])}",
        f"S{_safe_float_tag(case['entropy_scale'])}",
        f"sig{_safe_float_tag(case['stress_scale'])}",
        f"edot{_safe_float_tag(case['flow_epsdot_ref'])}",
        f"mg{_safe_float_tag(case['memory_gain'])}",
        f"M{_safe_float_tag(case['M_max'])}",
    ]
    if case.get("amp_max") is not None:
        parts.append(f"amp{_safe_float_tag(case['amp_max'])}")
    return "__".join(parts)


def _score_case(temp_rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Score a temperature series for peak-like or transition-like behavior.

    The score is not a physics law; it is a triage metric to rank cases for
    dense reruns.  It rewards a local K_selected peak and/or a clean
    brittle→mixed/soft→brittle/soft transition while penalizing invalid runaway.
    """
    if not temp_rows:
        return {
            "score": -1e9,
            "n_temps": 0,
            "n_invalid": 0,
            "n_soft": 0,
            "n_mixed": 0,
            "K_peak_MPa": np.nan,
            "K_peak_T": np.nan,
            "peak_excess_MPa": np.nan,
            "invalid_fraction": 1.0,
            "mode_sequence": "",
        }

    Ts = np.array([float(r.get("T", np.nan)) for r in temp_rows], dtype=float)
    Ks = np.array([float(r.get("KJ_selected_MPa_sqrt_m", np.nan)) for r in temp_rows], dtype=float)
    modes = [str(r.get("failure_mode", "")) for r in temp_rows]
    invalid = np.array([(m.startswith("invalid") or "invalid" in m) for m in modes], dtype=bool)
    soft = np.array([("soft_tearing" in m and not m.startswith("invalid")) for m in modes], dtype=bool)
    mixed = np.array([("mixed" in m and not m.startswith("invalid")) for m in modes], dtype=bool)
    brittle = np.array([("brittle" in m and not m.startswith("invalid")) for m in modes], dtype=bool)

    validK = np.isfinite(Ks) & (~invalid) & (Ks > 0)
    if np.any(validK):
        idx = int(np.nanargmax(np.where(validK, Ks, -np.inf)))
        Kpeak = float(Ks[idx]); Tpeak = float(Ts[idx])
        edge_vals = []
        if np.isfinite(Ks[0]) and not invalid[0]: edge_vals.append(float(Ks[0]))
        if np.isfinite(Ks[-1]) and not invalid[-1]: edge_vals.append(float(Ks[-1]))
        edge_ref = float(np.nanmean(edge_vals)) if edge_vals else float(np.nanmedian(Ks[validK]))
        peak_excess = Kpeak - edge_ref
    else:
        Kpeak = float("nan"); Tpeak = float("nan"); peak_excess = float("nan")

    n_soft = int(np.sum(soft)); n_mixed = int(np.sum(mixed)); n_invalid = int(np.sum(invalid)); n_brittle = int(np.sum(brittle))
    invalid_frac = n_invalid / max(len(temp_rows), 1)

    # Transition score: both brittle/reference points and soft/mixed points are useful.
    transition_bonus = 2.0 * min(n_brittle, n_soft + n_mixed)
    soft_bonus = 1.0 * n_soft + 0.4 * n_mixed
    peak_bonus = max(0.0, peak_excess if np.isfinite(peak_excess) else 0.0)
    invalid_penalty = 4.0 * n_invalid
    # Avoid overrewarding cases that are just flat brittle cracks.
    all_brittle_penalty = 3.0 if (n_brittle == len(temp_rows)) else 0.0
    score = peak_bonus + transition_bonus + soft_bonus - invalid_penalty - all_brittle_penalty

    return {
        "score": float(score),
        "n_temps": len(temp_rows),
        "n_invalid": n_invalid,
        "n_soft": n_soft,
        "n_mixed": n_mixed,
        "n_brittle": n_brittle,
        "K_peak_MPa": Kpeak,
        "K_peak_T": Tpeak,
        "peak_excess_MPa": float(peak_excess) if np.isfinite(peak_excess) else float("nan"),
        "invalid_fraction": float(invalid_frac),
        "mode_sequence": ";".join(modes),
    }


def _write_csv(path: str, rows: List[Dict[str, Any]]):
    if not rows:
        return
    os.makedirs(os.path.dirname(path), exist_ok=True)
    keys: List[str] = []
    for r in rows:
        for k in r.keys():
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        w.writeheader(); w.writerows(rows)


def _read_csv(path: str) -> List[Dict[str, Any]]:
    if not os.path.exists(path):
        return []
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def _to_float(x, default=float("nan")):
    try:
        if x is None or x == "":
            return default
        return float(x)
    except Exception:
        return default


def _planned_cases_path(output_dir: str) -> str:
    return os.path.join(output_dir, "planned_cases.csv")


def _load_planned_cases(output_dir: str) -> List[Dict[str, Any]]:
    """Load cases from a previous planned_cases.csv so resume preserves shuffled order."""
    rows = _read_csv(_planned_cases_path(output_dir))
    cases: List[Dict[str, Any]] = []
    for r in rows:
        cases.append({
            "system": r.get("system", "W[100]"),
            "Gc": _to_float(r.get("Gc")),
            "barrier_scale": _to_float(r.get("barrier_scale")),
            "entropy_scale": _to_float(r.get("entropy_scale")),
            "stress_scale": _to_float(r.get("stress_scale")),
            "flow_epsdot_ref": _to_float(r.get("flow_epsdot_ref")),
            "memory_gain": _to_float(r.get("memory_gain")),
            "M_max": _to_float(r.get("M_max")),
            "amp_max": _to_float(r.get("amp_max"), 5.0),
            "shield_max": _to_float(r.get("shield_max"), 0.85),
            "blunt_per_work": _to_float(r.get("blunt_per_work"), 0.35),
            "sharpen_per_damage": _to_float(r.get("sharpen_per_damage"), 0.35),
        })
    return cases


def _case_output_dir(args, case: Dict[str, Any], case_index: int) -> str:
    return os.path.join(args.output_dir, f"case_{case_index:04d}__{_case_tag(case)}")


def _existing_case_temperature_rows(args, case: Dict[str, Any], case_index: int) -> List[Dict[str, Any]]:
    """Return saved per-temperature rows for a completed case, decorated with case parameters."""
    case_dir = _case_output_dir(args, case, case_index)
    rows = _read_csv(os.path.join(case_dir, "summary_by_temperature.csv"))
    out: List[Dict[str, Any]] = []
    for r in rows:
        rr = dict(r)
        rr.update(case)
        rr["case_index"] = case_index
        rr["case_tag"] = _case_tag(case)
        # Normalize T key for scorer and downstream combined table.
        if "T" not in rr:
            if "T_K" in rr:
                rr["T"] = rr["T_K"]
            elif "temperature" in rr:
                rr["T"] = rr["temperature"]
        out.append(rr)
    return out


def _case_is_complete(args, case: Dict[str, Any], case_index: int) -> bool:
    rows = _existing_case_temperature_rows(args, case, case_index)
    if len(rows) < len(args.temperatures):
        return False
    have = set()
    for r in rows:
        t = _to_float(r.get("T", r.get("T_K", r.get("temperature"))))
        if np.isfinite(t):
            have.add(round(float(t), 6))
    want = {round(float(t), 6) for t in args.temperatures}
    return want.issubset(have)


def _case_row_from_existing(args, case: Dict[str, Any], case_index: int) -> Dict[str, Any]:
    temp_rows = _existing_case_temperature_rows(args, case, case_index)
    score = _score_case(temp_rows)
    case_row = dict(case)
    case_row.update(score)
    case_row["case_index"] = case_index
    case_row["case_tag"] = _case_tag(case)
    case_row["output_dir"] = _case_output_dir(args, case, case_index)
    case_row["resume_status"] = "completed_existing"
    return case_row


def build_cases(args) -> List[Dict[str, Any]]:
    systems = _expand_systems(args.plastic_barrier_json, args.barrier_systems)
    grid = itertools.product(
        systems,
        args.Gc_list,
        args.barrier_scales,
        args.entropy_scales,
        args.stress_scales,
        args.flow_epsdot_refs,
        args.memory_gains,
        args.M_max_values,
        args.amp_max_values,
        args.shield_max_values,
        args.blunt_work_values,
        args.sharpen_damage_values,
    )
    cases = []
    for (system, Gc, bscale, escale, sscale, edot, mgain, Mmax,
         ampmax, shieldmax, bluntwork, sharpdamage) in grid:
        cases.append({
            "system": system,
            "Gc": float(Gc),
            "barrier_scale": float(bscale),
            "entropy_scale": float(escale),
            "stress_scale": float(sscale),
            "flow_epsdot_ref": float(edot),
            "memory_gain": float(mgain),
            "M_max": float(Mmax),
            "amp_max": float(ampmax),
            "shield_max": float(shieldmax),
            "blunt_per_work": float(bluntwork),
            "sharpen_per_damage": float(sharpdamage),
        })
    if args.shuffle or args.max_cases is not None:
        rng = random.Random(args.seed)
        rng.shuffle(cases)
    if args.max_cases is not None and args.max_cases > 0:
        cases = cases[:args.max_cases]
    return cases


def run_case(args, case: Dict[str, Any], case_index: int, n_cases: int) -> Dict[str, Any]:
    cfg = make_emergent_config()
    cfg.T_list = list(args.temperatures)
    cfg.loading.n_steps = int(args.steps)
    cfg.loading.dU_top = float(args.dU_top)
    if args.nx is not None:
        cfg.mesh.nx = int(args.nx)
    if args.ny is not None:
        cfg.mesh.ny = int(args.ny)
    if args.mesh_jitter is not None:
        cfg.mesh.jitter = float(args.mesh_jitter)
    if args.ell_factor is not None:
        cfg.mesh.ell_factor = float(args.ell_factor)
    if args.ell is not None:
        cfg.mesh.ell_absolute_m = float(args.ell)
    cfg.diagnostics.save_every = int(args.save_every)
    cfg.diagnostics.make_plots = bool(args.make_plots)
    cfg.diagnostics.save_fields = bool(args.save_fields)
    cfg.diagnostics.save_field_pngs = bool(args.save_field_pngs)
    cfg.auto_stop.enabled = not args.no_auto_stop
    cfg.stop_on_invalid = bool(args.stop_on_invalid)
    cfg.invalid_min_step = int(args.invalid_min_step)
    cfg.invalid_wp_wext_pct = float(args.invalid_wp_wext_pct)
    if args.allow_invalid_soft_tearing:
        cfg.invalid_wp_wext_pct = np.inf

    cfg.dislocations.plastic_update_mode = "flow_stress"
    cfg.dislocations.flow_epsdot_ref = float(case["flow_epsdot_ref"])
    cfg.dislocations.freeze_rho = bool(args.freeze_rho)
    cfg.dislocations.enable_plasticity = not args.disable_plasticity
    cfg.dislocations.max_plastic_strain_increment = float(args.max_plastic_strain_increment)
    cfg.dislocations.max_rho_relative_increment = float(args.max_rho_relative_increment)
    cfg.dislocations.rho_cap = float(args.rho_cap)
    cfg.dislocations.dot_ep_max = float(args.dot_ep_max)

    cfg.phase_field.Gc0_athermal = float(case["Gc"])
    if args.disable_wp_gc_coupling:
        cfg.phase_field.wp_gc_coupling_mode = 'off'
        cfg.phase_field.plastic_work_to_Gc_efficiency = 0.0
        cfg.phase_field.Gc_local_cap_factor = 1.0
    else:
        cfg.phase_field.wp_gc_coupling_mode = str(args.wp_gc_coupling_mode)
        cfg.phase_field.plastic_work_to_Gc_efficiency = float(args.wp_gc_efficiency)
        cfg.phase_field.Gc_local_cap_factor = float(args.gc_local_cap_factor)
        cfg.phase_field.toughening_storage_coeff = float(args.toughening_storage_coeff)
        cfg.phase_field.toughening_dissipation_coeff = float(args.toughening_dissipation_coeff)

    cfg.tip_memory.enabled = not args.memory_off
    cfg.tip_memory.mode = "off" if args.memory_off else "stage1"
    cfg.tip_memory.state_gain = float(case["memory_gain"])
    cfg.tip_memory.M_max = float(case["M_max"])
    cfg.tip_memory.amp_max = float(case["amp_max"])
    cfg.tip_memory.shield_max = float(case["shield_max"])
    cfg.tip_memory.blunt_per_work = float(case["blunt_per_work"])
    cfg.tip_memory.sharpen_per_damage = float(case["sharpen_per_damage"])
    cfg.tip_memory.drive_exponent = float(args.tip_drive_exponent)
    cfg.tip_memory.couple_to_damage_drive = not args.no_tip_drive_coupling

    _load_exp_floor_plastic_barrier(
        cfg,
        args.plastic_barrier_json,
        case["system"],
        energy_scale=float(case["barrier_scale"]),
        entropy_scale=float(case["entropy_scale"]),
        stress_scale=float(case["stress_scale"]),
    )
    cfg.plasticity_barrier.exp_v_max_b3 = float(args.plastic_exp_vmax_b3)

    tag = _case_tag(case)
    cfg.output_dir = os.path.join(args.output_dir, f"case_{case_index:04d}__{tag}")

    print("\n" + "=" * 90)
    print(f"PARAM SWEEP CASE {case_index}/{n_cases}: {tag}")
    print("=" * 90)
    results = run_simulation(cfg)

    temp_rows: List[Dict[str, Any]] = []
    for T in sorted(results.keys()):
        row = history_summary(results[T])
        row.update(case)
        row["case_index"] = case_index
        row["case_tag"] = tag
        row["T"] = float(T)
        temp_rows.append(row)

    score = _score_case(temp_rows)
    case_row = dict(case)
    case_row.update(score)
    case_row["case_index"] = case_index
    case_row["case_tag"] = tag
    case_row["output_dir"] = cfg.output_dir
    return {"case_row": case_row, "temp_rows": temp_rows}


def main():
    p = argparse.ArgumentParser(description="Overnight EXP_floor barrier/memory diagnostic sweep")
    p.add_argument("--name", default="param_sweep", help="Sweep name, used in metadata only")
    p.add_argument("--output-dir", default="results_param_sweep", help="Output root directory")
    p.add_argument("--plastic-barrier-json", default="BarrierModel_Export.json")
    p.add_argument("--barrier-systems", nargs="+", default=["W[100]"],
                   help="EXP_floor systems; supports all, all_non_si")
    p.add_argument("--temperatures", nargs="+", type=float, default=[425, 450, 475, 500, 525, 550])
    p.add_argument("--Gc-list", nargs="+", type=float, default=[63.0])
    p.add_argument("--barrier-scales", nargs="+", type=float, default=[0.1, 0.2, 0.3])
    p.add_argument("--entropy-scales", nargs="+", type=float, default=[0.5, 1.0, 1.5])
    p.add_argument("--stress-scales", nargs="+", type=float, default=[0.8, 1.0, 1.25])
    p.add_argument("--flow-epsdot-refs", nargs="+", type=float, default=[1e-5, 1e-4, 1e-3])
    p.add_argument("--memory-gains", nargs="+", type=float, default=[0.25, 1.0])
    p.add_argument("--M-max-values", nargs="+", type=float, default=[3.0, 4.0, 5.0])
    p.add_argument("--amp-max-values", nargs="+", type=float, default=[5.0])
    p.add_argument("--shield-max-values", nargs="+", type=float, default=[0.85])
    p.add_argument("--blunt-work-values", nargs="+", type=float, default=[0.35])
    p.add_argument("--sharpen-damage-values", nargs="+", type=float, default=[0.35])
    p.add_argument("--tip-drive-exponent", type=float, default=2.0)
    p.add_argument("--no-tip-drive-coupling", action="store_true")
    p.add_argument("--memory-off", action="store_true")

    p.add_argument("--steps", type=int, default=140)
    p.add_argument("--dU-top", type=float, default=1e-6)
    p.add_argument("--nx", type=int, default=None, help="Mesh divisions in x; nodes = nx+1")
    p.add_argument("--ny", type=int, default=None, help="Mesh divisions in y; nodes = ny+1")
    p.add_argument("--mesh-jitter", type=float, default=None, help="Interior node jitter fraction; 0 disables jitter")
    p.add_argument("--ell-factor", type=float, default=None, help="Use ell = ell_factor*hbar unless --ell is supplied")
    p.add_argument("--ell", type=float, default=None, help="Fixed physical phase-field length ell [m]; overrides --ell-factor")
    p.add_argument("--save-every", type=int, default=999999,
                   help="Snapshot interval. Default effectively disables intermediate snapshots.")
    p.add_argument("--save-fields", action="store_true",
                   help="Save raw field snapshots for every case. Off by default for large sweeps.")
    p.add_argument("--save-field-pngs", action="store_true",
                   help="Save field PNGs for every case. Off by default for large sweeps.")
    p.add_argument("--make-plots", action="store_true",
                   help="Make per-case plots. Off by default for large sweeps.")
    p.add_argument("--no-auto-stop", action="store_true")
    p.add_argument("--stop-on-invalid", action="store_true", default=True)
    p.add_argument("--allow-invalid-soft-tearing", action="store_true")
    p.add_argument("--invalid-min-step", type=int, default=5)
    p.add_argument("--invalid-wp-wext-pct", type=float, default=2000.0)
    p.add_argument("--freeze-rho", action="store_true", default=True)
    p.add_argument("--disable-plasticity", action="store_true")
    p.add_argument("--disable-wp-gc-coupling", action="store_true", default=True)
    p.add_argument("--wp-gc-coupling-mode", choices=['off','direct','state'], default='state')
    p.add_argument("--wp-gc-efficiency", type=float, default=0.05)
    p.add_argument("--gc-local-cap-factor", type=float, default=100.0)
    p.add_argument("--toughening-storage-coeff", type=float, default=0.05)
    p.add_argument("--toughening-dissipation-coeff", type=float, default=0.02)
    p.add_argument("--rho-cap", type=float, default=1e16)
    p.add_argument("--dot-ep-max", type=float, default=1e3)
    p.add_argument("--max-plastic-strain-increment", type=float, default=2.5e-4)
    p.add_argument("--max-rho-relative-increment", type=float, default=0.25)
    p.add_argument("--plastic-exp-vmax-b3", type=float, default=1e4)

    p.add_argument("--max-cases", type=int, default=None, help="Random subset size after grid construction")
    p.add_argument("--shuffle", action="store_true", help="Shuffle case order")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--dry-run", action="store_true", help="Only write the planned cases, do not run")
    p.add_argument("--resume", action="store_true",
                   help="Resume an interrupted sweep: load existing planned_cases.csv when present, skip complete cases, and rerun incomplete/missing cases.")
    p.add_argument("--rebuild-summary-only", action="store_true",
                   help="Do not run cases; rebuild aggregate CSV summaries from completed case directories.")
    p.add_argument("--start-case", type=int, default=1,
                   help="First 1-based case index to consider. Useful for manual restarts.")
    p.add_argument("--end-case", type=int, default=None,
                   help="Last 1-based case index to consider. Useful for chunked/manual restarts.")
    args = p.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    planned_path = _planned_cases_path(args.output_dir)
    if args.resume and os.path.exists(planned_path):
        cases = _load_planned_cases(args.output_dir)
        print(f"Resume mode: loaded {len(cases)} planned cases from {planned_path}")
    else:
        cases = build_cases(args)
        _write_csv(os.path.join(args.output_dir, "planned_cases.csv"), [dict(c, case_index=i+1, case_tag=_case_tag(c)) for i, c in enumerate(cases)])
        print(f"Planned {len(cases)} cases. Wrote {os.path.join(args.output_dir, 'planned_cases.csv')}")

    with open(os.path.join(args.output_dir, "sweep_metadata.json"), "w") as f:
        json.dump({"name": args.name, "n_cases": len(cases), "args": vars(args)}, f, indent=2)
    if args.dry_run:
        return

    start_case = max(1, int(args.start_case))
    end_case = int(args.end_case) if args.end_case is not None else len(cases)
    end_case = min(end_case, len(cases))

    all_case_rows: List[Dict[str, Any]] = []
    all_temp_rows: List[Dict[str, Any]] = []
    t0 = time.time()
    n_ran = 0
    n_skipped = 0
    n_incomplete_before_start = 0

    for i, case in enumerate(cases, start=1):
        # Always include already-completed cases in the aggregate summary if possible.
        complete = _case_is_complete(args, case, i)
        if complete:
            all_case_rows.append(_case_row_from_existing(args, case, i))
            all_temp_rows.extend(_existing_case_temperature_rows(args, case, i))

        if i < start_case or i > end_case:
            if not complete:
                n_incomplete_before_start += 1
            continue

        if args.rebuild_summary_only:
            continue

        if args.resume and complete:
            n_skipped += 1
            print(f"Resume mode: skipping complete case {i}/{len(cases)}: {_case_tag(case)}")
            ranked = sorted(all_case_rows, key=lambda r: float(r.get("score", -1e9)), reverse=True)
            _write_csv(os.path.join(args.output_dir, "param_sweep_case_summary.csv"), ranked)
            _write_csv(os.path.join(args.output_dir, "param_sweep_temperature_summary.csv"), all_temp_rows)
            continue

        try:
            out = run_case(args, case, i, len(cases))
            # If this was an incomplete rerun, remove any previous partial/error row for this index.
            all_case_rows = [r for r in all_case_rows if int(_to_float(r.get("case_index"), -1)) != i]
            all_temp_rows = [r for r in all_temp_rows if int(_to_float(r.get("case_index"), -1)) != i]
            all_case_rows.append(out["case_row"])
            all_temp_rows.extend(out["temp_rows"])
            n_ran += 1
        except Exception as exc:
            err = dict(case)
            err.update({"case_index": i, "case_tag": _case_tag(case), "error": repr(exc), "score": -1e9, "resume_status": "error"})
            all_case_rows = [r for r in all_case_rows if int(_to_float(r.get("case_index"), -1)) != i]
            all_case_rows.append(err)
            print(f"ERROR in case {i}: {exc!r}")
        # write incremental summaries so an interrupted overnight run is still useful
        ranked = sorted(all_case_rows, key=lambda r: float(r.get("score", -1e9)), reverse=True)
        _write_csv(os.path.join(args.output_dir, "param_sweep_case_summary.csv"), ranked)
        _write_csv(os.path.join(args.output_dir, "param_sweep_temperature_summary.csv"), all_temp_rows)

    ranked = sorted(all_case_rows, key=lambda r: float(r.get("score", -1e9)), reverse=True)
    _write_csv(os.path.join(args.output_dir, "param_sweep_case_summary.csv"), ranked)
    _write_csv(os.path.join(args.output_dir, "param_sweep_temperature_summary.csv"), all_temp_rows)
    print(f"Resume/rebuild complete: ran {n_ran}, skipped complete {n_skipped}, aggregate rows {len(all_case_rows)} cases / {len(all_temp_rows)} temperatures in {(time.time()-t0)/3600:.2f} h")
    if n_incomplete_before_start:
        print(f"Note: {n_incomplete_before_start} incomplete cases were outside --start-case/--end-case and were not run.")


if __name__ == "__main__":
    main()
