#!/usr/bin/env python3
"""Search one common developed moving-process-zone closure for three classes.

Inputs are the three intrinsic EXP-floor rows and the single strict common
Peierls--Taylor closure selected from the v9.4 constitutive search.  The PT
closure and intrinsic cleavage/emission surfaces remain fixed.  Only active
moving-process-zone state parameters are sampled, and the same sampled closure
is applied to ceramic, weak-T, and DBTT-precursor rows.

This stage uses the reduced one-dimensional moving-front engine.  It evaluates
R-curve and state convergence over a finite crack extension before any 2-D
FEM/CZM calculation.  The obsolete fixed MPZ glide and detrap barriers are not
sampled because production transport uses the emission-derived PT model.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import math
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pandas as pd
from scipy.stats import qmc

from fit_mpz_three_classes import (
    class_constraint_loss,
    point_target_loss,
    simulate,
    state_regularization_loss,
)


REGION_TO_CLASS = {
    "ceramic_intrinsic": "ceramic",
    "weakT_intrinsic": "weakT",
    "DBTT_precursor": "DBTT",
}

SHARED_COLUMNS = [
    "mpz_source_sites_per_system",
    "mpz_trap_barrier_eV",
    "mpz_retained_recovery_barrier_eV",
    "c_blunt",
    "mpz_source_recovery_rate_s",
    "mpz_source_refresh_length_m",
    "mpz_length_m",
]


def _log_scale(u: np.ndarray, lo: float, hi: float) -> np.ndarray:
    return 10.0 ** (
        math.log10(lo) + u * (math.log10(hi) - math.log10(lo))
    )


def sample_shared_parameters(n: int, seed: int, base: pd.Series) -> pd.DataFrame:
    """Sobol sample of active common MPZ state parameters."""
    sampler = qmc.Sobol(7, scramble=True, seed=seed)
    if n > 0 and (n & (n - 1)) == 0:
        u = sampler.random_base2(int(math.log2(n)))
    else:
        u = sampler.random(n)
    sampled = pd.DataFrame({
        "mpz_source_sites_per_system": _log_scale(u[:, 0], 2.0, 300.0),
        "mpz_trap_barrier_eV": 0.25 + 1.55 * u[:, 1],
        "mpz_retained_recovery_barrier_eV": 0.80 + 2.20 * u[:, 2],
        "c_blunt": _log_scale(u[:, 3], 0.05, 3.0),
        "mpz_source_recovery_rate_s": _log_scale(u[:, 4], 1.0e-10, 1.0e-2),
        "mpz_source_refresh_length_m": _log_scale(u[:, 5], 5.0e-7, 2.0e-5),
        "mpz_length_m": _log_scale(u[:, 6], 1.0e-5, 2.0e-4),
    })

    # Preserve the inherited row and a small set of interpretable nearby
    # anchors.  They are common closures, not class-specific fits.
    anchor0 = {c: float(base[c]) for c in SHARED_COLUMNS}
    anchors = [anchor0]
    for site_scale, blunt_scale, recovery_scale in (
        (0.5, 0.5, 1.0),
        (1.0, 1.0, 0.25),
        (2.0, 1.5, 1.0),
        (1.0, 2.0, 4.0),
    ):
        rec = dict(anchor0)
        rec["mpz_source_sites_per_system"] = np.clip(
            rec["mpz_source_sites_per_system"] * site_scale, 2.0, 300.0
        )
        rec["c_blunt"] = np.clip(rec["c_blunt"] * blunt_scale, 0.05, 3.0)
        rec["mpz_source_recovery_rate_s"] = np.clip(
            max(rec["mpz_source_recovery_rate_s"], 1.0e-10) * recovery_scale,
            1.0e-10,
            1.0e-2,
        )
        anchors.append(rec)
    return pd.concat([pd.DataFrame(anchors), sampled], ignore_index=True)


def parse_temperatures(text: str) -> dict[str, list[float]]:
    """Parse ``class:T,T;class:T,T`` temperature specification."""
    result: dict[str, list[float]] = {}
    for block in str(text).split(";"):
        block = block.strip()
        if not block:
            continue
        name, values = block.split(":", 1)
        result[name.strip()] = [
            float(x) for x in values.replace(",", " ").split() if x
        ]
    missing = [x for x in ("ceramic", "weakT", "DBTT") if x not in result]
    if missing:
        raise ValueError(f"temperature specification lacks classes: {missing}")
    return result


def apply_shared(base: pd.Series, shared: dict[str, float]) -> pd.Series:
    row = base.copy()
    for key, value in shared.items():
        row[key] = float(value)
    # These legacy knobs remain in serialized rows for compatibility but are
    # inactive under emission-derived PT transport and are never searched here.
    row["mpz_pair_annihilation_rate_per_count_s"] = float(
        row.get("mpz_pair_annihilation_rate_per_count_s", 0.0)
    )
    return row


def convergence_length(
    events: list[dict[str, Any]],
    column: str,
    relative_tolerance: float = 0.12,
    absolute_tolerance: float = 1.0e-8,
) -> float:
    """First extension after which a smoothed state remains near its final level."""
    if len(events) < 8:
        return float("nan")
    df = pd.DataFrame(events)
    if column not in df or "a_um" not in df:
        return float("nan")
    y = pd.Series(df[column].to_numpy(float)).rolling(
        5, center=True, min_periods=1
    ).median().to_numpy()
    a = df["a_um"].to_numpy(float)
    n_tail = max(5, len(y) // 5)
    final = float(np.median(y[-n_tail:]))
    tol = max(abs(final) * relative_tolerance, absolute_tolerance)
    for i in range(len(y)):
        if np.all(np.abs(y[i:] - final) <= tol):
            return float(a[i])
    return float(a[-1])


def _evaluate_candidate(payload: dict[str, Any]) -> dict[str, Any]:
    candidate_index = int(payload["candidate_index"])
    shared = {k: float(v) for k, v in payload["shared"].items()}
    base_rows = {
        k: pd.Series(v) for k, v in payload["base_rows"].items()
    }
    targets = pd.DataFrame(payload["targets"])
    temperature_map = payload["temperature_map"]
    opt = SimpleNamespace(**payload["opt"])

    total_point = 0.0
    total_class = 0.0
    total_state = 0.0
    metrics_rows: list[dict[str, Any]] = []
    failed = False
    error = ""

    try:
        for klass in ("ceramic", "weakT", "DBTT"):
            row = apply_shared(base_rows[klass], shared)
            temperatures = [float(x) for x in temperature_map[klass]]
            sims = [simulate(row, T, opt) for T in temperatures]
            class_targets = targets[
                (targets["target_class"] == klass)
                & targets["T_K"].astype(float).isin(temperatures)
            ]
            point, _ = point_target_loss(class_targets, sims, "rcurve")
            class_loss, class_parts = class_constraint_loss(
                klass, sims, "rcurve"
            )
            state_loss, state_parts = state_regularization_loss(sims, row)
            total_point += float(point)
            total_class += float(class_loss)
            total_state += float(state_loss)

            for sim in sims:
                events = sim.get("events", [])
                record = {
                    "candidate_index": candidate_index,
                    "target_class": klass,
                    **shared,
                    **{k: v for k, v in sim.items() if k != "events"},
                    "K_convergence_length_um": convergence_length(
                        events, "K_MPa_sqrt_m", 0.08, 0.15
                    ),
                    "shield_convergence_length_um": convergence_length(
                        events, "shield_fraction", 0.12, 0.01
                    ),
                    "blunt_convergence_length_um": convergence_length(
                        events, "blunt_ratio", 0.12, 0.02
                    ),
                    "retained_convergence_length_um": convergence_length(
                        events, "retained_count", 0.15, 0.25
                    ),
                }
                for name, value in class_parts.items():
                    record[f"class_loss_{name}"] = float(value)
                for name, value in state_parts.items():
                    record[f"state_loss_{name}"] = float(value)
                metrics_rows.append(record)
    except Exception as exc:  # retain failed points for diagnosis and resume
        failed = True
        error = f"{type(exc).__name__}: {exc}"

    total = total_point + 2.0 * total_class + total_state
    if failed or not np.isfinite(total):
        total = 1.0e12
    return {
        "candidate_index": candidate_index,
        **shared,
        "score": float(total),
        "point_loss": float(total_point),
        "class_loss": float(total_class),
        "state_loss": float(total_state),
        "failed": bool(failed),
        "error": error,
        "metrics": metrics_rows,
    }


def _write_checkpoint(path: Path, result: dict[str, Any]) -> None:
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(result, allow_nan=True))
    tmp.replace(path)


def _load_checkpoint(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def save_best_events(
    shortlist: pd.DataFrame,
    base_rows: dict[str, pd.Series],
    targets: pd.DataFrame,
    temperature_map: dict[str, list[float]],
    opt: SimpleNamespace,
    out: Path,
    top_count: int,
) -> None:
    all_events: list[dict[str, Any]] = []
    all_metrics: list[dict[str, Any]] = []
    for _, cand in shortlist.head(top_count).iterrows():
        shared = {c: float(cand[c]) for c in SHARED_COLUMNS}
        for klass in ("ceramic", "weakT", "DBTT"):
            row = apply_shared(base_rows[klass], shared)
            for T in temperature_map[klass]:
                sim = simulate(row, float(T), opt)
                all_metrics.append({
                    "candidate_index": int(cand.candidate_index),
                    "target_class": klass,
                    **shared,
                    **{k: v for k, v in sim.items() if k != "events"},
                    "K_convergence_length_um": convergence_length(
                        sim.get("events", []), "K_MPa_sqrt_m", 0.08, 0.15
                    ),
                    "shield_convergence_length_um": convergence_length(
                        sim.get("events", []), "shield_fraction", 0.12, 0.01
                    ),
                    "blunt_convergence_length_um": convergence_length(
                        sim.get("events", []), "blunt_ratio", 0.12, 0.02
                    ),
                    "retained_convergence_length_um": convergence_length(
                        sim.get("events", []), "retained_count", 0.15, 0.25
                    ),
                })
                for ev in sim.get("events", []):
                    all_events.append({
                        "candidate_index": int(cand.candidate_index),
                        "target_class": klass,
                        "T_K": float(T),
                        **shared,
                        **ev,
                    })
    pd.DataFrame(all_metrics).to_csv(
        out / "developed_state_shortlist_metrics.csv", index=False
    )
    pd.DataFrame(all_events).to_csv(
        out / "developed_state_shortlist_events.csv", index=False
    )


def save_plots(out: Path) -> None:
    events_path = out / "developed_state_shortlist_events.csv"
    if not events_path.exists():
        return
    events = pd.read_csv(events_path)
    if events.empty:
        return
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return

    for candidate_index, group in events.groupby("candidate_index"):
        fig, axes = plt.subplots(2, 2, figsize=(11.5, 8.5))
        for (klass, T), g in group.groupby(["target_class", "T_K"]):
            g = g.sort_values("a_um")
            label = f"{klass} {T:g} K"
            axes[0, 0].plot(g.a_um, g.K_MPa_sqrt_m, label=label)
            axes[0, 1].plot(g.a_um, g.shield_fraction, label=label)
            axes[1, 0].plot(g.a_um, g.blunt_ratio, label=label)
            axes[1, 1].plot(g.a_um, g.retained_count, label=label)
        axes[0, 0].set_ylabel("K [MPa sqrt(m)]")
        axes[0, 1].set_ylabel("Shield fraction")
        axes[1, 0].set_ylabel("Blunting ratio")
        axes[1, 1].set_ylabel("Retained count")
        for ax in axes.flat:
            ax.set_xlabel("Crack extension [um]")
            ax.grid(True, alpha=0.25)
        axes[0, 0].legend(fontsize=7, ncol=2)
        fig.suptitle(f"Common MPZ candidate {int(candidate_index)}")
        fig.tight_layout()
        fig.savefig(
            out / f"developed_state_candidate_{int(candidate_index):05d}.png",
            dpi=180,
        )
        plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--selected-rows", type=Path,
        default=Path(
            "runs/mpz_v9_4_peierls_taylor_search_v1/strict_common_selection/"
            "pt_v9_4_recommended_intrinsic_rows.csv"
        ),
    )
    ap.add_argument(
        "--targets", type=Path,
        default=Path("mpz_three_class_design_targets.csv"),
    )
    ap.add_argument("--samples", type=int, default=128)
    ap.add_argument("--seed", type=int, default=94131)
    ap.add_argument("--max-workers", type=int, default=1)
    ap.add_argument(
        "--temperatures",
        default="ceramic:300,900,1200;weakT:300,700,1200;DBTT:300,700,900,1200",
    )
    ap.add_argument("--target-extension-um", type=float, default=400.0)
    ap.add_argument("--da-um", type=float, default=5.0)
    ap.add_argument("--dK", type=float, default=0.25)
    ap.add_argument("--Kdot", type=float, default=0.005)
    ap.add_argument("--Kmax", type=float, default=65.0)
    ap.add_argument("--top-count", type=int, default=12)
    ap.add_argument("--event-top-count", type=int, default=3)
    ap.add_argument("--resume", action="store_true", default=True)
    ap.add_argument(
        "--out", type=Path,
        default=Path("runs/mpz_v9_4_developed_state_search_v1"),
    )
    a = ap.parse_args()

    if not a.selected_rows.exists():
        raise SystemExit(f"selected intrinsic rows not found: {a.selected_rows}")
    rows = pd.read_csv(a.selected_rows)
    rows["target_class"] = rows["region"].map(REGION_TO_CLASS)
    if rows["target_class"].isna().any():
        raise SystemExit(
            "unrecognized selected region(s): "
            + ", ".join(rows.loc[rows.target_class.isna(), "region"].astype(str))
        )
    base_rows = {
        klass: group.sort_values("pt_score").iloc[0].copy()
        for klass, group in rows.groupby("target_class")
    }
    missing = [x for x in ("ceramic", "weakT", "DBTT") if x not in base_rows]
    if missing:
        raise SystemExit(f"selected rows lack classes: {missing}")

    targets = pd.read_csv(a.targets)
    temperature_map = parse_temperatures(a.temperatures)
    n_advances = int(round(a.target_extension_um / a.da_um)) + 1
    opt = SimpleNamespace(
        dK=float(a.dK),
        Kdot=float(a.Kdot),
        n_advances=n_advances,
        Kmax=float(a.Kmax),
        da_um=float(a.da_um),
        early_window_um=(20.0, min(180.0, 0.45 * a.target_extension_um)),
        plateau_window_um=(0.70 * a.target_extension_um, a.target_extension_um),
        target_dB_substep=0.25,
        target_emission_hazard_substep=1.0,
        source_active_fraction_min=1.0e-4,
        min_substep_fraction=1.0e-8,
        max_substeps=2_000_000,
        objective_mode="rcurve",
    )

    out = a.out.resolve()
    checkpoints = out / "checkpoints"
    checkpoints.mkdir(parents=True, exist_ok=True)
    samples = sample_shared_parameters(a.samples, a.seed, base_rows["weakT"])

    base_payload = {
        "base_rows": {k: v.to_dict() for k, v in base_rows.items()},
        "targets": targets.to_dict(orient="records"),
        "temperature_map": temperature_map,
        "opt": vars(opt),
    }
    results: list[dict[str, Any]] = []
    pending: list[dict[str, Any]] = []
    for idx, shared_row in samples.iterrows():
        checkpoint = checkpoints / f"candidate_{int(idx):05d}.json"
        if a.resume and checkpoint.exists():
            results.append(_load_checkpoint(checkpoint))
        else:
            pending.append({
                **base_payload,
                "candidate_index": int(idx),
                "shared": {
                    c: float(shared_row[c]) for c in SHARED_COLUMNS
                },
            })

    completed = len(results)
    total = len(samples)
    if pending:
        if a.max_workers <= 1:
            for payload in pending:
                result = _evaluate_candidate(payload)
                _write_checkpoint(
                    checkpoints / f"candidate_{result['candidate_index']:05d}.json",
                    result,
                )
                results.append(result)
                completed += 1
                print(
                    f"evaluated {completed}/{total} candidate={result['candidate_index']} "
                    f"score={result['score']:.6g} failed={result['failed']}",
                    flush=True,
                )
        else:
            with ProcessPoolExecutor(max_workers=a.max_workers) as pool:
                futures = {pool.submit(_evaluate_candidate, p): p for p in pending}
                for future in as_completed(futures):
                    result = future.result()
                    _write_checkpoint(
                        checkpoints / f"candidate_{result['candidate_index']:05d}.json",
                        result,
                    )
                    results.append(result)
                    completed += 1
                    print(
                        f"evaluated {completed}/{total} candidate={result['candidate_index']} "
                        f"score={result['score']:.6g} failed={result['failed']}",
                        flush=True,
                    )

    summary_rows = []
    metric_rows = []
    for result in results:
        summary_rows.append({k: v for k, v in result.items() if k != "metrics"})
        metric_rows.extend(result.get("metrics", []))
    all_df = pd.DataFrame(summary_rows).sort_values("score").reset_index(drop=True)
    metrics_df = pd.DataFrame(metric_rows)
    shortlist = all_df[~all_df.failed.astype(bool)].head(a.top_count).copy()
    shortlist.insert(0, "developed_rank", np.arange(1, len(shortlist) + 1))

    all_df.to_csv(out / "developed_state_search_all.csv", index=False)
    metrics_df.to_csv(out / "developed_state_search_metrics.csv", index=False)
    shortlist.to_csv(out / "developed_state_search_shortlist.csv", index=False)
    save_best_events(
        shortlist,
        base_rows,
        targets,
        temperature_map,
        opt,
        out,
        min(a.event_top_count, len(shortlist)),
    )
    save_plots(out)

    config = vars(a).copy()
    for key, value in list(config.items()):
        if isinstance(value, Path):
            config[key] = str(value)
    config.update({
        "temperature_map": temperature_map,
        "n_advances": n_advances,
        "shared_columns": SHARED_COLUMNS,
        "inactive_not_searched": [
            "mpz_glide_barrier_eV",
            "mpz_detrap_barrier_eV",
        ],
        "pt_closure_fixed": True,
        "intrinsic_surfaces_fixed": True,
        "status": "REDUCED_DEVELOPED_MPZ_SEARCH_NOT_2D_VALIDATED",
    })
    (out / "developed_state_search_config.json").write_text(
        json.dumps(config, indent=2)
    )
    print(
        json.dumps({
            "n_candidates": int(len(all_df)),
            "n_failed": int(all_df.failed.astype(bool).sum()),
            "best_score": float(shortlist.iloc[0].score) if len(shortlist) else None,
            "output": str(out),
        }, indent=2),
        flush=True,
    )


if __name__ == "__main__":
    main()
