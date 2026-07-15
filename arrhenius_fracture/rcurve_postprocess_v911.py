"""Cascade-aware resistance-curve postprocessing for v9.11.

A cohesive backend may serialize one unstable fixed-displacement jump into many
one-edge topology updates. Those rows are not independent resistance points.
This module preserves the raw events but clusters consecutive events that occur
with negligible change in remote displacement and reports each cluster as one
load event with an associated unstable jump span.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

A0_M = 0.5e-3


def find_steps_file(case_dir: str | Path, T_K: float) -> Path | None:
    root = Path(case_dir)
    exact = root / f"steps_{int(round(float(T_K))):04d}K.csv"
    if exact.exists():
        return exact
    matches = sorted(root.glob("steps_*K.csv"))
    return matches[0] if matches else None


def extract_raw_growth_events(case_dir: str | Path, T_K: float) -> pd.DataFrame:
    path = find_steps_file(case_dir, T_K)
    if path is None:
        return pd.DataFrame()
    st = pd.read_csv(path)
    required = {"KJ_Pa_sqrtm", "a_tip_m", "Uapp_m"}
    if not required.issubset(st.columns):
        return pd.DataFrame()

    if "crack_extension_m" in st.columns:
        ext_m = pd.to_numeric(st["crack_extension_m"], errors="coerce").to_numpy(float)
    else:
        ext_m = pd.to_numeric(st["a_tip_m"], errors="coerce").to_numpy(float) - A0_M
    ext_m = np.maximum(ext_m, 0.0)
    if "da_block_m" in st.columns:
        da_m = pd.to_numeric(st["da_block_m"], errors="coerce").fillna(0.0).to_numpy(float)
    else:
        da_m = np.r_[0.0, np.maximum(np.diff(ext_m), 0.0)]
    n_fire = (
        pd.to_numeric(st["n_fire"], errors="coerce").fillna(0.0).to_numpy(float)
        if "n_fire" in st.columns
        else np.zeros(len(st))
    )
    idx = np.flatnonzero((da_m > 1.0e-12) | (n_fire > 0.0))
    if idx.size == 0:
        return pd.DataFrame()

    out = pd.DataFrame({
        "raw_event_id": np.arange(1, idx.size + 1, dtype=int),
        "step": pd.to_numeric(st.iloc[idx].get("step", pd.Series(idx)), errors="coerce").to_numpy(),
        "Uapp_m": pd.to_numeric(st.iloc[idx]["Uapp_m"], errors="coerce").to_numpy(float),
        "KJ_MPa_sqrt_m": pd.to_numeric(st.iloc[idx]["KJ_Pa_sqrtm"], errors="coerce").to_numpy(float) / 1.0e6,
        "crack_extension_after_um": ext_m[idx] * 1.0e6,
        "da_block_um": da_m[idx] * 1.0e6,
        "n_fire": n_fire[idx],
    })
    out["crack_extension_before_um"] = np.maximum(
        out["crack_extension_after_um"] - out["da_block_um"], 0.0
    )

    optional = {
        "B_residual": "B",
        "N_em": "N_em",
        "sigma_tip_GPa": "sigma_tip_Pa",
        "K_shield_MPa_sqrt_m": "mpz_K_shield_Pa_sqrt_m",
        "adaptive_frac": "adaptive_frac",
        "dt_cur_s": "dt_cur_s",
        "B_target": "B_target",
        "B_fraction_of_target": "B_fraction_of_target",
        "stochastic_event_index": "stochastic_event_index",
    }
    for name, source in optional.items():
        if source in st.columns:
            vals = pd.to_numeric(st.iloc[idx][source], errors="coerce").to_numpy(float)
            if source == "sigma_tip_Pa":
                vals *= 1.0e-9
            if source == "mpz_K_shield_Pa_sqrt_m":
                vals *= 1.0e-6
            out[name] = vals
    return out.replace([np.inf, -np.inf], np.nan).dropna(
        subset=["Uapp_m", "KJ_MPa_sqrt_m", "crack_extension_after_um"]
    ).reset_index(drop=True)


def cluster_same_load_events(
    raw: pd.DataFrame,
    relative_load_tolerance: float = 1.0e-4,
    absolute_load_tolerance_m: float = 1.0e-12,
) -> pd.DataFrame:
    if raw.empty:
        return pd.DataFrame()
    rel_tol = max(float(relative_load_tolerance), 0.0)
    abs_tol = max(float(absolute_load_tolerance_m), 0.0)
    groups: list[list[int]] = [[0]]
    U = raw["Uapp_m"].to_numpy(float)
    for i in range(1, len(raw)):
        dU = abs(float(U[i] - U[i - 1]))
        scale = max(abs(float(U[i])), abs(float(U[i - 1])), 1.0e-30)
        same = dU <= max(abs_tol, rel_tol * scale)
        if same:
            groups[-1].append(i)
        else:
            groups.append([i])

    rows: list[dict[str, Any]] = []
    for cid, ids in enumerate(groups, start=1):
        block = raw.iloc[ids]
        first = block.iloc[0]
        last = block.iloc[-1]
        start = float(first["crack_extension_before_um"])
        end = float(last["crack_extension_after_um"])
        u0 = float(first["Uapp_m"])
        u1 = float(last["Uapp_m"])
        event_count = int(len(block))
        row = {
            "load_event_id": cid,
            "classification": (
                "unstable_same_load_cascade" if event_count > 1 else "single_topology_event"
            ),
            "raw_event_start": int(first["raw_event_id"]),
            "raw_event_end": int(last["raw_event_id"]),
            "topology_event_count": event_count,
            "step_start": float(first["step"]),
            "step_end": float(last["step"]),
            "crack_extension_um": start,
            "crack_extension_start_um": start,
            "crack_extension_end_um": end,
            "jump_span_um": max(end - start, 0.0),
            "Uapp_onset_m": u0,
            "Uapp_end_m": u1,
            "relative_load_change": abs(u1 - u0) / max(abs(u0), 1.0e-30),
            "KJ_MPa_sqrt_m": float(first["KJ_MPa_sqrt_m"]),
            "KJ_onset_MPa_sqrt_m": float(first["KJ_MPa_sqrt_m"]),
            "KJ_end_MPa_sqrt_m": float(last["KJ_MPa_sqrt_m"]),
            "KJ_min_MPa_sqrt_m": float(block["KJ_MPa_sqrt_m"].min()),
            "KJ_max_MPa_sqrt_m": float(block["KJ_MPa_sqrt_m"].max()),
        }
        for key in (
            "B_target", "B_fraction_of_target", "stochastic_event_index",
            "N_em", "K_shield_MPa_sqrt_m",
        ):
            if key in block.columns:
                row[f"{key}_onset"] = float(first[key])
                row[f"{key}_end"] = float(last[key])
        rows.append(row)
    return pd.DataFrame(rows)


def cascade_metrics(raw: pd.DataFrame, clustered: pd.DataFrame) -> dict[str, Any]:
    if raw.empty or clustered.empty:
        return {
            "n_raw_topology_events": 0,
            "n_independent_load_events": 0,
            "n_unstable_same_load_cascades": 0,
            "largest_same_load_jump_um": np.nan,
            "fraction_topology_events_in_cascades": np.nan,
            "rcurve_interpretation": "no_growth_events",
        }
    unstable = clustered[clustered["topology_event_count"] > 1]
    events_in = int(unstable["topology_event_count"].sum()) if not unstable.empty else 0
    return {
        "n_raw_topology_events": int(len(raw)),
        "n_independent_load_events": int(len(clustered)),
        "n_unstable_same_load_cascades": int(len(unstable)),
        "largest_same_load_jump_um": (
            float(unstable["jump_span_um"].max()) if not unstable.empty else 0.0
        ),
        "fraction_topology_events_in_cascades": float(events_in / max(len(raw), 1)),
        "rcurve_interpretation": (
            "contains_unstable_fixed_displacement_cascades"
            if not unstable.empty
            else "independent_reload_events"
        ),
    }


def write_cascade_aware_outputs(
    case_dir: str | Path,
    T_K: float,
    relative_load_tolerance: float = 1.0e-4,
    absolute_load_tolerance_m: float = 1.0e-12,
) -> dict[str, Any]:
    root = Path(case_dir)
    raw = extract_raw_growth_events(root, T_K)
    clustered = cluster_same_load_events(
        raw,
        relative_load_tolerance=relative_load_tolerance,
        absolute_load_tolerance_m=absolute_load_tolerance_m,
    )
    raw.to_csv(root / "R_curve_topology_events_raw.csv", index=False)
    clustered.to_csv(root / "R_curve_load_events_clustered.csv", index=False)
    # Compatibility filename now contains independent load events, not every
    # serialized cohesive-edge insertion.
    clustered.to_csv(root / "R_curve_event_sampled.csv", index=False)
    metrics = cascade_metrics(raw, clustered)
    pd.DataFrame([metrics]).to_csv(root / "R_curve_cascade_metrics.csv", index=False)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        if not clustered.empty:
            fig, ax = plt.subplots(figsize=(7.4, 4.9))
            singles = clustered[clustered["topology_event_count"] == 1]
            casc = clustered[clustered["topology_event_count"] > 1]
            if not singles.empty:
                ax.plot(
                    singles["crack_extension_um"], singles["KJ_MPa_sqrt_m"],
                    marker="o", linewidth=1.0, markersize=3, label="independent event",
                )
            for j, (_, row) in enumerate(casc.iterrows()):
                ax.hlines(
                    row["KJ_onset_MPa_sqrt_m"],
                    row["crack_extension_start_um"],
                    row["crack_extension_end_um"],
                    linewidth=2.0,
                    label="unstable same-load jump" if j == 0 else None,
                )
                ax.plot(row["crack_extension_start_um"], row["KJ_onset_MPa_sqrt_m"], "o")
            ax.set_xlabel(r"Projected crack extension $\Delta a_x$ ($\mu$m)")
            ax.set_ylabel(r"$K_J$ at load-event onset (MPa$\sqrt{m}$)")
            ax.set_title(f"Cascade-aware propagation events, {float(T_K):g} K")
            ax.grid(alpha=0.25)
            ax.legend(frameon=False)
            fig.tight_layout()
            fig.savefig(root / "R_curve_cascade_aware.png", dpi=220)
            plt.close(fig)
    except Exception:
        pass
    return metrics


__all__ = [
    "cascade_metrics",
    "cluster_same_load_events",
    "extract_raw_growth_events",
    "find_steps_file",
    "write_cascade_aware_outputs",
]
