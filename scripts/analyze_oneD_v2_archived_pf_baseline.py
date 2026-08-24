#!/usr/bin/env python3
"""Freeze lifecycle provenance and classify authoritative archived PF events."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_outputs/oneD_v2_mechanics_maps_and_baselines"
PF_REPO = Path("/Volumes/Data/Data/Nanopillar_calculation/PF-fracture-fatigue_v10_2_21_persistent_sites_top1")
FEM_REPO = Path("/Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM_MPZ_theta0_pf_parity_claude")
RUNROOT = PF_REPO / "runs/1_final_v10_2_30_theta0_rate_extremes_four_class_1000um_base3621_v1/rate1x"
CASES = {
    "Peak": RUNROOT / "v913_paper_peak01_0242980_persistent_sites/T1000K_th0_seed8666",
    "DBTT": RUNROOT / "v913_paper_dbtt01_0202500_persistent_sites/T1000K_th0_seed1008666",
}
PF_SOURCE = "ab6279331d050919f78e2e7f278bc332466e34f8"
PF_OBSERVER = "c695b44"
FEM_SOURCE = "931bed66913afc970117bce900805ecc9b6225f8"
ANALYSIS_PARENT = "3fd33a3d10537ed36a0ca7220d5796ac4f2477d3"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fingerprint(row: pd.Series, fields: list[str]) -> str:
    payload = {field: None if pd.isna(row[field]) else float(row[field]) for field in fields}
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def freeze_lifecycle() -> None:
    cap = ROOT / "analysis_outputs/oneD_v2_capability_v2"
    manifest = {
        "schema": "oneD_v2_backend_lifecycle_manifest_v1",
        "analysis_parent_commit": ANALYSIS_PARENT,
        "state_fingerprint_schema": {
            "PF": "oneD_v2_pf_observer_state_v1/canonical_sorted_json_sha256",
            "FEMCZM": "bp5_joint_late_veto_v1043/recursive_value_sha256",
        },
        "PF": {
            "source_repository": str(PF_REPO), "source_commit": PF_SOURCE,
            "analysis_observer_commit": PF_OBSERVER,
            "test_artifact": str(cap / "oneD_v2_pf_subdivision_fixture_audit.json"),
            "test_artifact_sha256": sha(cap / "oneD_v2_pf_subdivision_fixture_audit.json"),
            "normal_interval_evolution": "QUALIFIED",
            "rejected_trial_state_preservation": "QUALIFIED",
            "accepted_event_commit": "QUALIFIED", "post_event_renewal": "QUALIFIED",
            "threshold_RNG_fields_not_present": "NOT_APPLICABLE",
            "late_veto_fail_closed_termination": "QUALIFIED",
            "late_veto_rollback_and_continue": "UNSUPPORTED_BY_PRODUCTION",
        },
        "FEMCZM": {
            "source_repository": str(FEM_REPO), "source_commit": FEM_SOURCE,
            "analysis_commit": ANALYSIS_PARENT,
            "test_artifact": str(cap / "oneD_v2_fem_transaction_lineage.json"),
            "test_artifact_sha256": sha(cap / "oneD_v2_fem_transaction_lineage.json"),
            "normal_interval_evolution": "QUALIFIED",
            "rejected_trial_state_preservation": "QUALIFIED",
            "accepted_event_commit": "QUALIFIED", "post_event_renewal": "QUALIFIED",
            "threshold_RNG_snapshot_restore": "QUALIFIED",
            "late_veto_rollback_and_continue": "QUALIFIED",
        },
        "reopen_policy": "REOPEN_ONLY_IF_RECORDED_SOURCE_HASH_DIFFERS",
    }
    (OUT / "oneD_v2_backend_lifecycle_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


def classify_case(material: str, directory: Path):
    steps = pd.read_csv(directory / "steps_1000K.csv")
    steps["physical_time_s"] = steps["dt_cur_s"].cumsum()
    event_positions = np.flatnonzero(steps["n_fire"].to_numpy() > 0)
    pz_fields = [
        "N_em", "N_em_pre_renewal", "N_em_shed_to_wake", "sigma_back_Pa",
        "mpz_mobile_count", "mpz_retained_count", "mpz_escaped_total",
        "mpz_recovered_total", "mpz_wake_retained_total",
    ]
    rows = []
    avalanche_id = 0
    for transaction_index, pos in enumerate(event_positions):
        event = steps.iloc[pos]
        before = steps.iloc[max(pos - 1, 0)]
        next_pos = event_positions[transaction_index + 1] if transaction_index + 1 < len(event_positions) else None
        reload_slice = steps.iloc[pos + 1 : next_pos + 1] if next_pos is not None else steps.iloc[0:0]
        certified_reload = bool(next_pos is not None and (reload_slice["adaptive_frac"] >= 1.0 - 1e-12).any())
        if transaction_index > 0:
            prior_pos = event_positions[transaction_index - 1]
            prior_reload = steps.iloc[prior_pos + 1 : pos + 1]
            if (prior_reload["adaptive_frac"] >= 1.0 - 1e-12).any():
                avalanche_id += 1
        rows.append({
            "material_class": material,
            "event_transaction_index": transaction_index,
            "physical_avalanche_index": avalanche_id,
            "pre_event_step": int(event["step"]),
            "pre_event_time_s": float(event["physical_time_s"]),
            "pre_event_opening_m": float(event["Uapp_m"]),
            "pre_event_reaction_N_per_m": float(event["Ftop_N"]),
            "pre_event_native_J_J_per_m2": float(event["J_effective_direct_J_per_m2"]),
            "pre_event_native_KJ_MPa_sqrt_m": float(event["KJ_Pa_sqrtm"] / 1e6),
            "pre_event_projected_extension_um": float(before["crack_extension_m"] * 1e6),
            "post_event_projected_extension_um": float(event["crack_extension_m"] * 1e6),
            "event_extension_um": float(event["da_block_m"] * 1e6),
            "post_event_tip_x_m": float(event["a_tip_m"]),
            "post_event_tip_y_m": float(pd.read_csv(directory / "fronts_1000K.csv").iloc[pos]["y_m"]),
            "post_event_process_zone_fingerprint": fingerprint(event, pz_fields),
            "reload_time_to_next_event_s": np.nan if next_pos is None else float(steps.iloc[next_pos]["physical_time_s"] - event["physical_time_s"]),
            "reload_opening_to_next_event_m": np.nan if next_pos is None else float(steps.iloc[next_pos]["Uapp_m"] - event["Uapp_m"]),
            "certified_reload_before_next_event": certified_reload,
            "right_censored_at_target": bool(transaction_index == len(event_positions) - 1),
            "source_steps_sha256": sha(directory / "steps_1000K.csv"),
        })
    tx = pd.DataFrame(rows)
    avalanches = []
    for aid, group in tx.groupby("physical_avalanche_index", sort=True):
        avalanches.append({
            "material_class": material, "physical_avalanche_index": int(aid),
            "first_event_transaction_index": int(group.event_transaction_index.min()),
            "last_event_transaction_index": int(group.event_transaction_index.max()),
            "event_transaction_count": int(len(group)),
            "start_extension_um": float(group.pre_event_projected_extension_um.iloc[0]),
            "end_extension_um": float(group.post_event_projected_extension_um.iloc[-1]),
            "avalanche_extension_um": float(group.event_extension_um.sum()),
            "onset_opening_m": float(group.pre_event_opening_m.iloc[0]),
            "target_right_censored": bool(group.right_censored_at_target.iloc[-1]),
            "grouping_rule": "CONTIGUOUS_EVENTS_WITHOUT_INTERVENING_FULL_ACCEPTED_LOADING_INTERVAL",
        })
    ava = pd.DataFrame(avalanches)
    onset = tx.groupby("physical_avalanche_index", sort=True).head(1).copy()
    onset["onset_role"] = np.where(onset.physical_avalanche_index == 0,
                                    "INITIAL_ONSET_PRE_EVENT", "REINITIATION_ONSET_PRE_EVENT")
    interior = tx[~tx.index.isin(onset.index)].copy()
    interior["semantic_role"] = "INTERIOR_PRE_EVENT_NOT_RESISTANCE_POINT"
    return tx, ava, onset, interior


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    freeze_lifecycle()
    tables = [classify_case(material, directory) for material, directory in CASES.items()]
    names = ["pf_2d_event_transactions_v2.csv", "pf_2d_physical_avalanches_v2.csv",
             "pf_2d_onset_candidates_v2.csv", "pf_2d_in_avalanche_drive_v2.csv"]
    combined = []
    for index, name in enumerate(names):
        frame = pd.concat([item[index] for item in tables], ignore_index=True)
        frame.to_csv(OUT / name, index=False)
        combined.append(frame)
    tx, ava, onset, interior = combined
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)
    for ax, material in zip(axes, CASES):
        subset = tx[tx.material_class == material]
        for aid, group in subset.groupby("physical_avalanche_index"):
            ax.plot(group.pre_event_projected_extension_um, group.pre_event_native_KJ_MPa_sqrt_m,
                    marker=".", ms=3, lw=1, label=f"avalanche {aid}")
        ax.set(title=material, xlabel="Pre-event projected extension (µm)",
               ylabel="PF native $K_J$ (MPa√m)")
        ax.grid(alpha=.25); ax.legend(fontsize=8)
    fig.suptitle("PF 2-D numerical events grouped into physical avalanches")
    fig.savefig(OUT / "PF_2D_EVENT_AND_AVALANCHE_BASELINE.png", dpi=180)
    plt.close(fig)
    summary = ava.groupby("material_class").agg(
        physical_avalanche_count=("physical_avalanche_index", "count"),
        numerical_event_count=("event_transaction_count", "sum"),
        largest_avalanche_events=("event_transaction_count", "max"),
    ).reset_index()
    lines = ["# PF 2-D event and avalanche baseline", "",
             "Archived rate-1×, theta-0, 1000 K trajectories were read without modification.", "",
             "| Class | Numerical events | Physical avalanches | Largest avalanche | Target status |",
             "|---|---:|---:|---:|---|"]
    for row in summary.itertuples():
        lines.append(f"| {row.material_class} | {row.numerical_event_count} | {row.physical_avalanche_count} | {row.largest_avalanche_events} events | right-censored |")
    lines += ["", "A physical-avalanche boundary requires an intervening full accepted loading interval (`adaptive_frac = 1`). This uses the production driver's own accepted-interval identity rather than an undocumented time threshold. Interior event states are model-native drive diagnostics, not resistance points.", ""]
    (OUT / "PF_2D_EVENT_AND_AVALANCHE_BASELINE.md").write_text("\n".join(lines))


if __name__ == "__main__":
    main()
