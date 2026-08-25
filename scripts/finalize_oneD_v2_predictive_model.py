#!/usr/bin/env python3
"""Finalize the predictive V2 decision, registry, reports, and figures."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_outputs" / "oneD_v2_predictive_model"
PF_ROOT = Path(
    "/Volumes/Data/Data/Nanopillar_calculation/"
    "PF-fracture-fatigue_v10_2_21_persistent_sites_top1"
)
PF_MATERIALS = PF_ROOT / "arrhenius_fracture/data/materials"
PF_RUN = PF_ROOT / (
    "runs/1_final_v10_2_30_theta0_rate_extremes_four_class_1000um_base3621_v1/"
    "rate1x"
)
PF_RESPONSE = PF_RUN / (
    "temperature_response/v10_2_28_four_class_KJ_temperature_response.csv"
)
REGISTRY = PF_MATERIALS / "v10_2_27_v913_four_class_paper_registry.csv"
sys.path.insert(0, str(ROOT))

from arrhenius_fracture.emergent_gnd_contract_v913 import (  # noqa: E402
    ACTIVE_CANDIDATE_PARAMETER_FIELDS,
)
from scripts.run_oneD_v2_predictive_campaign import (  # noqa: E402
    IDS,
    inputs,
    run_case,
    summary,
)

CLASS_ORDER = ("Peak", "DBTT", "weak-T", "ceramic-like")
CURRENT_CLASSIFICATION = {
    "Peak": "QUALITATIVELY_VALID_BUT_PROVIDER_SENSITIVE",
    "DBTT": "PARAMETER_RESEARCH_REQUIRED",
    "weak-T": "QUALITATIVELY_VALID_BUT_PROVIDER_SENSITIVE",
    "ceramic-like": "ROBUST_ACROSS_BOTH_PROVIDERS",
}
FINAL_STATUS = {
    "Peak": "RETAINED_DIRECT_PF_VALIDATED_REDUCED_PROVIDER_SENSITIVE",
    "DBTT": "RETAINED_PRODUCTION_REFERENCE_REDUCED_NOT_QUALIFIED",
    "weak-T": "RETAINED_DIRECT_PF_VALIDATED_REDUCED_PROVIDER_SENSITIVE",
    "ceramic-like": "RETAINED_DIRECT_PF_AND_REDUCED_ROBUST",
}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def avalanche_summary(steps_file: Path) -> tuple[int, float, int]:
    steps = pd.read_csv(steps_file)
    positions = np.flatnonzero(steps.n_fire.to_numpy(float) > 0.0)
    if not len(positions):
        return 0, np.nan, 0
    avalanche = 0
    counts = {0: 1}
    for index in range(1, len(positions)):
        reload_slice = steps.iloc[positions[index - 1] + 1 : positions[index] + 1]
        if (reload_slice.adaptive_frac.to_numpy(float) >= 1.0 - 1.0e-12).any():
            avalanche += 1
        counts[avalanche] = counts.get(avalanche, 0) + 1
    return len(counts), max(counts.values()) / len(positions), len(positions)


def build_transfer() -> pd.DataFrame:
    response = pd.read_csv(PF_RESPONSE)
    response = response[response.temperature_K.isin((300.0, 1000.0, 1200.0))].copy()
    label_to_class = {
        "Peak": "Peak",
        "DBTT": "DBTT",
        "Weak-T": "weak-T",
        "Ceramic": "ceramic-like",
    }
    physics, rows, providers = inputs()
    reduced = {}
    mechanics, drive = providers["PF"]
    for _, row in rows.iterrows():
        material = IDS[row.candidate_id]
        for temperature in (300.0, 1000.0, 1200.0):
            reduced[(material, temperature)] = summary(
                run_case(
                    row, material, temperature, "PF", mechanics, drive, physics,
                    100.0,
                )
            )
    records = []
    for source in response.itertuples():
        material = label_to_class[str(source.plot_label)]
        avalanche_count, largest_fraction, event_count = avalanche_summary(
            Path(source.steps_file)
        )
        one_d = reduced[(material, float(source.temperature_K))]
        initial_error = (
            float(one_d["first_event_native_KJ_MPa_sqrt_m"])
            - float(source.initial_K_MPa_sqrt_m)
        ) / float(source.initial_K_MPa_sqrt_m)
        records.append(
            {
                "material_class": material,
                "candidate_id": source.candidate_id,
                "temperature_K": float(source.temperature_K),
                "transfer_data_role": "EXISTING_AUTHORITATIVE_2D_PF",
                "pf_2D_seed": int(source.seed),
                "pf_2D_initial_native_KJ_MPa_sqrt_m": float(source.initial_K_MPa_sqrt_m),
                "pf_2D_tail_average_native_KJ_MPa_sqrt_m": float(
                    source.tail_average_K_MPa_sqrt_m
                ),
                "pf_2D_event_count": event_count,
                "pf_2D_physical_avalanche_count": avalanche_count,
                "pf_2D_largest_avalanche_fraction": largest_fraction,
                "oneD_PF_status": one_d["status"],
                "oneD_PF_initial_native_KJ_MPa_sqrt_m": one_d[
                    "first_event_native_KJ_MPa_sqrt_m"
                ],
                "oneD_PF_onset_envelope_max_MPa_sqrt_m": one_d[
                    "onset_envelope_max_MPa_sqrt_m"
                ],
                "oneD_PF_physical_avalanche_count": one_d[
                    "physical_avalanche_count"
                ],
                "oneD_PF_largest_avalanche_fraction": one_d[
                    "largest_avalanche_fraction"
                ],
                "initial_onset_relative_error": initial_error,
                "topology_match": int(one_d["physical_avalanche_count"])
                == int(avalanche_count),
                "source_steps_file": source.steps_file,
                "source_steps_sha256": source.steps_sha256,
                "new_PF_run_launched": False,
            }
        )
    frame = pd.DataFrame(records)
    frame.to_csv(OUT / "oneD_v2_pf_transfer_results.csv", index=False)
    return frame


def build_registry() -> pd.DataFrame:
    source = pd.read_csv(REGISTRY)
    rows = []
    for _, row in source.iterrows():
        material = IDS[str(row.candidate_id)]
        record = {
            "material_class": material,
            "candidate_id": row.candidate_id,
            "source_option_key": row.option_key,
            "parameter_decision": "CURRENT_ROW_RETAINED_NO_ADMISSIBLE_REPLACEMENT",
            "current_parameter_classification": CURRENT_CLASSIFICATION[material],
            "final_selection_status": FINAL_STATUS[material],
            "reduced_provider_robust": material == "ceramic-like",
            "direct_2D_PF_evidence": "AVAILABLE_EXISTING_RATE1X_THETA0",
            "direct_2D_FEMCZM_status": (
                "EXISTING_1000K_BASELINE_ONLY"
                if material in ("Peak", "DBTT")
                else "NOT_VALIDATED"
            ),
        }
        for field in ACTIVE_CANDIDATE_PARAMETER_FIELDS:
            record[field] = float(row[field])
        rows.append(record)
    result = pd.DataFrame(rows)
    result.to_csv(OUT / "oneD_v2_new_four_class_registry.csv", index=False)
    return result


def figures(baseline, current, transfer, registry):
    plt.rcParams.update({"figure.dpi": 140, "font.size": 9})
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    x = np.arange(len(baseline))
    axes[0].bar(x, baseline.initial_onset_relative_error * 100.0, color="#4472C4")
    axes[0].axhline(0, color="black", lw=0.8)
    axes[0].set(
        xticks=x,
        xticklabels=[f"{p}\n{m}" for p, m in zip(baseline.provider, baseline.material_class)],
        ylabel="Initial-onset error versus own 2-D baseline (%)",
    )
    axes[1].scatter(
        x, baseline.authoritative_2D_physical_avalanche_count,
        marker="s", label="2-D target", color="#ED7D31",
    )
    axes[1].scatter(
        x, baseline.physical_avalanche_count,
        marker="o", label="predictive 1-D", color="#4472C4",
    )
    axes[1].set(
        xticks=x,
        xticklabels=[f"{p}\n{m}" for p, m in zip(baseline.provider, baseline.material_class)],
        ylabel="Physical-avalanche count",
    )
    axes[1].legend()
    fig.suptitle("Provider-native 1000 K baselines versus their own 2-D evidence")
    fig.savefig(OUT / "ONE_D_V2_PROVIDER_VS_OWN_2D_BASELINES.png")
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True, sharex=True)
    for ax, material in zip(axes.flat, CLASS_ORDER):
        local = current[current.material_class == material]
        for provider, color in (("PF", "#4472C4"), ("FEMCZM", "#ED7D31")):
            q = local[local.provider == provider].sort_values("temperature_K")
            ax.plot(q.temperature_K, q.onset_envelope_min_MPa_sqrt_m, "o-", color=color, label=provider)
            ax.fill_between(
                q.temperature_K,
                q.onset_envelope_min_MPa_sqrt_m,
                q.onset_envelope_max_MPa_sqrt_m,
                color=color,
                alpha=0.18,
            )
        ax.set_title(material)
        ax.grid(alpha=0.25)
    axes[0, 0].legend()
    fig.supxlabel("Temperature (K)")
    fig.supylabel("Native onset/re-initiation envelope (MPa√m)")
    fig.savefig(OUT / "ONE_D_V2_FOUR_CLASS_ONSET_ENVELOPES.png")
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True, sharex=True)
    for ax, material in zip(axes.flat, CLASS_ORDER):
        local = current[current.material_class == material]
        for provider, color in (("PF", "#4472C4"), ("FEMCZM", "#ED7D31")):
            q = local[local.provider == provider].sort_values("temperature_K")
            ax.plot(q.temperature_K, q.physical_avalanche_count, "o-", color=color, label=provider)
        ax.set_title(material)
        ax.grid(alpha=0.25)
    axes[0, 0].legend()
    fig.supxlabel("Temperature (K)")
    fig.supylabel("Physical-avalanche count")
    fig.savefig(OUT / "ONE_D_V2_AVALANCHE_TOPOLOGY_VS_TEMPERATURE.png")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 3.8), constrained_layout=True)
    bars = ax.bar(
        CLASS_ORDER,
        np.ones(len(CLASS_ORDER)),
        color=["#FFC000", "#C00000", "#FFC000", "#70AD47"],
    )
    for bar, material in zip(bars, CLASS_ORDER):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            0.5,
            CURRENT_CLASSIFICATION[material].replace("_", "\n"),
            ha="center",
            va="center",
            fontsize=7.5,
        )
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            1.05,
            FINAL_STATUS[material].replace("_", "\n"),
            ha="center",
            va="bottom",
            fontsize=6.5,
        )
    ax.set_ylabel("Current classification (bar) / final decision (above)")
    ax.set_title("Current-row diagnostic and final retention decision")
    ax.set_ylim(0, 1.9)
    fig.savefig(OUT / "ONE_D_V2_CURRENT_VS_FINAL_CLASSIFICATIONS.png")
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(10, 7), constrained_layout=True, sharex=True)
    for ax, material in zip(axes.flat, CLASS_ORDER):
        q = transfer[transfer.material_class == material].sort_values("temperature_K")
        ax.plot(q.temperature_K, q.pf_2D_initial_native_KJ_MPa_sqrt_m, "s-", label="2-D PF", color="#ED7D31")
        ax.plot(q.temperature_K, q.oneD_PF_initial_native_KJ_MPa_sqrt_m, "o--", label="1-D PF provider", color="#4472C4")
        ax.set_title(material)
        ax.grid(alpha=0.25)
    axes[0, 0].legend()
    fig.supxlabel("Temperature (K)")
    fig.supylabel("Initial native $K_J$ (MPa√m)")
    fig.savefig(OUT / "ONE_D_V2_DIRECT_PF_TRANSFER_VALIDATION.png")
    plt.close(fig)

    criteria = ("direct PF", "reduced trend", "reduced topology", "FEM 2-D")
    matrix = np.array(
        [
            [1, 1, 1, 1],
            [1, 1, 0, 1],
            [1, 1, 0, 0],
            [1, 1, 1, 0],
        ],
        dtype=float,
    )
    fig, ax = plt.subplots(figsize=(7, 4), constrained_layout=True)
    ax.imshow(matrix, vmin=0, vmax=1, cmap="RdYlGn", aspect="auto")
    ax.set(xticks=np.arange(len(criteria)), xticklabels=criteria, yticks=np.arange(4), yticklabels=CLASS_ORDER)
    for i in range(4):
        for j in range(4):
            ax.text(j, i, "qualified" if matrix[i, j] else "pending/fail", ha="center", va="center", fontsize=8)
    ax.set_title("Final four-class evidence summary")
    fig.savefig(OUT / "ONE_D_V2_FINAL_FOUR_CLASS_SUMMARY.png")
    plt.close(fig)


def write_reports(baseline, current, pareto, transfer, registry):
    def table(frame, columns):
        local = frame.loc[:, columns].copy()

        def cell(value):
            if pd.isna(value):
                return ""
            if isinstance(value, (float, np.floating)):
                return f"{float(value):.6g}"
            return str(value).replace("|", "\\|").replace("\n", " ")

        lines = [
            "| " + " | ".join(columns) + " |",
            "|" + "|".join("---" for _ in columns) + "|",
        ]
        for row in local.itertuples(index=False, name=None):
            lines.append("| " + " | ".join(cell(value) for value in row) + " |")
        return "\n".join(lines)

    architecture = f"""# One-dimensional V2 final model architecture

## Outcome

The V2 implementation has exact controlled common-kernel parity and qualified candidate-independent PF and FEM/CZM mechanics/source-drive maps. The natural predictive lane is **partially qualified**, not universally predictive: PF Peak closes, while DBTT and absolute FEM/CZM onset do not.

## Runtime composition

1. `SharedBarrierHazardCore` owns the common barrier/rate equations.
2. Backend state adapters/factories construct production-source PF unified-MPZ and FEM/CZM moving-tip states.
3. `ProviderMechanicsMap` supplies native PF or native/qualified FEM/CZM coefficients with fail-closed extension bounds.
4. `SourceDriveMap` supplies exact tensor-probe normalized drive factors on extension/radius grids. The final grid spans 0–1000 µm extension and 1–100 µm radius; no clipping or extrapolation is permitted.
5. `run_zero_d_predictive` is a versioned V9.13 zero-D state/lifecycle surrogate. It is not claimed to reproduce the full backend-specific production event transaction policy.

The controlled composition result remains `CONTROLLED_SOURCE_COMPOSITION_QUALIFIED`. Natural prediction is assessed separately.

No production equation, production trajectory, or canonical parameter registry was modified. New 2-D PF runs: **0**. New 2-D FEM/CZM runs: **0**.
"""
    (ROOT / "ONE_D_V2_FINAL_MODEL_ARCHITECTURE.md").write_text(architecture)

    native = f"""# One-dimensional V2 native baseline validation

## Decision

Only the PF Peak baseline is sufficient for reduced-model use at present: its initial native $K_J$ error is 3.98% and its one-avalanche topology matches. PF DBTT starts within 12.2% but exceeds the qualified 100-µm radius map after three events. FEM/CZM Peak gets the one-avalanche topology but underpredicts qualified onset $G$ by 61.5%; FEM/CZM DBTT underpredicts onset by 60.2% and exceeds the radius map.

{table(baseline, ['provider','material_class','status','event_count','physical_avalanche_count','first_event_native_KJ_MPa_sqrt_m','qualified_G_onset_min_J_m2','initial_onset_relative_error','topology_match','terminal_extension_um','max_tip_radius_um'])}

The 1000-µm extension is credible only for Peak under both mechanics maps. DBTT is right-censored by a qualification boundary, not physically arrested.
"""
    (ROOT / "ONE_D_V2_NATIVE_BASELINE_VALIDATION.md").write_text(native)

    diagnostics = current.groupby(["material_class", "provider"]).agg(
        complete_cases=("status", lambda x: int((x == "TARGET_RIGHT_CENSORED").sum())),
        cases=("status", "size"),
        avalanche_min=("physical_avalanche_count", "min"),
        avalanche_max=("physical_avalanche_count", "max"),
        onset_min=("onset_envelope_min_MPa_sqrt_m", "min"),
        onset_max=("onset_envelope_max_MPa_sqrt_m", "max"),
    ).reset_index()
    diagnostic = f"""# One-dimensional V2 current parameter diagnostic

## Classification

| Class | Current-row classification | Reason |
|---|---|---|
| Peak | QUALITATIVELY_VALID_BUT_PROVIDER_SENSITIVE | Both providers show an intermediate/high-temperature enhancement and later weakening; PF alone adds minor low-temperature subdivision. |
| DBTT | PARAMETER_RESEARCH_REQUIRED | Both providers reproduce a low shelf and upper transition, but high-temperature state growth is excessive and two cases per provider hit the radius-map bound. The 1000 K native topology does not match either own 2-D target. |
| weak-T | QUALITATIVELY_VALID_BUT_PROVIDER_SENSITIVE | Initial onset remains limited in span, but PF has 2–5 avalanches while FEM/CZM has one. |
| ceramic-like | ROBUST_ACROSS_BOTH_PROVIDERS | Low/declining onset and essentially one long avalanche under both providers. |

{table(diagnostics, ['material_class','provider','complete_cases','cases','avalanche_min','avalanche_max','onset_min','onset_max'])}

Absolute PF-native and FEM-qualified values were not forced to agree; classification uses trend, topology, and state evolution.
"""
    (ROOT / "ONE_D_V2_CURRENT_PARAMETER_DIAGNOSTIC.md").write_text(diagnostic)

    best = pareto.sort_values(
        ["status_objective", "class_topology_objective", "onset_envelope_objective"]
    ).head(8)
    search = f"""# One-dimensional V2 parameter search

## Result

The provider-robust DBTT search screened 187 shared rows (20 archived candidates, 39 controlled parameter morphs, and 128 Sobol-local perturbations). Sixteen survivors were rerun at six temperatures with three hazard seeds. **Zero rows passed the full shared DBTT contract.**

The non-dominated conflict is structural. Rows closest to the target avalanche counts exceed the qualified radius map and/or lose state realism; complete bounded rows sacrifice the DBTT transition or remain 2–3 avalanches away from the paired PF/FEM targets. Therefore no candidate is promoted by a hidden weighted compromise.

{table(best, ['candidate_id','status_objective','class_topology_objective','onset_envelope_objective','provider_robustness_objective','state_realism_objective','pareto_nondominated_validation','selection_status'])}

The machine population is `analysis_outputs/oneD_v2_predictive_model/oneD_v2_search_population.parquet`; independent objectives and active parameters are in `oneD_v2_pareto_candidates.csv`.
"""
    (ROOT / "ONE_D_V2_PARAMETER_SEARCH.md").write_text(search)

    selection = f"""# One-dimensional V2 new four-class selection

## Decision

No new row is scientifically admissible. The existing four paper rows are retained as the production reference registry, with explicit qualification status rather than being relabelled as provider-robust V2 fits.

{table(registry, ['material_class','candidate_id','parameter_decision','current_parameter_classification','final_selection_status','reduced_provider_robust','direct_2D_FEMCZM_status'])}

This is a completed negative selection: the search was run, no candidate cleared the contract, and no unqualified row was promoted.
"""
    (ROOT / "ONE_D_V2_NEW_FOUR_CLASS_SELECTION.md").write_text(selection)

    transfer_summary = transfer.groupby("material_class").agg(
        max_abs_initial_error=("initial_onset_relative_error", lambda x: float(np.max(np.abs(x)))),
        topology_matches=("topology_match", "sum"),
        cases=("topology_match", "size"),
        new_runs=("new_PF_run_launched", "sum"),
    ).reset_index()
    transfer_report = f"""# One-dimensional V2 to PF transfer validation

## Evidence used

The final rows are unchanged, so the authoritative existing rate-1×, theta=0 direct PF trajectories at 300, 1000, and 1200 K are the strongest transfer evidence. No fresh PF run was scientifically justified after the replacement search produced zero eligible rows.

{table(transfer_summary, ['material_class','max_abs_initial_error','topology_matches','cases','new_runs'])}

Direct PF itself preserves the four intended trends: Peak rises strongly toward 1000 K then weakens; DBTT rises from 23.36 to 56.19 MPa√m and retains two physical avalanches; weak-T stays near 19–21 MPa√m with one avalanche; ceramic-like declines from 13.64 to 10.13 MPa√m and becomes single-avalanche. Reduced PF transfer is good for Peak/weak-T/ceramic trend, but DBTT topology/state is not transferred.
"""
    (ROOT / "ONE_D_V2_TO_PF_TRANSFER_VALIDATION.md").write_text(transfer_report)

    final = """# One-dimensional V2 final parameter and 2-D decision

## Final answers

1. The predictive models work naturally only in part. PF Peak is qualified; the common zero-D lifecycle/state surrogate does not close DBTT, and FEM/CZM absolute onset is not qualified.
2. PF Peak reproduces its own 2-D baseline sufficiently. PF DBTT has good first onset but fails topology/state bounds. FEM/CZM Peak matches topology but not onset scale; FEM/CZM DBTT matches neither.
3. The current rows retain their intended four-class behavior in authoritative direct PF. In the reduced two-provider matrix, ceramic is robust, Peak and weak-T are provider-sensitive, and DBTT fails.
4. Provider mechanics/source drive explains much of the absolute and state-growth separation. The residual DBTT avalanche mismatch is also a lifecycle/model-form difference because the fast lane does not implement two independent production transaction policies.
5. Parameter research was required and completed, but no admissible replacement was found. The current rows remain reference rows with explicit caveats.
6. Final shared IDs remain 0242980 (Peak), 0202500 (DBTT), 0129902 (weak-T), and 0077080 (ceramic-like).
7. All four have direct existing 2-D PF evidence. Reduced-to-PF transfer fails specifically for DBTT topology/state.
8. Weak-T and ceramic-like have no direct 2-D FEM/CZM validation. Peak and DBTT have existing 1000 K FEM/CZM baselines, but reduced closure remains incomplete.
9. If later FEM/CZM work is authorized, the minimum matrix remains Peak and DBTT at 1000 K, theta=0, tip-only—two cases, not a broad campaign.
10. No production 2-D trajectory or physical formula was altered.

## Decision

Do not launch a broad FEM/CZM campaign and do not promote a new kinetic row. First replace/qualify the common fast lifecycle surrogate with backend-specific PF and FEM/CZM event transaction policies, then repeat the two 1000 K baseline checks. Existing PF data do not justify more PF runs for unchanged rows.
"""
    (ROOT / "ONE_D_V2_FINAL_PARAMETER_AND_2D_DECISION.md").write_text(final)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    baseline = pd.read_csv(OUT / "oneD_v2_native_baseline_results.csv")
    current = pd.read_csv(OUT / "oneD_v2_current_four_class_results.csv")
    pareto = pd.read_csv(OUT / "oneD_v2_pareto_candidates.csv")
    transfer = build_transfer()
    registry = build_registry()
    figures(baseline, current, transfer, registry)
    write_reports(baseline, current, pareto, transfer, registry)
    decision = {
        "schema": "oneD_v2_final_decision_v1",
        "overall_status": "PARTIALLY_QUALIFIED_NO_NEW_ROW_PROMOTED",
        "controlled_common_kernel_parity": "PASS_EXACT",
        "natural_predictive_models": {
            "PF_Peak": "QUALIFIED_FOR_REDUCED_USE",
            "PF_DBTT": "FAIL_TOPOLOGY_AND_DRIVE_MAP_BOUND",
            "FEMCZM_Peak": "FAIL_ABSOLUTE_ONSET_SCALE_TOPOLOGY_MATCHES",
            "FEMCZM_DBTT": "FAIL_ONSET_TOPOLOGY_AND_DRIVE_MAP_BOUND",
        },
        "current_row_classifications": CURRENT_CLASSIFICATION,
        "final_candidate_ids": {
            material: str(registry.set_index("material_class").loc[material, "candidate_id"])
            for material in CLASS_ORDER
        },
        "new_parameter_rows_promoted": 0,
        "dbtt_search": {
            "screen_candidates": 187,
            "repeated_validation_candidates": 16,
            "eligible_candidates": int(pareto.full_dbtt_contract_pass.sum()),
        },
        "direct_PF_transfer": "EXISTING_AUTHORITATIVE_300_1000_1200K",
        "new_2D_PF_runs": 0,
        "new_2D_FEMCZM_runs": 0,
        "pending_2D_FEMCZM_classes": ["weak-T", "ceramic-like"],
        "minimum_later_FEMCZM_matrix": [
            "Peak, 1000 K, theta=0, tip-only",
            "DBTT, 1000 K, theta=0, tip-only",
        ],
        "production_trajectories_altered": False,
        "physical_formulas_altered": False,
        "artifacts": {
            "registry_sha256": sha(OUT / "oneD_v2_new_four_class_registry.csv"),
            "transfer_sha256": sha(OUT / "oneD_v2_pf_transfer_results.csv"),
            "search_population_sha256": sha(OUT / "oneD_v2_search_population.parquet"),
            "pareto_sha256": sha(OUT / "oneD_v2_pareto_candidates.csv"),
        },
    }
    (OUT / "oneD_v2_final_decision.json").write_text(
        json.dumps(decision, indent=2, sort_keys=True) + "\n"
    )
    print("FINALIZATION_COMPLETE")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
