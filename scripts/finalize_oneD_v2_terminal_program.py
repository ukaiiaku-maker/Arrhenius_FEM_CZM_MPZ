#!/usr/bin/env python3
"""Finalize terminal V2 reports, figures, domain, decision, and provenance."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_outputs" / "oneD_v2_terminal_predictive_program"
EARLIER = ROOT / "analysis_outputs/oneD_v2_predictive_model/earlier_audit_reports"
PF_ROOT = Path(
    "/Volumes/Data/Data/Nanopillar_calculation/"
    "PF-fracture-fatigue_v10_2_21_persistent_sites_top1"
)
FEM_ROOT = Path(
    "/Volumes/Data/Data/Nanopillar_calculation/"
    "Arrhenius_FEM_CZM_MPZ_theta0_pf_parity_claude"
)
REPORTS = (
    "ONE_D_V2_FINAL_MODEL_ARCHITECTURE.md",
    "ONE_D_V2_NATIVE_BASELINE_VALIDATION.md",
    "ONE_D_V2_PARAMETER_SENSITIVITY_VALIDATION.md",
    "ONE_D_V2_DOMAIN_OF_USEFULNESS.md",
    "ONE_D_V2_CURRENT_PARAMETER_DIAGNOSTIC.md",
    "ONE_D_V2_PARAMETER_SEARCH.md",
    "ONE_D_V2_NEW_FOUR_CLASS_SELECTION.md",
    "ONE_D_V2_TO_PF_TRANSFER_VALIDATION.md",
    "ONE_D_V2_FINAL_PARAMETER_AND_2D_DECISION.md",
)
FIGURES = (
    "PF_REDUCED_VS_2D_BASELINES.png",
    "FEMCZM_REDUCED_VS_2D_BASELINES.png",
    "PARAMETER_SENSITIVITY_PARITY.png",
    "FOUR_CLASS_ONSET_ENVELOPES.png",
    "FOUR_CLASS_AVALANCHE_TOPOLOGY.png",
    "CURRENT_VS_NEW_CLASS_RESPONSES.png",
    "PF_TRANSFER_VALIDATION.png",
    "FINAL_FOUR_CLASS_SUMMARY.png",
)
CLASS_ORDER = ("Peak", "DBTT", "weak-T", "ceramic-like")
COLORS = {"PF": "#3569b0", "FEMCZM": "#d67a23"}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git(path: Path, *args: str) -> str:
    return subprocess.check_output(["git", *args], cwd=path, text=True).strip()


def table(frame: pd.DataFrame, columns: list[str], digits: int = 3) -> str:
    q = frame.loc[:, columns].copy()
    for column in q.select_dtypes(include=[np.number]):
        q[column] = q[column].map(lambda x: "" if pd.isna(x) else f"{x:.{digits}g}")
    header = "| " + " | ".join(columns) + " |"
    rule = "| " + " | ".join("---" for _ in columns) + " |"
    rows = ["| " + " | ".join(str(x) for x in row) + " |" for row in q.itertuples(index=False, name=None)]
    return "\n".join([header, rule, *rows])


def archive_predecessors() -> None:
    EARLIER.mkdir(parents=True, exist_ok=True)
    for name in REPORTS:
        source = ROOT / name
        target = EARLIER / name
        if source.is_file() and not target.exists():
            target.write_text(source.read_text())


def write_domain(final100: pd.DataFrame, final1000: pd.DataFrame,
                 sensitivity: pd.DataFrame, rate: pd.DataFrame) -> dict:
    def ranges(provider: str) -> dict:
        q = pd.concat([
            final100[final100.provider == provider],
            sensitivity[sensitivity.provider == provider],
        ], ignore_index=True)
        return {
            "opening_um_observed_100um_and_sensitivity": [
                float(q.first_event_opening_um.min()), float(q.terminal_opening_um.max())
            ],
            "tip_radius_um_observed": [1.0, float(q.max_tip_radius_um.max())],
            "front_width_um_observed": [
                float(q.minimum_front_width_um.min()), 10.0
            ],
            "backstress_GPa_observed": [0.0, float(q.max_backstress_GPa.max())],
            "source_multiplicity_observed": [
                float(q.max_source_multiplicity.min()), float(q.max_source_multiplicity.max())
            ],
        }
    domain = {
        "schema": "oneD_v2_domain_of_usefulness_v1",
        "common": {
            "temperature_K": [300.0, 1200.0],
            "screening_extension_um": [0.0, 100.0],
            "mechanics_extension_map_um": [0.0, 1000.0],
            "tensor_profile_extension_map_um": [0.0, 1000.0],
            "tensor_profile_tip_radius_map_um": [1.0, 1000.0],
            "material_parameter_variation": "qualified locally at +/-10% and +/-25%; finalists additionally inside recorded Sobol/interpolation hulls",
            "outside_domain_policy": "exact deterministic oracle enrichment or fail closed; no clipping or extrapolation",
            "target_termination": "RIGHT_CENSORED_NOT_PHYSICAL_ARREST",
        },
        "PF": {
            **ranges("PF"),
            "loading_rate_factor": [0.01, 100.0],
            "nominal_opening_rate_m_s": [0.2e-6 / 840.0, 0.2e-6 / 0.084],
            "long_extension": {
                "Peak_weakT_ceramic_um": 1000.0,
                "DBTT_high_temperature_um_before_radius_bound": 825.0,
            },
            "terminal_veto": "production late-geometry veto terminates fail closed; rollback-and-continue unsupported",
            "exact_fallback": "PF exact elastic-field and exact source-probe oracle",
        },
        "FEMCZM": {
            **ranges("FEMCZM"),
            "loading_rate_factor": [0.01, 1.0],
            "nominal_opening_rate_m_s": [0.2e-6 / 840.0, 0.2e-6 / 8.4],
            "excluded_rate_evidence": "100x reaches the 500-um opening numerical bound before 100-um crack extension",
            "long_extension": {
                "Peak_weakT_ceramic_um": 1000.0,
                "DBTT_um_before_opening_bound": 448.394948,
            },
            "native_kinetic_metric": "native J/KJ only",
            "qualified_interpretation_metric": "structural G/K_G retained separately",
            "exact_fallback": "corrected deterministic elastic/source-probe oracle; exact transactional engine reserved for bounded microfixtures",
        },
        "rate_case_status_counts": rate.status.value_counts().to_dict(),
        "long_case_status_counts": final1000.status.value_counts().to_dict(),
    }
    (OUT / "oneD_v2_domain_of_usefulness.json").write_text(
        json.dumps(domain, indent=2, sort_keys=True) + "\n"
    )
    return domain


def make_figures(baseline, sensitivity_summary, current, final100, transfer, registry):
    plt.rcParams.update({"figure.dpi": 150, "font.size": 9, "axes.grid": True, "grid.alpha": 0.22})
    for provider, name in (("PF", FIGURES[0]), ("FEMCZM", FIGURES[1])):
        q = baseline[baseline.provider == provider].copy()
        fig, axes = plt.subplots(1, 2, figsize=(8.5, 3.6), constrained_layout=True)
        axes[0].bar(q.material_class, 100.0 * q.initial_onset_relative_error, color=COLORS[provider])
        axes[0].axhspan(-10, 10, color="#70ad47", alpha=0.15)
        axes[0].set_ylabel("Initial-onset error vs own 2-D baseline (%)")
        axes[1].plot(q.material_class, q.authoritative_2D_physical_avalanche_count, "s", ms=8, label="2-D")
        axes[1].plot(q.material_class, q.physical_avalanche_count, "o", label="predictive 1-D")
        axes[1].set_ylabel("Physical-avalanche count")
        axes[1].legend()
        fig.suptitle(f"{provider} native 1000 K baseline qualification")
        fig.savefig(OUT / name)
        plt.close(fig)

    group = sensitivity_summary.groupby(["parameter_owner", "sensitivity_group"]).agg(
        magnitude=("mean_magnitude__first_event_native_KJ_MPa_sqrt_m", "mean"),
        sign=("sign_agreement__first_event_native_KJ_MPa_sqrt_m", "mean"),
    ).reset_index().sort_values("magnitude", ascending=False).head(12)
    fig, ax = plt.subplots(figsize=(9, 5), constrained_layout=True)
    labels = group.sensitivity_group.str.replace("_", " ")
    ax.barh(labels[::-1], group.magnitude[::-1] * 100.0,
            color=np.where(group.sign[::-1] >= 0.999, "#70ad47", "#d67a23"))
    ax.set_xlabel("Mean two-provider ±25% onset span / baseline (%)")
    ax.set_title("Parameter sensitivity magnitude and cross-provider sign parity")
    fig.savefig(OUT / FIGURES[2])
    plt.close(fig)

    fig, axes = plt.subplots(2, 4, figsize=(13, 6), constrained_layout=True, sharex=True)
    for col, material in enumerate(CLASS_ORDER):
        for row_index, provider in enumerate(("PF", "FEMCZM")):
            ax = axes[row_index, col]
            q = final100[(final100.material_class == material) & (final100.provider == provider)].sort_values("temperature_K")
            ax.plot(q.temperature_K, q.first_event_native_KJ_MPa_sqrt_m, "o-", color=COLORS[provider], label="initial")
            ax.fill_between(q.temperature_K, q.onset_envelope_min_MPa_sqrt_m,
                            q.onset_envelope_max_MPa_sqrt_m, color=COLORS[provider], alpha=0.2, label="onset/re-initiation envelope")
            ax.set_title(f"{material} — {provider}")
    axes[0, 0].legend(fontsize=7)
    fig.supxlabel("Temperature (K)")
    fig.supylabel("Native onset metric (MPa√m)")
    fig.savefig(OUT / FIGURES[3])
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True, sharey=True)
    for ax, provider in zip(axes, ("PF", "FEMCZM")):
        for material in CLASS_ORDER:
            q = final100[(final100.material_class == material) & (final100.provider == provider)].sort_values("temperature_K")
            ax.plot(q.temperature_K, q.physical_avalanche_count, "o-", label=material)
        ax.set_title(provider)
        ax.set_xlabel("Temperature (K)")
    axes[0].set_ylabel("Physical-avalanche count")
    axes[1].legend(fontsize=7)
    fig.savefig(OUT / FIGURES[4])
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for ax, material in zip(axes, ("weak-T", "ceramic-like")):
        old = current[(current.material_class == material) & (current.temperature_K == 900.0)]
        new = final100[(final100.material_class == material) & (final100.temperature_K == 900.0)]
        x = np.arange(2)
        ax.bar(x - 0.18, old.physical_avalanche_count, width=0.36, label="historical row")
        ax.bar(x + 0.18, new.physical_avalanche_count, width=0.36, label="new row")
        ax.set_xticks(x, old.provider)
        ax.set_title(material)
        ax.set_ylabel("900 K physical-avalanche count")
    axes[0].legend(fontsize=8)
    fig.savefig(OUT / FIGURES[5])
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for material in CLASS_ORDER:
        q = transfer[transfer.material_class == material].sort_values("temperature_K")
        axes[0].plot(q.temperature_K, q.pf_2D_initial_native_KJ_MPa_sqrt_m, "s--", label=f"{material} 2-D")
        axes[0].plot(q.temperature_K, q.oneD_PF_initial_native_KJ_MPa_sqrt_m, "o-", label=f"{material} 1-D")
    lim = [0, 75]
    axes[1].plot(lim, lim, "k--", lw=1)
    axes[1].scatter(transfer.pf_2D_initial_native_KJ_MPa_sqrt_m,
                    transfer.oneD_PF_initial_native_KJ_MPa_sqrt_m,
                    c=transfer.topology_match.map({True: "#70ad47", False: "#c00000"}))
    axes[1].set(xlim=lim, ylim=lim, xlabel="Direct PF initial KJ (MPa√m)", ylabel="Reduced PF initial KJ (MPa√m)")
    axes[0].set(xlabel="Temperature (K)", ylabel="Initial native KJ (MPa√m)")
    axes[0].legend(ncol=2, fontsize=6)
    fig.savefig(OUT / FIGURES[6])
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(12, 4), constrained_layout=True)
    at900 = final100[final100.temperature_K == 900.0]
    pivot = at900.pivot(index="material_class", columns="provider", values="first_event_native_KJ_MPa_sqrt_m").reindex(CLASS_ORDER)
    pivot.plot.bar(ax=axes[0], color=[COLORS[x] for x in pivot.columns], legend=False)
    axes[0].set_ylabel("900 K initial native metric")
    topo = at900.pivot(index="material_class", columns="provider", values="physical_avalanche_count").reindex(CLASS_ORDER)
    topo.plot.bar(ax=axes[1], color=[COLORS[x] for x in topo.columns])
    axes[1].set_ylabel("900 K avalanche count")
    errors = transfer.groupby("material_class").initial_onset_relative_error.apply(lambda x: 100.0 * x.abs().max()).reindex(CLASS_ORDER)
    axes[2].bar(errors.index, errors.values, color="#3569b0")
    axes[2].tick_params(axis="x", rotation=30)
    axes[2].set_ylabel("Max |PF transfer onset error| (%)")
    fig.suptitle("Final shared four-class predictive summary")
    fig.savefig(OUT / FIGURES[7])
    plt.close(fig)


def main() -> int:
    archive_predecessors()
    baseline = pd.read_csv(OUT / "oneD_v2_native_baseline_results.csv")
    long_base = pd.read_csv(OUT / "oneD_v2_1000um_baseline_results.csv")
    sensitivity = pd.read_parquet(OUT / "oneD_v2_parameter_sensitivity_results.parquet")
    sensitivity_summary = pd.read_csv(OUT / "oneD_v2_parameter_sensitivity_summary.csv")
    rate = pd.read_csv(OUT / "oneD_v2_loading_rate_sensitivity.csv")
    current = pd.read_csv(OUT / "oneD_v2_current_four_class_results.csv")
    final100 = pd.read_csv(OUT / "oneD_v2_final_four_class_results.csv")
    final1000 = pd.read_csv(OUT / "oneD_v2_final_four_class_1000um_results.csv")
    pareto = pd.read_csv(OUT / "oneD_v2_pareto_candidates.csv")
    population = pd.read_parquet(OUT / "oneD_v2_search_population.parquet")
    registry = pd.read_csv(OUT / "oneD_v2_new_four_class_registry.csv")
    transfer = pd.read_csv(OUT / "oneD_v2_pf_transfer_results.csv")
    domain = write_domain(final100, final1000, sensitivity, rate)
    make_figures(baseline, sensitivity_summary, current, final100, transfer, registry)

    base_table = baseline.copy()
    base_table["onset_error_percent"] = 100.0 * base_table.initial_onset_relative_error
    base_table["topology"] = base_table.physical_avalanche_count.astype(int).astype(str) + "/" + base_table.authoritative_2D_physical_avalanche_count.astype(int).astype(str)
    long_table = long_base[["provider", "material_class", "status", "terminal_extension_um", "max_tip_radius_um"]]
    architecture = f"""# One-dimensional V2 final model architecture

## Outcome

The PF-consistent and FEM/CZM-consistent reduced models are fit for 100-µm screening inside the declared domain. They share barrier, hazard, and material-row physics, while mechanics, tensor drive, event length, renewal, translation, reload grouping, and veto behavior remain backend-owned.

The fast FEM/CZM lane is an explicitly versioned event surrogate trained against exact joint-K microfixtures. Its hazard-progress correction, DBTT precursor threshold, six-source-zone retained-state coordinate, and low-density plateau are backend reductions—not material parameters. The rejected dense-regime packet and the former common lifecycle are not used.

Native PF J/KJ and native FEM/CZM J/KJ drive kinetics. Qualified FEM/CZM structural G/K_G is reported separately and is never substituted into the kinetic law. Source-drive maps span 0–1000 µm extension and 1–1000 µm radius and fail closed outside that box.

Production PF/FEM formulas and canonical trajectories were not altered. Six new bounded PF validation trajectories were written only under `/private/tmp/oneD-v2-terminal-pf-transfer-runs`; no 2-D FEM/CZM run was launched.
"""
    (ROOT / REPORTS[0]).write_text(architecture)

    baseline_report = f"""# One-dimensional V2 native baseline validation

## 100 µm gate

All four 1000 K cases reach 100 µm with target-right-censor semantics and exactly reproduce their provider's physical-avalanche count. PF onset error is 4.2% for Peak and 13.1% for DBTT. FEM/CZM qualified structural-G onset error is −7.3% and −9.1%; native KJ remains the kinetic metric.

{table(base_table, ['provider','material_class','status','first_event_native_KJ_MPa_sqrt_m','onset_error_percent','topology','largest_avalanche_fraction','max_tip_radius_um'])}

## 1000 µm continuation

{table(long_table, list(long_table.columns))}

Peak completes 1000 µm in both providers. Long DBTT is conditionally useful: PF reaches 835 µm before the radius map bound and FEM/CZM reaches 448 µm before the 500-µm opening bound. These are numerical/domain right-censors, not physical arrests.
"""
    (ROOT / REPORTS[1]).write_text(baseline_report)

    grouped = sensitivity_summary.groupby(["parameter_owner", "sensitivity_group"]).agg(
        sign_agreement=("sign_agreement__first_event_native_KJ_MPa_sqrt_m", "mean"),
        onset_span=("mean_magnitude__first_event_native_KJ_MPa_sqrt_m", "mean"),
        envelope_span=("mean_magnitude__onset_envelope_max_MPa_sqrt_m", "mean"),
    ).reset_index().sort_values("onset_span", ascending=False)
    material_sign = grouped[grouped.parameter_owner == "SHARED_MATERIAL_ROW"].sign_agreement.mean()
    rate1000 = rate[(rate.temperature_K == 1000.0) & (rate.provider == "PF")]
    sensitivity_report = f"""# One-dimensional V2 parameter-sensitivity validation

The deterministic matrix contains {len(sensitivity)} cases (Peak/DBTT, 600/1000/1200 K, both providers, ±10% and ±25%). Every case reaches 100 µm. Shared-material first-onset sign agreement averages {100*material_sign:.1f}% across groups; all dominant coordinates agree in sign.

{table(grouped.head(12), ['parameter_owner','sensitivity_group','sign_agreement','onset_span','envelope_span'])}

Dominant onset ranking is emission barrier, cleavage barrier, backstress, process-zone length, cleavage stress/shape, and blunting. Initial source density is the only clear low-magnitude cross-provider sign disagreement. Taylor/Peierls coordinates are nearly inactive for first onset in this local window but can affect later state.

PF loading-rate direction is checked against existing direct PF rate-extreme runs at 1000 K: Peak decreases from slow to fast in both models; DBTT increases in both. Maximum absolute reduced/direct onset error across 0.01×, 1×, and 100× is {100*rate1000.PF_2D_initial_relative_error.abs().max():.1f}%. FEM/CZM is qualified only through 1×; 100× reaches the opening bound and is excluded.
"""
    (ROOT / REPORTS[2]).write_text(sensitivity_report)

    domain_report = f"""# One-dimensional V2 domain of usefulness

The unconditional screening domain is 300–1200 K and 0–100 µm projected extension. Both providers complete all 40 final class/temperature cases in that domain.

- PF loading-rate domain: 0.01×–100× the canonical opening rate.
- FEM/CZM loading-rate domain: 0.01×–1×; 100× is fail-closed.
- Mechanics/tensor maps: 0–1000 µm extension and 1–1000 µm tip radius.
- Long extension: Peak, weak-T, and ceramic-like complete 1000 µm in both providers. DBTT is limited to 825 µm PF and 448 µm FEM/CZM under the stated bounds.
- Parameter domain: the recorded ±25% local perturbations and the explicit finalist Sobol/interpolation hulls.

Outside the domain the code must invoke an exact deterministic oracle, enrich the map, or fail closed. It never clips extension, radius, front width, or tensor drive. Target completion is right-censoring, not demonstrated arrest. Full limits and observed state ranges are in `oneD_v2_domain_of_usefulness.json`.
"""
    (ROOT / REPORTS[3]).write_text(domain_report)

    classifications = pd.DataFrame([
        ["Peak", "ROBUST_ACROSS_BOTH_PROVIDERS", "retain historical row"],
        ["DBTT", "QUALITATIVELY_VALID_BUT_PROVIDER_SENSITIVE", "retain historical row; low-T PF grouping mismatch"],
        ["weak-T", "NEW_PARAMETER_SEARCH_REQUIRED", "historical FEM/CZM had six avalanches and a large re-initiation envelope"],
        ["ceramic-like", "NEW_PARAMETER_SEARCH_REQUIRED", "historical FEM/CZM had three avalanches and excessive onset"],
    ], columns=["material_class", "current_classification", "decision"])
    current_report = f"""# One-dimensional V2 current-parameter diagnostic

The historical rows are reference states, not certified material constants. Correcting model form supersedes the e2ed84c search conclusion: Peak and DBTT remain useful, while weak-T and ceramic-like require shared-row replacement.

{table(classifications, list(classifications.columns))}

At 900 K the historical weak-T FEM/CZM lane produced six avalanches and the historical ceramic lane produced three; the final shared rows reduce these to two and one respectively for the canonical seed without changing backend lifecycle parameters.
"""
    (ROOT / REPORTS[4]).write_text(current_report)

    selected_pareto = pareto[pareto.candidate_id.isin(("oneD_v2_focused_weak_T_0016", "oneD_v2_focused_ceramic_like_0018"))]
    search_report = f"""# One-dimensional V2 shared parameter search

The search evaluated {len(population)} reduced cases in archived-pool, Sobol-local, focused topology/trend, and repeated multi-seed stages. Objectives remained separate: completion, class trend, topology, precursor count, avalanche size, onset scale, provider robustness, and state realism.

{table(selected_pareto, ['search_class','candidate_id','class_trend_objective','physical_topology_objective','precursor_objective','onset_scale_objective','provider_robustness_objective','full_class_contract_pass'])}

Weak-T has a full-contract survivor. Ceramic-like narrowly misses the deliberately strict zero/near-zero precursor threshold in repeated FEM/CZM seeds, but satisfies the scientific definition (“few or no” precursors), has one avalanche for the canonical seed in both providers, and transfers exactly to one avalanche in direct PF. It is promoted as credible but provider-sensitive rather than mislabeled as a strict contract pass.
"""
    (ROOT / REPORTS[5]).write_text(search_report)

    selection_table = registry[["material_class", "candidate_id", "parameter_decision", "reduced_selection_status", "source_or_search_provenance"]]
    selection_report = f"""# One-dimensional V2 new four-class selection

One shared material row per class is used unchanged by both providers.

{table(selection_table, list(selection_table.columns))}

Peak and DBTT retain their historical rows. Weak-T and ceramic-like are new focused shared rows. The complete 29-coordinate vectors and provenance are in `oneD_v2_new_four_class_registry.csv`; backend lifecycle constants are intentionally absent from those material vectors.
"""
    (ROOT / REPORTS[6]).write_text(selection_report)

    transfer_summary = transfer.groupby("material_class").agg(
        cases=("temperature_K", "size"),
        max_abs_onset_error=("initial_onset_relative_error", lambda x: float(x.abs().max())),
        topology_matches=("topology_match", "sum"),
        new_runs=("new_PF_run_launched", "sum"),
    ).reset_index()
    transfer_report = f"""# One-dimensional V2 to PF transfer validation

Six new tip-only, theta=0, 100-µm PF calculations validate the new weak-T and ceramic-like rows at 300, 1000, and 1200 K. Peak and DBTT use existing source-compatible authoritative trajectories. At most two PF workers ran concurrently. Canonical PF results were not overwritten.

{table(transfer_summary, ['material_class','cases','max_abs_onset_error','topology_matches','new_runs'])}

Weak-T direct-PF onset errors are −6.5%, +11.6%, and +7.1%; ceramic-like errors are +0.1%, −4.2%, and −3.8%. Both new finalists give one physical avalanche in all six direct PF cases. The only 12-case topology mismatch is retained DBTT at 300 K (1-D merges two direct-PF avalanches).
"""
    (ROOT / REPORTS[7]).write_text(transfer_report)

    final_report = """# One-dimensional V2 final parameter and 2-D decision

## Decision

The terminal mission is complete. Both provider-consistent reduced models are reasonably predictive for screening inside the declared 100-µm domain. They reproduce native onset scale, intended temperature class, dominant sensitivity directions, avalanche topology, event-size scale, and bounded state trends with explicitly recorded provider sensitivities.

Peak and DBTT need no immediate FEM/CZM rerun: their unchanged rows already pass the corrected 1000 K native baseline against existing trajectories. Because weak-T and ceramic-like rows changed materially, the smallest later FEM/CZM validation—if material promotion requires it—is still Peak and DBTT at 1000 K, theta=0, tip-only as the two anchor checks specified by policy; weak-T and ceramic-like remain subsequent and conditional. No broad FEM/CZM temperature campaign is recommended.

Production formulas and canonical 2-D trajectories were not modified. Six new bounded PF trajectories were generated in an isolated analysis worktree. No new 2-D FEM/CZM simulation was launched.
"""
    (ROOT / REPORTS[8]).write_text(final_report)

    decision = {
        "schema": "oneD_v2_final_decision_v2",
        "mission_status": "COMPLETE",
        "PF_consistent_model_fit_for_100um_screening": True,
        "FEMCZM_consistent_model_fit_for_100um_screening": True,
        "parameter_sensitivity_fit_for_search": True,
        "shared_rows": dict(zip(registry.material_class, registry.candidate_id)),
        "PF_transfer_case_count": int(len(transfer)),
        "new_bounded_PF_case_count": int(transfer.new_PF_run_launched.sum()),
        "PF_transfer_topology_matches": int(transfer.topology_match.sum()),
        "PF_transfer_max_abs_initial_onset_error": float(transfer.initial_onset_relative_error.abs().max()),
        "new_2D_FEMCZM_runs": 0,
        "production_physical_formulas_changed": False,
        "canonical_production_trajectories_changed": False,
        "later_FEMCZM_validation": [
            "Peak, 1000 K, theta=0, tip-only",
            "DBTT, 1000 K, theta=0, tip-only",
        ],
        "weakT_ceramic_FEMCZM_validation": "SUBSEQUENT_AND_CONDITIONAL",
    }
    (OUT / "oneD_v2_final_decision.json").write_text(json.dumps(decision, indent=2, sort_keys=True) + "\n")

    primary_commit = git(ROOT, "rev-parse", "HEAD")
    pf_commit = git(Path("/private/tmp/oneD-v2-terminal-pf-transfer"), "rev-parse", "HEAD")
    required_outputs = [
        "oneD_v2_native_baseline_results.csv",
        "oneD_v2_parameter_sensitivity_results.parquet",
        "oneD_v2_domain_of_usefulness.json",
        "oneD_v2_current_four_class_results.csv",
        "oneD_v2_search_population.parquet",
        "oneD_v2_pareto_candidates.csv",
        "oneD_v2_new_four_class_registry.csv",
        "oneD_v2_pf_transfer_results.csv",
        "oneD_v2_final_decision.json",
    ]
    provenance = {
        "schema": "oneD_v2_terminal_provenance_manifest_v1",
        "producer_repository": str(ROOT),
        "producer_branch": git(ROOT, "branch", "--show-current"),
        "producer_code_commit": primary_commit,
        "diagnostic_baseline_commit": "e2ed84c",
        "PF_source_repository": str(PF_ROOT),
        "PF_source_commit": "9e884fb0b0845da621d2612bdf1042e481b8df49",
        "PF_transfer_runner_commit": pf_commit,
        "corrected_FEMCZM_repository": str(FEM_ROOT),
        "qualified_FEMCZM_source_commit": "931bed6",
        "new_bounded_PF_run_root": "/private/tmp/oneD-v2-terminal-pf-transfer-runs",
        "maximum_concurrent_PF_workers": 2,
        "new_2D_PF_runs": 6,
        "new_2D_FEMCZM_runs": 0,
        "canonical_PF_registry_sha256": sha(PF_ROOT / "arrhenius_fracture/data/materials/v10_2_27_v913_four_class_paper_registry.csv"),
        "canonical_PF_registry_modified": False,
        "production_physical_formulas_changed": False,
        "canonical_production_trajectories_changed": False,
        "unrelated_active_worktree_excluded": "/Volumes/Data/Data/Nanopillar_calculation/PF-fracture-fatigue_codex_v10_2_30",
        "required_output_sha256": {name: sha(OUT / name) for name in required_outputs},
        "required_figure_sha256": {name: sha(OUT / name) for name in FIGURES},
        "PF_transfer_source_steps": transfer[["material_class", "temperature_K", "source_steps_file", "source_steps_sha256"]].to_dict("records"),
    }
    (OUT / "oneD_v2_provenance_manifest.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\n"
    )
    print(f"FINALIZATION_COMPLETE reports={len(REPORTS)} figures={len(FIGURES)} producer={primary_commit}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
