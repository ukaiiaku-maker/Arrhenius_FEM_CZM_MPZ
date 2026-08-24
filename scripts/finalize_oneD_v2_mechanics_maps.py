#!/usr/bin/env python3
"""Finalize map manifest and straight-versus-kinked path sensitivity."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_outputs/oneD_v2_mechanics_maps_and_baselines"
QUAL = Path("/private/tmp/oneD-v2-qualified-fem-analysis")
PF_REPO = "/Volumes/Data/Data/Nanopillar_calculation/PF-fracture-fatigue_v10_2_21_persistent_sites_top1"
FEM_REPO = "/Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM_MPZ_theta0_pf_parity_claude"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    pf = pd.read_csv(OUT / "oneD_v2_pf_mechanics_map.csv")
    native = pd.read_csv(OUT / "oneD_v2_fem_native_mechanics_map.csv")
    structural = pd.read_csv(OUT / "oneD_v2_fem_qualified_G_map.csv")
    validation = pd.read_csv(OUT / "oneD_v2_map_interpolation_validation.csv")
    peak = json.loads((QUAL / "pf_fem_arbitrary_path_qualification_summary.json").read_text())
    dbtt = json.loads((QUAL / "dbtt500_arbitrary_path_qualification_summary.json").read_text())
    straight_g = float(structural.loc[structural.actual_extension_um == 500, "qualified_G_J_per_m2"].iloc[0])
    straight_j = float(native.loc[native.actual_extension_um == 500, "native_J_aggregate_J_per_m2"].iloc[0])
    rows = []
    for material, anchor in (("Peak", peak), ("DBTT", dbtt)):
        opening = float(anchor["linearity_check"]["openings"][-1])
        factor = (12.5e-6 / opening) ** 2
        arbitrary_g_same_opening = float(anchor["G_reference_J_per_m2"] * factor)
        arbitrary_j_same_opening = float(anchor["J_finest_median_J_per_m2"] * factor)
        rows.append({
            "material_class": material, "comparison_extension_um": 500,
            "comparison_opening_m": 12.5e-6, "straight_path_G_J_per_m2": straight_g,
            "arbitrary_path_G_same_opening_J_per_m2": arbitrary_g_same_opening,
            "straight_minus_arbitrary_G_fraction": (straight_g-arbitrary_g_same_opening)/arbitrary_g_same_opening,
            "straight_path_native_J_J_per_m2": straight_j,
            "arbitrary_path_native_J_same_opening_J_per_m2": arbitrary_j_same_opening,
            "straight_minus_arbitrary_native_J_fraction": (straight_j-arbitrary_j_same_opening)/arbitrary_j_same_opening,
            "arbitrary_path_G_uncertainty_same_opening_J_per_m2": float(anchor["G_uncertainty_J_per_m2"] * factor),
            "interpretation": "PATH_REDUCTION_SENSITIVITY_NOT_PARAMETER_TUNING",
        })
    sensitivity = pd.DataFrame(rows)
    sensitivity.to_csv(OUT / "oneD_v2_path_reduction_sensitivity.csv", index=False)
    manifest = {
        "schema": "oneD_v2_mechanics_map_manifest_v2",
        "interpolation_method": "BOUNDED_PIECEWISE_LINEAR",
        "PCHIP_status": "DIAGNOSTIC_NOT_SELECTED",
        "extrapolation_policy": "FAIL_CLOSED_NO_EXTRAPOLATION",
        "campaign_outputs_populated": False,
        "maps": {
            "PF": {
                "provider": "PFMechanicsProvider", "source_repository": PF_REPO,
                "source_commit": "ab6279331d050919f78e2e7f278bc332466e34f8",
                "bulk_constraint": "PLANE_STRAIN", "crack_representation": "FINITE_WIDTH_STIFFNESS_KILLED_WAKE",
                "failed_region_state": "BINARY_STIFFNESS_KILL", "local_process_zone": "FINITE_SIZE_REDUCED_STATE",
                "native_driving_metric": "PF_NATIVE_DOMAIN_J", "continuum_energy_release_status": "UNQUALIFIED",
                "veto_policy": "FAIL_CLOSED_TERMINATE", "map_nodes_um": pf.actual_extension_um.tolist(),
                "raw_file": "oneD_v2_pf_mechanics_map.csv", "raw_file_sha256": sha(OUT / "oneD_v2_pf_mechanics_map.csv"),
                "valid_interval_um": [0, 1000], "candidate_independence": "NO_CANDIDATE_SEED_OR_TEMPERATURE_INPUT",
                "temperature_policy": "TEMPERATURE_INDEPENDENT_FIXED_ELASTICITY",
                "opening_scaling_max_relative_error": float(pf.filter(like="opening_scaling_").max().max()),
                "maximum_LOO_relative_error": float(validation[validation.provider == "PF"].relative_error.max()),
                "qualification_status": "PF_PRODUCTION_DISCRETE_MAP_QUALIFIED",
            },
            "FEMCZM_NATIVE": {
                "provider": "FEMCZMMechanicsProvider", "source_repository": FEM_REPO,
                "source_commit": "931bed66913afc970117bce900805ecc9b6225f8",
                "qualification_operator_commit": "15bfd51", "bulk_constraint": "PLANE_STRAIN",
                "crack_representation": "CODIMENSION_ONE_DISPLACEMENT_DISCONTINUITY",
                "failed_interface_state": "TRACTION_FREE", "local_process_zone": "FINITE_SIZE_REDUCED_STATE",
                "native_driving_metric": "FEM_NATIVE_DOMAIN_J", "map_nodes_um": native.actual_extension_um.tolist(),
                "raw_file": "oneD_v2_fem_native_mechanics_map.csv", "raw_file_sha256": sha(OUT / "oneD_v2_fem_native_mechanics_map.csv"),
                "valid_interval_um": [0, 1000], "candidate_independence": "NO_PROCESS_ZONE_CANDIDATE_SEED_OR_TEMPERATURE_INPUT",
                "temperature_policy": "TEMPERATURE_INDEPENDENT_FIXED_ELASTICITY",
                "maximum_LOO_relative_error": float(validation[(validation.provider == "FEMCZM") & validation.quantity.str.contains("native|F_over")].relative_error.max()),
                "qualification_status": "FEMCZM_PRODUCTION_NATIVE_MAP_QUALIFIED",
            },
            "FEMCZM_STRUCTURAL_G": {
                "provider": "FEMCZMMechanicsProvider", "source_repository": FEM_REPO,
                "source_commit": "931bed66913afc970117bce900805ecc9b6225f8",
                "qualification_operator_commits": ["15bfd51", "dc1146e", "44c2d95"],
                "bulk_constraint": "PLANE_STRAIN", "crack_representation": "CODIMENSION_ONE_DISPLACEMENT_DISCONTINUITY",
                "failed_interface_state": "TRACTION_FREE", "local_process_zone": "FINITE_SIZE_REDUCED_STATE",
                "qualified_structural_metric": "G_ENERGY_COMPLIANCE_VCCT", "map_nodes_um": structural.actual_extension_um.tolist(),
                "raw_file": "oneD_v2_fem_qualified_G_map.csv", "raw_file_sha256": sha(OUT / "oneD_v2_fem_qualified_G_map.csv"),
                "valid_interval_um": [0, 1000], "candidate_independence": "NO_PROCESS_ZONE_CANDIDATE_SEED_OR_TEMPERATURE_INPUT",
                "temperature_policy": "TEMPERATURE_INDEPENDENT_FIXED_ELASTICITY",
                "qualified_node_count": int(structural.node_qualification.eq("QUALIFIED").sum()),
                "maximum_LOO_relative_error": float(validation[(validation.provider == "FEMCZM") & validation.quantity.str.contains("G_over|KG_over")].relative_error.max()),
                "qualification_status": "FEMCZM_STRUCTURAL_G_MAP_QUALIFIED",
            },
        },
    }
    (OUT / "oneD_v2_mechanics_map_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    lines = ["# One-dimensional V2 path-reduction sensitivity", "",
             "The standardized straight map was not tuned to either archived path.", "",
             "| Class | Straight G bias | Straight native-J bias | Interpretation |", "|---|---:|---:|---|"]
    for row in sensitivity.itertuples():
        lines.append(f"| {row.material_class} | {100*row.straight_minus_arbitrary_G_fraction:.2f}% | {100*row.straight_minus_arbitrary_native_J_fraction:.2f}% | path-reduction sensitivity |")
    lines += ["", "Peak and DBTT arbitrary-path anchors are separately qualified at their archived openings and rescaled only by verified linear-elastic opening² scaling. These differences are reported, not absorbed into material parameters.", ""]
    (OUT / "ONE_D_V2_PATH_REDUCTION_SENSITIVITY.md").write_text("\n".join(lines))
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.4), constrained_layout=True)
    axes[0].bar(sensitivity.material_class, 100*sensitivity.straight_minus_arbitrary_G_fraction)
    axes[0].set(ylabel="Straight minus kinked G (%)", title="Qualified structural G")
    axes[1].bar(sensitivity.material_class, 100*sensitivity.straight_minus_arbitrary_native_J_fraction)
    axes[1].set(ylabel="Straight minus kinked native J (%)", title="Native domain J")
    for ax in axes: ax.axhline(0, color="k", lw=.8); ax.grid(axis="y", alpha=.25)
    fig.savefig(OUT / "STRAIGHT_VS_KINKED_PATH_SENSITIVITY.png", dpi=180); plt.close(fig)


if __name__ == "__main__":
    main()
