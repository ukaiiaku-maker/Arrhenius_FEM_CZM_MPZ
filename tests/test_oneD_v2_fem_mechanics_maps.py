from pathlib import Path

import pandas as pd


OUT = Path("analysis_outputs/oneD_v2_mechanics_maps_and_baselines")


def native():
    return pd.read_csv(OUT / "oneD_v2_fem_native_mechanics_map.csv")


def structural():
    return pd.read_csv(OUT / "oneD_v2_fem_qualified_G_map.csv")


def test_native_and_structural_maps_are_distinct():
    assert "J_native_over_U2" in native()
    assert "G_over_U2" in structural()
    assert "G_over_U2" not in native()


def test_structural_G_requires_three_method_agreement():
    data = structural()
    qualified = data[data.node_qualification == "QUALIFIED"]
    values = qualified[["G_energy_J_per_m2", "G_compliance_J_per_m2", "G_VCCT_J_per_m2"]]
    relative_spread = (values.max(axis=1) - values.min(axis=1)) / values.max(axis=1)
    assert (relative_spread <= .10).all()
    assert qualified.topology_invariants_pass.all()


def test_all_required_structural_nodes_qualified():
    data = structural()
    assert len(data) == 18
    assert data.node_qualification.eq("QUALIFIED").all()
    assert data.qualified_G_J_per_m2.notna().all()


def test_plane_strain_does_not_encode_crack_thickness():
    text = (OUT / "ONE_D_V2_FEMCZM_NATIVE_MAP_QUALIFICATION.md").read_text()
    assert "plane-strain" in text
    assert "codimension-one displacement-discontinuity" in text


def test_global_reactions_close():
    for data in (native(), structural()):
        assert data.reaction_residual_relative_error.max() < 1e-9
        assert data.reaction_energy_relative_error.max() < 1e-9


def test_selected_opening_scaling_checks_pass():
    data = structural().dropna(subset=["opening_scaling_G_ratio_error"])
    assert set(data.actual_extension_um) == {0, 500, 1000}
    columns = [column for column in data if column.startswith("opening_scaling_")]
    assert data[columns].to_numpy().max() < 1e-9


def test_native_domain_spread_does_not_erase_qualified_G():
    assert (native().J_domain_spread_J_per_m2 >= 0).all()
    assert structural().node_qualification.eq("QUALIFIED").all()


def test_manifest_keeps_metadata_axes_separate_and_bounds_closed():
    import json
    manifest = json.loads((OUT / "oneD_v2_mechanics_map_manifest.json").read_text())
    assert manifest["extrapolation_policy"] == "FAIL_CLOSED_NO_EXTRAPOLATION"
    for item in manifest["maps"].values():
        assert item["bulk_constraint"] == "PLANE_STRAIN"
        assert item["crack_representation"] != item["bulk_constraint"]
        assert item["local_process_zone"] == "FINITE_SIZE_REDUCED_STATE"


def test_path_reduction_error_is_reported_not_tuned():
    sensitivity = pd.read_csv(OUT / "oneD_v2_path_reduction_sensitivity.csv")
    assert set(sensitivity.material_class) == {"Peak", "DBTT"}
    assert sensitivity.interpretation.eq("PATH_REDUCTION_SENSITIVITY_NOT_PARAMETER_TUNING").all()
    assert sensitivity.straight_minus_arbitrary_G_fraction.abs().max() < .03
