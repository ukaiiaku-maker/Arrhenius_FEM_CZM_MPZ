from pathlib import Path

import pandas as pd


OUT = Path("analysis_outputs/oneD_v2_mechanics_maps_and_baselines")


def frame():
    return pd.read_csv(OUT / "oneD_v2_pf_mechanics_map.csv")


def test_pf_map_uses_sequential_production_sharp_wake():
    data = frame()
    assert data.geometry_event_count.is_monotonic_increasing
    # The finite-width backend changes only when a new mesh node enters the
    # killed band, so several physical segments may intentionally share a
    # damage fingerprint.  Those plateaus are a discrete feature, not a replay
    # failure.
    assert 1 < data.damage_wake_fingerprint.nunique() <= data.actual_extension_um.nunique()
    assert (data.wake_width_m == 2 * data.mesh_scale_m).all()


def test_pf_map_is_candidate_independent():
    forbidden = {"candidate", "material_class", "seed", "temperature_K"}
    assert not forbidden.intersection(frame().columns)


def test_pf_map_never_exposes_continuum_G():
    assert not any("continuum_G" in column or "qualified_G" in column for column in frame().columns)
    text = (OUT / "ONE_D_V2_PF_MECHANICS_MAP_QUALIFICATION_V3.md").read_text()
    assert "not continuum G" in text


def test_pf_opening_scaling_is_exact():
    data = frame()
    columns = [column for column in data if column.startswith("opening_scaling_")]
    assert data[columns].to_numpy().max() < 1e-9


def test_pf_interpolation_is_bounded_and_extrapolation_closed():
    validation = pd.read_csv(OUT / "oneD_v2_map_interpolation_validation.csv")
    assert validation.left_out_extension_um.between(0, 1000).all()
    text = (OUT / "ONE_D_V2_PF_MECHANICS_MAP_QUALIFICATION_V3.md").read_text()
    assert "extrapolation fails closed" in text


def test_nominal_two_micron_node_does_not_truncate_segment():
    row = frame().query("target_extension_um == 2").iloc[0]
    assert row.actual_extension_um == 5
    assert row.geometry_event_count == 1
