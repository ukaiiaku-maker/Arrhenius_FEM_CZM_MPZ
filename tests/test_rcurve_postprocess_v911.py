from __future__ import annotations

import numpy as np
import pandas as pd

from arrhenius_fracture.rcurve_postprocess_v911 import (
    cascade_metrics,
    cluster_same_load_events,
    write_cascade_aware_outputs,
)
from run_mpz_v9_11_stochastic_continuation_700K import continuation_metrics


def test_same_load_topology_events_are_one_unstable_jump():
    raw = pd.DataFrame({
        "raw_event_id": [1, 2, 3, 4],
        "step": [10, 11, 12, 20],
        "Uapp_m": [1.0e-4, 1.00001e-4, 1.00002e-4, 1.2e-4],
        "KJ_MPa_sqrt_m": [20.0, 21.0, 22.0, 30.0],
        "crack_extension_before_um": [0.0, 5.0, 10.0, 15.0],
        "crack_extension_after_um": [5.0, 10.0, 15.0, 20.0],
        "da_block_um": [5.0, 5.0, 5.0, 5.0],
        "n_fire": [1, 1, 1, 1],
    })
    out = cluster_same_load_events(raw, relative_load_tolerance=1.0e-4)
    assert len(out) == 2
    first = out.iloc[0]
    assert first["classification"] == "unstable_same_load_cascade"
    assert first["topology_event_count"] == 3
    assert np.isclose(first["jump_span_um"], 15.0)
    assert np.isclose(first["KJ_onset_MPa_sqrt_m"], 20.0)
    metrics = cascade_metrics(raw, out)
    assert metrics["n_raw_topology_events"] == 4
    assert metrics["n_independent_load_events"] == 2
    assert metrics["n_unstable_same_load_cascades"] == 1
    assert np.isclose(metrics["fraction_topology_events_in_cascades"], 0.75)


def test_postprocessor_preserves_raw_and_writes_clustered_compatibility_file(tmp_path):
    rows = pd.DataFrame({
        "step": [1, 2, 3, 4],
        "Uapp_m": [1.0e-4, 1.0e-4, 1.0e-4, 1.2e-4],
        "KJ_Pa_sqrtm": [20e6, 21e6, 22e6, 30e6],
        "a_tip_m": [0.505e-3, 0.510e-3, 0.515e-3, 0.520e-3],
        "crack_extension_m": [5e-6, 10e-6, 15e-6, 20e-6],
        "da_block_m": [5e-6, 5e-6, 5e-6, 5e-6],
        "n_fire": [1, 1, 1, 1],
        "B": [0.1, 0.2, 0.3, 0.4],
        "sigma_tip_Pa": [5e9, 6e9, 7e9, 8e9],
        "mpz_K_shield_Pa_sqrt_m": [0.1e6, 0.2e6, 0.3e6, 0.4e6],
    })
    rows.to_csv(tmp_path / "steps_0700K.csv", index=False)
    metrics = write_cascade_aware_outputs(tmp_path, 700.0)
    raw = pd.read_csv(tmp_path / "R_curve_topology_events_raw.csv")
    clustered = pd.read_csv(tmp_path / "R_curve_load_events_clustered.csv")
    compat = pd.read_csv(tmp_path / "R_curve_event_sampled.csv")
    assert len(raw) == 4
    assert np.allclose(raw["sigma_tip_GPa"], [5.0, 6.0, 7.0, 8.0])
    assert np.allclose(raw["K_shield_MPa_sqrt_m"], [0.1, 0.2, 0.3, 0.4])
    assert len(clustered) == 2
    pd.testing.assert_frame_equal(clustered, compat)
    assert metrics["largest_same_load_jump_um"] == 15.0


def test_no_growth_case_writes_parseable_schema_and_nan_metrics(tmp_path):
    rows = pd.DataFrame({
        "step": [1, 2, 3],
        "Uapp_m": [1.0e-5, 2.0e-5, 3.0e-5],
        "KJ_Pa_sqrtm": [1.0e6, 1.5e6, 2.0e6],
        "a_tip_m": [0.5e-3, 0.5e-3, 0.5e-3],
        "crack_extension_m": [0.0, 0.0, 0.0],
        "da_block_m": [0.0, 0.0, 0.0],
        "n_fire": [0, 0, 0],
    })
    rows.to_csv(tmp_path / "steps_0700K.csv", index=False)
    metrics = write_cascade_aware_outputs(tmp_path, 700.0)
    raw = pd.read_csv(tmp_path / "R_curve_topology_events_raw.csv")
    clustered = pd.read_csv(tmp_path / "R_curve_load_events_clustered.csv")
    assert raw.empty
    assert clustered.empty
    assert "KJ_onset_MPa_sqrt_m" in clustered.columns
    assert metrics["rcurve_interpretation"] == "no_growth_events"
    load = continuation_metrics(tmp_path, None)
    assert np.isnan(load["Kload_200_500um_median"])
    assert np.isnan(load["delta_Kload_median_minus_init"])


def test_zero_byte_cluster_file_is_treated_as_no_growth(tmp_path):
    (tmp_path / "R_curve_load_events_clustered.csv").write_bytes(b"")
    load = continuation_metrics(tmp_path, 12.0)
    assert np.isnan(load["Kload_200_500um_mean"])
    assert np.isnan(load["delta_Kload_median_minus_init"])
