from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from arrhenius_fracture import pf_theta0_driving_force_trajectory_v10051840 as pf_trajectory

REAL_PF_STEPS_CSV = Path(
    "/Volumes/Data/Data/Nanopillar_calculation/PF-fracture-fatigue_v10_2_21_persistent_sites_top1"
    "/runs/v10_4_1_theta0_rate1x_bulk_PT_four_class_1000um_selective_reuse_base3621_v1"
    "/v913_paper_peak01_0242980_persistent_sites/T1000K_th0_seed8666/steps_1000K.csv"
)


def _write_steps(
    path: Path,
    *,
    n: int,
    dt_s: float = 8.4,
    first_passage_row: int | None = None,
) -> None:
    step = np.arange(1, n + 1, dtype=float)
    KJ = np.linspace(2.0e6, 9.0e6, n)
    J = np.linspace(50.0, 6700.0, n)
    crack_extension_m = np.zeros(n)
    B = np.zeros(n)
    N_em = np.zeros(n)
    if first_passage_row is not None:
        crack_extension_m[first_passage_row:] = 5.0e-6
        B[first_passage_row:] = 0.4
        N_em[first_passage_row:] = 500.0
    frame = pd.DataFrame(
        {
            "step": step,
            "dt_cur_s": np.full(n, dt_s),
            "J_effective_direct_J_per_m2": J,
            "J_signed_direct_J_per_m2": J,
            "KJ_Pa_sqrtm": KJ,
            "B": B,
            "N_em": N_em,
            "crack_extension_m": crack_extension_m,
            "W_bulk_plastic_cumulative_J_per_m": np.linspace(0.0, 1.0e-4, n),
        }
    )
    frame.to_csv(path, index=False)


def test_load_reconstructs_physical_time_from_dt_and_records_provenance(tmp_path: Path):
    csv_path = tmp_path / "steps_1000K.csv"
    _write_steps(csv_path, n=10, dt_s=8.4)

    trajectory = pf_trajectory.load_pf_theta0_driving_force_trajectory(csv_path)

    assert trajectory.n_rows == 10
    np.testing.assert_allclose(
        trajectory.physical_time_s, 8.4 * np.arange(1, 11, dtype=float)
    )
    assert trajectory.source_path == str(csv_path.resolve())
    assert trajectory.source_sha256 == hashlib.sha256(csv_path.read_bytes()).hexdigest()
    assert trajectory.first_passage_row_index is None


def test_first_passage_row_detected_and_excluded_from_prefracture_slice(tmp_path: Path):
    csv_path = tmp_path / "steps_1000K.csv"
    _write_steps(csv_path, n=20, first_passage_row=12)

    trajectory = pf_trajectory.load_pf_theta0_driving_force_trajectory(csv_path)
    assert trajectory.first_passage_row_index == 12
    assert trajectory.first_passage_step == 13.0
    assert trajectory.first_passage_KJ_Pa_sqrtm == pytest.approx(trajectory.KJ_target_Pa_sqrtm[12])

    pre = pf_trajectory.prefracture_slice(trajectory)
    assert pre.n_rows == 12
    assert np.all(pre.crack_extension_m == 0.0)
    # the crossing row itself must never be folded into the prefracture slice
    assert pre.physical_time_s[-1] < trajectory.first_passage_time_s


def test_prefracture_slice_is_full_trajectory_when_no_growth_recorded(tmp_path: Path):
    csv_path = tmp_path / "steps_1000K.csv"
    _write_steps(csv_path, n=5)

    trajectory = pf_trajectory.load_pf_theta0_driving_force_trajectory(csv_path)
    pre = pf_trajectory.prefracture_slice(trajectory)
    assert pre.n_rows == trajectory.n_rows == 5


def test_interpolate_target_is_bounded_and_rejects_extrapolation(tmp_path: Path):
    csv_path = tmp_path / "steps_1000K.csv"
    _write_steps(csv_path, n=10, dt_s=8.4)
    trajectory = pf_trajectory.load_pf_theta0_driving_force_trajectory(csv_path)

    mid_time = 0.5 * (trajectory.physical_time_s[0] + trajectory.physical_time_s[-1])
    kj_mid = pf_trajectory.interpolate_target(trajectory, mid_time, channel="KJ")
    assert trajectory.KJ_target_Pa_sqrtm[0] < kj_mid < trajectory.KJ_target_Pa_sqrtm[-1]

    with pytest.raises(ValueError, match="extrapolation is not permitted"):
        pf_trajectory.interpolate_target(trajectory, trajectory.physical_time_s[-1] + 1.0, channel="KJ")
    with pytest.raises(ValueError, match="extrapolation is not permitted"):
        pf_trajectory.interpolate_target(trajectory, trajectory.physical_time_s[0] - 1.0, channel="J")


def test_missing_required_column_raises(tmp_path: Path):
    csv_path = tmp_path / "steps_1000K.csv"
    frame = pd.DataFrame({"step": [1.0, 2.0], "dt_cur_s": [8.4, 8.4]})
    frame.to_csv(csv_path, index=False)

    with pytest.raises(ValueError, match="lacks required columns"):
        pf_trajectory.load_pf_theta0_driving_force_trajectory(csv_path)


def test_nonpositive_dt_raises(tmp_path: Path):
    csv_path = tmp_path / "steps_1000K.csv"
    _write_steps(csv_path, n=5, dt_s=8.4)
    frame = pd.read_csv(csv_path)
    frame.loc[2, "dt_cur_s"] = 0.0
    frame.to_csv(csv_path, index=False)

    with pytest.raises(ValueError, match="non-positive dt_cur_s"):
        pf_trajectory.load_pf_theta0_driving_force_trajectory(csv_path)


@pytest.mark.skipif(
    not REAL_PF_STEPS_CSV.is_file(),
    reason="real PF reference case not present on this machine",
)
def test_real_pf_reference_prefracture_first_passage_matches_handoff():
    trajectory = pf_trajectory.load_pf_theta0_driving_force_trajectory(REAL_PF_STEPS_CSV)

    assert trajectory.source_sha256 == (
        "666d839ae8ef4267ae4ac7ab52a4220de676525bfb70917744e25a119c291b0c"
    )
    assert trajectory.first_passage_row_index == 176
    assert trajectory.first_passage_KJ_Pa_sqrtm / 1.0e6 == pytest.approx(53.9019, rel=1.0e-4)
    assert trajectory.first_passage_J_J_per_m2 == pytest.approx(6531.0, rel=0.05)

    pre = pf_trajectory.prefracture_slice(trajectory)
    assert pre.n_rows == 176
    assert np.all(pre.crack_extension_m == 0.0)
