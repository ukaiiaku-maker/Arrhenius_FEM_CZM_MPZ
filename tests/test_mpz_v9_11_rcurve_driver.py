from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

from run_mpz_v9_11_mode_i_rcurve_3T import build_command


def test_long_growth_command_uses_target_extension_without_first_fire_stop(tmp_path):
    args = SimpleNamespace(
        mpz_length_um=100.0,
        mpz_n_bins=200,
        mpz_profile_sector_half_angle_deg=45.0,
        mpz_profile_damage_cutoff=0.85,
        nx=36,
        ny=72,
        tip_h_fine=1.0e-6,
        tip_ratio=1.2,
        dU=2.0e-7,
        dt=8.4,
        steps=6000,
        n_stagger=2,
        print_every=25,
        target_extension_um=500.0,
        adaptive_event_target=0.15,
        da_phys_um=5.0,
        rJ_cluster_um=20.0,
        rJ_outer_um=25.0,
        crystal_theta_deg=45.0,
        save_snapshots=12,
        snapshot_cols=4,
        snapshot_by_extension_um=50.0,
        make_solver_plots=False,
    )
    cmd = build_command(
        "python",
        args,
        "DBTT",
        Path("mpz_v9_11_parameters/DBTT/spatial_promotion_manifest.csv"),
        700,
        tmp_path,
    )
    assert "arrhenius_fracture.mode_i_first_passage_v9_11" in cmd
    assert "--target-crack-extension-um" in cmd
    i = cmd.index("--target-crack-extension-um")
    assert float(cmd[i + 1]) == 500.0
    assert "--stop-after-first-fire" not in cmd
    assert cmd[cmd.index("--max-fronts") + 1] == "1"
    assert cmd[cmd.index("--temperatures") + 1] == "700"
    assert "--no-plots" in cmd
