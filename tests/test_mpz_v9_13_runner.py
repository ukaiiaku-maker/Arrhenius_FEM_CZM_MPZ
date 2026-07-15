from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pandas as pd

import run_mpz_v9_13_deterministic_material_transfer as runner


def _args():
    return SimpleNamespace(
        parameter_root=Path("mpz_v9_11_parameters"), T_K=700.0,
        target_extension_um=100.0, steps=6000, nx=36, ny=72,
        tip_h_fine=1e-6, tip_ratio=1.2, dU=2e-7, dt=8.4,
        n_stagger=2, print_every=25, adaptive_event_target=0.15,
        da_phys_um=5.0, mpz_length_um=100.0, mpz_n_bins=200,
        crystal_theta_deg=45.0, save_snapshots=5, snapshot_cols=5,
        snapshot_by_extension_um=25.0,
    )


def test_build_command_uses_v913_direct_driver_and_full_fields(tmp_path):
    cmd = runner.build_command(_args(), "weakT", tmp_path, force_rerun=True)
    assert cmd[1] == "run_mpz_v9_13_mode_i_rcurve.py"
    assert "--make-solver-plots" in cmd
    assert "--no-skip-existing" in cmd
    assert cmd[cmd.index("--bulk-plasticity-mode") + 1] == "tip_only"


def test_temperature_summary_is_read_from_run_root_and_copied(tmp_path):
    run_root = tmp_path / "run"
    case_dir = run_root / "weakT" / "T700_th45"
    case_dir.mkdir(parents=True)
    pd.DataFrame([{
        "class": "weakT", "T_K": 700, "status": "complete",
        "final_extension_um": 101.0, "K_init_MPa_sqrt_m": 17.0,
    }]).to_csv(run_root / "rcurve_temperature_summary.csv", index=False)
    row, copied = runner._copy_temperature_summary(run_root, case_dir, "weakT", 700.0)
    assert row["status"] == "complete"
    assert row["final_extension_um"] == 101.0
    assert copied is not None
    assert (case_dir / "rcurve_temperature_summary_v913.csv").exists()


def test_effective_seed_common_is_paired_and_independent_is_offset():
    assert runner.effective_seed(1, "ceramic", "common") == 1
    assert runner.effective_seed(1, "DBTT", "common") == 1
    assert runner.effective_seed(1, "ceramic", "independent") != runner.effective_seed(1, "DBTT", "independent")
