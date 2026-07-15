from __future__ import annotations

from argparse import Namespace
from pathlib import Path

from run_mpz_v9_12_tip_only_material_rcurve import (
    build_command,
    effective_seed,
)


def _args():
    return Namespace(
        parameter_root=Path("mpz_v9_11_parameters"),
        T_K=700.0,
        target_extension_um=100.0,
        steps=6000,
        nx=36,
        ny=72,
        tip_h_fine=1e-6,
        tip_ratio=1.2,
        dU=2e-7,
        dt=8.4,
        n_stagger=2,
        print_every=25,
        adaptive_event_target=0.15,
        da_phys_um=5.0,
        mpz_length_um=100.0,
        mpz_n_bins=200,
        crystal_theta_deg=45.0,
        save_snapshots=5,
        snapshot_cols=5,
        snapshot_by_extension_um=25.0,
    )


def test_independent_rng_streams_differ_by_material():
    vals = {
        effective_seed(1, "ceramic", "independent"),
        effective_seed(1, "weakT", "independent"),
        effective_seed(1, "DBTT", "independent"),
    }
    assert len(vals) == 3
    assert effective_seed(7, "ceramic", "common") == 7
    assert effective_seed(7, "DBTT", "common") == 7


def test_command_requires_tip_only_and_full_field_rendering(tmp_path):
    cmd = build_command(_args(), "weakT", tmp_path, force_rerun=True)
    joined = " ".join(cmd)
    assert "--bulk-plasticity-mode tip_only" in joined
    assert "--make-solver-plots" in cmd
    assert "--no-skip-existing" in cmd
    assert "--save-snapshots 5" in joined
    assert "--snapshot-by-extension-um 25.0" in joined
