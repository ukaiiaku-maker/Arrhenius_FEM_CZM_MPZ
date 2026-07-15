from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import run_mpz_v9_14_event_remesh_gate as runner


def _args():
    return SimpleNamespace(
        parameter_root=Path("mpz_v9_11_parameters"),
        T_K=700.0,
        target_extension_um=50.0,
        steps=12000,
        nx=36,
        ny=72,
        tip_h_fine=1e-6,
        tip_ratio=1.2,
        dU=2e-7,
        dt=8.4,
        n_stagger=2,
        print_every=25,
        adaptive_event_target=0.05,
        da_phys_um=5.0,
        mpz_length_um=100.0,
        mpz_n_bins=200,
        crystal_theta_deg=45.0,
        save_snapshots=5,
        snapshot_cols=5,
        snapshot_by_extension_um=10.0,
    )


def test_build_command_selects_v914_tip_only_full_field_runner(tmp_path):
    cmd = runner.build_command(_args(), "DBTT", tmp_path, force_rerun=True)
    assert cmd[1] == "run_mpz_v9_14_mode_i_rcurve.py"
    assert cmd[cmd.index("--bulk-plasticity-mode") + 1] == "tip_only"
    assert cmd[cmd.index("--adaptive-event-target") + 1] == "0.05"
    assert "--make-solver-plots" in cmd
    assert "--no-skip-existing" in cmd


def test_v914_defaults_use_short_deterministic_remesh_gate(monkeypatch, tmp_path):
    # Parser defaults are exercised indirectly by capturing the first run_case call.
    captured = {}

    def fake_run_case(args, seed, cls, root):
        captured.update({
            "target": args.target_extension_um,
            "statistics": args.event_statistics,
            "stochastic_emission": args.stochastic_emission,
            "target_h": args.event_remesh_target_h_m,
            "patch_radius_um": args.event_remesh_patch_radius_um,
        })
        return {"class": cls}

    monkeypatch.setattr(runner, "run_case", fake_run_case)
    monkeypatch.setattr(runner, "audit_campaign", lambda *a, **k: {
        "numerical_event_remesh_gate_passed": False,
        "material_transfer_gate_passed_v914": False,
        "interpretation": "test",
        "failed_numerical_remesh_cases": [],
    })
    monkeypatch.setattr(
        "sys.argv",
        ["run_mpz_v9_14_event_remesh_gate.py", "--outroot", str(tmp_path), "--classes", "ceramic"],
    )
    runner.main()
    assert captured["target"] == 50.0
    assert captured["statistics"] == "deterministic"
    assert captured["stochastic_emission"] is False
    assert captured["target_h"] == 1.0e-6
    assert captured["patch_radius_um"] == 25.0
