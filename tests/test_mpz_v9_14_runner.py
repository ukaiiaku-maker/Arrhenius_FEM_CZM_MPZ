from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np

import arrhenius_fracture.mode_i_first_passage_v9_14 as mode_v914
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
        adaptive_event_target=0.01,
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
    assert cmd[cmd.index("--adaptive-event-target") + 1] == "0.01"
    assert "--make-solver-plots" in cmd
    assert "--no-skip-existing" in cmd


def test_absolute_action_predictor_does_not_divide_by_remaining_threshold(monkeypatch):
    raw_dB = 0.017

    def raw_predict(self, K_cleave, K_emit, T, dt):
        return raw_dB

    monkeypatch.setattr(
        mode_v914._engine_v911._BaseEngine,
        "predict_clock_increment_drives",
        raw_predict,
    )

    class Dummy:
        B = 0.999999
        B_target = 1.0

        @staticmethod
        def _reload_gate_active(_K):
            return False

    dummy = Dummy()
    predicted = mode_v914._absolute_action_predictor(
        dummy, 20.0e6, 20.0e6, 700.0, 1.0
    )
    assert np.isclose(predicted, raw_dB)
    assert dummy.adaptive_prediction_coordinate_v914 == (
        "absolute_integrated_hazard_action"
    )
    assert np.isclose(dummy.adaptive_predicted_absolute_dB_v914, raw_dB)


def test_v914_defaults_use_short_deterministic_remesh_gate(monkeypatch, tmp_path):
    # Parser defaults are exercised indirectly by capturing the first run_case call.
    captured = {}

    def fake_run_case(args, seed, cls, root):
        captured.update({
            "target": args.target_extension_um,
            "statistics": args.event_statistics,
            "stochastic_emission": args.stochastic_emission,
            "adaptive_event_target": args.adaptive_event_target,
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
    assert captured["adaptive_event_target"] == 0.01
    assert captured["target_h"] == 1.0e-6
    assert captured["patch_radius_um"] == 25.0
