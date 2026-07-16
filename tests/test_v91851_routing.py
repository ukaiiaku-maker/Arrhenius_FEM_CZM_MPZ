from __future__ import annotations

from types import SimpleNamespace

import run_mpz_v9_18_5_1_mode_i_rcurve as mode_runner
import run_mpz_v9_18_5_1_persistent_plastic_wake as campaign_runner


def test_mode_i_command_routes_to_v91851(monkeypatch, tmp_path):
    def fake_build(py, args, class_name, manifest, T_K, case_dir):
        return [py, "-m", "arrhenius_fracture.mode_i_first_passage_v9_11", "--out", str(case_dir)]

    monkeypatch.setattr(mode_runner, "_original_build_command", fake_build)
    cmd = mode_runner._build_command_v91851(
        "python", SimpleNamespace(), "ceramic", {}, 700.0, tmp_path
    )
    assert "arrhenius_fracture.mode_i_first_passage_v9_18_5_1" in cmd
    assert "arrhenius_fracture.mode_i_first_passage_v9_11" not in cmd


def test_campaign_command_routes_to_v91851(monkeypatch, tmp_path):
    def fake_original(args, class_name, run_root, force_rerun):
        return ["python", "run_mpz_v9_18_mode_i_rcurve.py", "--outroot", str(run_root)]

    monkeypatch.setattr(campaign_runner._v9181, "_original_build", fake_original)
    cmd = campaign_runner._build_command_v91851(
        SimpleNamespace(), "ceramic", tmp_path, False
    )
    assert "run_mpz_v9_18_5_1_mode_i_rcurve.py" in cmd
    assert "run_mpz_v9_18_mode_i_rcurve.py" not in cmd
