from __future__ import annotations

import run_mpz_v9_18_5_8_mode_i_rcurve as mode
import run_mpz_v9_18_5_8_persistent_plastic_wake as wake


def test_mode_i_route_replaces_v911_module(monkeypatch):
    monkeypatch.setattr(
        mode,
        "_original_build_command",
        lambda py, args, class_name, manifest, T_K, case_dir: [
            py,
            "-m",
            "arrhenius_fracture.mode_i_first_passage_v9_11",
        ],
    )
    cmd = mode._build_command_v91858("python", None, "weakT", "m.csv", 700, "out")
    assert "arrhenius_fracture.mode_i_first_passage_v9_18_5_6" in cmd
    assert "arrhenius_fracture.mode_i_first_passage_v9_11" not in cmd


def test_persistent_wake_route_replaces_driver(monkeypatch):
    monkeypatch.setattr(
        wake._v9181,
        "_original_build",
        lambda args, class_name, run_root, force_rerun: [
            "python",
            "run_mpz_v9_18_mode_i_rcurve.py",
        ],
    )
    cmd = wake._build_command_v91858(None, "DBTT", "out", False)
    assert "run_mpz_v9_18_5_8_mode_i_rcurve.py" in cmd
    assert "run_mpz_v9_18_mode_i_rcurve.py" not in cmd
