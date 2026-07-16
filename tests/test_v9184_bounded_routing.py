from __future__ import annotations

from types import SimpleNamespace

import run_mpz_v9_18_4_mode_i_rcurve as runner
from arrhenius_fracture import mode_i_first_passage_v9_18_4_bounded as bounded


def test_mode_i_runner_routes_to_bounded_runtime(monkeypatch):
    def fake_build(py, args, class_name, manifest, T_K, case_dir):
        return [py, "-m", "arrhenius_fracture.mode_i_first_passage_v9_11", "--out", str(case_dir)]

    monkeypatch.setattr(runner, "_original_build_command", fake_build)
    cmd = runner._build_command_v9184(
        "python", SimpleNamespace(), "ceramic", "manifest.csv", 700, "out"
    )
    assert "arrhenius_fracture.mode_i_first_passage_v9_18_4_bounded" in cmd
    assert "arrhenius_fracture.mode_i_first_passage_v9_11" not in cmd


def test_bounded_runtime_exports_transactional_guards():
    assert callable(bounded._bounded_insert)
    assert callable(bounded._bounded_advance)
    assert callable(bounded.main)
