from __future__ import annotations

from types import SimpleNamespace

from arrhenius_fracture import mode_i_first_passage_v9_17 as v917
from arrhenius_fracture import mode_i_first_passage_v9_17_1 as v9171


def test_internal_one_fire_initializer_sets_front_configuration():
    calls = []

    def original_init(obj, marker):
        calls.append(marker)
        obj.f = SimpleNamespace(max_advances_per_step=99.0)

    wrapped = v9171._one_fire_engine_initializer(original_init)
    obj = SimpleNamespace()
    wrapped(obj, "called")

    assert calls == ["called"]
    assert obj.f.max_advances_per_step == 1.0
    assert obj.v9171_one_fire_internal_routing is True


def test_main_suppresses_only_unsupported_cli_injection(monkeypatch):
    captured = {}

    def fake_v916_main(argv):
        captured["argv"] = list(argv)
        return ["ok"]

    monkeypatch.setattr(v917._v916, "main", fake_v916_main)
    result = v9171.main(["--out", "unused-test-output", "--crystal-aniso"])

    assert result == ["ok"]
    assert "--max-advances-per-step" not in captured["argv"]
    assert "--crystal-aniso" in captured["argv"]
