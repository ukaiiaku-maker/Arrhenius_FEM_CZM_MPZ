from pathlib import Path


def test_v1002_smoke_runner_uses_supported_straight_mode_i_cli():
    root = Path(__file__).resolve().parents[1]
    text = (root / "run_v10_0_2_progressive_one_segment_smoke.sh").read_text()

    assert "--da-phys 5e-6" in text
    assert "--da-phys-um" not in text
    assert "--no-crystal-branch" not in text
    assert "--crystal-aniso" not in text
    assert "--crystal-compete" not in text
    assert "--max-fronts 1" in text
    assert "mode_i_first_passage_v10_0_2_progressive" in text
