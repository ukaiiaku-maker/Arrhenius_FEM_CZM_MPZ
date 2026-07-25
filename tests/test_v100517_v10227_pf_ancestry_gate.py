from pathlib import Path


def test_runner_accepts_descendants_of_canonical_parameter_commit():
    text = Path("run_v10_0_5_17_v10227_four_class_campaign.sh").read_text()
    assert "PF_CANONICAL_COMMIT_REQUIRED" in text
    assert "merge-base --is-ancestor" in text
    assert "PF_COMMIT_PREFIX_REQUIRED" not in text
    assert "v10.2.27-paper-four-class-30deg-long-rcurves" in text


def test_smoke_exports_canonical_ancestor_not_moving_head():
    text = Path("run_v10_0_5_17_v10227_four_class_300_1000K_20um_smoke.sh").read_text()
    assert "PF_CANONICAL_COMMIT_REQUIRED" in text
    assert "73a97ff5fc15c3be3f513e6ce3f219dcf028580c" in text
    assert "PF_COMMIT_PREFIX_REQUIRED" not in text
