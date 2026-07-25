from pathlib import Path


RUNNER = Path("run_v10_0_5_19_frozen_four_class_temperature_sweep.sh")


def test_temperature_sweep_uses_v100519_standalone_path():
    text = RUNNER.read_text()
    assert "mode_i_first_passage_v10_0_5_19_four_class_standalone" in text
    assert "PFROOT" not in text
    assert "phase_field_solver_dependency=false" in text
    assert "phase_field_result_comparison_enabled=false" in text


def test_temperature_sweep_has_exact_four_parameterizations():
    text = RUNNER.read_text()
    for option in (
        "v913_paper_peak01_0242980_persistent_sites",
        "v913_paper_dbtt01_0202500_persistent_sites",
        "v913_paper_weakT01_0129902_persistent_sites",
        "v913_paper_ceramic01_0077080_persistent_sites",
    ):
        assert text.count(option) == 1
    assert "expected 56 cases" in text


def test_temperature_sweep_production_defaults():
    text = RUNNER.read_text()
    assert 'TEMPERATURES=${TEMPERATURES:-"300 400 500 600 700 800 850 900 950 1000 1050 1100 1150 1200"}' in text
    assert "TARGET_EXT_UM=${TARGET_EXT_UM:-1000}" in text
    assert "STEPS=${STEPS:-1000000}" in text
    assert "MAX_JOBS=${MAX_JOBS:-2}" in text
    assert "CRYSTAL_THETA_DEG=${CRYSTAL_THETA_DEG:-30}" in text
    assert "MIN_GLOBAL_FORWARD=${MIN_GLOBAL_FORWARD:-0.05}" in text
    assert "--max-fronts 1" in text
    assert "--crystal-aniso" in text
    assert "--crystal-compete" in text


def test_temperature_sweep_keeps_quality_floors_and_restart_policy():
    text = RUNNER.read_text()
    assert "triangle_quality_floor_relaxed=false" in text
    assert "child_area_ratio_floor_relaxed=false" in text
    assert "ARRHENIUS_QUALITY_PARTITION_MAX_DEPTH" in text
    assert "ARRHENIUS_QUALITY_PARTITION_MIN_SEGMENT_M" in text
    assert "SKIP_FINISHED=${SKIP_FINISHED:-1}" in text
    assert "RESTART_INCOMPLETE=${RESTART_INCOMPLETE:-1}" in text
