from scripts.validate_oneD_v2_peak_dbtt_R import (
    CONTROL_IDS, SEEDS, SURVIVORS, TEMPERATURES, planned_cases,
)


def test_survivor_set_retains_controls_and_diverse_variants():
    for material in ("Peak", "DBTT"):
        assert SURVIVORS[material][0] == CONTROL_IDS[material]
        assert len(SURVIVORS[material]) >= 4


def test_dense_validation_has_three_seeds_and_required_temperatures():
    assert TEMPERATURES == (300.0, 600.0, 800.0, 900.0, 1000.0, 1100.0, 1200.0)
    assert all(len(values) >= 3 for values in SEEDS.values())


def test_planned_cases_use_both_providers_and_only_reduced_targets():
    rows = planned_cases()
    assert {row[5] for row in rows} == {"PF", "FEMCZM"}
    assert {row[4] for row in rows} == {100.0, 300.0}
    for material, candidate, temperature, seed, target, provider in rows:
        if target == 300.0:
            assert candidate != CONTROL_IDS[material]
            assert seed == SEEDS[material][0]


def test_all_strict_canonical_dbtt_survivors_are_stress_tested():
    required = {
        "oneD_v2_dbtt_R_3e168d9381eafa7b",
        "oneD_v2_dbtt_R_41a3f2096b3c4c54",
        "oneD_v2_dbtt_R_7b3b5b45354e64ef",
        "oneD_v2_dbtt_R_c2dabb42101cf812",
        "oneD_v2_dbtt_R_c4859a34963f15af",
    }
    assert required <= set(SURVIVORS["DBTT"])
