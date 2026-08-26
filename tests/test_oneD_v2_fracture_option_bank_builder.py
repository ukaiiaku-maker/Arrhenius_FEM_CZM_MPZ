from scripts.finalize_oneD_v2_peak_dbtt_R_and_option_bank import (
    CLASS_ORDER, CONTROLS, FORBIDDEN_MATERIAL_FIELDS, SHORTLIST_COUNTS,
    TARGET_BANK_COUNTS, canonical_json, material_hash,
)


def test_option_bank_targets_all_four_classes():
    assert tuple(TARGET_BANK_COUNTS) == CLASS_ORDER
    assert all(12 <= count <= 24 for count in TARGET_BANK_COUNTS.values())
    assert all(2 <= SHORTLIST_COUNTS[name] <= 8 for name in CLASS_ORDER)


def test_controls_are_explicit_for_all_classes():
    assert set(CONTROLS) == set(CLASS_ORDER)


def test_backend_reduction_fields_are_forbidden_material_coordinates():
    assert "hazard_progress_scale" in FORBIDDEN_MATERIAL_FIELDS
    assert "reload_gap_threshold_m" in FORBIDDEN_MATERIAL_FIELDS
