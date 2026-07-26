from __future__ import annotations

from pathlib import Path

import pytest

from arrhenius_fracture.four_class_parameter_bridge_v100518 import (
    EXPECTED_OPTIONS,
    INSTALLED_FOUR_ROW_FINGERPRINT_SHA256,
    SOURCE_REGISTRY_SHA256,
    SOURCE_WEAKT_CERAMIC_FINGERPRINT_SHA256,
    active_fingerprint,
    load_four_class_parameter_option,
)

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "runtime_inputs" / "v10_0_5_18_four_class"


def test_exact_four_class_registry_loads_all_canonical_rows():
    observed = {}
    for option, expected in EXPECTED_OPTIONS.items():
        candidate, audit = load_four_class_parameter_option(SOURCE, option)
        observed[option] = (candidate.candidate_id, audit["material_class"])
        assert audit["source_registry_sha256"] == SOURCE_REGISTRY_SHA256
        assert (
            audit["source_weakT_ceramic_active_parameter_fingerprint_sha256"]
            == SOURCE_WEAKT_CERAMIC_FINGERPRINT_SHA256
        )
        assert (
            audit["installed_four_row_active_parameter_fingerprint_sha256"]
            == INSTALLED_FOUR_ROW_FINGERPRINT_SHA256
        )
        assert audit["PF_runtime_dependency"] is False
        assert audit["parameter_values_manually_reconstructed"] is False
        assert audit["mechanics_changed_by_parameter_transfer"] is False
        assert audit["source_closure_changed_by_parameter_transfer"] is False
        assert audit["persistent_sites"] is True
        assert audit["finite_source_inventory"] is False
        assert audit["source_refresh"] is False
        assert audit["explicit_recovery"] is False
        assert audit["dynamic_tip_radius"] is True
        assert audit["dynamic_front_width"] is True
    assert observed == EXPECTED_OPTIONS


def test_new_weakT_and_ceramic_rows_have_exact_selected_values():
    weak, weak_audit = load_four_class_parameter_option(
        SOURCE, "v913_paper_weakT01_0129902_persistent_sites"
    )
    ceramic, ceramic_audit = load_four_class_parameter_option(
        SOURCE, "v913_paper_ceramic01_0077080_persistent_sites"
    )
    assert weak.candidate_id == "v913_zeroD_sobol_0129902"
    assert weak.cleave_G00_eV == pytest.approx(2.683425398543477, rel=0, abs=0)
    assert weak.emit_G00_eV == pytest.approx(1.6985005037859082, rel=0, abs=0)
    assert weak.rho_source0_m2 == pytest.approx(18895562156334.47, rel=0, abs=0)
    assert weak.c_blunt == pytest.approx(1.9507201686501503, rel=0, abs=0)
    assert weak_audit["material_class"] == "weakT"

    assert ceramic.candidate_id == "v913_zeroD_sobol_0077080"
    assert ceramic.cleave_G00_eV == pytest.approx(3.673669350333512, rel=0, abs=0)
    assert ceramic.emit_G00_eV == pytest.approx(2.798553348518908, rel=0, abs=0)
    assert ceramic.rho_source0_m2 == pytest.approx(1174563748149.2874, rel=0, abs=0)
    assert ceramic.c_blunt == pytest.approx(0.3964464701712131, rel=0, abs=0)
    assert ceramic_audit["material_class"] == "ceramic"


def test_source_fingerprint_is_two_row_and_installed_fingerprint_is_four_row():
    rows = []
    for option in EXPECTED_OPTIONS:
        _, audit = load_four_class_parameter_option(SOURCE, option)
        rows.append(audit["selected_registry_row"])
    assert active_fingerprint(rows) == INSTALLED_FOUR_ROW_FINGERPRINT_SHA256
    transferred = [
        row for row in rows if row["material_class"] in {"weakT", "ceramic"}
    ]
    assert (
        active_fingerprint(transferred)
        == SOURCE_WEAKT_CERAMIC_FINGERPRINT_SHA256
    )


def test_backup_options_are_not_exposed_as_production_rows():
    with pytest.raises(KeyError):
        load_four_class_parameter_option(
            SOURCE, "v913_paper_weakT02_0008816_persistent_sites"
        )
    with pytest.raises(KeyError):
        load_four_class_parameter_option(
            SOURCE, "v913_paper_ceramic02_0008536_persistent_sites"
        )
