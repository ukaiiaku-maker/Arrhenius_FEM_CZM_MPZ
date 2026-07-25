from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from arrhenius_fracture.audited_pf_parameter_bridge_v100517 import (
    load_audited_parameter_option,
)


def _row(option: str, candidate: str) -> dict[str, object]:
    return {
        "option_key": option,
        "candidate_id": candidate,
        "material_class": "DBTT",
        "role": "paper primary",
        "mechanism_summary": "synthetic audited test row",
        "validation_status": "test",
        "Tref_K": 481.33,
        "n_slip_channels": 2,
        "rho_forest_floor_m2": 5.0e12,
        "peierls_stress_fraction": 0.5773502691896258,
        "taylor_stress_fraction": 0.5773502691896258,
        "mobile_shield_fraction": 0.0,
        "source_recovery_rate_s": 0.0,
        "L_pz_um_recommended": 50.0,
        "n_bins_recommended": 80,
        "cleave_G00_eV": 2.1,
        "cleave_gT_eV_per_K": 0.003,
        "cleave_sigc0_GPa": 6.4,
        "cleave_sT_GPa_per_K": -0.001,
        "cleave_exp_a": 0.7,
        "cleave_exp_n": 1.5,
        "cleave_floor_frac": 0.02,
        "emit_G00_eV": 3.2,
        "emit_gT_eV_per_K": -0.002,
        "emit_sigc0_GPa": 7.1,
        "emit_sT_GPa_per_K": 0.001,
        "emit_exp_a": 0.4,
        "emit_exp_n": 2.2,
        "emit_floor_frac": 0.03,
        "peierls_H0_eV": 4.3,
        "peierls_activation_entropy_kB": -8.0,
        "peierls_exp_a": 0.8,
        "peierls_exp_n": 1.2,
        "taylor_H0_eV": 1.4,
        "taylor_activation_entropy_kB": 18.0,
        "taylor_exp_a": 0.6,
        "taylor_exp_n": 1.1,
        "taylor_corr_rho_c_m2": 2.0e13,
        "taylor_corr_scale": 1.7,
        "source_sites_per_system": 141.0590567476921,
        "encounter_efficiency": 9.160246308716648,
        "retained_recovery_rate_s": 0.0,
        "source_refresh_length_um": 0.0,
        "c_blunt": 2.5,
        "peierls_nu0_s": 1.0e12,
        "taylor_nu0_s": 1.0e11,
        "rho_source0_m2": 3.0e15,
        "recovery_nu0_s": 0.0,
        "recovery_H0_eV": 0.0,
        "recovery_activation_entropy_kB": 0.0,
        "reference_source_area_um2": 25.0,
        "reference_front_width_um": 10.0,
        "source_zone_length_um": 2.0,
        "legacy_source_sites_active": 0,
        "legacy_source_refresh_active": 0,
        "explicit_recovery_active": 0,
    }


def _install_tree(root: Path, entry: str, option: str, candidate: str, row: dict[str, object]) -> None:
    package = root / "arrhenius_fracture"
    materials = package / "data" / "materials"
    materials.mkdir(parents=True)
    if entry == "v10.2.25":
        stem = "sharp_front_v10_2_25"
        registry = "v10_2_25_v913_paper_campaign_registry.csv"
        selection = "v10_2_25_v913_paper_campaign_selection.json"
    else:
        stem = "sharp_front_v10_2_26"
        registry = "v10_2_26_v913_weakT_ceramic_registry.csv"
        selection = "v10_2_26_v913_weakT_ceramic_selection.json"
    (package / f"{stem}_audited.py").write_text(
        f"from . import {stem} as _entry\n"
        "from .persistent_site_audited_engine_v10221 import AuditedPersistentSiteStateResolvedTipEngine\n"
        "from .persistent_site_bracket_fix_v10221 import install_backstress_complementarity_fix\n"
        "from .persistent_site_physical_width_v10222 import install_physical_front_width\n"
    )
    (package / f"{stem}.py").write_text(
        "from pathlib import Path\n"
        f'DEFAULT_REGISTRY = Path(__file__).parent / "data" / "materials" / "{registry}"\n'
        f'SELECTION_RECORD = Path(__file__).parent / "data" / "materials" / "{selection}"\n'
        f"VALID_OPTIONS = {{{option!r}: {candidate!r}}}\n"
    )
    with (materials / registry).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row))
        writer.writeheader()
        writer.writerow(row)
    (materials / selection).write_text(
        json.dumps(
            {
                "primary_candidates": [
                    {
                        "option_key": option,
                        "candidate_id": candidate,
                        "response_class": "synthetic",
                    }
                ]
            }
        )
    )


@pytest.mark.parametrize("entry", ["v10.2.25", "v10.2.26"])
def test_bridge_loads_exact_row_through_audited_option_map(tmp_path, entry):
    option = "v913_paper_test_persistent_sites"
    candidate_id = "v913_zeroD_sobol_test"
    row = _row(option, candidate_id)
    _install_tree(tmp_path, entry, option, candidate_id, row)
    candidate, audit = load_audited_parameter_option(tmp_path, entry, option)
    assert candidate.option_key == option
    assert candidate.candidate_id == candidate_id
    assert candidate.cleave_G00_eV == pytest.approx(float(row["cleave_G00_eV"]))
    assert candidate.emit_sigc0_GPa == pytest.approx(float(row["emit_sigc0_GPa"]))
    assert candidate.peierls_H0_eV == pytest.approx(float(row["peierls_H0_eV"]))
    assert candidate.taylor_corr_scale == pytest.approx(float(row["taylor_corr_scale"]))
    assert candidate.rho_source0_m2 == pytest.approx(float(row["rho_source0_m2"]))
    assert audit["selected_through_audited_entry_option_map"] is True
    assert audit["parameter_values_manually_reconstructed"] is False
    assert audit["selection_metadata"]["candidate_id"] == candidate_id
    assert audit["finite_source_inventory"] is False
    assert audit["source_refresh"] is False
    assert audit["explicit_recovery"] is False


def test_bridge_rejects_nonzero_refresh_or_recovery(tmp_path):
    option = "v913_paper_bad_persistent_sites"
    candidate_id = "v913_zeroD_sobol_bad"
    row = _row(option, candidate_id)
    row["source_refresh_length_um"] = 5.0
    _install_tree(tmp_path, "v10.2.25", option, candidate_id, row)
    with pytest.raises(ValueError, match="violates the audited persistent-site closure"):
        load_audited_parameter_option(tmp_path, "v10.2.25", option)


def test_bridge_rejects_mapping_registry_disagreement(tmp_path):
    option = "v913_paper_badmap_persistent_sites"
    row = _row(option, "candidate_in_registry")
    _install_tree(tmp_path, "v10.2.26", option, "candidate_in_mapping", row)
    with pytest.raises(ValueError, match="stable option mapping disagrees"):
        load_audited_parameter_option(tmp_path, "v10.2.26", option)
