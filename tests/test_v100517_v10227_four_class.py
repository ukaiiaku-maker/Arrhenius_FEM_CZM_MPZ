from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import pytest

import arrhenius_fracture.audited_pf_parameter_bridge_v100517_v10227 as bridge


RUNNER = Path("run_v10_0_5_17_v10227_four_class_campaign.sh")
SMOKE = Path("run_v10_0_5_17_v10227_four_class_300_1000K_20um_smoke.sh")


def _base_row(option: str, candidate: str, material_class: str) -> dict[str, object]:
    return {
        "option_key": option,
        "candidate_id": candidate,
        "material_class": material_class,
        "role": f"paper primary {material_class}",
        "mechanism_summary": "synthetic exact-transfer test row",
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


def _install_tree(root: Path, rows: list[dict[str, object]]) -> tuple[Path, str, str]:
    package = root / "arrhenius_fracture"
    materials = package / "data" / "materials"
    materials.mkdir(parents=True)
    (package / bridge.MAPPING_MODULE).write_text(
        "from pathlib import Path\n"
        f'REGISTRY = "{bridge.REGISTRY_NAME}"\n'
        f'SELECTION = "{bridge.SELECTION_NAME}"\n'
        "def load_valid_options(): return {}\n"
        "VALID_OPTIONS = load_valid_options()\n"
    )
    (package / bridge.AUDITED_MODULE).write_text(
        "from . import sharp_front_v10_2_27 as _entry\n"
        "from .persistent_site_audited_engine_v10221 import AuditedPersistentSiteStateResolvedTipEngine\n"
        "from .persistent_site_bracket_fix_v10221 import install_backstress_complementarity_fix\n"
        "from .persistent_site_physical_width_v10222 import install_physical_front_width\n"
    )
    registry = materials / bridge.REGISTRY_NAME
    with registry.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    fingerprint = bridge.active_fingerprint(rows)
    selection = {
        "candidate_count": 4,
        "source_active_parameter_fingerprint_sha256": fingerprint,
        "primary_candidates": [
            {
                "option_key": row["option_key"],
                "candidate_id": row["candidate_id"],
                "paper_material_class": row["material_class"],
            }
            for row in rows
        ],
    }
    (materials / bridge.SELECTION_NAME).write_text(json.dumps(selection))
    registry_hash = hashlib.sha256(registry.read_bytes()).hexdigest()
    return root, registry_hash, fingerprint


def _four_rows() -> list[dict[str, object]]:
    rows = []
    for index, (option, (candidate, material_class)) in enumerate(bridge.EXPECTED_OPTIONS.items()):
        row = _base_row(option, candidate, material_class)
        row["cleave_G00_eV"] = 2.1 + 0.1 * index
        row["rho_source0_m2"] = 3.0e15 + index * 1.0e14
        rows.append(row)
    return rows


def test_bridge_loads_exact_v10227_row_and_fingerprint(tmp_path, monkeypatch):
    root, registry_hash, fingerprint = _install_tree(tmp_path, _four_rows())
    monkeypatch.setattr(bridge, "EXPECTED_REGISTRY_SHA256", registry_hash)
    monkeypatch.setattr(bridge, "EXPECTED_ACTIVE_FINGERPRINT_SHA256", fingerprint)
    option = "v913_paper_weakT01_0129902_persistent_sites"
    candidate, audit = bridge.load_audited_parameter_option(root, option)
    assert candidate.option_key == option
    assert candidate.candidate_id == "v913_zeroD_sobol_0129902"
    assert audit["material_class"] == "weakT"
    assert audit["active_parameter_fingerprint_sha256"] == fingerprint
    assert audit["mechanics_changed"] is False
    assert audit["source_closure_changed"] is False
    assert audit["stochastic_cleavage_law_changed"] is False


def test_bridge_rejects_nonzero_recovery(tmp_path, monkeypatch):
    rows = _four_rows()
    rows[2]["recovery_H0_eV"] = 0.2
    root, registry_hash, fingerprint = _install_tree(tmp_path, rows)
    monkeypatch.setattr(bridge, "EXPECTED_REGISTRY_SHA256", registry_hash)
    monkeypatch.setattr(bridge, "EXPECTED_ACTIVE_FINGERPRINT_SHA256", fingerprint)
    with pytest.raises(ValueError, match="persistent-site closure requires recovery_H0_eV=0"):
        bridge.load_audited_parameter_option(
            root, "v913_paper_weakT01_0129902_persistent_sites"
        )


def test_runner_selects_only_new_four_rows_and_preserves_mechanics():
    text = RUNNER.read_text()
    for option in bridge.EXPECTED_OPTIONS:
        assert option in text
    assert "v913_paper_weakT01_0257068_persistent_sites" not in text
    assert "v913_paper_ceramic01_0189364_persistent_sites" not in text
    assert "TARGET_EXT_UM=${TARGET_EXT_UM:-1000}" in text
    assert "CRYSTAL_THETA_DEG=${CRYSTAL_THETA_DEG:-30}" in text
    assert "MIN_GLOBAL_FORWARD=${MIN_GLOBAL_FORWARD:-0.05}" in text
    assert "--plane-gate-global" in text
    assert "--max-fronts 1" in text
    assert "--crack-backend adaptive_czm" in text
    assert "diagnose_v100517_j_scale.py" in text
    assert "mode_i_first_passage_v10_0_5_17_v10227_four_class" in text


def test_smoke_wrapper_is_eight_cases_at_20um():
    text = SMOKE.read_text()
    assert 'TEMPERATURES=${TEMPERATURES:-"300 1000"}' in text
    assert "TARGET_EXT_UM=${TARGET_EXT_UM:-20}" in text
    assert "MAX_JOBS=${MAX_JOBS:-2}" in text
    assert "RUN_J_SCALE_DIAGNOSTIC=${RUN_J_SCALE_DIAGNOSTIC:-1}" in text
