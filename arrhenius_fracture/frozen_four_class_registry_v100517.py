"""Standalone frozen four-class parameter registry for FEM/CZM v10.0.5.17.

This module contains the exact four numerical rows selected for the paper
campaign.  It has no runtime dependency on a phase-field checkout, branch,
commit, package, or result.  The rows are validated by an embedded registry
SHA-256, exact option/candidate/class identities, persistent-site closure, and
a normalized active-parameter fingerprint before conversion to the existing
FEM/CZM persistent-site state contract.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
from dataclasses import asdict
from typing import Any

from . import audited_pf_parameter_bridge_v100517_v10227 as _legacy

BRIDGE_SCHEMA = "frozen_four_class_registry_v10_0_5_17"
FROZEN_REGISTRY_SHA256 = "4ba723c80abcfdd7101cee0afaa9b4104eebf2a8e847892528128664829c7760"
EXPECTED_ACTIVE_FINGERPRINT_SHA256 = "923463f936ddbdb5a997de7d1859b82774df2257b1faf1dde9fef44283816010"
EXPECTED_OPTIONS = _legacy.EXPECTED_OPTIONS

FROZEN_REGISTRY_CSV = """option_key,candidate_id,material_class,role,mechanism_summary,validation_status,Tref_K,n_slip_channels,rho_forest_floor_m2,peierls_stress_fraction,taylor_stress_fraction,mobile_shield_fraction,source_recovery_rate_s,L_pz_um_recommended,n_bins_recommended,cleave_G00_eV,cleave_gT_eV_per_K,cleave_sigc0_GPa,cleave_sT_GPa_per_K,cleave_exp_a,cleave_exp_n,cleave_floor_frac,emit_G00_eV,emit_gT_eV_per_K,emit_sigc0_GPa,emit_sT_GPa_per_K,emit_exp_a,emit_exp_n,emit_floor_frac,peierls_H0_eV,peierls_activation_entropy_kB,peierls_exp_a,peierls_exp_n,taylor_H0_eV,taylor_activation_entropy_kB,taylor_exp_a,taylor_exp_n,taylor_corr_rho_c_m2,taylor_corr_scale,source_sites_per_system,encounter_efficiency,retained_recovery_rate_s,source_refresh_length_um,c_blunt,peierls_nu0_s,taylor_nu0_s,rho_source0_m2,recovery_nu0_s,recovery_H0_eV,recovery_activation_entropy_kB,reference_source_area_um2,reference_front_width_um,source_zone_length_um,legacy_source_sites_active,legacy_source_refresh_active,explicit_recovery_active
v913_paper_peak01_0242980_persistent_sites,v913_zeroD_sobol_0242980,peak,paper primary peak 1,Delayed abrupt peak centered at 1100 K with strong high-temperature softening.,"v10.2.25 paper campaign selection; exact active parameters transferred from v10.2.23; persistent sites, no finite source inventory, no source refresh, and no explicit recovery.",481.33,2,5000000000000.0,0.5773502691896258,0.5773502691896258,0.0,0.0,50.0,80,4.011803912930191,0.0081468212408944,7.349638969637454,-0.0012099820682778,1.3227705595549195,1.2558047282509506,0.0184513317135146,2.0446315124630927,0.000941995373927,7.205215461552143,0.0013325903080403,0.0858719444833695,1.1423492641188204,0.0058420097608974,6.368386985268444,12.2265643812716,1.9097672308795155,1.1372957159765065,1.6986502355337143,38.45378890633583,0.4561179580539465,1.2506696581840515,58780073970578.01,0.8638954413677038,141.0590567476921,9.160246308716648,0.0,0.0,2.687256144359708,1000000000000.0,100000000000.0,1.4102520844479966e+16,0.0,0.0,0.0,25.0,10.0,2.0,0,0,0
v913_paper_dbtt01_0202500_persistent_sites,v913_zeroD_sobol_0202500,DBTT,paper primary DBTT 1,"Preferred classical DBTT example: low-temperature shelf, sharp transition near 1000 K, and stable upper shelf.","v10.2.25 paper campaign selection; exact active parameters transferred from v10.2.24; persistent sites, no finite source inventory, no source refresh, and no explicit recovery.",481.33,2,5000000000000.0,0.5773502691896258,0.5773502691896258,0.0,0.0,50.0,80,2.357721022795886,0.0056085130451247,6.7501407796517015,-0.0012041387008503,0.4968366474844515,2.618804107885808,0.026239358201278,3.344896922819316,0.0061668124636635,2.7143102660775185,-0.0052233074344694,0.1572539726179093,0.5156166724395007,0.0103682254353579,3.525720003992319,-36.33932862430811,1.448395521938801,1.871259921230376,1.4998139148950578,31.365934647619724,0.7048910272121429,0.5389550719410181,22042481205294.297,1.2690231726199324,141.0590567476921,9.160246308716648,0.0,0.0,2.929164256900549,1000000000000.0,100000000000.0,523948859759023.8,0.0,0.0,0.0,25.0,10.0,2.0,0,0,0
v913_paper_weakT01_0129902_persistent_sites,v913_zeroD_sobol_0129902,weakT,paper primary weak-temperature/FCC-like,Strict 100 um 1-D weak-temperature/FCC-like primary selected by minimum class score.,v10.2.27 exact active-parameter transfer from the final v9.13 five-temperature 100 um 1-D selection; 2-D validation required.,481.33,2,5000000000000.0,0.5773502691896258,0.5773502691896258,0.0,0.0,50.0,80,2.683425398543477,0.0060197069197893,4.890972775872797,-0.0002585870958864,1.1411482973955571,1.1264505988918243,0.0019271603769193,1.6985005037859082,0.0025603877846151,0.7007976029999554,-0.000113367350772,0.6522800833545626,0.8178969756700099,0.001963405800003,7.117342075519264,-47.75738567113876,0.2264547564554959,1.9374389592558143,3.2621826685033737,39.933723136782646,1.055419295541942,1.0061173520516604,237130686553740.03,0.0927493991876699,141.0590567476921,9.160246308716648,0.0,0.0,1.9507201686501503,1000000000000.0,100000000000.0,18895562156334.47,0.0,0.0,0.0,25.0,10.0,2.0,0,0,0
v913_paper_ceramic01_0077080_persistent_sites,v913_zeroD_sobol_0077080,ceramic,paper primary ceramic-like,Strict 100 um 1-D ceramic-like primary selected by minimum class score.,v10.2.27 exact active-parameter transfer from the final v9.13 five-temperature 100 um 1-D selection; 2-D validation required.,481.33,2,5000000000000.0,0.5773502691896258,0.5773502691896258,0.0,0.0,50.0,80,3.673669350333512,0.00892692707479,4.497151743154973,-0.0022528592869639,1.4499619576428089,1.0335969480220228,0.0436070518879531,2.798553348518908,-0.002825030600652,9.278617694508284,0.0099400414247065,0.4376551472768187,0.791314263548702,0.0228240380649044,6.716475480981171,-43.48882682621479,1.178178071556613,1.0056940987706184,0.7904905830975623,-25.948952957987785,1.0212089663371444,0.599093577126041,4.043330255464039e+16,1.8013635924189044,141.0590567476921,9.160246308716648,0.0,0.0,0.3964464701712131,1000000000000.0,100000000000.0,1174563748149.2874,0.0,0.0,0.0,25.0,10.0,2.0,0,0,0
"""

FROZEN_SELECTION = {
    "schema": "v10.0.5.17_frozen_four_class_selection_v1",
    "candidate_count": 4,
    "canonical_option_order": list(EXPECTED_OPTIONS),
    "primary_candidates": [
        {
            "option_key": option,
            "candidate_id": candidate,
            "paper_material_class": material_class,
        }
        for option, (candidate, material_class) in EXPECTED_OPTIONS.items()
    ],
    "phase_field_solver_dependency": False,
    "phase_field_result_comparison_enabled": False,
    "transfer_policy": (
        "Exact frozen numerical rows only; no fitting, transformation, rounding, "
        "rescaling, temperature shifting, or substitution."
    ),
}


def _rows() -> list[dict[str, str]]:
    observed = hashlib.sha256(FROZEN_REGISTRY_CSV.encode("utf-8")).hexdigest()
    if observed != FROZEN_REGISTRY_SHA256:
        raise RuntimeError(
            "embedded frozen registry SHA-256 mismatch: "
            f"expected={FROZEN_REGISTRY_SHA256} observed={observed}"
        )
    rows = list(csv.DictReader(io.StringIO(FROZEN_REGISTRY_CSV)))
    _legacy._validate_rows(rows, FROZEN_SELECTION)
    fingerprint = _legacy.active_fingerprint(rows)
    if fingerprint != EXPECTED_ACTIVE_FINGERPRINT_SHA256:
        raise RuntimeError(
            "embedded frozen active-parameter fingerprint mismatch: "
            f"expected={EXPECTED_ACTIVE_FINGERPRINT_SHA256} observed={fingerprint}"
        )
    return rows


def load_frozen_parameter_option(parameter_option: str):
    if parameter_option not in EXPECTED_OPTIONS:
        raise KeyError(
            f"unknown frozen four-class option {parameter_option!r}; "
            f"allowed={sorted(EXPECTED_OPTIONS)}"
        )
    rows = _rows()
    selected = next(row for row in rows if row["option_key"] == parameter_option)
    candidate = _legacy._candidate_from_row(selected)
    audit: dict[str, Any] = {
        "schema": BRIDGE_SCHEMA,
        "parameter_entry": "v10.0.5.17-frozen-four-class",
        "parameter_option": parameter_option,
        "candidate_id": candidate.candidate_id,
        "material_class": selected.get("material_class"),
        "paper_role": selected.get("role"),
        "mechanism_summary": selected.get("mechanism_summary"),
        "validation_status": selected.get("validation_status"),
        "registry_path": "embedded://arrhenius_fracture/frozen_four_class_registry_v100517.py",
        "selection_path": "embedded://arrhenius_fracture/frozen_four_class_registry_v100517.py",
        "registry_sha256": FROZEN_REGISTRY_SHA256,
        "selection_sha256": hashlib.sha256(
            json.dumps(FROZEN_SELECTION, sort_keys=True, separators=(",", ":")).encode("utf-8")
        ).hexdigest(),
        "active_parameter_fingerprint_sha256": EXPECTED_ACTIVE_FINGERPRINT_SHA256,
        "selected_registry_row": dict(selected),
        "converted_candidate": asdict(candidate),
        "selected_row_sha256": hashlib.sha256(
            json.dumps(selected, sort_keys=True, separators=(",", ":")).encode("utf-8")
        ).hexdigest(),
        "parameter_values_manually_reconstructed": False,
        "selected_through_frozen_FEM_registry": True,
        "phase_field_solver_dependency": False,
        "phase_field_result_comparison_enabled": False,
        "mechanics_changed": False,
        "source_closure_changed": False,
        "stochastic_cleavage_law_changed": False,
        "persistent_sites": True,
        "finite_source_inventory": False,
        "source_refresh": False,
        "explicit_recovery": False,
    }
    return candidate, audit


def validate_frozen_registry() -> dict[str, Any]:
    rows = _rows()
    return {
        "schema": "v10.0.5.17_frozen_four_class_preflight",
        "n_options": len(rows),
        "options_in_production_order": list(EXPECTED_OPTIONS),
        "classes": [row["material_class"] for row in rows],
        "registry_sha256": FROZEN_REGISTRY_SHA256,
        "active_parameter_fingerprint_sha256": EXPECTED_ACTIVE_FINGERPRINT_SHA256,
        "phase_field_solver_dependency": False,
        "phase_field_result_comparison_enabled": False,
        "persistent_sites": True,
        "finite_source_inventory": False,
        "source_refresh": False,
        "explicit_recovery": False,
    }


__all__ = [
    "BRIDGE_SCHEMA",
    "EXPECTED_ACTIVE_FINGERPRINT_SHA256",
    "EXPECTED_OPTIONS",
    "FROZEN_REGISTRY_SHA256",
    "load_frozen_parameter_option",
    "validate_frozen_registry",
]
