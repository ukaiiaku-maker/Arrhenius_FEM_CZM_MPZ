from __future__ import annotations

import hashlib
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from arrhenius_fracture import emission_derived_plasticity as installed_pt
from arrhenius_fracture import mode_i_first_passage_v10_0_5_13_5_barrier_only as v135
from arrhenius_fracture import mode_i_first_passage_v10_0_5_14_persistent_site as v14
from arrhenius_fracture import mode_i_first_passage_v10_0_5_18_3_9_atomic_path_corridor as base
from arrhenius_fracture import mode_i_first_passage_v10_0_5_18_4_0_theta0_pf_full_field_parity as full
from arrhenius_fracture import mode_i_first_passage_v10_0_5_18_4_0_theta0_pf_parity as control
from arrhenius_fracture import refinement_audit_v100516 as refinement_audit
from arrhenius_fracture.active_only_kernel_family_compat_v10051840 import (
    _normalized_payload,
)
from arrhenius_fracture.bulk_pt_detailed_balance_v1041_exact import (
    EmissionDerivedPeierlsTaylorModel,
    ExpFloorSurface,
    config_from_dislocation_config,
)
from arrhenius_fracture.signed_kernel_family_v1005141 import FAMILY_SCHEMA


def _row() -> dict[str, float]:
    return {
        "Tref_K": 481.33,
        "rho_forest_floor_m2": 5.0e12,
        "emit_G00_eV": 2.0446315124630927,
        "emit_gT_eV_per_K": 0.000941995373927,
        "emit_sigc0_GPa": 7.205215461552143,
        "emit_sT_GPa_per_K": 0.0013325903080403,
        "emit_exp_a": 0.0858719444833695,
        "emit_exp_n": 1.1423492641188204,
        "emit_floor_frac": 0.0058420097608974,
        "peierls_H0_eV": 6.368386985268444,
        "peierls_activation_entropy_kB": 12.2265643812716,
        "peierls_exp_a": 1.9097672308795155,
        "peierls_exp_n": 1.1372957159765065,
        "peierls_stress_fraction": 0.5773502691896258,
        "peierls_nu0_s": 1.0e12,
        "taylor_H0_eV": 1.6986502355337143,
        "taylor_activation_entropy_kB": 38.45378890633583,
        "taylor_exp_a": 0.4561179580539465,
        "taylor_exp_n": 1.2506696581840515,
        "taylor_stress_fraction": 0.5773502691896258,
        "taylor_nu0_s": 1.0e11,
        "taylor_corr_rho_c_m2": 5.878007397057801e13,
        "taylor_corr_scale": 0.8638954413677038,
    }


def _candidate() -> SimpleNamespace:
    row = _row()
    return SimpleNamespace(
        L_pz_um_recommended=50.0,
        n_bins_recommended=80,
        source_sites_per_system_provenance=141.0590567476921,
        source_refresh_length_um_provenance=0.0,
        source_zone_length_um=2.0,
        peierls_stress_fraction=row["peierls_stress_fraction"],
        taylor_stress_fraction=row["taylor_stress_fraction"],
        taylor_corr_rho_c_m2=row["taylor_corr_rho_c_m2"],
        taylor_corr_scale=row["taylor_corr_scale"],
        encounter_efficiency=9.160246308716648,
        rho_forest_floor_m2=row["rho_forest_floor_m2"],
    )


def test_detailed_balance_is_exactly_zero_at_zero_stress():
    surface = ExpFloorSurface(
        G00_eV=2.0,
        gT_eV_per_K=1.0e-3,
        sigc0_Pa=5.0e9,
        sT_Pa_per_K=0.0,
        alpha=0.5,
        exponent=1.2,
        floor_fraction=0.01,
        floor_min_eV=1.0e-4,
        floor_max_fraction=0.95,
        Tref_K=481.33,
        attempt_frequency_s=1.0e12,
    )
    rates = surface.rate_s(np.zeros(4), 1000.0)
    assert np.array_equal(rates, np.zeros(4))
    assert surface.rate_s(np.array([1.0e9]), 1000.0)[0] > 0.0


def test_exact_bulk_mapping_uses_direct_selected_row_fields():
    cfg = SimpleNamespace()
    full._exact_bulk_parameters(cfg, _row())
    model_cfg = config_from_dislocation_config(cfg)
    row = _row()
    assert model_cfg.peierls.G00_eV == row["peierls_H0_eV"]
    assert model_cfg.taylor.G00_eV == row["taylor_H0_eV"]
    assert model_cfg.peierls.sigc0_Pa == row["peierls_stress_fraction"] * row["emit_sigc0_GPa"] * 1.0e9
    assert model_cfg.taylor.sigc0_Pa == row["taylor_stress_fraction"] * row["emit_sigc0_GPa"] * 1.0e9
    assert model_cfg.taylor_renewal_time_s == 1.0e-9
    assert model_cfg.mobile_density_floor_m2 == row["rho_forest_floor_m2"]


def test_exact_bulk_rate_uses_pf_hit_order_and_mobile_law():
    cfg = SimpleNamespace()
    full._exact_bulk_parameters(cfg, _row())
    model = EmissionDerivedPeierlsTaylorModel(config_from_dislocation_config(cfg))
    stress = np.array([2.0e9, 4.0e9])
    rho = np.array([5.0e12, 2.0e14])
    out = model.rates(stress, rho, 1000.0, 2.5e-10)
    expected_order = 1.0 + _row()["taylor_corr_scale"] * (
        rho / _row()["taylor_corr_rho_c_m2"]
    )
    expected_mobile = np.minimum(0.01 * rho, 1.0e14)
    assert np.allclose(out["taylor_m_eff"], expected_order, rtol=1.0e-13)
    assert np.allclose(out["rho_mobile_m2"], expected_mobile, rtol=1.0e-13)
    assert np.all(out["equivalent_plastic_rate_s"] >= 0.0)


def test_active_only_empty_wake_schema_is_zero_only_and_does_not_touch_active_data():
    active_i = [[1.0, 2.0], [3.0, 4.0]]
    active_ii = [[0.1, 0.2], [0.3, 0.4]]
    payload = {
        "schema": FAMILY_SCHEMA,
        "wake_kernel_forced_zero": True,
        "wake_shielding_supported": False,
        "constitutive_K_shield_cap": False,
        "wake_x_m": [],
        "activation_to_line_content_by_system": [2.0, 3.0],
        "states": [
            {
                "active_kernel_I_Pa_sqrt_m_per_signed_line": active_i,
                "active_kernel_II_Pa_sqrt_m_per_signed_line": active_ii,
                "wake_kernel_I_Pa_sqrt_m_per_signed_line": [[], []],
                "wake_kernel_II_Pa_sqrt_m_per_signed_line": [[], []],
            },
            {
                "active_kernel_I_Pa_sqrt_m_per_signed_line": active_i,
                "active_kernel_II_Pa_sqrt_m_per_signed_line": active_ii,
            },
        ],
    }

    normalized, applied = _normalized_payload(payload)

    assert applied is True
    assert normalized["wake_x_m"] == [0.0]
    assert normalized["states"][0][
        "active_kernel_I_Pa_sqrt_m_per_signed_line"
    ] == active_i
    assert normalized["states"][1][
        "active_kernel_II_Pa_sqrt_m_per_signed_line"
    ] == active_ii
    for state in normalized["states"]:
        assert state["wake_kernel_I_Pa_sqrt_m_per_signed_line"] == [[0.0], [0.0]]
        assert state["wake_kernel_II_Pa_sqrt_m_per_signed_line"] == [[0.0], [0.0]]
    contract = normalized["fem_empty_wake_schema_compatibility"]
    assert contract["runtime_wake_shielding_enabled"] is False
    assert contract["active_kernel_modified"] is False
    assert contract["physics_modified"] is False


def test_long_corridor_constructor_is_captured_for_post_run_audit(monkeypatch):
    mesh = SimpleNamespace(
        production_refinement_radius_m=330.0e-6,
        production_refinement_policy="test_long_corridor",
    )
    monkeypatch.setattr(
        v135,
        "make_physical_refinement_mesh_v1005135",
        lambda *args, **kwargs: mesh,
    )

    with refinement_audit.installed_refinement_audit_v100516():
        observed = v135.make_physical_refinement_mesh_v1005135()
        payload = refinement_audit.refinement_audit_payload_v100516()

    assert observed is mesh
    assert payload["captured_mesh_available"] is True
    assert payload["captured_mesh_has_refinement_metadata"] is True
    assert payload["production_refinement_radius_m"] == 330.0e-6
    assert payload["long_corridor_constructor_capture_active"] is True


def test_full_field_entry_installs_and_restores_overlay(monkeypatch, tmp_path: Path):
    kernel = tmp_path / "family.json"
    kernel.write_text("reference kernel\n")
    digest = hashlib.sha256(kernel.read_bytes()).hexdigest()
    monkeypatch.setattr(control, "REFERENCE_KERNEL_SHA256", digest)
    monkeypatch.setenv("CLEAVAGE_HAZARD_SEED", "8666")

    original_model = installed_pt.EmissionDerivedPeierlsTaylorModel
    original_policy = v14.persistent_site_policy
    original_fields = base._fields
    observed = {}

    def fake_base_main(args):
        observed["args"] = list(args)
        observed["model"] = installed_pt.EmissionDerivedPeierlsTaylorModel
        observed["policy"] = v14.persistent_site_policy(_candidate())
        observed["fields"] = base._fields()
        return "ok"

    monkeypatch.setattr(base, "main", fake_base_main)
    out = tmp_path / "out"
    args = [
        "--parameter-option", control.REFERENCE_OPTION,
        "--temperatures", "1000",
        "--bulk-plasticity-mode", "full_field",
        "--j-decomposition", "cluster",
        "--nx", "36", "--ny", "72",
        "--n-stagger", "2", "--max-fronts", "1",
        "--crystal-theta-deg", "0",
        "--dU", "2e-7", "--dt", "8.4",
        "--tip-h-fine", "1e-6", "--tip-ratio", "1.20",
        "--da-phys", "5e-6", "--adaptive-event-target", "0.15",
        "--signed-kernel-family", str(kernel),
        "--out", str(out),
    ]
    assert full.main(args) == "ok"
    assert observed["model"] is EmissionDerivedPeierlsTaylorModel
    assert observed["policy"]["bulk_plasticity_mode"] == full.SOLVER_BULK_MODE
    assert observed["policy"]["bulk_plasticity_semantic_mode"] == full.SEMANTIC_BULK_MODE
    mode_index = observed["args"].index("--bulk-plasticity-mode")
    assert observed["args"][mode_index + 1] == full.SOLVER_BULK_MODE
    assert observed["fields"]["PF_v10_4_1_bulk_parity_active"] is True
    assert installed_pt.EmissionDerivedPeierlsTaylorModel is original_model
    assert v14.persistent_site_policy is original_policy
    assert base._fields is original_fields
