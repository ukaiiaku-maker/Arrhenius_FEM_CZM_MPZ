"""Audited v10.0.3 progressive Mode-I kinetic-CZM entry point.

This point release replaces the eager copied-global transform used by v10.0.2
with a delayed transform constructed after the live v9.11 engine, mechanics, J,
and plasticity wrappers are installed.
"""
from __future__ import annotations

import csv
import json
import math
import os
from pathlib import Path
import sys
from typing import Any

from . import crack_backend as _cb
from . import fem as _fem
from . import mesh as _mesh
from . import mixed_mode_first_passage_v9_11 as v911
from . import mode_i_first_passage_v9_11 as modei911
from . import mode_i_first_passage_v9_18_3 as _v9183
from . import mode_i_first_passage_v9_18_5 as _v9185
from . import mode_i_first_passage_v9_18_5_6 as _v91856
from . import sharp_front as _sharp_front
from .kinetic_campaign_czm import KineticCampaignCZMConfig
from .kinetic_campaign_czm_v1003 import engine_factory_v1003
from .kinetic_progressive_2d_v1003 import (
    build_delayed_progressive_run_2d_v1003,
    progressive_runtime_payload_v1003,
    reset_progressive_runtime_v1003,
    write_progressive_runtime_audit_v1003,
)
from .mode_i_first_passage_v10_0 import (
    _option_value,
    _replace_option,
    parser,
)
from .pf_equivalent_material_manifest import (
    PF_SOURCE,
    load_material_manifest,
    pf_manifest_path,
)

MODEL_ID = (
    "FEM_CZM_Mode_I_kinetic_campaign_czm_v10_0_3_"
    "delayed_live_binding_progressive_event_lifecycle"
)


def _set_env(name: str, value: str, saved: dict[str, str | None]) -> None:
    if name not in saved:
        saved[name] = os.environ.get(name)
    os.environ[name] = str(value)


def _restore_env(saved: dict[str, str | None]) -> None:
    for name, value in saved.items():
        if value is None:
            os.environ.pop(name, None)
        else:
            os.environ[name] = value


def _last_step(root: Path, temperature_K: float) -> dict[str, str]:
    path = root / f"steps_{int(round(float(temperature_K))):04d}K.csv"
    if not path.exists():
        return {}
    with path.open(newline="") as fp:
        rows = list(csv.DictReader(fp))
    return rows[-1] if rows else {}


def _as_float(row: dict[str, str], key: str) -> float | None:
    try:
        value = float(row.get(key, ""))
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def _max_column(root: Path, temperature_K: float, key: str) -> float | None:
    path = root / f"steps_{int(round(float(temperature_K))):04d}K.csv"
    if not path.exists():
        return None
    values = []
    with path.open(newline="") as fp:
        for row in csv.DictReader(fp):
            value = _as_float(row, key)
            if value is not None:
                values.append(value)
    return max(values) if values else None


def source_population_bound(
    capacity: float,
    extension_m: float,
    refresh_length_m: float,
) -> float:
    """Conservative finite-source population bound including advance refresh.

    At zero extension, at most one complete source inventory can have emitted.
    During advance, the maximum possible additional refreshed inventory is
    bounded by ``capacity * extension / refresh_length``. This deliberately
    assumes instantaneous depletion after every refresh, making it conservative
    for retained/mobile population checks without weakening the finite-source law.
    """

    cap = max(float(capacity), 0.0)
    distance = max(float(extension_m), 0.0)
    length = max(float(refresh_length_m), 1.0e-30)
    return cap * (1.0 + distance / length)


def main(argv: list[str] | None = None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    opts, remaining = parser().parse_known_args(user_args)
    if opts.czm_opening_coupling != "clock_linear":
        raise SystemExit(
            "v10.0.3 progressive validation requires "
            "--czm-opening-coupling clock_linear"
        )
    if opts.v10_material_source != PF_SOURCE:
        raise SystemExit(
            "v10.0.3 progressive equivalence requires PF v10.1.7.1 rows"
        )

    material_class = (
        opts.v10_material_class
        or _option_value(remaining, "--mpz-material-class")
        or "ceramic"
    )
    manifest = load_material_manifest(material_class, parameter_source=PF_SOURCE)
    _replace_option(
        remaining,
        "--mpz-material-manifest",
        str(pf_manifest_path(material_class)),
    )
    _replace_option(remaining, "--mpz-material-class", manifest.name)
    _replace_option(remaining, "--crack-backend", "adaptive_czm")
    if _option_value(remaining, "--mpz-length-um") is None:
        remaining.extend(["--mpz-length-um", "100"])
    if _option_value(remaining, "--mpz-n-bins") is None:
        remaining.extend(["--mpz-n-bins", "200"])
    if _option_value(remaining, "--crystal-theta-deg") is None:
        remaining.extend(["--crystal-theta-deg", "45"])

    kinetic_cfg = KineticCampaignCZMConfig(
        max_action_substep=opts.max_action_substep,
        max_translation_substep_m=opts.max_translation_substep_m,
        min_substep_s=opts.min_kinetic_substep_s,
        max_internal_steps=opts.max_internal_steps,
        coupling_scheme="strang",
        wake_shielding=False,
        active_shielding=True,
        signed_active_shielding=True,
        mobile_shield_fraction=1.0,
        backstress_scale=1.0,
        source_refresh_scale=1.0,
    ).validate()

    original_factory = v911._engine_factory
    original_run_2d = _sharp_front.run_2d
    original_base_build = _sharp_front.build_engine
    original_make_mesh = _mesh.make_tri_mesh
    original_assemble = _fem.assemble_mechanics
    original_solve = _fem.solve_dirichlet
    original_insert = _cb.AdaptiveCZMBackend._insert_target_in_incident_triangle
    original_advance = _cb.AdaptiveCZMBackend.advance
    engines: list[Any] = []
    env_saved: dict[str, str | None] = {}
    error: BaseException | None = None

    def patched_factory(original_build, context, mm, row):
        return engine_factory_v1003(
            original_build,
            context,
            mm,
            row,
            manifest,
            kinetic_cfg,
            engines,
        )

    reset_progressive_runtime_v1003()
    _v9185._RUNTIME["mesh"] = None
    _v9185._RUNTIME["corridor_centers"] = []
    _v9185._RUNTIME["component_anchor_history"] = []
    _v9185._RUNTIME["quality_vetoes"] = []
    _v91856._AUDIT["accepted_events"] = []
    _v91856._AUDIT["resolution_warnings"] = []
    _v91856._AUDIT["quality_vetoes"] = []
    _v91856._AUDIT["consecutive_veto_abort"] = None

    _v9185._make_corridor_mesh._original = original_make_mesh
    _v9185._assemble_with_mesh_capture._original = original_assemble
    _v91856._strict_quality_advance_v91856._original = original_advance

    _set_env("ARRHENIUS_CZM_OPENING_COUPLING", "clock_linear", env_saved)
    _set_env(
        "ARRHENIUS_MAX_TRIAL_DAMAGE_CHANGE",
        os.environ.get("ARRHENIUS_MAX_TRIAL_DAMAGE_CHANGE", "0.05"),
        env_saved,
    )
    target = _option_value(remaining, "--target-crack-extension-um")
    if target is not None and "ARRHENIUS_COMMITTED_TARGET_EXTENSION_UM" not in os.environ:
        _set_env("ARRHENIUS_COMMITTED_TARGET_EXTENSION_UM", target, env_saved)

    v911._engine_factory = patched_factory
    _mesh.make_tri_mesh = _v9185._make_corridor_mesh
    _fem.assemble_mechanics = _v9185._assemble_with_mesh_capture
    _fem.solve_dirichlet = _v9185._component_anchored_solve
    _cb.AdaptiveCZMBackend._insert_target_in_incident_triangle = (
        _v9183._edge_aware_insert_target_in_incident_triangle
    )
    _cb.AdaptiveCZMBackend.advance = _v91856._strict_quality_advance_v91856
    _sharp_front.run_2d = build_delayed_progressive_run_2d_v1003(
        original_run_2d,
        original_base_build,
    )

    results = None
    try:
        results = modei911.main(remaining)
    except BaseException as exc:
        error = exc
        raise
    finally:
        v911._engine_factory = original_factory
        _sharp_front.run_2d = original_run_2d
        _mesh.make_tri_mesh = original_make_mesh
        _fem.assemble_mechanics = original_assemble
        _fem.solve_dirichlet = original_solve
        _cb.AdaptiveCZMBackend._insert_target_in_incident_triangle = original_insert
        _cb.AdaptiveCZMBackend.advance = original_advance
        _restore_env(env_saved)
        _v91856._write_audit(remaining, error)

    out_value = _option_value(remaining, "--out")
    if out_value is None:
        raise RuntimeError("v10.0.3 progressive validation requires --out")
    out = Path(out_value)
    progressive_path = write_progressive_runtime_audit_v1003(out)
    progressive = progressive_runtime_payload_v1003()

    if not progressive.get("full_progressive_trial_loop_active", False):
        raise RuntimeError(
            "v10.0.3 run returned without an active audited trial lifecycle"
        )
    if len(engines) != 1:
        raise RuntimeError(
            f"v10.0.3 expected exactly one campaign engine; found {len(engines)}"
        )

    temperatures = [
        float(row.get("T_K", row.get("T", 0.0))) for row in (results or [])
    ]
    if not temperatures or any(T <= 0.0 for T in temperatures):
        raise RuntimeError(
            f"v10.0.3 could not recover valid temperatures from results: {temperatures}"
        )

    checks = []
    for T in temperatures:
        final = _last_step(out, T)
        B_final = _as_float(final, "B")
        extension = _as_float(final, "crack_extension_m")
        max_N_em = _max_column(out, T, "N_em")
        budget = progressive.get("source_budget_total")
        if extension is None or max_N_em is None or budget is None:
            raise RuntimeError(
                f"v10.0.3 missing source/extension diagnostics at {T:g} K"
            )
        population_bound = source_population_bound(
            float(budget),
            extension,
            float(manifest.source_refresh_length_m),
        )
        tolerance = max(1.0e-8, 1.0e-10 * max(population_bound, 1.0))
        if max_N_em > population_bound + tolerance:
            raise RuntimeError(
                f"v10.0.3 finite-source violation at {T:g} K: "
                f"max N_em={max_N_em:.16g} > refresh-aware bound="
                f"{population_bound:.16g}"
            )
        checks.append({
            "T_K": T,
            "B_final": B_final,
            "crack_extension_m": extension,
            "max_N_em": max_N_em,
            "source_budget_total": float(budget),
            "source_refresh_length_m": float(manifest.source_refresh_length_m),
            "source_population_bound": population_bound,
        })

    for row in results or []:
        row.update({
            "model": MODEL_ID,
            "point_release": "10.0.3",
            "front_state_model": "kinetic_campaign_czm",
            "front_state_model_detail": (
                "pf_v10_1_7_1_campaign_calibrated_continuous_tip_reset_safe_v1003"
            ),
            "progressive_runtime_audit": progressive_path.name,
            "live_binding_capture_verified": True,
            "delayed_transform_active": True,
        })
        T = float(row.get("T_K", row.get("T", 0.0)))
        match = next((x for x in checks if x["T_K"] == T), None)
        if match is not None:
            row["B_final"] = match["B_final"]
            row["crack_extension_final_m"] = match["crack_extension_m"]
            row["max_N_em"] = match["max_N_em"]
            row["source_budget_total"] = match["source_budget_total"]
            row["source_refresh_length_m"] = match["source_refresh_length_m"]
            row["source_population_bound"] = match["source_population_bound"]

    payload = {
        "model": MODEL_ID,
        "point_release": "10.0.3",
        "front_state_model": "kinetic_campaign_czm",
        "material_parameter_source": PF_SOURCE,
        "material": manifest.as_dict(),
        "kinetic_config": vars(kinetic_cfg),
        "opening_coupling": "clock_linear",
        "full_progressive_trial_loop_active": True,
        "delayed_transform_active": True,
        "live_binding_capture_verified": True,
        "campaign_dispatch_active": True,
        "anisotropic_elasticity_active": True,
        "anisotropic_J_active": True,
        "anisotropic_path_selection_active": False,
        "straight_single_front_mode_I_checkpoint": True,
        "stress_channels_separated": True,
        "wake_shielding_active": False,
        "stored_energy_cleavage_active": False,
        "artificial_sigma_cap_active": False,
        "artificial_emission_cap_active": False,
        "artificial_N_sat_active": False,
        "progressive_runtime_audit": progressive_path.name,
        "runtime": progressive,
        "result_checks": checks,
        "engine_audits": [eng.audit_payload() for eng in engines],
        "long_progressive_runs_authorized": False,
        "penalty_convergence_authorized": False,
    }
    (out / "kinetic_campaign_czm_v10_0_3_audit.json").write_text(
        json.dumps(payload, indent=2, default=str)
    )
    (out / "mode_i_v10_0_3_results.json").write_text(
        json.dumps(results, indent=2, default=str)
    )
    print("V10.0.3 DELAYED LIVE-BINDING PROGRESSIVE INTEGRATION COMPLETE")
    return results


if __name__ == "__main__":
    main()
