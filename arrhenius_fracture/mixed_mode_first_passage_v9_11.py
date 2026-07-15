"""Mixed-mode 2-D FEM/CZM validation using the selected v9.10.2/v9.10.3 MPZ law.

Version 9.11 is an integration branch, not a re-fit. It preserves the v8
interaction between domain-integral K and dimensionless anisotropic direction
factors, replaces the front-local scalar closure with the independently shaped
MPZ state, and makes the 2-D bulk plasticity use the same independent Peierls and
Taylor barriers.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from . import mixed_mode_first_passage_v8 as v8
from .bulk_plasticity_v9102 import independent_bulk_pt_active
from .mpz_front_engine_v911 import MovingProcessZone2DFrontEngine
from .mpz_parameterization_v911 import (
    apply_exact_barrier_args,
    apply_pt_dislocation_config,
    build_mpz_config,
    compact_audit,
    load_selected_row,
)
from .process_zone_2d_v911 import sample_process_zone_profile

MODEL_ID = "FEM_CZM_mixed_mode_MPZ_v9_11_independent_PT_non_double_counted_shielding"


class CalibratedV911TipEngineMixin(v8.CalibratedTipEngineMixin):
    """Apply the most recent 2-D profile before each front hazard update."""

    def _mm_drives(self, fallback_K):
        profile = getattr(self._mm, "mpz_profile_2d", None)
        if hasattr(self, "set_2d_process_zone_profile"):
            self.set_2d_process_zone_profile(profile)
        return super()._mm_drives(fallback_K)


def parser():
    p = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
    p.add_argument("--mpz-material-manifest", required=True)
    p.add_argument("--mpz-material-class", required=True)
    p.add_argument("--mpz-length-um", type=float, default=100.0)
    p.add_argument("--mpz-n-bins", type=int, default=200)
    p.add_argument("--mpz-profile-sector-half-angle-deg", type=float, default=45.0)
    p.add_argument("--mpz-profile-damage-cutoff", type=float, default=0.85)
    p.add_argument("--mpz-profile-min-elements", type=int, default=8)
    p.add_argument("--mixity-loading-angle-deg", type=float, default=None)
    p.add_argument("--mixity-open-coeff", type=float, default=None)
    p.add_argument("--mixity-shear-coeff", type=float, default=None)
    p.add_argument("--target-traction-phase-deg", type=float, required=True)
    p.add_argument("--traction-shear-sign", type=float, default=1.0)
    p.add_argument("--traction-probe-radius-m", type=float, default=10e-6)
    p.add_argument("--traction-annulus-half-width", type=float, default=0.45)
    p.add_argument("--traction-sector-half-angle-deg", type=float, default=40.0)
    p.add_argument("--traction-damage-cutoff", type=float, default=0.85)
    p.add_argument("--reference-cleavage-shape", type=float, required=True)
    p.add_argument("--reference-slip-shape", type=float, default=0.0)
    p.add_argument("--shear-emission-weight", type=float, default=1.0)
    p.add_argument("--directional-factor-max", type=float, default=5.0)
    p.add_argument("--solver-seed", type=int, default=1)
    return p


def _set_bulk_pt_namespace_defaults(args: Any, row: dict[str, Any]) -> None:
    """Set fields consumed while ``sharp_front.run_2d`` builds its bulk config."""
    emit0 = max(float(row["emit_G00_eV"]), 1.0e-30)
    values = {
        "bulk_kinetics_model": "emission_derived_peierls_taylor_multihit",
        "peierls_energy_scale": float(row["peierls_H0_eV"]) / emit0,
        "peierls_entropy_scale": float(row["peierls_activation_entropy_kB"]),
        "peierls_stress_scale": 1.0,
        "nu0_peierls": 1.0e12,
        "taylor_energy_scale": float(row["taylor_H0_eV"]) / emit0,
        "taylor_entropy_scale": float(row["taylor_activation_entropy_kB"]),
        "taylor_stress_scale": 1.0,
        "nu0_taylor": 1.0e11,
        "pt_taylor_corr_rho_c": float(row["taylor_corr_rho_c_m2"]),
        "pt_taylor_renewal_time_s": 1.0,
        "pt_taylor_m_exponent": 1.0,
        "pt_taylor_m_scale": float(row["taylor_corr_scale"]),
        "pt_taylor_m_cap": math.inf,
        "pt_mobile_fraction": 0.01,
        "pt_mobile_saturation_density_m2": math.inf,
        "pt_mobile_density_floor_m2": 0.0,
        "pt_jump_fraction": 1.0,
        "pt_jump_length_min_m": 0.0,
        "pt_taylor_phi_max": math.inf,
        "front_state_model": "moving_pz",
        "pz_store_to_rho_scale": 0.0,
        "tip_source_rho_per_emit": 0.0,
        "exhaustion": False,
    }
    for name, value in values.items():
        setattr(args, name, value)


def _plasticity_capture_factory(original, context, row):
    def wrapped(ep_gp, rho_gp, sigma_gp, mat, T, dt, plast_model, disl_cfg, *a, **kw):
        apply_pt_dislocation_config(disl_cfg, row)
        result = original(
            ep_gp, rho_gp, sigma_gp, mat, T, dt, plast_model, disl_cfg, *a, **kw
        )
        try:
            context.bulk_rho_gp = np.asarray(result[1], float).copy()
        except Exception:
            context.bulk_rho_gp = None
        return result
    return wrapped


def _j_profile_wrapper_factory(original_compute, context, mpz_args):
    base = v8._j_wrapper_factory(original_compute, context)

    def wrapped(mesh, u, sigma_gp, psi_e_gp, d, crack_tip, crack_direction,
                mat, ell, cfg=None, crack_segments=None, exclude_radius=0.0):
        J, KJ, info = base(
            mesh, u, sigma_gp, psi_e_gp, d, crack_tip, crack_direction, mat, ell,
            cfg=cfg, crack_segments=crack_segments,
            exclude_radius=exclude_radius,
        )
        try:
            profile = sample_process_zone_profile(
                mesh,
                sigma_gp,
                getattr(context, "bulk_rho_gp", None),
                d,
                crack_tip,
                crack_direction,
                length_m=float(mpz_args.mpz_length_um) * 1.0e-6,
                n_bins=int(mpz_args.mpz_n_bins),
                sector_half_angle_deg=float(mpz_args.mpz_profile_sector_half_angle_deg),
                damage_cutoff=float(mpz_args.mpz_profile_damage_cutoff),
                min_elements=int(mpz_args.mpz_profile_min_elements),
                poisson=float(getattr(mat, "nu", 0.28)),
            )
            context.mpz_profile_2d = profile
            diag = profile.diagnostics()
        except Exception as exc:
            context.mpz_profile_2d = None
            diag = {
                "mpz_2d_profile_reliable": False,
                "mpz_2d_profile_reason": f"profile_exception:{type(exc).__name__}",
                "bulk_scalar_rho_used_for_signed_shielding": False,
            }
        info.update(diag)
        context.latest.update(diag)
        if context.records:
            context.records[-1].update(diag)
        return J, KJ, info
    return wrapped


def _engine_factory(original_build, context, args, row):
    def build(parsed_args, material):
        base = original_build(parsed_args, material)
        cfg = build_mpz_config(args, row)
        base.f.r0 = 1.0e-6
        base.f.L_pz = cfg.length_m
        base.f.c_blunt = float(row["c_blunt"])
        base.f.nu0_c = 1.0e12
        base.f.nu0_e = 1.0e11
        base.f.m_hits = 3.0
        base.f.tau_c = 1.0e-6
        base.f.sigma_cap = 0.0
        base.f.dN_cap = math.inf
        base.f.N_sat = math.inf
        base.f.recover_k = 0.0
        base.f.k_shield = 0.0
        base.f.chi_shield = 0.0
        base.f.v_emb_b3 = 0.0

        Engine = type(
            "CalibratedMixedModeMPZV911Engine",
            (CalibratedV911TipEngineMixin, MovingProcessZone2DFrontEngine),
            {},
        )
        eng = Engine(base.f, base.cb, base.eb, base.G, base.nu, base.b, cfg)
        eng._mm_init(context)
        return eng
    return build


def _write_integration_audit(out: Path, row: dict[str, Any]) -> None:
    payload = {
        "model": MODEL_ID,
        **compact_audit(row),
        "bulk_PT": {
            "form": "independent_v9.10.2_EXP_floor_Peierls_then_Taylor",
            "Peierls_attempt_frequency_s-1": 1.0e12,
            "Taylor_attempt_frequency_s-1": 1.0e11,
            "Taylor_order_cap_active": False,
            "mobile_density_saturation_active": False,
            "minimum_jump_length_active": False,
            "plastic_rate_cap_active": False,
            "scalar_density_evolution": "existing_2D_storage_recovery_law_not_part_of_v9.10.2_fit",
            "scalar_density_role": "forest_density_for_Taylor_and_bulk_carrier_fallback",
        },
        "process_zone_coupling": {
            "absolute_tip_stress": "calibrated_K_drive/sqrt(2*pi*r_eff)",
            "FEM_finite_radius_stress_role": "dimensionless_spatial_and_directional_shape_only",
            "bulk_scalar_rho_role": "Taylor_forest_density_baseline_only",
            "bulk_scalar_rho_signed_shielding": False,
            "bulk_plastic_shielding": "already_in_domain_integral_J_via_FEM_redistribution",
            "unresolved_MPZ_shielding": "retained_line_K_integral_subtracted_once",
            "explicit_GND_backstress": False,
            "reason_GND_backstress_inactive": "2-D state is scalar rho, not signed slip-system GND/Nye tensor",
        },
    }
    (out / "mpz_v9_11_integration_audit.json").write_text(
        json.dumps(payload, indent=2, default=str)
    )


def main(argv=None):
    from . import fem as femmod
    from . import j_integral as jimod
    from . import plasticity as plastmod
    from . import sharp_front as sf

    mm, remaining = parser().parse_known_args(argv)
    row = load_selected_row(mm.mpz_material_manifest, mm.mpz_material_class)
    args = sf._build_parser().parse_args(remaining)
    apply_exact_barrier_args(args, row)
    _set_bulk_pt_namespace_defaults(args, row)

    if args.mode != "2d":
        raise SystemExit("v9.11 MPZ validation requires --mode 2d")
    if not bool(getattr(args, "crystal_aniso", False)):
        raise SystemExit("v9.11 requires --crystal-aniso")
    if not bool(getattr(args, "crystal_compete", False)):
        raise SystemExit("v9.11 requires --crystal-compete")
    if bool(getattr(args, "crystal_branch", False)) or int(getattr(args, "max_fronts", 1)) != 1:
        raise SystemExit("v9.11 first validation requires branching off and --max-fronts 1")

    if mm.mixity_open_coeff is not None or mm.mixity_shear_coeff is not None:
        if mm.mixity_open_coeff is None or mm.mixity_shear_coeff is None:
            raise SystemExit("provide both --mixity-open-coeff and --mixity-shear-coeff")
        qo, qs = v8.normalize_loading_coefficients(mm.mixity_open_coeff, mm.mixity_shear_coeff)
    else:
        alpha = float(0.0 if mm.mixity_loading_angle_deg is None else mm.mixity_loading_angle_deg)
        qo, qs = math.cos(math.radians(alpha)), math.sin(math.radians(alpha))

    context = v8.ProductionBackendControlContext(
        qo, qs, mm.target_traction_phase_deg,
        float(getattr(args, "crystal_theta_deg", 45.0) or 45.0),
        float(0.3 if getattr(args, "cleave_gamma_aniso", None) is None else getattr(args, "cleave_gamma_aniso")),
        mm.traction_probe_radius_m, mm.traction_annulus_half_width,
        mm.traction_sector_half_angle_deg, mm.traction_damage_cutoff,
        mm.traction_shear_sign, mm.reference_cleavage_shape,
        mm.reference_slip_shape, mm.shear_emission_weight,
        mm.directional_factor_max, mm.solver_seed,
    )
    context.bulk_rho_gp = None
    context.mpz_profile_2d = None
    context.material_row = row

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "anisotropic_calibrated_tip_run_config.json").write_text(json.dumps({
        "model": MODEL_ID,
        **vars(mm),
        **compact_audit(row),
        "crystal_theta_deg": context.crystal_theta_deg,
        "cleavage_gamma_aniso": context.cleavage_gamma_aniso,
        "note": "v8 exact mixed-mode control + v9.10.2 independent PT + non-double-counted 2-D MPZ coupling",
    }, indent=2, default=str))
    _write_integration_audit(out, row)

    osolve = femmod.solve_dirichlet
    oJ = jimod.compute_J_integral
    obuild = sf.build_engine
    oplast = plastmod.update_plasticity
    try:
        femmod.solve_dirichlet = v8._mixed_solve_factory(osolve, context)
        jimod.compute_J_integral = _j_profile_wrapper_factory(oJ, context, mm)
        plastmod.update_plasticity = _plasticity_capture_factory(oplast, context, row)
        sf.build_engine = _engine_factory(obuild, context, mm, row)
        with independent_bulk_pt_active():
            base = sf.run_2d(args)
    finally:
        femmod.solve_dirichlet = osolve
        jimod.compute_J_integral = oJ
        plastmod.update_plasticity = oplast
        sf.build_engine = obuild

    v8._write_records(out, context)
    vals = [v8._summary(out, T, context, base) for T in args.temperatures]
    for payload in vals:
        payload.update({
            "model": MODEL_ID,
            **compact_audit(row),
            "front_state_model_detail": "moving_pz_v911_independent_shapes_2d_profile",
            "bulk_pt_model_v9102_active": True,
            "bulk_scalar_rho_used_for_signed_shielding": False,
            "explicit_GND_backstress_active": False,
        })
        summary_json = out / "anisotropic_calibrated_tip_first_passage_summary.json"
        summary_json.write_text(json.dumps(payload, indent=2, default=str))
        try:
            import csv
            with (out / "anisotropic_calibrated_tip_first_passage_summary.csv").open("w", newline="") as fp:
                writer = csv.DictWriter(fp, fieldnames=list(payload))
                writer.writeheader()
                writer.writerow(payload)
        except Exception:
            pass
    print("MIXED_MODE_MPZ_V9_11 complete:", json.dumps(vals, indent=2, default=str))
    return vals


if __name__ == "__main__":
    main()
