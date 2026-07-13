"""Anisotropic mixed-mode first passage with the calibrated sharp-tip magnitude.

V8 preserves the v5 calibrated-tip correction and adds exact production-backend loading control.  The existing sharp-front/CZM engine
continues to convert the domain-integral KJ amplitude into the calibrated
process-zone stress

    sigma_tip = Kdrive / sqrt(2*pi*r_eff)

including the existing blunting, shielding, and barrier parameters.  Cubic FEM
stresses sampled at a finite physical radius are used only to construct
DIMENSIONLESS directional factors for cleavage and emission, to select the
candidate crack direction, and to audit the achieved local mode mixity.

No absolute finite-radius FEM traction is inserted into a barrier calibrated to
the sharp-tip stress scale.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

MODEL_ID = "FEM_CZM_mixed_mode_first_passage_v8_exact_backend_full_circle_boolean_safe"


def _unit(v):
    v = np.asarray(v, float)
    n = float(np.linalg.norm(v))
    return np.array([1.0, 0.0]) if n < 1e-30 else v / n


def angle_error_deg(value, target):
    return (float(value) - float(target) + 180.0) % 360.0 - 180.0


def traction_phase_deg(sigma_nn, tau_tn, shear_sign=1.0):
    return math.degrees(math.atan2(float(shear_sign) * float(tau_tn), float(sigma_nn)))


def shear_sign_from_basis(response_raw, min_diagonal=1e-12):
    M = np.asarray(response_raw, float)
    if M.shape != (2, 2) or not np.all(np.isfinite(M)):
        raise ValueError("response_raw must be a finite 2x2 matrix")
    scale = max(float(np.max(np.abs(M))), 1.0)
    if M[0, 0] <= min_diagonal * scale:
        raise ValueError("opening basis does not produce positive opening traction")
    if abs(M[1, 1]) <= min_diagonal * scale:
        raise ValueError("sliding basis has negligible shear response")
    return 1.0 if M[1, 1] >= 0 else -1.0


def loading_angle_from_response_basis(response_matrix, target_phase_deg, max_abs_alpha_deg=89.9):
    M = np.asarray(response_matrix, float)
    if M.shape != (2, 2) or not np.all(np.isfinite(M)):
        raise ValueError("response_matrix must be a finite 2x2 matrix")
    cond = float(np.linalg.cond(M))
    if not np.isfinite(cond) or cond > 1e12:
        raise ValueError(f"response matrix is singular or ill-conditioned: cond={cond:.6g}")
    p = math.radians(float(target_phase_deg))
    q = np.linalg.solve(M, np.array([math.cos(p), math.sin(p)]))
    if not np.all(np.isfinite(q)) or np.linalg.norm(q) < 1e-30:
        raise ValueError("invalid boundary vector")
    alpha = math.degrees(math.atan2(float(q[1]), float(q[0])))
    if abs(alpha) > float(max_abs_alpha_deg):
        raise ValueError(f"required loading angle {alpha:.6g} exceeds {max_abs_alpha_deg:.6g}")
    return float(alpha)


def phase_derivative_deg_per_deg(response_matrix, alpha_deg):
    """Derivative d psi / d alpha for the linear traction response matrix."""
    M = np.asarray(response_matrix, float)
    a = math.radians(float(alpha_deg))
    q = np.array([math.cos(a), math.sin(a)])
    dq = np.array([-math.sin(a), math.cos(a)])
    r = M @ q
    dr = M @ dq
    den = float(r @ r)
    if den <= 1e-30:
        return float("nan")
    return float((r[0] * dr[1] - r[1] * dr[0]) / den)


def safeguarded_alpha_update(alpha_deg, achieved_phase_deg, target_phase_deg,
                             response_matrix, max_step_deg=12.0, max_abs_alpha_deg=89.9):
    err = angle_error_deg(achieved_phase_deg, target_phase_deg)
    deriv = phase_derivative_deg_per_deg(response_matrix, alpha_deg)
    if not np.isfinite(deriv) or abs(deriv) < 1e-4:
        step = -math.copysign(min(abs(err), float(max_step_deg)), err)
    else:
        step = float(np.clip(-err / deriv, -abs(max_step_deg), abs(max_step_deg)))
    return float(np.clip(float(alpha_deg) + step, -abs(max_abs_alpha_deg), abs(max_abs_alpha_deg)))


def normalize_loading_coefficients(open_coeff, shear_coeff):
    q = np.asarray([open_coeff, shear_coeff], float)
    if q.shape != (2,) or not np.all(np.isfinite(q)):
        raise ValueError("loading coefficients must be finite")
    n = float(np.linalg.norm(q))
    if n <= 1e-30:
        raise ValueError("loading coefficient vector is zero")
    q /= n
    return float(q[0]), float(q[1])


def loading_coefficients_from_response_basis(response_matrix, target_phase_deg):
    """Return the full-circle normalized boundary vector for a target phase.

    Unlike v6, the sign of the opening coefficient is retained.  Some mixed-mode
    targets in an anisotropic geometry require a small compressive remote opening
    component while the local crack-tip normal traction remains tensile.
    """
    M = np.asarray(response_matrix, float)
    if M.shape != (2, 2) or not np.all(np.isfinite(M)):
        raise ValueError("response_matrix must be a finite 2x2 matrix")
    cond = float(np.linalg.cond(M))
    if not np.isfinite(cond) or cond > 1e12:
        raise ValueError("response matrix is singular or ill-conditioned")
    p = math.radians(float(target_phase_deg))
    q = np.linalg.solve(M, np.array([math.cos(p), math.sin(p)]))
    return normalize_loading_coefficients(q[0], q[1])


def wrap_loading_angle_deg(alpha_deg):
    return (float(alpha_deg) + 180.0) % 360.0 - 180.0


def loading_coefficients_from_alpha_deg(alpha_deg):
    a = math.radians(float(alpha_deg))
    return float(math.cos(a)), float(math.sin(a))


def loading_alpha_deg_from_coefficients(open_coeff, shear_coeff):
    qo, qs = normalize_loading_coefficients(open_coeff, shear_coeff)
    return wrap_loading_angle_deg(math.degrees(math.atan2(qs, qo)))


def phase_from_response_alpha(response_matrix, alpha_deg, shear_sign=1.0):
    qo, qs = loading_coefficients_from_alpha_deg(alpha_deg)
    r = np.asarray(response_matrix, float) @ np.array([qo, qs])
    return traction_phase_deg(r[0], r[1], shear_sign)


def safeguarded_event_alpha_update(samples, target_phase_deg, response_matrix,
                                   max_step_deg=20.0):
    """Empirical event-state update on the full loading circle.

    Samples carry an *unwrapped* `loading_alpha_unwrapped_deg` coordinate.  The
    boundary coefficients are periodic, but retaining an unwrapped coordinate
    makes local secant and bracket operations well defined across +/-180 deg.
    """
    valid = []
    for row in samples:
        a = float(row.get("loading_alpha_unwrapped_deg", np.nan))
        ph = float(row.get("achieved_psi_deg", np.nan))
        if np.isfinite(a) and np.isfinite(ph):
            valid.append((a, angle_error_deg(ph, target_phase_deg)))
    if not valid:
        raise ValueError("no finite event-state phase samples")

    sv = sorted(valid)
    brackets = []
    for (a0, e0), (a1, e1) in zip(sv[:-1], sv[1:]):
        if e0 == 0.0:
            return float(a0)
        if e1 == 0.0:
            return float(a1)
        if e0 * e1 < 0.0:
            brackets.append((abs(a1-a0), a0, e0, a1, e1))
    if brackets:
        _, a0, e0, a1, e1 = min(brackets)
        anew = a1 - e1*(a1-a0)/(e1-e0)
        if not (min(a0, a1) < anew < max(a0, a1)):
            anew = 0.5*(a0+a1)
        return float(anew)

    if len(valid) >= 2:
        a0, e0 = valid[-2]
        a1, e1 = valid[-1]
        da = a1-a0
        de = e1-e0
        if abs(da) > 1e-10 and abs(de/da) > 1e-4:
            step = float(np.clip(-e1/(de/da), -abs(max_step_deg), abs(max_step_deg)))
            return float(a1+step)

    a1, e1 = valid[-1]
    h = 0.05
    pp = phase_from_response_alpha(response_matrix, a1+h)
    pm = phase_from_response_alpha(response_matrix, a1-h)
    deriv = angle_error_deg(pp, pm)/(2*h)
    if not np.isfinite(deriv) or abs(deriv) < 1e-4:
        step = -math.copysign(min(abs(max_step_deg), max(abs(e1), 1.0)), e1)
    else:
        step = float(np.clip(-e1/deriv, -abs(max_step_deg), abs(max_step_deg)))
    return float(a1+step)


def energy_matrix_from_basis(J_open, J_slide, J_equal, U_cal=1.0):
    u2 = max(float(U_cal) ** 2, 1e-300)
    g11 = float(J_open) / u2
    g22 = float(J_slide) / u2
    g12 = float(J_equal) / u2 - 0.5 * (g11 + g22)
    return np.array([[g11, g12], [g12, g22]], float)


def _element_centroids(mesh):
    return np.asarray(mesh.nodes, float)[np.asarray(mesh.elems, int)].mean(axis=1)


def process_zone_traction_probe(mesh, sigma_gp, d, tip, crack_direction, radius_m,
                                annulus_half_width=0.45, sector_half_angle_deg=40.0,
                                damage_cutoff=0.85, min_elements=4):
    tip = np.asarray(tip, float)
    t = _unit(crack_direction)
    n = np.array([-t[1], t[0]])
    c = _element_centroids(mesh)
    dx = c - tip
    x1 = dx @ t
    x2 = dx @ n
    rr = np.sqrt(x1*x1 + x2*x2)
    th = np.degrees(np.arctan2(x2, x1))
    r0 = max(float(radius_m), 1e-30)
    hw = float(np.clip(annulus_half_width, 0.05, 0.95))
    dg = np.asarray(d, float)[np.asarray(mesh.elems, int)].mean(axis=1)
    expansion = 1.0
    for expansion in (1.0, 1.5, 2.0, 3.0):
        sel = ((rr >= max((1-hw)*r0/expansion, 0.25*r0)) &
               (rr <= (1+hw)*r0*expansion) &
               (np.abs(th) <= min(float(sector_half_angle_deg)*expansion, 85.0)) &
               (dg < float(damage_cutoff)))
        if int(np.count_nonzero(sel)) >= int(min_elements):
            break
    idx = np.flatnonzero(sel)
    if len(idx) < int(min_elements):
        return {"reliable": False, "n_elements": int(len(idx)),
                "probe_radius_m": r0, "expansion": float(expansion)}
    area = np.maximum(np.asarray(mesh.area_e, float)[idx], 1e-30)
    w = area / area.sum()
    sig = np.asarray(sigma_gp, float)
    sxx = float(w @ sig[0, idx])
    syy = float(w @ sig[1, idx])
    sxy = float(w @ sig[2, idx])
    S = np.array([[sxx, sxy], [sxy, syy]])
    stt = float(t @ S @ t)
    snn = float(n @ S @ n)
    ttn = float(t @ S @ n)
    eig = np.linalg.eigvalsh(S)
    return {"reliable": True, "n_elements": int(len(idx)),
            "probe_radius_m": r0, "expansion": float(expansion),
            "sigma_tt_Pa": stt, "sigma_nn_Pa": snn, "tau_tn_Pa": ttn,
            "sigma1_Pa": float(eig[-1]), "sigma2_Pa": float(eig[0]),
            "sigma_xx_Pa": sxx, "sigma_yy_Pa": syy, "sigma_xy_Pa": sxy,
            "stress_tensor": S}


def directional_shape_metrics(S, crystal_theta_deg, gamma_aniso, forward,
                              min_forward=0.2):
    """Dimensionless directional measures from the local anisotropic stress shape."""
    from .crystal import cleave_direction_competition, bcc_slip_traces

    S = np.asarray(S, float)
    eig = np.linalg.eigvalsh(S)
    sigma1 = max(float(eig[-1]), 0.0)
    if not np.isfinite(sigma1) or sigma1 <= 1e-30:
        return {"reliable": False, "sigma1_Pa": sigma1,
                "cleavage_shape": 0.0, "slip_shape": 0.0}

    selected, all_candidates = cleave_direction_competition(
        S, float(crystal_theta_deg), _unit(forward), min_forward=float(min_forward),
        gamma_aniso=float(gamma_aniso), branch_ratio=0.9)
    if not selected:
        return {"reliable": False, "sigma1_Pa": sigma1,
                "cleavage_shape": 0.0, "slip_shape": 0.0}
    winner = selected[0]

    best_slip = None
    for p in bcc_slip_traces(float(crystal_theta_deg)):
        tau = float(np.asarray(p["t"]) @ S @ np.asarray(p["n"]))
        row = {"name": p["name"], "tau_Pa": tau, "abs_tau_Pa": abs(tau),
               "angle_deg": float(p["angle_deg"])}
        if best_slip is None or row["abs_tau_Pa"] > best_slip["abs_tau_Pa"]:
            best_slip = row
    best_slip = best_slip or {"name": "none", "tau_Pa": 0.0,
                              "abs_tau_Pa": 0.0, "angle_deg": float("nan")}

    return {
        "reliable": True,
        "sigma1_Pa": sigma1,
        "cleavage_shape": float(winner["overdrive"]) / sigma1,
        "slip_shape": float(best_slip["abs_tau_Pa"]) / sigma1,
        "candidate_angle_deg": float(winner["angle_deg"]),
        "candidate_sigma_nn_Pa": float(winner["sigma_nn"]),
        "candidate_overdrive_Pa": float(winner["overdrive"]),
        "candidate_gamma_rel": float(winner["gamma"]),
        "candidate_t": np.asarray(winner["t"], float),
        "candidate_n": np.asarray(winner["n"], float),
        "slip_system_name": best_slip["name"],
        "slip_tau_signed_Pa": float(best_slip["tau_Pa"]),
        "slip_tau_abs_Pa": float(best_slip["abs_tau_Pa"]),
        "n_directional_candidates": int(len(all_candidates)),
    }


def directional_drive_factors(cleavage_shape, slip_shape,
                              reference_cleavage_shape, reference_slip_shape,
                              shear_emission_weight=1.0, factor_max=5.0):
    """Map anisotropic stress SHAPE to scalar K-drive multipliers.

    The cleavage multiplier is normalized to one in the calibrated Mode-I
    reference state.  Emission retains the same Mode-I baseline and receives an
    additional shear contribution only above the reference-state shear content.
    """
    rc = float(reference_cleavage_shape)
    rs = max(float(reference_slip_shape), 0.0)
    if not np.isfinite(rc) or rc <= 1e-12:
        raise ValueError("reference_cleavage_shape must be positive and finite")
    cs = max(float(cleavage_shape), 0.0)
    ss = max(float(slip_shape), 0.0)
    fc_raw = cs / rc
    shear_excess2 = max(ss*ss - rs*rs, 0.0)
    fe_raw = math.sqrt(max(fc_raw, 0.0)**2 +
                       max(float(shear_emission_weight), 0.0) * shear_excess2)
    cap = max(float(factor_max), 1e-12)
    fc = float(np.clip(fc_raw, 0.0, cap))
    fe = float(np.clip(fe_raw, 0.0, cap))
    return {
        "cleavage_factor": fc,
        "emission_factor": fe,
        "cleavage_factor_raw": float(fc_raw),
        "emission_factor_raw": float(fe_raw),
        "shear_excess_shape": math.sqrt(shear_excess2),
        "directional_factor_cap_active": bool(fc != fc_raw or fe != fe_raw),
    }


@dataclass
class ProductionBackendControlContext:
    loading_open_coeff: float
    loading_shear_coeff: float
    target_traction_phase_deg: float
    crystal_theta_deg: float = 45.0
    cleavage_gamma_aniso: float = 0.3
    probe_radius_m: float = 10e-6
    annulus_half_width: float = 0.45
    sector_half_angle_deg: float = 40.0
    damage_cutoff: float = 0.85
    shear_sign: float = 1.0
    reference_cleavage_shape: float = 1.0
    reference_slip_shape: float = 0.0
    shear_emission_weight: float = 1.0
    directional_factor_max: float = 5.0
    solver_seed: int = 1
    records: list[dict[str, Any]] = field(default_factory=list)
    latest: dict[str, Any] = field(default_factory=dict)
    last_loading: dict[str, float] = field(default_factory=dict)

    def __post_init__(self):
        qo, qs = normalize_loading_coefficients(self.loading_open_coeff, self.loading_shear_coeff)
        self.loading_open_coeff = qo
        self.loading_shear_coeff = qs

    @property
    def loading_angle_deg(self):
        return math.degrees(math.atan2(self.loading_shear_coeff, self.loading_open_coeff))


class CalibratedTipEngineMixin:
    """Preserve the original tip-stress scale and apply anisotropic K factors."""

    def _mm_init(self, context):
        self._mm = context
        self._mm_threshold = 1.0
        self._mm_prev_Kcleave = None

    def _mm_drives(self, fallback_K):
        KJ = max(float(fallback_K), 0.0)
        r = self._mm.latest or {}
        fc = float(r.get("cleavage_factor", 1.0))
        fe = float(r.get("emission_factor", 1.0))
        return KJ*fc, KJ*fe, {"KJ": KJ, "fc": fc, "fe": fe, **r}

    def _sigma_from_drive(self, Kdrive):
        K_eff = max(float(Kdrive) - self.f.k_shield*self.N_em*self.G*self.b /
                    np.sqrt(2.0*np.pi*self.f.L_pz), 0.0)
        s = K_eff / np.sqrt(2.0*np.pi*self.r_eff())
        if self.f.sigma_cap > 0:
            s = min(s, self.f.sigma_cap)
        return float(s)

    def predict_clock_increment(self, K, T, dt):
        Kc, _Ke, _ = self._mm_drives(K)
        sig2 = self._sigma_from_drive(Kc)
        lam2, _, _ = self.lambda_cleave(sig2, T)
        if self._mm_prev_Kcleave is not None and self._mm_prev_Kcleave > 0 and Kc > 0:
            sig1 = self._sigma_from_drive(self._mm_prev_Kcleave)
            lam1, _, _ = self.lambda_cleave(sig1, T)
            lo, hi = sorted((max(lam1, 0.0), max(lam2, 0.0)))
            if lo <= 0:
                leff = 0.5*hi
            elif abs(hi-lo) <= 1e-12*max(hi, 1e-300):
                leff = hi
            else:
                leff = (hi-lo)/math.log(hi/lo)
        else:
            leff = max(lam2, 0.0)
        remaining = max(self._mm_threshold-self.B, 1e-12)
        return float(max(leff*float(dt)/remaining, 0.0))

    def step(self, K, T, dt):
        Kc, Ke, mm = self._mm_drives(K)
        sig_emit = self._sigma_from_drive(Ke)
        K_eff_audit = max(Ke - self.f.k_shield*self.N_em*self.G*self.b /
                          np.sqrt(2.0*np.pi*self.f.L_pz), 0.0)
        sigma_tip_uncapped = K_eff_audit/np.sqrt(2.0*np.pi*self.r_eff())
        sigma_cap_active = bool(self.f.sigma_cap > 0 and sigma_tip_uncapped > self.f.sigma_cap)

        lam_e, sig_em_eff, Ge = self.lambda_emit(sig_emit, T)
        prod_raw = lam_e*dt
        prod = min(prod_raw, self.f.dN_cap)
        dN_cap_active = bool(np.isfinite(self.f.dN_cap) and prod_raw > self.f.dN_cap)
        N_sat_factor = 1.0
        if np.isfinite(self.f.N_sat) and self.f.N_sat > 0:
            N_sat_factor = max(1.0-self.N_em/self.f.N_sat, 0.0)
            prod *= N_sat_factor
        ann = self.f.recover_k*self.N_em*dt
        self.W_emit += sig_em_eff*self.b*self.f.L_pz*prod
        self.N_em = max(self.N_em+prod-ann, 0.0)

        sig_cleave = self._sigma_from_drive(Kc)
        lam_c, lam_c_raw, Gc_eff = self.lambda_cleave(sig_cleave, T)
        if self.f.tau_B > 0 and dt > 0:
            self.B *= np.exp(-min(dt/self.f.tau_B, 80.0))
        if self._mm_prev_Kcleave is not None and self._mm_prev_Kcleave > 0 and Kc > 0:
            sig1 = self._sigma_from_drive(self._mm_prev_Kcleave)
            lam1, _, _ = self.lambda_cleave(sig1, T)
            lo, hi = sorted((max(lam1, 0.0), max(lam_c, 0.0)))
            if lo <= 0:
                leff = 0.5*hi
            elif abs(hi-lo) <= 1e-12*max(hi, 1e-300):
                leff = hi
            else:
                leff = (hi-lo)/math.log(hi/lo)
        else:
            leff = lam_c
        self.B += leff*dt
        self._mm_prev_Kcleave = Kc
        self.K_prev = Kc
        self.t += dt

        Npre = self.N_em
        spre = self.sigma_back()
        rpre = self.r_eff()
        dGpre = self.dG_emb()
        nfire = 0
        if not np.isfinite(self.B):
            self.B = 0.0
        while self.B >= self._mm_threshold and nfire < 100000:
            self.B -= self._mm_threshold
            nfire += 1
        fired = nfire > 0
        retained = Npre
        shed = 0.0
        if fired:
            retain = float(np.clip(self.f.wake_retain, 0, 1))**nfire
            retained = Npre*retain
            shed = Npre-retained
            self.N_em = retained
            self.a_adv += self.f.da*nfire
            self.n_adv += nfire

        info = {
            "fired": fired, "n_fire": nfire,
            "v_crack": self.f.da*nfire/dt if dt > 0 else 0.0,
            "sigma_tip": sig_cleave, "sigma_emit_tip": sig_emit,
            "sigma_back": self.sigma_back(), "lambda_e": lam_e,
            "lambda_c": lam_c, "lambda_c_raw": lam_c_raw, "B": self.B,
            "N_em": self.N_em, "r_eff": self.r_eff(),
            "dG_emb_eV": self.dG_emb()/1.602176634e-19,
            "G_cleave_eff_eV": Gc_eff/1.602176634e-19,
            **self.cleavage_diagnostics(sig_cleave, T),
            "G_emit_eV": Ge/1.602176634e-19, "W_emit": self.W_emit,
            "sigma_tip_uncapped": float(sigma_tip_uncapped),
            "sigma_cap_active": sigma_cap_active,
            "dN_emit_raw": float(prod_raw), "dN_cap_active": dN_cap_active,
            "N_sat_factor": float(N_sat_factor),
            "N_sat_active": bool(N_sat_factor < 0.999999),
            "N_em_pre_renewal": Npre, "N_em_retained": retained,
            "N_em_shed_to_wake": shed,
            "sigma_back_pre_renewal": spre, "r_eff_pre_renewal": rpre,
            "dG_emb_pre_renewal_eV": dGpre/1.602176634e-19,
            "anisotropic_KJ_Pa_sqrt_m": mm.get("KJ"),
            "anisotropic_Kcleave_Pa_sqrt_m": Kc,
            "anisotropic_Kemit_Pa_sqrt_m": Ke,
            "anisotropic_cleavage_factor": mm.get("fc"),
            "anisotropic_emission_factor": mm.get("fe"),
            "anisotropic_reference_phase_deg": mm.get("reference_traction_phase_deg"),
            "anisotropic_candidate_angle_deg": mm.get("candidate_angle_deg"),
            "anisotropic_candidate_sigma_nn_Pa": mm.get("candidate_sigma_nn_Pa"),
            "anisotropic_probe_sigma1_Pa": mm.get("probe_sigma1_Pa"),
            "anisotropic_probe_tau_slip_abs_Pa": mm.get("slip_tau_abs_Pa"),
            "anisotropic_probe_sigma_cleave_overdrive_Pa": mm.get("candidate_overdrive_Pa"),
            "anisotropic_slip_system": mm.get("slip_system_name"),
            "anisotropic_directional_factor_cap_active": mm.get("directional_factor_cap_active", False),
        }
        return info


def _mixed_solve_factory(original_solve, context):
    def solve_mixed(K, Rint, u, bnd, Uy_top, Uy_bot):
        from scipy.sparse.linalg import spsolve
        qo = float(context.loading_open_coeff)
        qs = float(context.loading_shear_coeff)
        alpha = math.atan2(qs, qo)
        Atotal = float(Uy_top-Uy_bot)
        Un = Atotal*qo
        Us = Atotal*qs
        u_open, _ = original_solve(K, Rint, u, bnd, 0.5*Un, -0.5*Un)
        Kc = K.tocsr()
        Ropen = Rint + Kc@(u_open-u)
        if abs(Us) <= 1e-30:
            unew = u_open
            Rfull = Ropen
        else:
            ndof = len(u)
            prescribed = np.zeros(ndof, bool)
            target = u_open.copy()
            tn, bn = bnd.top_nodes, bnd.bot_nodes
            prescribed[2*tn] = True
            prescribed[2*tn+1] = True
            prescribed[2*bn] = True
            prescribed[2*bn+1] = True
            target[2*tn] = u_open[2*tn] + 0.5*Us
            target[2*bn] = u_open[2*bn] - 0.5*Us
            target[2*tn+1] = 0.5*Un
            target[2*bn+1] = -0.5*Un
            free = ~prescribed
            dup = target[prescribed]-u_open[prescribed]
            rhs = -Ropen[free]-Kc[np.ix_(free, prescribed)]@dup
            unew = u_open.copy()
            unew[free] = u_open[free] + spsolve(Kc[np.ix_(free, free)], rhs)
            unew[prescribed] = target[prescribed]
            Rfull = Ropen + Kc@(unew-u_open)
        Fx = float(np.sum(Rfull[2*bnd.top_nodes]))
        Fy = float(np.sum(Rfull[2*bnd.top_nodes+1]))
        Fgen = Fx*math.sin(alpha)+Fy*math.cos(alpha)
        context.last_loading = {
            "U_total_m": Atotal, "U_open_m": Un, "U_shear_m": Us,
            "loading_angle_deg": context.loading_angle_deg,
            "loading_open_coeff": qo, "loading_shear_coeff": qs,
            "loading_alpha_deg": loading_alpha_deg_from_coefficients(qo, qs),
            "loading_open_is_tensile": bool(qo >= 0.0),
            "generalized_reaction_N": Fgen,
            "reaction_x_N": Fx, "reaction_y_N": Fy,
        }
        return unew, Fgen
    return solve_mixed


def _j_wrapper_factory(original_compute, context):
    def wrapped(mesh, u, sigma_gp, psi_e_gp, d, crack_tip, crack_direction,
                mat, ell, cfg=None, crack_segments=None, exclude_radius=0.0):
        J, KJ, info = original_compute(
            mesh, u, sigma_gp, psi_e_gp, d, crack_tip, crack_direction, mat, ell,
            cfg=cfg, crack_segments=crack_segments, exclude_radius=exclude_radius)
        probe = process_zone_traction_probe(
            mesh, sigma_gp, d, crack_tip, crack_direction,
            context.probe_radius_m, context.annulus_half_width,
            context.sector_half_angle_deg, context.damage_cutoff)
        phase_probe_reliable = bool(probe.get("reliable", False))
        if phase_probe_reliable:
            sn = float(probe.get("sigma_nn_Pa", np.nan))
            tt = float(probe.get("tau_tn_Pa", np.nan))
            phase_probe_reliable = bool(np.isfinite(sn) and np.isfinite(tt) and
                                        math.hypot(sn, tt) > 1e-12)
        if phase_probe_reliable:
            phase = traction_phase_deg(probe["sigma_nn_Pa"], probe["tau_tn_Pa"], context.shear_sign)
            metrics = directional_shape_metrics(
                probe["stress_tensor"], context.crystal_theta_deg,
                context.cleavage_gamma_aniso, crack_direction)
        else:
            phase = float("nan")
            metrics = {"reliable": False}
        directional_reliable = bool(metrics.get("reliable", False))
        if directional_reliable:
            factors = directional_drive_factors(
                metrics["cleavage_shape"], metrics["slip_shape"],
                context.reference_cleavage_shape, context.reference_slip_shape,
                context.shear_emission_weight, context.directional_factor_max)
        else:
            factors = {"cleavage_factor": 1.0, "emission_factor": 1.0,
                       "cleavage_factor_raw": float("nan"),
                       "emission_factor_raw": float("nan"),
                       "shear_excess_shape": float("nan"),
                       "directional_factor_cap_active": False}
        reliable = bool(phase_probe_reliable and directional_reliable)
        md = {
            "model": MODEL_ID,
            "J_J_per_m2": float(J),
            "KJ_reference_Pa_sqrt_m": float(KJ),
            "traction_phase_probe_reliable": phase_probe_reliable,
            "directional_metrics_reliable": directional_reliable,
            "traction_probe_reliable": reliable,
            "traction_probe_radius_m": context.probe_radius_m,
            "reference_sigma_nn_Pa": probe.get("sigma_nn_Pa"),
            "reference_tau_tn_Pa": probe.get("tau_tn_Pa"),
            "reference_traction_phase_deg": phase,
            "target_traction_phase_deg": context.target_traction_phase_deg,
            "reference_phase_error_deg": angle_error_deg(phase, context.target_traction_phase_deg) if np.isfinite(phase) else np.nan,
            "probe_sigma1_Pa": metrics.get("sigma1_Pa"),
            "cleavage_shape": metrics.get("cleavage_shape"),
            "slip_shape": metrics.get("slip_shape"),
            "candidate_angle_deg": metrics.get("candidate_angle_deg"),
            "candidate_sigma_nn_Pa": metrics.get("candidate_sigma_nn_Pa"),
            "candidate_overdrive_Pa": metrics.get("candidate_overdrive_Pa"),
            "candidate_gamma_rel": metrics.get("candidate_gamma_rel"),
            "slip_system_name": metrics.get("slip_system_name"),
            "slip_tau_signed_Pa": metrics.get("slip_tau_signed_Pa"),
            "slip_tau_abs_Pa": metrics.get("slip_tau_abs_Pa"),
            "reference_cleavage_shape": context.reference_cleavage_shape,
            "reference_slip_shape": context.reference_slip_shape,
            **factors,
            **context.last_loading,
            "tip_x_m": float(np.asarray(crack_tip)[0]),
            "tip_y_m": float(np.asarray(crack_tip)[1]),
        }
        info.update(md)
        context.latest = md.copy()
        context.records.append(md.copy())
        return J, KJ, info
    return wrapped


def _engine_factory(original_build, context, base_class):
    class Engine(CalibratedTipEngineMixin, base_class):
        def __init__(self, *a, **kw):
            base_class.__init__(self, *a, **kw)
            self._mm_init(context)
    def build(args, mat):
        base = original_build(args, mat)
        return Engine(base.f, base.cb, base.eb, base.G, base.nu, base.b)
    return build


def _write_records(out, context):
    if not context.records:
        return
    cols = sorted({k for r in context.records for k in r})
    with (out/"anisotropic_calibrated_tip_calls.csv").open("w", newline="") as fp:
        w = csv.DictWriter(fp, fieldnames=cols)
        w.writeheader()
        w.writerows(context.records)


def _summary(out, T, context, base_summary):
    steps = out/f"steps_{int(T):04d}K.csv"
    accepted = {}
    if steps.exists():
        a = np.genfromtxt(steps, delimiter=",", names=True)
        if np.size(a):
            a = np.atleast_1d(a)
            fire = np.flatnonzero(a["n_fire"] > 0)
            i = int(fire[0]) if len(fire) else len(a)-1
            accepted = {n: float(a[n][i]) for n in a.dtype.names}
    U = accepted.get("Uapp_m", np.nan)
    Kacc = accepted.get("KJ_Pa_sqrtm", np.nan)
    rec = context.records[-1] if context.records else {}
    if context.records and np.isfinite(U):
        du = np.array([abs(float(r.get("U_total_m", np.nan))-U) for r in context.records])
        if np.any(np.isfinite(du)):
            best = float(np.nanmin(du))
            tol = max(1e-14, 1e-10*max(abs(U), 1e-12))
            near = [r for r, q in zip(context.records, du) if np.isfinite(q) and q <= best+tol]
            if near and np.isfinite(Kacc):
                rec = min(near, key=lambda r: abs(float(r.get("KJ_reference_Pa_sqrt_m", np.nan))-Kacc)
                          if np.isfinite(float(r.get("KJ_reference_Pa_sqrt_m", np.nan))) else float("inf"))
            elif near:
                rec = near[-1]
    b = base_summary[0] if base_summary else {}
    KJ = float(rec.get("KJ_reference_Pa_sqrt_m", np.nan))
    fc = float(rec.get("cleavage_factor", np.nan))
    fe = float(rec.get("emission_factor", np.nan))
    payload = {
        "model": MODEL_ID,
        "T_K": float(T),
        "loading_angle_deg": context.loading_angle_deg,
        "loading_open_coeff": context.loading_open_coeff,
        "loading_shear_coeff": context.loading_shear_coeff,
        "loading_alpha_deg": loading_alpha_deg_from_coefficients(context.loading_open_coeff, context.loading_shear_coeff),
        "loading_open_is_tensile": bool(context.loading_open_coeff >= 0.0),
        "target_traction_phase_deg": context.target_traction_phase_deg,
        "traction_phase_first_deg": rec.get("reference_traction_phase_deg"),
        "traction_phase_error_first_deg": rec.get("reference_phase_error_deg"),
        "J_first_J_per_m2": rec.get("J_J_per_m2"),
        "KJ_reference_first_MPa_sqrt_m": KJ/1e6,
        "Kcleave_calibrated_first_MPa_sqrt_m": KJ*fc/1e6,
        "Kemit_calibrated_first_MPa_sqrt_m": KJ*fe/1e6,
        "cleavage_factor_first": fc,
        "emission_factor_first": fe,
        "cleavage_shape_first": rec.get("cleavage_shape"),
        "slip_shape_first": rec.get("slip_shape"),
        "reference_cleavage_shape": context.reference_cleavage_shape,
        "reference_slip_shape": context.reference_slip_shape,
        "reference_sigma_nn_first_GPa": float(rec.get("reference_sigma_nn_Pa", np.nan))/1e9,
        "reference_tau_tn_first_GPa": float(rec.get("reference_tau_tn_Pa", np.nan))/1e9,
        "probe_sigma1_first_GPa": float(rec.get("probe_sigma1_Pa", np.nan))/1e9,
        "candidate_sigma_nn_first_GPa": float(rec.get("candidate_sigma_nn_Pa", np.nan))/1e9,
        "candidate_overdrive_first_GPa": float(rec.get("candidate_overdrive_Pa", np.nan))/1e9,
        "slip_tau_abs_first_GPa": float(rec.get("slip_tau_abs_Pa", np.nan))/1e9,
        "candidate_angle_first_deg": rec.get("candidate_angle_deg"),
        "candidate_gamma_rel": rec.get("candidate_gamma_rel"),
        "slip_system_first": rec.get("slip_system_name"),
        "traction_phase_probe_reliable": rec.get("traction_phase_probe_reliable"),
        "directional_metrics_reliable": rec.get("directional_metrics_reliable"),
        "traction_probe_reliable": rec.get("traction_probe_reliable"),
        "directional_factor_cap_active": rec.get("directional_factor_cap_active"),
        "control_state": "first_passage" if accepted.get("n_fire", 0) > 0 else "right_censored_endpoint",
        "Kc_first_existing_MPa_sqrt_m": b.get("Kc_first_MPa_sqrt_m"),
        "N_em_final": b.get("N_em_final"),
        "B_final": b.get("B_final"),
        "mode_classification": b.get("mode"),
        "crystal_aniso": True,
        "crystal_compete": True,
        "crystal_theta_deg": context.crystal_theta_deg,
        "cleavage_gamma_aniso": context.cleavage_gamma_aniso,
        "absolute_stress_driver": "existing_sharp_front_sigma_tip_from_directionally_scaled_KJ",
        "finite_radius_traction_role": "dimensionless_direction_and_mixity_only",
        "event_phase_within_2deg": bool(np.isfinite(float(rec.get("reference_phase_error_deg", np.nan))) and
                                         abs(float(rec.get("reference_phase_error_deg", np.nan))) <= 2.0),
    }
    (out/"anisotropic_calibrated_tip_first_passage_summary.json").write_text(
        json.dumps(payload, indent=2, default=str))
    with (out/"anisotropic_calibrated_tip_first_passage_summary.csv").open("w", newline="") as fp:
        w = csv.DictWriter(fp, fieldnames=list(payload))
        w.writeheader()
        w.writerow(payload)
    return payload


def parser():
    p = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
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


def main(argv=None):
    from . import sharp_front as sf
    from . import fem as femmod
    from . import j_integral as jimod

    mm, remaining = parser().parse_known_args(argv)
    args = sf._build_parser().parse_args(remaining)
    if args.mode != "2d":
        raise SystemExit("v8 anisotropic mixed mode requires --mode 2d")
    if not bool(getattr(args, "crystal_aniso", False)):
        raise SystemExit("v8 requires --crystal-aniso")
    if not bool(getattr(args, "crystal_compete", False)):
        raise SystemExit("v8 requires --crystal-compete")
    if bool(getattr(args, "crystal_branch", False)) or int(getattr(args, "max_fronts", 1)) != 1:
        raise SystemExit("v8 first-passage screen requires branching off and --max-fronts 1")

    if mm.mixity_open_coeff is not None or mm.mixity_shear_coeff is not None:
        if mm.mixity_open_coeff is None or mm.mixity_shear_coeff is None:
            raise SystemExit("provide both --mixity-open-coeff and --mixity-shear-coeff")
        qo, qs = normalize_loading_coefficients(mm.mixity_open_coeff, mm.mixity_shear_coeff)
    else:
        alpha = float(0.0 if mm.mixity_loading_angle_deg is None else mm.mixity_loading_angle_deg)
        qo, qs = math.cos(math.radians(alpha)), math.sin(math.radians(alpha))
    context = ProductionBackendControlContext(
        qo, qs, mm.target_traction_phase_deg,
        float(getattr(args, "crystal_theta_deg", 45.0) or 45.0),
        float(0.3 if getattr(args, "cleave_gamma_aniso", None) is None else getattr(args, "cleave_gamma_aniso")),
        mm.traction_probe_radius_m, mm.traction_annulus_half_width,
        mm.traction_sector_half_angle_deg, mm.traction_damage_cutoff,
        mm.traction_shear_sign, mm.reference_cleavage_shape,
        mm.reference_slip_shape, mm.shear_emission_weight,
        mm.directional_factor_max, mm.solver_seed)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out/"anisotropic_calibrated_tip_run_config.json").write_text(json.dumps({
        "model": MODEL_ID, **vars(mm),
        "crystal_theta_deg": context.crystal_theta_deg,
        "cleavage_gamma_aniso": context.cleavage_gamma_aniso,
        "note": "exact production-backend loading coefficients + anisotropic direction partition + calibrated sharp-tip magnitude",
    }, indent=2))

    osolve = femmod.solve_dirichlet
    oJ = jimod.compute_J_integral
    obuild = sf.build_engine
    try:
        femmod.solve_dirichlet = _mixed_solve_factory(osolve, context)
        jimod.compute_J_integral = _j_wrapper_factory(oJ, context)
        sf.build_engine = _engine_factory(obuild, context, sf.FrontEngine)
        base = sf.run_2d(args)
    finally:
        femmod.solve_dirichlet = osolve
        jimod.compute_J_integral = oJ
        sf.build_engine = obuild
    _write_records(out, context)
    vals = [_summary(out, T, context, base) for T in args.temperatures]
    print("MIXED_MODE_V8 complete:", json.dumps(vals, indent=2, default=str))
    return vals


if __name__ == "__main__":
    main()
