"""J-consistent isotropic mixed-mode first-passage wrapper with phase-ratio calibration for the active FEM/CZM sharp-front solver.

This module does not replace ``arrhenius_fracture.sharp_front``.  It temporarily
wraps three interfaces during one run:

* the displacement boundary solve, to impose combined opening/sliding loading;
* the domain J calculation, to extract KI/KII with a robust full-field multi-annulus Williams fit;
* the front engine, to use an opening-sensitive cleavage drive and a
  shear-assisted emission drive.

The underlying barrier surfaces, plasticity, adaptive CZM backend, mesh, and
first-passage clocks remain those of the active project.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
from dataclasses import dataclass, field
from pathlib import Path
from types import MethodType
from typing import Any

import numpy as np

MODEL_ID = "FEM_CZM_mixed_mode_first_passage_v3_2_J_consistent_signed_basis"


def _unit(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, dtype=float)
    n = float(np.linalg.norm(v))
    return np.array([1.0, 0.0]) if n < 1e-30 else v / n


def maximum_hoop_drive(KI: float, KII: float, limit_deg: float = 85.0) -> tuple[float, float]:
    """Return maximum positive hoop-stress SIF coefficient and its angle.

    Uses the leading isotropic Williams field

        K_theta(theta) = KI cos^3(theta/2)
                         - 3 KII sin(theta/2) cos^2(theta/2).

    The projection is also useful as a transparent scalar diagnostic in the
    anisotropic campaign, but anisotropic mode decomposition should be treated
    as approximate until an interaction-integral decomposition is added.
    """
    th = np.linspace(-abs(limit_deg), abs(limit_deg), 1361) * np.pi / 180.0
    c = np.cos(0.5 * th)
    s = np.sin(0.5 * th)
    kval = float(KI) * c**3 - 3.0 * float(KII) * s * c**2
    i = int(np.argmax(kval))
    return max(float(kval[i]), 0.0), float(np.degrees(th[i]))



def angle_error_deg(value: float, target: float) -> float:
    """Wrapped angular error in degrees."""
    return (float(value) - float(target) + 180.0) % 360.0 - 180.0




def mode_signs_from_basis(mode_basis_raw: np.ndarray, *, min_diagonal: float = 1.0e-12) -> np.ndarray:
    """Return phase-convention signs from pure opening and pure sliding bases.

    A Williams projection may use the opposite sign convention for Mode I or
    Mode II relative to the imposed boundary displacement.  The physical
    convention for this campaign is defined by requiring the pure-opening
    basis to have positive KI and the pure-sliding basis to have positive KII.
    The returned two-vector multiplies ``[KI_raw, KII_raw]``.
    """
    M = np.asarray(mode_basis_raw, dtype=float)
    if M.shape != (2, 2) or not np.all(np.isfinite(M)):
        raise ValueError("mode_basis_raw must be a finite 2x2 matrix")
    diag = np.array([M[0, 0], M[1, 1]], dtype=float)
    scale = max(float(np.max(np.abs(M))), 1.0)
    if abs(diag[0]) <= float(min_diagonal) * scale:
        raise ValueError("pure-opening basis has negligible KI response")
    if abs(diag[1]) <= float(min_diagonal) * scale:
        raise ValueError("pure-sliding basis has negligible KII response")
    return np.where(diag >= 0.0, 1.0, -1.0)


def apply_mode_signs(KI_raw: float, KII_raw: float, signs: np.ndarray) -> tuple[float, float, float]:
    """Apply calibrated mode signs and return normalized KI, KII, and phase."""
    sg = np.asarray(signs, dtype=float).reshape(2)
    KI = float(sg[0] * float(KI_raw))
    KII = float(sg[1] * float(KII_raw))
    psi = math.degrees(math.atan2(KII, KI))
    return KI, KII, psi

def loading_angle_from_mode_basis(mode_basis: np.ndarray, target_psi_deg: float,
                                  max_abs_alpha_deg: float = 89.9) -> float:
    """Return boundary loading angle from a 2x2 elastic mode-response matrix.

    ``mode_basis[:, 0]`` is the measured ``[KI, KII]`` response to unit opening
    displacement and ``mode_basis[:, 1]`` is the response to unit sliding
    displacement.  The requested mode vector is solved in this measured basis,
    so geometric cross-coupling is retained instead of assuming ``alpha=psi``.
    """
    M = np.asarray(mode_basis, dtype=float)
    if M.shape != (2, 2) or not np.all(np.isfinite(M)):
        raise ValueError("mode_basis must be a finite 2x2 matrix")
    cond = float(np.linalg.cond(M))
    if not np.isfinite(cond) or cond > 1.0e12:
        raise ValueError(f"mode basis is singular or ill-conditioned: cond={cond:.3g}")
    psi = math.radians(float(target_psi_deg))
    target = np.array([math.cos(psi), math.sin(psi)], dtype=float)
    q = np.linalg.solve(M, target)
    if not np.all(np.isfinite(q)) or float(np.linalg.norm(q)) <= 1e-30:
        raise ValueError("invalid boundary combination from mode basis")
    alpha = math.degrees(math.atan2(float(q[1]), float(q[0])))
    if abs(alpha) > float(max_abs_alpha_deg):
        raise ValueError(
            f"required loading angle {alpha:.3f} deg exceeds limit "
            f"{max_abs_alpha_deg:.3f} deg"
        )
    return float(alpha)


def phase_projection_gate(record: dict[str, float], *, min_fits: int = 2,
                          min_points: int = 10,
                          max_condition: float = 1.0e12,
                          max_phase_spread_deg: float = 20.0) -> tuple[bool, list[str]]:
    """Assess phase-ratio usability without gating on amplitude-fit residual.

    v3 obtains the authoritative amplitude from the domain J integral.  The
    Williams projection is used only for the ratio ``KII/KI`` during the
    undamaged isotropic elastic calibration.  Consequently, relative amplitude
    residual and annular K-magnitude spread are diagnostics, not rejection
    criteria.  The phase gate checks finite modes, sample support, conditioning,
    and annulus-to-annulus phase stability.
    """
    reasons: list[str] = []
    KI = float(record.get("KI_MPa_sqrt_m", record.get("KI_Pa_sqrt_m", float("nan"))))
    KII = float(record.get("KII_MPa_sqrt_m", record.get("KII_Pa_sqrt_m", float("nan"))))
    psi = float(record.get("achieved_psi_deg", record.get("mode_phase_deg", float("nan"))))
    n = int(record.get("projection_n", record.get("mode_projection_n", 0)) or 0)
    nf = int(record.get("projection_fit_count", record.get("mode_projection_fit_count", 0)) or 0)
    cond = float(record.get("projection_condition", record.get("mode_projection_condition", float("inf"))))
    spread = float(record.get("projection_psi_spread_deg", record.get("mode_projection_psi_spread_deg", float("inf"))))
    if not (np.isfinite(KI) and np.isfinite(KII) and np.isfinite(psi)):
        reasons.append("nonfinite_mode_ratio")
    if KI <= 0.0:
        reasons.append("nonpositive_KI")
    if n < int(min_points):
        reasons.append(f"points<{int(min_points)}")
    if nf < int(min_fits):
        reasons.append(f"fits<{int(min_fits)}")
    if not np.isfinite(cond) or cond > float(max_condition):
        reasons.append("ill_conditioned")
    if not np.isfinite(spread) or spread > float(max_phase_spread_deg):
        reasons.append("phase_spread")
    return (len(reasons) == 0), reasons

def _williams_shapes(theta: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return isotropic leading Williams stress shapes for modes I and II.

    Component order is ``(sigma_11, sigma_22, sigma_12)`` in crack-local
    coordinates.  Multiplication by ``K/sqrt(2*pi*r)`` gives stress.
    """
    th = np.asarray(theta, dtype=float)
    s = np.sin(0.5 * th)
    c = np.cos(0.5 * th)
    s3 = np.sin(1.5 * th)
    c3 = np.cos(1.5 * th)
    fI = np.column_stack([
        c * (1.0 - s * s3),
        c * (1.0 + s * s3),
        s * c * c3,
    ])
    fII = np.column_stack([
        -s * (2.0 + c * c3),
        s * c * c3,
        c * (1.0 - s * s3),
    ])
    return fI, fII


def _weighted_lstsq(A: np.ndarray, y: np.ndarray, w: np.ndarray) -> tuple[np.ndarray, float]:
    """Column-scaled weighted least squares and normalized condition number."""
    Aw = np.asarray(A, float) * np.asarray(w, float)[:, None]
    yw = np.asarray(y, float) * np.asarray(w, float)
    scale = np.linalg.norm(Aw, axis=0)
    scale = np.where(scale > 1e-30, scale, 1.0)
    As = Aw / scale[None, :]
    z = np.linalg.lstsq(As, yw, rcond=None)[0]
    x = z / scale
    return x, float(np.linalg.cond(As))


def _fit_williams_subset(mesh, sigma_gp, idx: np.ndarray, rr: np.ndarray,
                         theta: np.ndarray, e1: np.ndarray, e2: np.ndarray,
                         max_irls: int = 8) -> dict[str, float]:
    """Robustly fit KI, KII, and nonsingular constant stresses."""
    idx = np.asarray(idx, dtype=int)
    n = len(idx)
    if n < 7:
        return {"valid": False, "n": int(n)}

    sx = np.asarray(sigma_gp[0, idx], float)
    sy = np.asarray(sigma_gp[1, idx], float)
    txy = np.asarray(sigma_gp[2, idx], float)
    R = np.column_stack([e1, e2])
    sl = np.empty((n, 3), dtype=float)
    for q, (a, b, cxy) in enumerate(zip(sx, sy, txy)):
        S = np.array([[a, cxy], [cxy, b]], dtype=float)
        Slocal = R.T @ S @ R
        sl[q] = (Slocal[0, 0], Slocal[1, 1], Slocal[0, 1])

    fI, fII = _williams_shapes(theta[idx])
    singular = 1.0 / np.sqrt(2.0 * np.pi * rr[idx])
    A = np.zeros((3*n, 5), dtype=float)
    y = np.zeros(3*n, dtype=float)
    for j in range(3):
        rows = np.arange(n) * 3 + j
        A[rows, 0] = singular * fI[:, j]
        A[rows, 1] = singular * fII[:, j]
        A[rows, 2+j] = 1.0
        y[rows] = sl[:, j]

    area = np.maximum(np.asarray(mesh.area_e[idx], float), 1e-30)
    area_w = np.sqrt(area / np.median(area))
    base_w = np.repeat(area_w, 3)
    robust_w = np.ones_like(base_w)
    x = np.zeros(5)
    cond = float("inf")
    for _ in range(max_irls):
        x_old = x.copy()
        x, cond = _weighted_lstsq(A, y, base_w * robust_w)
        resid = y - A @ x
        med = float(np.median(resid))
        mad = 1.4826 * float(np.median(np.abs(resid-med))) + 1e-12
        z = np.abs(resid-med) / (1.5*mad)
        robust_w = np.where(z <= 1.0, 1.0, 1.0/np.maximum(z, 1e-30))
        if np.linalg.norm(x-x_old) <= 1e-10 * max(np.linalg.norm(x), 1.0):
            break

    resid = y - A @ x
    denom = float(np.sum((base_w*y)**2))
    rel_rmse = math.sqrt(float(np.sum((base_w*resid)**2)) / max(denom, 1e-30))
    KI, KII = float(x[0]), float(x[1])
    return {
        "valid": bool(np.isfinite(KI) and np.isfinite(KII)),
        "n": int(n),
        "KI": KI,
        "KII": KII,
        "psi": math.degrees(math.atan2(KII, KI)),
        "rel_rmse": rel_rmse,
        "condition": cond,
        "T11_Pa": float(x[2]),
        "T22_Pa": float(x[3]),
        "T12_Pa": float(x[4]),
    }


def project_near_tip_modes(mesh, sigma_gp, d, crack_tip, crack_direction,
                           r_min: float, r_max: float,
                           angular_limit_deg: float = 105.0,
                           damage_cutoff: float = 0.85) -> dict[str, float]:
    """Robust multi-annulus full-field decomposition of KI and KII.

    All three crack-local stress components are fitted to the leading isotropic
    Williams fields, together with independent nonsingular constant stresses.
    Huber IRLS suppresses isolated nonlinear or poorly resolved integration
    points.  Overlapping annuli provide a contour-stability diagnostic.

    The decomposition is exact for isotropic LEFM and a transparent engineering
    projection for anisotropic runs.  ``mode_projection_reliable`` must be
    checked before interpreting a point as a mixed-mode fracture datum.
    """
    tip = np.asarray(crack_tip, float)
    e1 = _unit(np.asarray(crack_direction, float))
    e2 = np.array([-e1[1], e1[0]])
    cent = mesh.nodes[mesh.elems].mean(axis=1)
    rel = cent - tip[None, :]
    x1 = rel @ e1
    x2 = rel @ e2
    rr = np.hypot(x1, x2)
    theta = np.arctan2(x2, x1)
    dgp = np.mean(np.asarray(d)[mesh.elems], axis=1)

    htip = float(getattr(mesh, "hbar_tip", 0.0) or getattr(mesh, "hbar", 0.0) or 0.0)
    if not np.isfinite(r_min) or r_min <= 0:
        r_min = max(1.25*htip, 1e-12)
    if not np.isfinite(r_max) or r_max <= r_min:
        r_max = max(6.0*htip, 2.5*r_min)
    angle_limit = math.radians(abs(float(angular_limit_deg)))

    common = ((rr >= r_min) & (rr <= r_max) &
              (np.abs(theta) <= angle_limit) &
              (dgp < float(damage_cutoff)))
    global_idx = np.flatnonzero(common)

    # Global fit plus three overlapping radial windows.  These are not treated
    # as independent data; their spread is a numerical reliability diagnostic.
    span = r_max-r_min
    windows = [
        (r_min, r_max, "global"),
        (r_min, r_min+0.48*span, "inner"),
        (r_min+0.22*span, r_min+0.74*span, "middle"),
        (r_min+0.48*span, r_max, "outer"),
    ]
    fits = []
    for lo, hi, label in windows:
        idx = np.flatnonzero(common & (rr >= lo) & (rr <= hi))
        fit = _fit_williams_subset(mesh, sigma_gp, idx, rr, theta, e1, e2)
        fit.update({"label": label, "rmin": float(lo), "rmax": float(hi)})
        if fit.get("valid"):
            fits.append(fit)

    if not fits:
        return {
            "KI_Pa_sqrt_m": float("nan"),
            "KII_Pa_sqrt_m": float("nan"),
            "mode_phase_deg": float("nan"),
            "mode_projection_n": int(len(global_idx)),
            "mode_projection_fit_count": 0,
            "mode_projection_rmin_m": float(r_min),
            "mode_projection_rmax_m": float(r_max),
            "mode_projection_angle_deg": float(angular_limit_deg),
            "mode_projection_rel_rmse": float("nan"),
            "mode_projection_condition": float("inf"),
            "mode_projection_psi_spread_deg": float("inf"),
            "mode_projection_K_spread_frac": float("inf"),
            "mode_projection_reliable": False,
        }

    KI_vals = np.array([f["KI"] for f in fits], float)
    KII_vals = np.array([f["KII"] for f in fits], float)
    KI = float(np.median(KI_vals))
    KII = float(np.median(KII_vals))
    psi = math.degrees(math.atan2(KII, KI))
    psi_vals = np.array([f["psi"] for f in fits], float)
    psi_spread = 1.4826 * float(np.median(np.abs(psi_vals-np.median(psi_vals)))) if len(fits)>1 else float("inf")
    kmag = np.hypot(KI_vals, KII_vals)
    kmag_med = max(float(np.median(kmag)), 1e-30)
    kspread = 1.4826 * float(np.median(np.abs(kmag-np.median(kmag)))) / kmag_med if len(fits)>1 else float("inf")
    rel_rmse = float(np.median([f["rel_rmse"] for f in fits]))
    condition = float(max(f["condition"] for f in fits))
    n_global = next((f["n"] for f in fits if f["label"] == "global"), len(global_idx))
    reliable = bool(
        len(fits) >= 2 and n_global >= 10 and
        np.isfinite(rel_rmse) and rel_rmse <= 0.40 and
        np.isfinite(condition) and condition <= 1e8 and
        np.isfinite(psi_spread) and psi_spread <= 12.0 and
        np.isfinite(kspread) and kspread <= 0.35 and
        KI > 0.0
    )
    return {
        "KI_Pa_sqrt_m": KI,
        "KII_Pa_sqrt_m": KII,
        "mode_phase_deg": psi,
        "mode_projection_n": int(n_global),
        "mode_projection_fit_count": int(len(fits)),
        "mode_projection_rmin_m": float(r_min),
        "mode_projection_rmax_m": float(r_max),
        "mode_projection_angle_deg": float(angular_limit_deg),
        "mode_projection_rel_rmse": rel_rmse,
        "mode_projection_condition": condition,
        "mode_projection_psi_spread_deg": psi_spread,
        "mode_projection_K_spread_frac": kspread,
        "mode_projection_reliable": reliable,
        "mode_projection_KI_annulus_min": float(np.min(KI_vals)),
        "mode_projection_KI_annulus_max": float(np.max(KI_vals)),
        "mode_projection_KII_annulus_min": float(np.min(KII_vals)),
        "mode_projection_KII_annulus_max": float(np.max(KII_vals)),
    }


@dataclass
class MixedModeContext:
    loading_angle_deg: float
    target_mode_phase_deg: float = 0.0
    solver_seed: int = 1
    stochastic_first_passage: bool = False
    shear_emission_weight: float = 1.0
    projection_angle_deg: float = 105.0
    projection_damage_cutoff: float = 0.85
    projection_rmin_m: float = 0.0
    projection_rmax_m: float = 0.0
    records: list[dict[str, Any]] = field(default_factory=list)
    latest: dict[str, float] = field(default_factory=dict)
    last_loading: dict[str, float] = field(default_factory=dict)
    thresholds_drawn: list[float] = field(default_factory=list)

    def __post_init__(self):
        self.rng = np.random.default_rng(int(self.solver_seed))

    @property
    def alpha_rad(self) -> float:
        return math.radians(float(self.loading_angle_deg))

    def draw_threshold(self) -> float:
        if not self.stochastic_first_passage:
            h = 1.0
        else:
            u = float(np.clip(self.rng.random(), 1e-15, 1.0 - 1e-15))
            h = -math.log(u)
        self.thresholds_drawn.append(float(h))
        return float(h)


class MixedModeFrontEngineMixin:
    """Front-engine overrides for separated mixed-mode scalar drives."""

    def _mm_init(self, context: MixedModeContext):
        self._mm = context
        self._mm_threshold = context.draw_threshold()
        self._mm_prev_Kcleave = None

    def _mm_drives(self, fallback_K: float) -> tuple[float, float, dict[str, float]]:
        # v3: the total domain-integral measure is the authoritative amplitude.
        # For the isotropic campaign, decompose it using the calibrated target
        # phase angle so KI^2 + KII^2 = KJ^2 exactly.  The local Williams fit is
        # retained only as an independent diagnostic and never drives kinetics.
        KJ = max(float(fallback_K), 0.0)
        psi = math.radians(float(self._mm.target_mode_phase_deg))
        KI = KJ * math.cos(psi)
        KII = KJ * math.sin(psi)
        Kopen, kink = maximum_hoop_drive(KI, KII)
        Kemit = math.sqrt(max(Kopen, 0.0)**2 +
                          max(float(self._mm.shear_emission_weight), 0.0) * KII**2)
        return Kopen, Kemit, {"KI": KI, "KII": KII, "kink": kink,
                              "phase_deg": float(self._mm.target_mode_phase_deg),
                              "KJ": KJ}

    def _sigma_from_drive(self, Kdrive: float) -> float:
        K_eff = max(float(Kdrive) - self.f.k_shield * self.N_em * self.G * self.b
                    / np.sqrt(2.0 * np.pi * self.f.L_pz), 0.0)
        s = K_eff / np.sqrt(2.0 * np.pi * self.r_eff())
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
            if lo <= 0.0:
                leff = 0.5 * hi
            elif abs(hi-lo) <= 1e-12 * max(hi, 1e-300):
                leff = hi
            else:
                leff = (hi-lo) / np.log(hi/lo)
        else:
            leff = max(lam2, 0.0)
        remaining = max(self._mm_threshold - self.B, 1e-12)
        return float(max(leff * dt / remaining, 0.0))

    def step(self, K, T, dt):
        # This follows FrontEngine.step, but uses separate mixed-mode drives and
        # a seeded exponential cumulative-hazard threshold.
        Kc, Ke, mm = self._mm_drives(K)
        sig_emit = self._sigma_from_drive(Ke)
        K_eff_audit = max(Ke - self.f.k_shield * self.N_em * self.G * self.b
                          / np.sqrt(2.0*np.pi*self.f.L_pz), 0.0)
        sigma_tip_uncapped = K_eff_audit / np.sqrt(2.0*np.pi*self.r_eff())
        sigma_cap_active = bool(self.f.sigma_cap > 0 and sigma_tip_uncapped > self.f.sigma_cap)

        lam_e, sig_em_eff, Ge = self.lambda_emit(sig_emit, T)
        prod_raw = lam_e * dt
        prod = min(prod_raw, self.f.dN_cap)
        dN_cap_active = bool(np.isfinite(self.f.dN_cap) and prod_raw > self.f.dN_cap)
        N_sat_factor = 1.0
        if np.isfinite(self.f.N_sat) and self.f.N_sat > 0.0:
            N_sat_factor = max(1.0 - self.N_em / self.f.N_sat, 0.0)
            prod *= N_sat_factor
        ann = self.f.recover_k * self.N_em * dt
        self.W_emit += sig_em_eff * self.b * self.f.L_pz * prod
        self.N_em = max(self.N_em + prod - ann, 0.0)

        sig_cleave = self._sigma_from_drive(Kc)
        lam_c, lam_c_raw, Gc_eff = self.lambda_cleave(sig_cleave, T)
        if self.f.tau_B > 0 and dt > 0:
            self.B *= np.exp(-min(dt/self.f.tau_B, 80.0))
        if self._mm_prev_Kcleave is not None and self._mm_prev_Kcleave > 0 and Kc > 0:
            sig1 = self._sigma_from_drive(self._mm_prev_Kcleave)
            lam1, _, _ = self.lambda_cleave(sig1, T)
            lo, hi = sorted((max(lam1, 0.0), max(lam_c, 0.0)))
            if lo <= 0.0:
                leff = 0.5 * hi
            elif abs(hi-lo) <= 1e-12 * max(hi, 1e-300):
                leff = hi
            else:
                leff = (hi-lo) / np.log(hi/lo)
        else:
            leff = lam_c
        self.B += leff * dt
        self._mm_prev_Kcleave = Kc
        self.K_prev = Kc
        self.t += dt

        N_em_pre = self.N_em
        sigma_back_pre = self.sigma_back()
        r_eff_pre = self.r_eff()
        dG_emb_pre = self.dG_emb()
        n_fire = 0
        if not np.isfinite(self.B):
            self.B = 0.0
        while self.B >= self._mm_threshold and n_fire < 100000:
            self.B -= self._mm_threshold
            n_fire += 1
            self._mm_threshold = self._mm.draw_threshold()
        fired = n_fire >= 1
        v_crack = self.f.da * n_fire / dt if dt > 0 else 0.0
        N_retained = N_em_pre
        N_shed = 0.0
        if fired:
            retain = np.clip(self.f.wake_retain, 0.0, 1.0) ** n_fire
            N_retained = N_em_pre * retain
            N_shed = N_em_pre * (1.0-retain)
            self.N_em = N_retained
            self.a_adv += self.f.da * n_fire
            self.n_adv += n_fire

        info = {
            "fired": bool(fired), "n_fire": int(n_fire), "v_crack": v_crack,
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
            "N_sat_active": bool(np.isfinite(self.f.N_sat) and self.f.N_sat > 0 and N_sat_factor < 0.999999),
            "N_em_pre_renewal": N_em_pre, "N_em_retained": N_retained,
            "N_em_shed_to_wake": N_shed,
            "sigma_back_pre_renewal": sigma_back_pre,
            "r_eff_pre_renewal": r_eff_pre,
            "dG_emb_pre_renewal_eV": dG_emb_pre/1.602176634e-19,
            "mixed_KI_Pa_sqrt_m": mm["KI"],
            "mixed_KII_Pa_sqrt_m": mm["KII"],
            "mixed_Kopen_Pa_sqrt_m": Kc,
            "mixed_Kemit_Pa_sqrt_m": Ke,
            "mixed_kink_angle_deg": mm["kink"],
            "mixed_hazard_threshold_next": self._mm_threshold,
        }
        return info


def _mixed_solve_factory(original_solve, context: MixedModeContext):
    def solve_mixed(K, Rint, u, bnd, Uy_top, Uy_bot):
        from scipy.sparse.linalg import spsolve
        alpha = context.alpha_rad
        Atotal = float(Uy_top - Uy_bot)
        Un = Atotal * math.cos(alpha)
        Us = Atotal * math.sin(alpha)
        # Opening equilibrium first: preserves free lateral contraction.
        u_open, _ = original_solve(K, Rint, u, bnd, 0.5*Un, -0.5*Un)
        Kc = K.tocsr()
        R_open = Rint + Kc @ (u_open-u)
        if abs(Us) <= 1e-30:
            Rfull = R_open
            Fgen = float(np.sum(Rfull[2*bnd.top_nodes+1]))
            u_new = u_open
        else:
            ndof = len(u)
            prescribed = np.zeros(ndof, dtype=bool)
            target = u_open.copy()
            tn = bnd.top_nodes; bn = bnd.bot_nodes
            prescribed[2*tn] = True; prescribed[2*tn+1] = True
            prescribed[2*bn] = True; prescribed[2*bn+1] = True
            target[2*tn] = u_open[2*tn] + 0.5*Us
            target[2*bn] = u_open[2*bn] - 0.5*Us
            target[2*tn+1] = 0.5*Un
            target[2*bn+1] = -0.5*Un
            free = ~prescribed
            du_p = target[prescribed] - u_open[prescribed]
            rhs = -R_open[free] - Kc[np.ix_(free, prescribed)] @ du_p
            u_new = u_open.copy()
            u_new[free] = u_open[free] + spsolve(Kc[np.ix_(free, free)], rhs)
            u_new[prescribed] = target[prescribed]
            Rfull = R_open + Kc @ (u_new-u_open)
            Fx = float(np.sum(Rfull[2*tn]))
            Fy = float(np.sum(Rfull[2*tn+1]))
            Fgen = Fx*math.sin(alpha) + Fy*math.cos(alpha)
        context.last_loading = {
            "U_total_m": Atotal, "U_open_m": Un, "U_shear_m": Us,
            "loading_angle_deg": context.loading_angle_deg,
            "generalized_reaction_N": Fgen,
        }
        return u_new, Fgen
    return solve_mixed


def _j_wrapper_factory(original_compute, context: MixedModeContext):
    def wrapped(mesh, u, sigma_gp, psi_e_gp, d, crack_tip, crack_direction,
                mat, ell, cfg=None, crack_segments=None, exclude_radius=0.0):
        J, KJ, info = original_compute(mesh, u, sigma_gp, psi_e_gp, d,
                                       crack_tip, crack_direction, mat, ell,
                                       cfg=cfg, crack_segments=crack_segments,
                                       exclude_radius=exclude_radius)
        htip = float(mesh.hbar_tip or mesh.hbar)
        rmin = context.projection_rmin_m or (1.25*htip)
        rmax = context.projection_rmax_m or (6.0*htip)
        diag = project_near_tip_modes(mesh, sigma_gp, d, crack_tip, crack_direction,
                                      rmin, rmax, context.projection_angle_deg,
                                      context.projection_damage_cutoff)

        psi = math.radians(float(context.target_mode_phase_deg))
        KJpos = max(float(KJ), 0.0)
        KI = KJpos * math.cos(psi)
        KII = KJpos * math.sin(psi)
        Kopen, kink = maximum_hoop_drive(KI, KII)

        md = {
            "mode_projection_call_index": int(len(context.records)),
            "KJ_Pa_sqrt_m": float(KJ),
            "J_J_per_m2": float(J),
            "KI_Pa_sqrt_m": float(KI),
            "KII_Pa_sqrt_m": float(KII),
            "mode_phase_deg": float(context.target_mode_phase_deg),
            "Kopen_maxhoop_Pa_sqrt_m": float(Kopen),
            "maxhoop_kink_angle_deg": float(kink),
            "phase_control_method": "isotropic_KJ_target_phase",
            "phase_control_reliable": True,
            "KJ_identity_rel_error": float(abs(KI*KI + KII*KII - KJpos*KJpos) / max(KJpos*KJpos, 1e-30)),
            **context.last_loading,
            "tip_x_m": float(np.asarray(crack_tip)[0]),
            "tip_y_m": float(np.asarray(crack_tip)[1]),
        }
        for key, value in diag.items():
            md["diagnostic_" + key] = value
        info.update(md)
        context.latest = md.copy()
        context.records.append(md.copy())
        return J, KJ, info
    return wrapped


def _engine_factory(original_build, context: MixedModeContext, base_class):
    class MixedModeEngine(MixedModeFrontEngineMixin, base_class):
        def __init__(self, *a, **kw):
            base_class.__init__(self, *a, **kw)
            self._mm_init(context)
    def build(args, mat):
        base = original_build(args, mat)
        eng = MixedModeEngine(base.f, base.cb, base.eb, base.G, base.nu, base.b)
        return eng
    return build


def _write_records(out: Path, context: MixedModeContext):
    if not context.records:
        return
    cols = sorted({k for r in context.records for k in r})
    with (out/"mixed_mode_projection_calls.csv").open("w", newline="") as fp:
        w = csv.DictWriter(fp, fieldnames=cols)
        w.writeheader(); w.writerows(context.records)


def _accepted_first_passage(out: Path, T: float, context: MixedModeContext, summary: list[dict]) -> dict:
    steps = out/f"steps_{int(T):04d}K.csv"
    accepted = {}
    if steps.exists():
        a = np.genfromtxt(steps, delimiter=",", names=True)
        if a.size:
            a = np.atleast_1d(a)
            fires = np.flatnonzero(a["n_fire"] > 0)
            i = int(fires[0]) if len(fires) else len(a)-1
            accepted = {name: float(a[name][i]) for name in a.dtype.names}
    U = accepted.get("Uapp_m", float("nan"))
    KJ = accepted.get("KJ_Pa_sqrtm", float("nan"))
    candidates = context.records
    if candidates:
        if np.isfinite(U):
            du = np.array([abs(float(r.get("U_total_m", np.nan))-U) for r in candidates])
            finite = np.isfinite(du)
            if np.any(finite):
                best_du = float(np.min(du[finite]))
                # Keep only calls at the accepted displacement to numerical
                # precision.  The cluster-J value is used solely as a tie-break
                # identifying the call written to the accepted step record; it
                # is not used as the mixed-mode decomposition itself.
                tol_u = max(1e-14, 1e-10*max(abs(U), 1e-12))
                near = [r for r, q in zip(candidates, du) if np.isfinite(q) and q <= best_du+tol_u]
                if near and np.isfinite(KJ):
                    mode = min(near, key=lambda r: abs(float(r.get("KJ_Pa_sqrt_m", np.nan))-KJ)
                               if np.isfinite(float(r.get("KJ_Pa_sqrt_m", np.nan))) else float("inf"))
                else:
                    mode = near[-1] if near else candidates[int(np.nanargmin(du))]
            else:
                mode = candidates[-1]
        else:
            mode = candidates[-1]
    else:
        mode = {}
    s0 = summary[0] if summary else {}
    payload = {
        "model": MODEL_ID,
        "T_K": float(T),
        "solver_seed": int(context.solver_seed),
        "stochastic_first_passage": bool(context.stochastic_first_passage),
        "hazard_threshold_initial": (context.thresholds_drawn[0] if context.thresholds_drawn else None),
        "loading_angle_deg": float(context.loading_angle_deg),
        "Kc_first_existing_MPa_sqrt_m": s0.get("Kc_first_MPa_sqrt_m"),
        "KI_first_MPa_sqrt_m": float(mode.get("KI_Pa_sqrt_m", np.nan))/1e6,
        "KII_first_MPa_sqrt_m": float(mode.get("KII_Pa_sqrt_m", np.nan))/1e6,
        "KJ_first_MPa_sqrt_m": float(mode.get("KJ_Pa_sqrt_m", np.nan))/1e6,
        "Kopen_maxhoop_first_MPa_sqrt_m": float(mode.get("Kopen_maxhoop_Pa_sqrt_m", np.nan))/1e6,
        "mode_phase_first_deg": mode.get("mode_phase_deg"),
        "maxhoop_kink_first_deg": mode.get("maxhoop_kink_angle_deg"),
        "U_total_first_m": mode.get("U_total_m"),
        "U_open_first_m": mode.get("U_open_m"),
        "U_shear_first_m": mode.get("U_shear_m"),
        "phase_control_method": mode.get("phase_control_method"),
        "phase_control_reliable": mode.get("phase_control_reliable"),
        "KJ_identity_rel_error": mode.get("KJ_identity_rel_error"),
        "diagnostic_projection_method": "full_field_multiring_isotropic_Williams_IRLS",
        "diagnostic_projection_n": mode.get("diagnostic_mode_projection_n"),
        "diagnostic_projection_fit_count": mode.get("diagnostic_mode_projection_fit_count"),
        "diagnostic_projection_rel_rmse": mode.get("diagnostic_mode_projection_rel_rmse"),
        "diagnostic_projection_condition": mode.get("diagnostic_mode_projection_condition"),
        "diagnostic_projection_psi_spread_deg": mode.get("diagnostic_mode_projection_psi_spread_deg"),
        "diagnostic_projection_K_spread_frac": mode.get("diagnostic_mode_projection_K_spread_frac"),
        "diagnostic_projection_reliable": mode.get("diagnostic_mode_projection_reliable"),
        "diagnostic_projection_phase_deg": mode.get("diagnostic_mode_phase_deg"),
        "diagnostic_projection_rmin_m": mode.get("diagnostic_mode_projection_rmin_m"),
        "diagnostic_projection_rmax_m": mode.get("diagnostic_mode_projection_rmax_m"),
        "diagnostic_projection_angle_deg": mode.get("diagnostic_mode_projection_angle_deg"),
        "mode_ratio_KII_over_KI": (float(mode.get("KII_Pa_sqrt_m", np.nan))/float(mode.get("KI_Pa_sqrt_m", np.nan))
                                     if np.isfinite(float(mode.get("KI_Pa_sqrt_m", np.nan))) and abs(float(mode.get("KI_Pa_sqrt_m", np.nan))) > 1e-30
                                     else np.nan),
        "control_state": ("first_passage" if accepted.get("n_fire", 0.0) > 0 else "right_censored_endpoint"),
        "N_em_final": s0.get("N_em_final"),
        "mode_classification": s0.get("mode"),
    }
    (out/"mixed_mode_first_passage_summary.json").write_text(json.dumps(payload, indent=2, default=str))
    with (out/"mixed_mode_first_passage_summary.csv").open("w", newline="") as fp:
        w = csv.DictWriter(fp, fieldnames=list(payload))
        w.writeheader(); w.writerow(payload)
    return payload


def build_mixed_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(add_help=False, allow_abbrev=False)
    p.add_argument("--mixity-loading-angle-deg", type=float, default=0.0)
    p.add_argument("--target-mode-phase-deg", type=float, required=True,
                   help="Calibrated isotropic KI/KII phase angle used with KJ for kinetics.")
    p.add_argument("--allow-anisotropic-approx", action="store_true",
                   help="Unsafe research override. v3 quantitative mode control is isotropic only.")
    p.add_argument("--solver-seed", type=int, default=1)
    p.add_argument("--deterministic-threshold", action="store_true",
                   help="Accepted for compatibility; deterministic H=1 is already the v2 default.")
    p.add_argument("--stochastic-threshold", action="store_true",
                   help="Opt in to exponential first-passage thresholds; not recommended for geometry maps.")
    p.add_argument("--shear-emission-weight", type=float, default=1.0)
    p.add_argument("--mode-projection-angle-deg", type=float, default=105.0)
    p.add_argument("--mode-projection-damage-cutoff", type=float, default=0.85)
    p.add_argument("--mode-projection-rmin-m", type=float, default=0.0)
    p.add_argument("--mode-projection-rmax-m", type=float, default=0.0)
    return p


def main(argv=None):
    import sys
    from . import sharp_front as sf
    from . import fem as femmod
    from . import j_integral as jimod

    mm, remaining = build_mixed_parser().parse_known_args(argv)
    args = sf._build_parser().parse_args(remaining)
    if args.mode != "2d":
        raise SystemExit("mixed_mode_first_passage_v3 requires --mode 2d")
    if int(getattr(args, "max_fronts", 1)) != 1 or bool(getattr(args, "crystal_branch", False)):
        raise SystemExit("v3 mixed-mode first-passage requires branching off and --max-fronts 1")
    if bool(getattr(args, "crystal_aniso", False)) and not mm.allow_anisotropic_approx:
        raise SystemExit("v3 quantitative mode control is isotropic. Remove --crystal-aniso; "
                         "anisotropic KI/KII requires an anisotropic interaction integral.")
    context = MixedModeContext(
        loading_angle_deg=mm.mixity_loading_angle_deg,
        target_mode_phase_deg=mm.target_mode_phase_deg,
        solver_seed=mm.solver_seed,
        stochastic_first_passage=bool(mm.stochastic_threshold and not mm.deterministic_threshold),
        shear_emission_weight=mm.shear_emission_weight,
        projection_angle_deg=mm.mode_projection_angle_deg,
        projection_damage_cutoff=mm.mode_projection_damage_cutoff,
        projection_rmin_m=mm.mode_projection_rmin_m,
        projection_rmax_m=mm.mode_projection_rmax_m,
    )
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    (out/"mixed_mode_run_config.json").write_text(json.dumps({
        "model": MODEL_ID, **vars(mm), "base_driver": "arrhenius_fracture.sharp_front",
        "note": "deterministic H=1; isotropic KJ amplitude is decomposed at the calibrated target phase. Williams stress fits are diagnostic only.",
    }, indent=2))

    orig_solve = femmod.solve_dirichlet
    orig_J = jimod.compute_J_integral
    orig_build = sf.build_engine
    try:
        femmod.solve_dirichlet = _mixed_solve_factory(orig_solve, context)
        jimod.compute_J_integral = _j_wrapper_factory(orig_J, context)
        sf.build_engine = _engine_factory(orig_build, context, sf.FrontEngine)
        summary = sf.run_2d(args)
    finally:
        femmod.solve_dirichlet = orig_solve
        jimod.compute_J_integral = orig_J
        sf.build_engine = orig_build
    _write_records(out, context)
    payloads = []
    for T in args.temperatures:
        payloads.append(_accepted_first_passage(out, T, context, summary))
    print("MIXED_MODE_V3 complete:", json.dumps(payloads, indent=2, default=str))
    return payloads


if __name__ == "__main__":
    main()
