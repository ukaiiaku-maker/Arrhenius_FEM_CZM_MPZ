"""Mixed-mode first-passage wrapper for the active FEM/CZM sharp-front solver.

This module does not replace ``arrhenius_fracture.sharp_front``.  It temporarily
wraps three interfaces during one run:

* the displacement boundary solve, to impose combined opening/sliding loading;
* the domain J calculation, to extract local KI/KII from the near-tip FEM stress;
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

MODEL_ID = "FEM_CZM_mixed_mode_first_passage_v1"


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


def project_near_tip_modes(mesh, sigma_gp, d, crack_tip, crack_direction,
                           r_min: float, r_max: float,
                           wedge_deg: float = 18.0) -> dict[str, float]:
    """Project KI and KII from stresses directly ahead of the crack tip.

    In crack-local coordinates, the leading isotropic fields at theta=0 are

        sigma_22 = KI / sqrt(2*pi*r),
        sigma_12 = KII / sqrt(2*pi*r).

    We fit each stress component to ``K/sqrt(2*pi*r) + constant`` over an
    annular forward wedge.  The constant absorbs the leading nonsingular stress.
    """
    tip = np.asarray(crack_tip, float)
    e1 = _unit(np.asarray(crack_direction, float))
    e2 = np.array([-e1[1], e1[0]])
    cent = mesh.nodes[mesh.elems].mean(axis=1)
    rel = cent - tip[None, :]
    x1 = rel @ e1
    x2 = rel @ e2
    rr = np.hypot(x1, x2)
    ang = np.degrees(np.arctan2(x2, np.maximum(x1, 1e-300)))
    dgp = np.mean(np.asarray(d)[mesh.elems], axis=1)
    # Use the requested forward wedge first.  On a coarse graded mesh the
    # annulus can contain fewer than four centroids, so expand in controlled
    # stages rather than returning a silent NaN.
    trials = [
        (float(r_min), float(r_max), abs(float(wedge_deg))),
        (max(0.5*float(getattr(mesh, "hbar_tip", 0.0) or getattr(mesh, "hbar", 0.0)), 0.25*r_min),
         max(1.5*r_max, 4.0*float(getattr(mesh, "hbar_tip", 0.0) or getattr(mesh, "hbar", 0.0))),
         max(abs(float(wedge_deg)), 30.0)),
        (max(0.25*r_min, 1e-12), max(2.5*r_max, 6.0*float(getattr(mesh, "hbar_tip", 0.0) or getattr(mesh, "hbar", 0.0))), 45.0),
    ]
    idx = np.array([], dtype=int)
    used = trials[0]
    for rlo, rhi, aw in trials:
        mask = ((x1 > 0.0) & (rr >= rlo) & (rr <= rhi) &
                (np.abs(ang) <= aw) & (dgp < 0.90))
        idx = np.flatnonzero(mask)
        used = (rlo, rhi, aw)
        if len(idx) >= 4:
            break
    if len(idx) < 4:
        return {
            "KI_Pa_sqrt_m": float("nan"),
            "KII_Pa_sqrt_m": float("nan"),
            "mode_phase_deg": float("nan"),
            "mode_projection_n": int(len(idx)),
            "mode_projection_rmin_m": float(used[0]),
            "mode_projection_rmax_m": float(used[1]),
            "mode_projection_wedge_deg": float(used[2]),
        }

    # Transform stress to crack-local coordinates.
    sx = np.asarray(sigma_gp[0, idx], float)
    sy = np.asarray(sigma_gp[1, idx], float)
    txy = np.asarray(sigma_gp[2, idx], float)
    R = np.column_stack([e1, e2])
    s22 = np.empty(len(idx))
    s12 = np.empty(len(idx))
    for q, (a, b, c) in enumerate(zip(sx, sy, txy)):
        S = np.array([[a, c], [c, b]])
        Sl = R.T @ S @ R
        s22[q] = Sl[1, 1]
        s12[q] = Sl[0, 1]

    basis = 1.0 / np.sqrt(2.0 * np.pi * rr[idx])
    A = np.column_stack([basis, np.ones_like(basis)])
    # Area weighting reduces sensitivity to irregular point density.
    w = np.sqrt(np.maximum(mesh.area_e[idx], 1e-30))
    Aw = A * w[:, None]
    KI = float(np.linalg.lstsq(Aw, s22 * w, rcond=None)[0][0])
    KII = float(np.linalg.lstsq(Aw, s12 * w, rcond=None)[0][0])
    psi = math.degrees(math.atan2(KII, KI))
    return {
        "KI_Pa_sqrt_m": KI,
        "KII_Pa_sqrt_m": KII,
        "mode_phase_deg": psi,
        "mode_projection_n": int(len(idx)),
        "mode_projection_rmin_m": float(used[0]),
        "mode_projection_rmax_m": float(used[1]),
        "mode_projection_wedge_deg": float(used[2]),
    }


@dataclass
class MixedModeContext:
    loading_angle_deg: float
    solver_seed: int = 1
    stochastic_first_passage: bool = True
    shear_emission_weight: float = 1.0
    projection_wedge_deg: float = 18.0
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
        d = self._mm.latest
        KI = float(d.get("KI_Pa_sqrt_m", float("nan")))
        KII = float(d.get("KII_Pa_sqrt_m", float("nan")))
        if not np.isfinite(KI) or not np.isfinite(KII):
            KI = max(float(fallback_K), 0.0)
            KII = 0.0
        Kopen, kink = maximum_hoop_drive(KI, KII)
        Kemit = math.sqrt(max(Kopen, 0.0)**2 +
                          max(float(self._mm.shear_emission_weight), 0.0) * KII**2)
        return Kopen, Kemit, {"KI": KI, "KII": KII, "kink": kink}

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
        outer = float(info.get("r_outer", 8.0*ell))
        # Keep the mode projection tied to the tip resolution rather than the
        # cluster-J contour.  The primary cluster contour can be much larger
        # than the asymptotic K-dominant zone and otherwise drives psi toward
        # the far-field loading state.  Explicit CLI radii still override.
        htip = float(mesh.hbar_tip or mesh.hbar)
        rmin = context.projection_rmin_m or (1.5*htip)
        rmax = context.projection_rmax_m or (2.25*htip)
        md = project_near_tip_modes(mesh, sigma_gp, d, crack_tip, crack_direction,
                                    rmin, rmax, context.projection_wedge_deg)
        Kopen, kink = maximum_hoop_drive(md["KI_Pa_sqrt_m"], md["KII_Pa_sqrt_m"]) \
            if np.isfinite(md["KI_Pa_sqrt_m"]) else (float("nan"), float("nan"))
        md.update({
            "KJ_Pa_sqrt_m": float(KJ), "J_J_per_m2": float(J),
            "Kopen_maxhoop_Pa_sqrt_m": float(Kopen),
            "maxhoop_kink_angle_deg": float(kink),
            **context.last_loading,
            "tip_x_m": float(np.asarray(crack_tip)[0]),
            "tip_y_m": float(np.asarray(crack_tip)[1]),
        })
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
        def score(r):
            su = abs(float(r.get("U_total_m", 0.0))-U) / max(abs(U), 1e-20) if np.isfinite(U) else 0.0
            sk = abs(float(r.get("KJ_Pa_sqrt_m", 0.0))-KJ) / max(abs(KJ), 1.0) if np.isfinite(KJ) else 0.0
            return su + sk
        mode = min(candidates, key=score)
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
        "projection_method": "near_tip_forward_wedge_isotropic_Williams_fit",
        "projection_n": mode.get("mode_projection_n"),
        "projection_rmin_m": mode.get("mode_projection_rmin_m"),
        "projection_rmax_m": mode.get("mode_projection_rmax_m"),
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
    p.add_argument("--solver-seed", type=int, default=1)
    p.add_argument("--deterministic-threshold", action="store_true")
    p.add_argument("--shear-emission-weight", type=float, default=1.0)
    p.add_argument("--mode-projection-wedge-deg", type=float, default=18.0)
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
        raise SystemExit("mixed_mode_first_passage_v1 requires --mode 2d")
    if int(getattr(args, "max_fronts", 1)) != 1 or bool(getattr(args, "crystal_branch", False)):
        raise SystemExit("v1 mixed-mode first-passage requires branching off and --max-fronts 1")
    context = MixedModeContext(
        loading_angle_deg=mm.mixity_loading_angle_deg,
        solver_seed=mm.solver_seed,
        stochastic_first_passage=not mm.deterministic_threshold,
        shear_emission_weight=mm.shear_emission_weight,
        projection_wedge_deg=mm.mode_projection_wedge_deg,
        projection_rmin_m=mm.mode_projection_rmin_m,
        projection_rmax_m=mm.mode_projection_rmax_m,
    )
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    (out/"mixed_mode_run_config.json").write_text(json.dumps({
        "model": MODEL_ID, **vars(mm), "base_driver": "arrhenius_fracture.sharp_front",
        "note": "loading angle is boundary displacement angle; achieved SIF phase angle is measured from FEM",
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
    print("MIXED_MODE_V1 complete:", json.dumps(payloads, indent=2, default=str))
    return payloads


if __name__ == "__main__":
    main()
