#!/usr/bin/env python3
"""V5.7 incremental extension of an existing V5.6 barrier atlas (corrected).

This driver deliberately does NOT recompute the strength map or the blunt-notch
S-N calculations already stored in a V5.6 atlas.  It reads those raw tables and
adds two projections:

1. monotonic K-ramp fracture response Kc(T) under several existing crack-opening
   / shielding contexts (Option B: several fracture contexts per plastic surface);
2. ACTUAL simulated 1-D notch crack-growth da/dN(DeltaK,T) per (surface, fracture
   context, T), with rate-defined thresholds DeltaK_th at da/dN = 1e-10 m/cycle
   (primary) and 1e-12 m/cycle (sensitivity), found by adaptive DeltaK bracketing
   of the simulated crack-growth curve.

CORRECTION vs the first V5.7 draft: the life-defined "notch-equivalent DeltaK"
geometric conversion has been REMOVED.  It was not the same observable as the
existing rate-defined crack-growth threshold and would have created an
inconsistency between the atlas and the fatigue section of the paper.  The
fatigue-growth stage below calls the same established 1-D block-stepping path
used by the refined two-barrier study (run_v1_two_barrier_dbtt_fatigue_map
_corrected.py::_map_cycle_step) and the same adaptive bracketing numerics
(run_adaptive_two_barrier_threshold_study.py::locate_bracket/crossing_estimate/
threshold_record), imported at run time so the physics and threshold logic are
never duplicated.

Restartability:
  monotonic task key: (surface_id, context_id, T_K)              -> skipped under --resume
  fatigue raw key:    (surface_id, context_id, T_K, DeltaK)      -> NEVER recomputed
The adaptive threshold controller reads all existing DeltaK points for a state,
checks whether every criterion is bracketed to tolerance, proposes only the next
missing refinement point, runs it through the V1 driver, appends it to the raw
CSV, and repeats until converged.  An interrupted run therefore resumes even in
the middle of a threshold refinement.

The four observables of the atlas remain SEPARATE, as in the refined two-barrier
study: sigma_y(T) (reused from V5.6), blunt-feature S-N initiation (reused from
V5.6), monotonic Kc(T) (stage 1), and rate-defined DeltaK_th(T) (stage 2).
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
from pathlib import Path
from types import SimpleNamespace
from typing import Dict, List, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.special import gammainc
from scipy.stats import chi2_contingency, spearmanr

from arrhenius_fracture.config import ElasticProperties, FractureBarrier, KB, EV_TO_J
from arrhenius_fracture.sharp_front import FrontConfig, FrontEngine
from arrhenius_fracture.sn_arrhenius_chain import build_chain_from_namespace
from arrhenius_fracture.fatigue_v1 import (
    FatigueWaveform,
    FatigueControllerConfig,
    FatigueCycleHazardController,
)

KBEV = 8.617333262145e-5

# W-like elastic constants used by the existing V1 sharp-front model.
_MAT = ElasticProperties()
E_PA = _MAT.E
NU = _MAT.nu
G_PA = _MAT.G
B_M = _MAT.b


# ----------------------------------------------------------------------------
# run-time imports of the established machinery (no physics duplication)
# ----------------------------------------------------------------------------

def _load_module(path: Path, name: str, required: Sequence[str]):
    spec = importlib.util.spec_from_file_location(name, str(path))
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot import {path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    missing = [r for r in required if not hasattr(mod, r)]
    if missing:
        raise AttributeError(f"{path} does not define required members {missing}")
    return mod


def load_map_runner(path: Path):
    """Corrected two-barrier map runner: the established 1-D fatigue block path."""
    return _load_module(path, "two_barrier_map_runner_v57", [
        "_map_cycle_step", "_cycle_cleave_hazard", "_renew_from_clock",
        "AnchoredCleavageBarrier", "_compat_barrier_diagnostics",
        "_parse_float_or_inf",
    ])


def load_adaptive_study(path: Path):
    """Adaptive threshold study: the established bracketing/refinement numerics."""
    return _load_module(path, "adaptive_threshold_study_v57", [
        "effective_rate", "locate_bracket", "crossing_estimate", "threshold_record",
    ])


# ----------------------------------------------------------------------------
# shared small utilities (unchanged from the first V5.7 draft)
# ----------------------------------------------------------------------------

def find_one(directory: Path, patterns: Sequence[str]) -> Path:
    for pat in patterns:
        hits = sorted(directory.glob(pat))
        if hits:
            return hits[-1]
    raise FileNotFoundError(f"Could not find any of {patterns} in {directory}")


def parse_inf(x) -> float:
    if isinstance(x, str) and x.strip().lower() in {"inf", "+inf", "infinity", "+infinity"}:
        return float("inf")
    try:
        return float(x)
    except Exception:
        return float("inf")


def append_rows(path: Path, rows: Sequence[dict]) -> None:
    rows = list(rows)
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    # Stable field order: existing header first if present, otherwise union.
    if path.exists() and path.stat().st_size > 0:
        with path.open("r", newline="") as f:
            header = next(csv.reader(f))
        extra = [k for r in rows for k in r if k not in header]
        if extra:
            # Schema changed: rewrite through pandas to preserve existing data.
            old = pd.read_csv(path)
            new = pd.DataFrame(rows)
            pd.concat([old, new], ignore_index=True, sort=False).to_csv(path, index=False)
            return
        fields = header
        mode, write_header = "a", False
    else:
        fields = []
        for r in rows:
            for k in r:
                if k not in fields:
                    fields.append(k)
        mode, write_header = "w", True
    with path.open(mode, newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        if write_header:
            w.writeheader()
        w.writerows(rows)


def _emit_surface_terms(design: pd.DataFrame, T: float):
    # Same emission scaling used by V5.6 ArrheniusPlasticChain.
    G0 = np.maximum(
        0.75 * design["exp_G00_eV"].to_numpy(float)
        + 0.75 * design["exp_gT_eV_per_K"].to_numpy(float)
        * (float(T) - design["exp_Tref_K"].to_numpy(float)),
        1.0e-10,
    )
    sigc = np.maximum(
        design["exp_sigc0_GPa"].to_numpy(float) * 1.0e9
        + design["exp_sT_MPa_per_K"].to_numpy(float) * 1.0e6
        * (float(T) - design["exp_Tref_K"].to_numpy(float)),
        1.0e6,
    )
    floor_frac = design["exp_floor_frac"].to_numpy(float)
    Gfloor = np.minimum(0.95 * G0, np.maximum(1.0e-4 * 0.75, floor_frac * G0))
    a = design["exp_a"].to_numpy(float)
    n = design["exp_n"].to_numpy(float)
    return G0, sigc, Gfloor, a, n


def _exp_floor(G0, Gfloor, sigc, a, n, sigma):
    x = np.maximum(np.abs(np.asarray(sigma, float)), 0.0) / np.maximum(sigc, 1.0e6)
    return np.maximum(Gfloor + (G0 - Gfloor) * np.exp(-a * np.power(x, n)), 0.0)


def _cleavage_G_eV(context: dict, sigma_Pa: np.ndarray, T_K: float) -> np.ndarray:
    G00 = float(context["cleave_G00_eV"])
    sigc = float(context["cleave_sigc0_GPa"]) * 1.0e9
    a = float(context["cleave_exp_a"])
    n = float(context["cleave_exp_n"])
    floor = max(1.0e-4, float(context["cleave_floor_frac"]) * G00)
    x = np.maximum(np.asarray(sigma_Pa, float), 0.0) / max(sigc, 1.0e6)
    G300 = floor + (G00 - floor) * np.exp(-a * np.power(x, n))
    shift = -(float(T_K) - 300.0) * float(context.get("cleave_S_kB", 0.0)) * KBEV
    return np.maximum(G300 + shift, 0.0)


def _logmean_rate(lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    lo = np.maximum(np.asarray(lo, float), 0.0)
    hi = np.maximum(np.asarray(hi, float), 0.0)
    a = np.minimum(lo, hi)
    b = np.maximum(lo, hi)
    out = np.empty_like(a)
    zero = a <= 0.0
    close = (~zero) & (np.abs(b - a) <= 1.0e-12 * np.maximum(b, 1.0))
    reg = (~zero) & (~close)
    out[zero] = 0.5 * b[zero]
    out[close] = b[close]
    out[reg] = (b[reg] - a[reg]) / np.log(b[reg] / a[reg])
    return out


# ----------------------------------------------------------------------------
# stage 1: vectorized monotonic K-ramp Kc(T) (UNCHANGED from the first draft)
# ----------------------------------------------------------------------------

def run_monotonic_batch(design: pd.DataFrame, context: dict, T_K: float,
                        Kmax_MPa_sqrtm: float, dK_MPa_sqrtm: float,
                        Kdot_MPa_sqrtm_per_s: float) -> List[dict]:
    """Vectorized first-passage K-ramp, exact in the same V1 state equations.

    The loop is only over K increments; all surfaces in the block are advanced
    simultaneously as arrays.  This is the speed-critical path that makes
    several fracture contexts per V5.6 surface practical.
    """
    if design.empty:
        return []
    N = len(design)
    G0e, sigce, Gfe, ae, ne = _emit_surface_terms(design, T_K)

    # FrontConfig-equivalent context parameters.
    r0 = 1.0e-6
    sigma_cap = 30.0e9
    nu0_c = 1.0e12
    nu0_e = 1.0e11
    m_hits = 3.0
    tau_c = 1.0e-6
    dN_cap = 50.0
    beta_back = 1.0
    L_pz = 1.0e-6
    rho0 = 5.0e12

    chi = float(context["chi_shield"])
    Nsat = parse_inf(context["N_sat"])
    emb_sat_frac = float(context["emb_sat_frac"])
    c_blunt = float(context["c_blunt"])
    v_emb_b3 = float(context["v_emb_b3"])
    recover_k = float(context["recover_k"])

    dK = float(dK_MPa_sqrtm) * 1.0e6
    Kmax = float(Kmax_MPa_sqrtm) * 1.0e6
    Kdot = float(Kdot_MPa_sqrtm_per_s) * 1.0e6
    dt = dK / max(Kdot, 1.0e-300)
    nsteps = int(math.ceil(Kmax / dK))

    Nem = np.zeros(N, float)
    Bclock = np.zeros(N, float)
    Kprev = np.zeros(N, float)
    active = np.ones(N, bool)
    Kc = np.full(N, np.nan)
    Nem_fire = np.full(N, np.nan)
    shield_fire = np.full(N, np.nan)
    blunt_fire = np.full(N, np.nan)
    dGemb_fire = np.full(N, np.nan)

    back_coeff = beta_back * G_PA * B_M / (2.0 * np.pi * (1.0 - NU) * L_pz)
    stored_coeff_eV = (0.5 * G_PA * B_M**2) * (v_emb_b3 * B_M**3) / EV_TO_J

    for istep in range(1, nsteps + 1):
        if not np.any(active):
            break
        idx = np.flatnonzero(active)
        K = min(istep * dK, Kmax)
        Ni = Nem[idx]

        r_eff = r0 + c_blunt * B_M * Ni
        sig_tip = np.minimum(K / np.sqrt(2.0 * np.pi * r_eff), sigma_cap)
        sig_back = back_coeff * Ni
        sig_emit = np.maximum(sig_tip - sig_back, 0.0)
        Ge = _exp_floor(G0e[idx], Gfe[idx], sigce[idx], ae[idx], ne[idx], sig_emit)
        lam_e = nu0_e * np.exp(np.clip(-Ge / max(KBEV * T_K, 1.0e-30), -700.0, 0.0))
        prod = np.minimum(lam_e * dt, dN_cap)
        if np.isfinite(Nsat) and Nsat > 0.0:
            prod *= np.maximum(1.0 - Ni / Nsat, 0.0)
        ann = recover_k * Ni * dt
        Ni_new = np.maximum(Ni + prod - ann, 0.0)
        Nem[idx] = Ni_new

        # Post-emission cleavage state, matching FrontEngine.step ordering.
        r_eff2 = r0 + c_blunt * B_M * Ni_new
        sig_tip2 = np.minimum(K / np.sqrt(2.0 * np.pi * r_eff2), sigma_cap)
        sig_back2 = back_coeff * Ni_new
        sig_c = np.maximum(sig_tip2 - chi * sig_back2, 0.0)
        Gc = _cleavage_G_eV(context, sig_c, T_K)
        rho = rho0 + Ni_new / (L_pz**2)
        dGemb = stored_coeff_eV * rho
        dGemb_eff = np.minimum(dGemb, emb_sat_frac * Gc)
        Geff = np.maximum(Gc - dGemb_eff, 0.0)
        lam_raw = nu0_c * np.exp(np.clip(-Geff / max(KBEV * T_K, 1.0e-30), -700.0, 0.0))
        lam_c = gammainc(m_hits, np.minimum(lam_raw * tau_c, 1.0e12)) / tau_c

        # Previous-K endpoint evaluated at the current post-emission state.
        kp = Kprev[idx]
        sig_prev = np.zeros_like(sig_tip2)
        pos = kp > 0.0
        sig_prev[pos] = np.minimum(kp[pos] / np.sqrt(2.0 * np.pi * r_eff2[pos]), sigma_cap)
        sig_c_prev = np.maximum(sig_prev - chi * sig_back2, 0.0)
        Gc_prev = _cleavage_G_eV(context, sig_c_prev, T_K)
        dGemb_prev = np.minimum(dGemb, emb_sat_frac * Gc_prev)
        Geff_prev = np.maximum(Gc_prev - dGemb_prev, 0.0)
        lam_raw_prev = nu0_c * np.exp(np.clip(-Geff_prev / max(KBEV * T_K, 1.0e-30), -700.0, 0.0))
        lam_prev = gammainc(m_hits, np.minimum(lam_raw_prev * tau_c, 1.0e12)) / tau_c
        lam_eff = lam_c.copy()
        lam_eff[pos] = _logmean_rate(lam_prev[pos], lam_c[pos])

        Bnew = Bclock[idx] + lam_eff * dt
        Bclock[idx] = Bnew
        Kprev[idx] = K
        fired_local = Bnew >= 1.0
        if np.any(fired_local):
            fi = idx[fired_local]
            Kc[fi] = K / 1.0e6
            Nem_fire[fi] = Nem[fi]
            sb = back_coeff * Nem[fi]
            rt = r0 + c_blunt * B_M * Nem[fi]
            st = np.minimum(K / np.sqrt(2.0 * np.pi * rt), sigma_cap)
            shield_fire[fi] = sb / np.maximum(st, 1.0)
            blunt_fire[fi] = rt / r0
            rho_f = rho0 + Nem[fi] / (L_pz**2)
            dGemb_fire[fi] = stored_coeff_eV * rho_f
            active[fi] = False

    rows = []
    for j, (_, r) in enumerate(design.iterrows()):
        if np.isfinite(Kc[j]):
            mode = "ductile" if ((shield_fire[j] > 0.25) or (blunt_fire[j] > 1.25)) else "brittle"
        else:
            mode = "no_fracture_in_window"
        rows.append({
            "surface_id": str(r["surface_id"]),
            "surface_index": int(r["surface_index"]),
            "context_id": str(context["context_id"]),
            "source_case_label": str(context["source_case_label"]),
            "source_response_regime": str(context["source_response_regime"]),
            "T_K": float(T_K),
            "Kc_first_MPa_sqrtm": float(Kc[j]) if np.isfinite(Kc[j]) else np.nan,
            "reached_monotonic_Kmax": bool(not np.isfinite(Kc[j])),
            "N_em_at_fire": float(Nem_fire[j]) if np.isfinite(Nem_fire[j]) else float(Nem[j]),
            "shield_frac_at_fire": float(shield_fire[j]) if np.isfinite(shield_fire[j]) else np.nan,
            "blunt_ratio_at_fire": float(blunt_fire[j]) if np.isfinite(blunt_fire[j]) else np.nan,
            "dG_emb_eV_at_fire": float(dGemb_fire[j]) if np.isfinite(dGemb_fire[j]) else np.nan,
            "mode_at_fire": mode,
            "monotonic_Kmax_MPa_sqrtm": float(Kmax_MPa_sqrtm),
            "monotonic_dK_MPa_sqrtm": float(dK_MPa_sqrtm),
            "Kdot_MPa_sqrtm_per_s": float(Kdot_MPa_sqrtm_per_s),
        })
    return rows


# ----------------------------------------------------------------------------
# temperature-curve descriptors (shared by Kc(T) and DeltaK_th(T))
# ----------------------------------------------------------------------------

def interp_at(x, y, xq):
    x = np.asarray(x, float); y = np.asarray(y, float)
    good = np.isfinite(x) & np.isfinite(y)
    if good.sum() < 2 or xq < np.min(x[good]) or xq > np.max(x[good]):
        return np.nan
    return float(np.interp(xq, x[good], y[good]))


def _classify_temperature_curve(Tg: np.ndarray, Kg: np.ndarray, prefix: str) -> dict:
    """Shared temperature-response taxonomy for Kc(T) and DeltaK_th(T).

    The classifier is IDENTICAL for both observables so the three-way
    class-association analysis compares like with like.
    """
    order = np.argsort(Tg); Tg, Kg = Tg[order], Kg[order]
    Kref = interp_at(Tg, Kg, 300.0)
    if not np.isfinite(Kref) or Kref <= 0:
        Kref = float(np.median(Kg))
    low = float(np.mean(Kg[:min(2, len(Kg))]))
    high = float(np.mean(Kg[-min(2, len(Kg)):]))
    ipeak = int(np.argmax(Kg))
    interior = 0 < ipeak < len(Kg) - 1
    endpoint_max = max(float(Kg[0]), float(Kg[-1]))
    peak_prom = (float(Kg[ipeak]) - endpoint_max) / max(Kref, 1e-12)
    post_drop = (float(Kg[ipeak]) - float(Kg[-1])) / max(Kref, 1e-12)
    total_range = (float(np.max(Kg)) - float(np.min(Kg))) / max(Kref, 1e-12)
    slope = np.diff(Kg) / np.diff(Tg) / max(Kref, 1e-12) * 100.0
    pos_frac = float(np.mean(slope > 0.01)) if len(slope) else np.nan
    neg_frac = float(np.mean(slope < -0.01)) if len(slope) else np.nan
    high_low_ratio = high / max(low, 1e-12)

    if interior and peak_prom >= 0.12 and post_drop >= 0.08:
        cls = "peak_shaped"
    elif high_low_ratio >= 1.50 and pos_frac >= 0.35:
        cls = "DBTT_like"
    elif total_range <= 0.25:
        cls = "weak_temperature"
    elif high_low_ratio <= 0.80 and neg_frac >= 0.50:
        cls = "ceramic_like"
    else:
        cls = "mixed_competition"

    # transition temperature from maximum positive normalized slope.
    if len(slope) and np.any(np.isfinite(slope)):
        imax = int(np.nanargmax(slope))
        Ttrans = 0.5 * (Tg[imax] + Tg[imax + 1])
    else:
        Ttrans = np.nan
    return {
        f"{prefix}_response_class": cls,
        f"{prefix}_ref_MPa_sqrtm": float(Kref),
        f"{prefix}_lowT_MPa_sqrtm": low,
        f"{prefix}_highT_MPa_sqrtm": high,
        f"{prefix}_high_over_low": high_low_ratio,
        f"{prefix}_total_range_norm": total_range,
        f"{prefix}_peak_T_K": float(Tg[ipeak]),
        f"{prefix}_peak_prominence_norm": float(peak_prom),
        f"{prefix}_post_peak_drop_norm": float(post_drop),
        f"{prefix}_max_positive_slope_norm_per_100K": float(np.max(slope)) if len(slope) else np.nan,
        f"{prefix}_max_negative_slope_norm_per_100K": float(np.min(slope)) if len(slope) else np.nan,
        f"{prefix}_transition_T_K": float(Ttrans),
        f"{prefix}_positive_slope_fraction": pos_frac,
        f"{prefix}_negative_slope_fraction": neg_frac,
    }


def fracture_curve_descriptors(points: pd.DataFrame) -> pd.DataFrame:
    rows = []
    if points.empty:
        return pd.DataFrame()
    for (sid, cid), g in points.groupby(["surface_id", "context_id"], sort=False):
        g = g.sort_values("T_K")
        T = g["T_K"].to_numpy(float)
        K = g["Kc_first_MPa_sqrtm"].to_numpy(float)
        good = np.isfinite(T) & np.isfinite(K)
        Tg, Kg = T[good], K[good]
        base = {
            "surface_id": sid,
            "context_id": cid,
            "source_case_label": str(g["source_case_label"].iloc[0]),
            "source_response_regime": str(g["source_response_regime"].iloc[0]),
            "n_temperature_points": int(len(T)),
            "n_resolved_Kc": int(good.sum()),
        }
        if good.sum() < 5:
            rows.append({**base, "fracture_response_class": "insufficient_or_unresolved"})
            continue
        d = _classify_temperature_curve(Tg, Kg, "Kc")
        d["fracture_response_class"] = d.pop("Kc_response_class")
        d["fracture_transition_T_K"] = d.pop("Kc_transition_T_K")
        d["fracture_positive_slope_fraction"] = d.pop("Kc_positive_slope_fraction")
        d["fracture_negative_slope_fraction"] = d.pop("Kc_negative_slope_fraction")
        rows.append({**base, **d})
    return pd.DataFrame(rows)


def fatigue_threshold_curve_descriptors(thr: pd.DataFrame, criterion: float) -> pd.DataFrame:
    """DeltaK_th(T) temperature-curve descriptors per (surface, context) at one criterion.

    Only bracketed thresholds contribute the point value (the log-interpolated
    estimate).  below/above_search_range states are censored bounds, retained in
    the threshold table but excluded from the curve fit rather than being
    assigned a value.
    """
    rows = []
    if thr.empty:
        return pd.DataFrame()
    sel = thr[np.isclose(thr["rate_criterion_m_per_cycle"].astype(float), criterion,
                         rtol=1e-6, atol=0.0)]
    for (sid, cid), g in sel.groupby(["surface_id", "fracture_context"], sort=False):
        g = g.sort_values("temperature_K")
        T = g["temperature_K"].to_numpy(float)
        K = g["DeltaK_th_MPa_sqrtm"].to_numpy(float)
        good = np.isfinite(T) & np.isfinite(K)
        Tg, Kg = T[good], K[good]
        n_cens_low = int((g["threshold_status"] == "below_search_range").sum())
        n_cens_high = int((g["threshold_status"] == "above_search_range").sum())
        base = {
            "surface_id": sid,
            "context_id": cid,
            "rate_criterion_m_per_cycle": float(criterion),
            "n_temperature_points": int(len(T)),
            "n_bracketed_thresholds": int(good.sum()),
            "n_below_search_range": n_cens_low,
            "n_above_search_range": n_cens_high,
        }
        if good.sum() < 5:
            rows.append({**base, "DKth_response_class": "insufficient_or_unresolved"})
            continue
        rows.append({**base, **_classify_temperature_curve(Tg, Kg, "DKth")})
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------
# stage 2: genuine V1 crack-growth da/dN(DeltaK,T) and rate-defined thresholds
# ----------------------------------------------------------------------------

class SurfaceEmitBarrier:
    """FrontEngine-facing adapter around the V5.6 per-surface emission barrier.

    The wrapped ScaledExpFloorBarrier (built through the SAME
    build_chain_from_namespace path used by the V5.6 S-N atlas, with
    energy/entropy scales 0.75/0.75 and the design-row EXP-floor overrides)
    is evaluated at the ACTUAL temperature: the design surface carries its own
    gT temperature slope, so no anchored-entropy wrapper is applied here.
    """

    def __init__(self, base, map_mod):
        self.base = base
        self._map_mod = map_mod
        self.mechanism = getattr(base, "mechanism", "surface_dislocation_emission")
        self.rate_prefactor = float(getattr(base, "rate_prefactor", 1.0e11))

    def deltaG_eV(self, sigma_Pa, T_K: float):
        return np.asarray(self.base.deltaG_eV(sigma_Pa, T_K), dtype=float)

    def rate(self, sigma_Pa, T_K: float):
        return self.base.rate(sigma_Pa, T_K)

    def G_barrier(self, sigma, T: float = 0.0, b: float = 2.74e-10):
        return self.deltaG_eV(sigma, T) * EV_TO_J

    def G_J(self, sigma, T: float = 0.0):
        return self.G_barrier(sigma, T)

    def S(self, sigma, T: float = 0.0):
        # -dG/dT at fixed sigma; for the ScaledExpFloorBarrier this is
        # -entropy_scale*gT (plus floor-branch terms) -- use the barrier's own
        # numerical diagnostic to stay faithful to any branch switching.
        sig = np.asarray(sigma, dtype=float)
        out = np.empty_like(sig)
        for idx, val in np.ndenumerate(sig):
            out[idx] = self.base.entropy_over_kB_numeric(float(val), float(T)) * KB
        return out

    def entropy_over_kB_numeric(self, sigma_Pa: float, T_K: float, dT: float = 1.0) -> float:
        return float(self.base.entropy_over_kB_numeric(sigma_Pa, T_K, dT))

    def diagnostics(self, sigma_Pa, T_K: float, b: float = 2.74e-10):
        return self._map_mod._compat_barrier_diagnostics(self, sigma_Pa, T_K, b)

    def __getattr__(self, name):
        return getattr(self.base, name)


def _design_namespace(row: pd.Series) -> SimpleNamespace:
    """Design-row EXP-floor overrides in the exact form build_chain_from_namespace consumes."""
    return SimpleNamespace(
        exp_system="W[100]",
        exp_G00_eV=float(row["exp_G00_eV"]),
        exp_gT_eV_per_K=float(row["exp_gT_eV_per_K"]),
        exp_sigc0_GPa=float(row["exp_sigc0_GPa"]),
        exp_sT_MPa_per_K=float(row["exp_sT_MPa_per_K"]),
        exp_Tref_K=float(row["exp_Tref_K"]),
        exp_a=float(row["exp_a"]),
        exp_n=float(row["exp_n"]),
        exp_floor_frac=float(row["exp_floor_frac"]),
        # V5.6 chain defaults (case-64-M1 fatigue scaling) left explicit:
        emit_energy_scale=0.75, emit_entropy_scale=0.75, emit_stress_scale=1.0,
        nu0_emit_pz=1.0e11,
    )


def _context_front_config(context: dict, da_m: float, map_mod) -> FrontConfig:
    """FrontConfig from a fracture-context row.

    Constants (r0, sigma_cap, nu0, m_hits, tau_c, beta_back, L_pz, rho0) match
    run_monotonic_batch and the corrected map runner's _front_config exactly, so
    stage-1 Kc(T) and stage-2 da/dN share one physics parameterization.
    """
    f = FrontConfig()
    f.r0 = 1.0e-6
    f.sigma_cap = 30.0e9
    f.m_hits = 3.0
    f.tau_c = 1.0e-6
    f.nu0_c = 1.0e12
    f.nu0_e = 1.0e11
    f.beta_back = 1.0
    f.L_pz = 1.0e-6
    f.rho0 = 5.0e12
    f.c_blunt = float(context.get("c_blunt", 1.0))
    f.v_emb_b3 = float(context.get("v_emb_b3", 500.0))
    f.wake_retain = float(context.get("wake_retain", 0.3))
    f.chi_shield = float(context.get("chi_shield", 0.0))
    f.emb_sat_frac = float(context.get("emb_sat_frac", 1.0))
    f.N_sat = map_mod._parse_float_or_inf(context.get("N_sat", float("inf")))
    f.recover_k = float(context.get("recover_k", 0.0))
    f.da = float(da_m)
    return f


def _context_cleavage_barrier(context: dict, map_mod):
    """Anchored EXP-floor cleavage barrier from a fracture-context row.

    Same construction as the corrected map runner's make_cleavage_barrier, with
    parameters from the context row and the anchored entropy taken from the
    context's cleave_S_kB at T_anchor = 300 K.  This reproduces _cleavage_G_eV
    (the stage-1 vectorized form) exactly.
    """
    cb = FractureBarrier()
    cb.barrier_kind = "exp_floor"
    cb.ef_G00_eV = float(context["cleave_G00_eV"])
    cb.ef_gT_eV_per_K = 0.0
    cb.ef_sigc0_Pa = float(context["cleave_sigc0_GPa"]) * 1e9
    cb.ef_sT_Pa_per_K = 0.0
    cb.ef_a = float(context["cleave_exp_a"])
    cb.ef_n = float(context["cleave_exp_n"])
    cb.ef_Tref_K = 300.0
    cb.ef_floor_frac = float(context["cleave_floor_frac"])
    cb.ef_floor_min_eV = 1e-4
    cb.ef_floor_max_frac = 0.95
    cb.ef_S_hs_kB = 0.0
    return map_mod.AnchoredCleavageBarrier(
        cb, S_kB=float(context.get("cleave_S_kB", 0.0)), T_anchor_K=300.0)


def make_surface_context_front(row: pd.Series, context: dict, args, map_mod):
    chain = build_chain_from_namespace(_design_namespace(row), b_m=B_M)
    emit = SurfaceEmitBarrier(chain.emit, map_mod)
    cbA = _context_cleavage_barrier(context, map_mod)
    fcfg = _context_front_config(context, args.da_m, map_mod)
    front = FrontEngine(fcfg, cbA, emit, G_PA, NU, B_M)
    cfg = FatigueControllerConfig(
        n_phase=int(args.n_phase),
        block_cycles=1e5,
        adaptive_cycles=True,
        max_block_cycles=float("inf"),
        min_block_cycles=1e-6,
        target_dB=float(args.target_dB),
        target_dN_store=float(args.target_dN_store),
        recovery_per_s=float(context.get("recover_k", 0.0)),
        N_sat=map_mod._parse_float_or_inf(context.get("N_sat", float("inf"))),
        storage_model=str(context.get("map_storage_model", "all_retained")),
        fixed_retained_fraction=1.0,
    )
    for name, value in [
        ("target_dN_emit", 0.20),
        ("target_dN_mobile", 0.20),
        ("target_dN_escape", float("inf")),
        ("target_dN_peierls", float("inf")),
        ("target_dN_taylor", float("inf")),
    ]:
        if hasattr(cfg, name):
            setattr(cfg, name, value)
    controller = FatigueCycleHazardController(cfg, emit, chain.peierls, chain.taylor)
    return front, controller


def run_fatigue_point(row: pd.Series, context: dict, T_K: float, DeltaK_MPa: float,
                      args, map_mod) -> Dict:
    """One (surface, context, T, DeltaK) crack-growth point via the established
    corrected-map block-stepping path (mod._map_cycle_step)."""
    front, controller = make_surface_context_front(row, context, args, map_mod)
    Kmax_MPa = float(DeltaK_MPa) / max(1.0 - args.R, 1e-12)
    wave = FatigueWaveform(
        Kmax=Kmax_MPa * 1e6,
        R=float(args.R),
        frequency_Hz=float(args.frequency_Hz),
        closure_clip=True,
    )
    ns = SimpleNamespace(
        max_block_cycles=float(args.max_block_cycles),
        min_block_cycles=float(args.min_block_cycles),
        target_dB=float(args.target_dB),
        target_state_fraction=float(args.target_state_fraction),
        saturation_tol_fraction=float(args.saturation_tol_fraction),
        target_dN_store_unbounded=float(args.target_dN_store_unbounded),
    )
    pred0 = controller.integrate_one_cycle(front, wave, T_K)
    log_ratio0 = math.log10(max(pred0.mu_emit, 1e-300) / max(pred0.mu_cleave, 1e-300))

    cycles_done = 0.0
    cycles_first = float("nan")
    blocks = 0
    for ib in range(int(args.max_blocks)):
        if cycles_done >= args.cycles_max or front.n_adv >= args.n_advances:
            break
        remaining = args.cycles_max - cycles_done
        last = map_mod._map_cycle_step(front, controller, wave, T_K, remaining, ns)
        cycles = float(last.get("cycles", 0.0))
        if not math.isfinite(cycles) or cycles <= 0.0:
            break
        cycles_done += cycles
        blocks = ib + 1
        if math.isnan(cycles_first) and int(last.get("n_fire", 0)) > 0:
            cycles_first = cycles_done

    measured = front.n_adv > 0 and cycles_done > 0.0
    da_dN = front.a_adv / cycles_done if measured else float("nan")
    ub = args.da_m / cycles_done if (not measured and cycles_done > 0.0) else float("nan")
    if measured:
        status = "measured"
    elif cycles_done >= 0.999 * args.cycles_max:
        status = "censored_cycle_horizon"
    else:
        status = "unresolved_block_limited"

    return {
        "surface_id": str(row["surface_id"]),
        "surface_index": int(row["surface_index"]),
        "fracture_context": str(context["context_id"]),
        "source_case_label": str(context.get("source_case_label", "")),
        "temperature_K": float(T_K),
        "DeltaK_MPa_sqrtm": float(DeltaK_MPa),
        "Kmax_MPa_sqrtm": float(Kmax_MPa),
        "R": float(args.R),
        "cycles_total": float(cycles_done),
        "cycles_to_first_fire": float(cycles_first),
        "a_adv_m": float(front.a_adv),
        "n_adv": int(front.n_adv),
        "da_dN_m_per_cycle": float(da_dN),
        "da_dN_upper_bound_m_per_cycle": float(ub),
        "status": status,
        "direct_lt_1_cycle": bool(math.isfinite(cycles_first) and cycles_first < 1.0),
        "N_em_final": float(front.N_em),
        "sigma_back_Pa": float(front.sigma_back()),
        "dG_emb_eV": float(front.dG_emb() / EV_TO_J),
        "B_final": float(front.B),
        "blocks_completed": int(blocks),
        "mu_emit_initial_per_cycle": float(pred0.mu_emit),
        "mu_cleave_initial_per_cycle": float(pred0.mu_cleave),
        "log10_mu_emit_over_cleave_initial": float(log_ratio0),
    }


def fatigue_growth_and_thresholds(design: pd.DataFrame, contexts: pd.DataFrame,
                                  fatigue_temps: Sequence[float], args,
                                  map_mod, adapt_mod, out: Path,
                                  analysis_only: bool) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Adaptive rate-defined DeltaK thresholds per (surface, context, T).

    Restartable at the RAW-point level: the persistent calculation key is
    (surface_id, fracture_context, temperature_K, DeltaK).  For each state the
    controller (1) reads all existing DeltaK points, (2) checks whether the
    criteria are bracketed to tolerance, (3) proposes only the next missing
    refinement point, (4) runs it through the V1 driver, (5) appends it to the
    raw CSV, and (6) repeats until converged.  Threshold summaries are always
    rebuilt from the complete raw table, so they stay consistent after resume.
    """
    points_path = out / "fatigue_growth_points_v5_7.csv"
    thr_path = out / "fatigue_thresholds_v5_7.csv"

    points = pd.read_csv(points_path).to_dict("records") if points_path.exists() else []

    def pkey(sid, cid, T, dk):
        return (str(sid), str(cid), round(float(T), 9), round(float(dk), 7))

    cache: Dict[Tuple, Dict] = {
        pkey(r["surface_id"], r["fracture_context"], r["temperature_K"], r["DeltaK_MPa_sqrtm"]): r
        for r in points
    }
    if cache:
        print(f"fatigue resume database: {len(cache)} completed (surface, context, T, DeltaK) keys")

    criteria = sorted(float(c) for c in args.rate_criteria)
    seeds = sorted(set(float(s) for s in args.DeltaK_seeds
                       if args.DeltaK_min <= float(s) <= args.DeltaK_max))

    def evaluate(row, context, T, dk) -> Dict:
        dk = min(max(float(dk), args.DeltaK_min), args.DeltaK_max)
        k = pkey(row["surface_id"], context["context_id"], T, dk)
        if k in cache:
            return cache[k]
        rec = run_fatigue_point(row, context, T, dk, args, map_mod)
        cache[k] = rec
        points.append(rec)
        append_rows(points_path, [rec])
        eff, src = adapt_mod.effective_rate(rec)
        print(f"  F {row['surface_id']} ctx={context['context_id']} T={T:g} "
              f"DK={dk:.5g} rate={eff:.3e} ({src})")
        return rec

    thresholds: List[Dict] = []
    if not analysis_only:
        n_groups = len(design) * len(contexts) * len(set(fatigue_temps))
        ig = 0
        for _, cser in contexts.iterrows():
            context = cser.to_dict()
            cid = str(context["context_id"])
            for _, row in design.iterrows():
                for T in sorted(set(float(x) for x in fatigue_temps)):
                    ig += 1
                    group_points = [r for r in points
                                    if str(r["surface_id"]) == str(row["surface_id"])
                                    and str(r["fracture_context"]) == cid
                                    and round(float(r["temperature_K"]), 9) == round(T, 9)]
                    need_seed = [dk for dk in seeds
                                 if pkey(row["surface_id"], cid, T, dk) not in cache]
                    if need_seed or not group_points:
                        print(f"\n=== fatigue group {ig}/{n_groups}: "
                              f"{row['surface_id']} ctx={cid} T={T:g} K ===")
                    for dk in seeds:
                        rec = evaluate(row, context, T, dk)
                        if rec not in group_points:
                            group_points.append(rec)

                    for crit in criteria:
                        for _ in range(int(args.max_refine_iters)):
                            a, b, vals = adapt_mod.locate_bracket(group_points, crit)
                            if a is None:
                                if vals and vals[0][1] >= crit and vals[0][0] > args.DeltaK_min * 1.001:
                                    new_dk = max(args.DeltaK_min, vals[0][0] / 2.0)
                                elif vals and vals[-1][1] < crit and vals[-1][0] < args.DeltaK_max / 1.001:
                                    new_dk = min(args.DeltaK_max, vals[-1][0] * 1.5)
                                else:
                                    break
                            else:
                                est = adapt_mod.crossing_estimate(a, b, crit)
                                width = b[0] - a[0]
                                if width <= max(args.threshold_abs_tol,
                                                args.threshold_rel_tol * max(est, 1e-12)):
                                    break
                                new_dk = est
                            existing = np.array([float(r["DeltaK_MPa_sqrtm"]) for r in group_points])
                            if len(existing) and np.min(np.abs(existing - new_dk)) < 1e-7:
                                break
                            group_points.append(evaluate(row, context, T, new_dk))

    # Threshold summary is ALWAYS rebuilt from the full raw table (including
    # under --analysis-only), one row per (surface, context, T, criterion).
    all_points = pd.DataFrame(points)
    if not all_points.empty:
        for (sid, cid, T), g in all_points.groupby(
                ["surface_id", "fracture_context", "temperature_K"], sort=False):
            grp = g.to_dict("records")
            for crit in criteria:
                tr = adapt_mod.threshold_record(grp, crit, args.threshold_abs_tol,
                                                args.threshold_rel_tol)
                thresholds.append({
                    "surface_id": str(sid),
                    "fracture_context": str(cid),
                    "temperature_K": float(T),
                    "rate_criterion_m_per_cycle": float(crit),
                    "DeltaK_th_MPa_sqrtm": tr.get("DeltaK_threshold_estimate_MPa_sqrtm", np.nan),
                    "DeltaK_lo_MPa_sqrtm": tr.get("DeltaK_threshold_lower_MPa_sqrtm", np.nan),
                    "DeltaK_hi_MPa_sqrtm": tr.get("DeltaK_threshold_upper_MPa_sqrtm", np.nan),
                    "threshold_status": tr.get("threshold_class", "unresolved"),
                    "threshold_width_MPa_sqrtm": tr.get("threshold_width_MPa_sqrtm", np.nan),
                    "converged": bool(tr.get("converged", False)),
                    "lower_rate_m_per_cycle": tr.get("lower_rate_m_per_cycle", np.nan),
                    "upper_rate_m_per_cycle": tr.get("upper_rate_m_per_cycle", np.nan),
                    "lower_source": tr.get("lower_source", ""),
                    "upper_source": tr.get("upper_source", ""),
                    "n_DeltaK_evaluations": int(tr.get("n_evaluated_DeltaK", 0)),
                    "threshold_definition": "rate_defined_da_dN_crossing_of_simulated_V1_growth_curve",
                })
    thr = pd.DataFrame(thresholds)
    thr.to_csv(thr_path, index=False)
    return all_points, thr


# ----------------------------------------------------------------------------
# association / correlation analysis
# ----------------------------------------------------------------------------

def cramers_v_table(df: pd.DataFrame, row_col: str, col_col: str,
                    context_id: str, family: str) -> Tuple[dict, List[dict]]:
    tab = pd.crosstab(df[row_col], df[col_col])
    if tab.shape[0] < 2 or tab.shape[1] < 2:
        return ({"analysis_family": family, "context_id": context_id,
                 "cramers_V": np.nan, "chi2_p": np.nan, "n": int(tab.values.sum())}, [])
    chi2, p, _, exp = chi2_contingency(tab, correction=False)
    n = tab.values.sum(); k = min(tab.shape)
    V = math.sqrt(chi2 / max(n * (k - 1), 1.0))
    cells = []
    obs = tab.to_numpy(float)
    z = (obs - exp) / np.sqrt(np.maximum(exp, 1e-30))
    for i, rv in enumerate(tab.index):
        for j, cv in enumerate(tab.columns):
            cells.append({"analysis_family": family, "context_id": context_id,
                          "row_class": str(rv), "column_class": str(cv),
                          "observed": int(obs[i, j]), "expected": float(exp[i, j]),
                          "standardized_residual": float(z[i, j])})
    return ({"analysis_family": family, "context_id": context_id,
             "cramers_V": float(V), "chi2_p": float(p), "n": int(n),
             "n_row_classes": tab.shape[0], "n_column_classes": tab.shape[1]}, cells)


def build_fourway(frac_desc: pd.DataFrame, thr_desc: pd.DataFrame,
                  phenotype: pd.DataFrame, primary_criterion: float
                  ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if frac_desc.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    merged = frac_desc.merge(phenotype, on="surface_id", how="left")
    if not thr_desc.empty:
        td = thr_desc[np.isclose(thr_desc["rate_criterion_m_per_cycle"].astype(float),
                                 primary_criterion, rtol=1e-6, atol=0.0)]
        keep = ["surface_id", "context_id", "DKth_response_class",
                "DKth_ref_MPa_sqrtm", "DKth_high_over_low", "DKth_transition_T_K",
                "n_bracketed_thresholds"]
        keep = [c for c in keep if c in td.columns]
        merged = merged.merge(td[keep], on=["surface_id", "context_id"], how="left")
    stats, cells = [], []
    for cid, g in merged.groupby("context_id", sort=False):
        pairs = [
            ("strength_response_class", "fracture_response_class", "strength_vs_fracture_class"),
            ("fracture_response_class", "fatigue_temperature_pattern", "fracture_vs_fatigue_temperature_phenotype"),
            ("strength_response_class", "fatigue_temperature_pattern", "strength_vs_fatigue_temperature_phenotype"),
            ("fracture_response_class", "DKth_response_class", "fracture_vs_growth_threshold_class"),
            ("strength_response_class", "DKth_response_class", "strength_vs_growth_threshold_class"),
            ("DKth_response_class", "fatigue_temperature_pattern", "growth_threshold_vs_SN_temperature_phenotype"),
        ]
        for rcol, ccol, fam in pairs:
            if rcol not in g or ccol not in g:
                continue
            st, ce = cramers_v_table(g.dropna(subset=[rcol, ccol]), rcol, ccol, cid, fam)
            stats.append(st); cells.extend(ce)
    return merged, pd.DataFrame(stats), pd.DataFrame(cells)


def scalar_threshold_correlations(thr: pd.DataFrame, mono: pd.DataFrame,
                                  surface_desc: pd.DataFrame,
                                  primary_criterion: float) -> pd.DataFrame:
    """Spearman correlations linking DeltaK_th to Kc, S-N endurance level, and
    the strength anomaly, per fracture context (bracketed thresholds only)."""
    rows: List[dict] = []
    if thr.empty:
        return pd.DataFrame()
    t = thr[np.isclose(thr["rate_criterion_m_per_cycle"].astype(float),
                       primary_criterion, rtol=1e-6, atol=0.0)].copy()
    t = t[t["threshold_status"] == "bracketed"]
    if t.empty:
        return pd.DataFrame()

    def add(cid, family, x, y, n_note=""):
        x = np.asarray(x, float); y = np.asarray(y, float)
        good = np.isfinite(x) & np.isfinite(y)
        if good.sum() < 5:
            rows.append({"context_id": cid, "analysis_family": family,
                         "spearman_rho": np.nan, "p_value": np.nan,
                         "n": int(good.sum()), "note": "insufficient_pairs"})
            return
        rho, p = spearmanr(x[good], y[good])
        rows.append({"context_id": cid, "analysis_family": family,
                     "spearman_rho": float(rho), "p_value": float(p),
                     "n": int(good.sum()), "note": n_note})

    # 1) DeltaK_th vs monotonic Kc at matched (surface, context, T).
    if not mono.empty:
        m = mono.rename(columns={"context_id": "fracture_context", "T_K": "temperature_K"})
        j = t.merge(m[["surface_id", "fracture_context", "temperature_K",
                       "Kc_first_MPa_sqrtm"]],
                    on=["surface_id", "fracture_context", "temperature_K"], how="inner")
        for cid, g in j.groupby("fracture_context", sort=False):
            add(cid, "DKth_vs_Kc_matched_T",
                g["DeltaK_th_MPa_sqrtm"], g["Kc_first_MPa_sqrtm"])

    # 2) DeltaK_th vs blunt-feature S-N endurance level (N=1e12 bracket midpoint)
    #    at matched fatigue temperature, both shielded and no-shield S-N cases.
    if surface_desc is not None and not surface_desc.empty:
        sd = surface_desc.copy()
        sd["temperature_K"] = sd["fatigue_T_K"].astype(float)
        for case_prefix in ["shielded", "no_shield"]:
            mid_col = f"{case_prefix}_sn_sigma_N12_mid_norm"
            if mid_col not in sd.columns:
                continue
            sn = sd[["surface_id", "temperature_K", mid_col, "sigma_ref_nominal_MPa"]].copy()
            sn["sigma_N12_mid_MPa"] = sn[mid_col].astype(float) * sn["sigma_ref_nominal_MPa"].astype(float)
            j = t.merge(sn[["surface_id", "temperature_K", "sigma_N12_mid_MPa"]],
                        on=["surface_id", "temperature_K"], how="inner")
            for cid, g in j.groupby("fracture_context", sort=False):
                add(cid, f"DKth_vs_SN_endurance_level_{case_prefix}",
                    g["DeltaK_th_MPa_sqrtm"], g["sigma_N12_mid_MPa"])

        # 3) DeltaK_th near 300 K vs strength-anomaly gain (one value per surface).
        if "strength_anomaly_gain_frac" in sd.columns:
            gain = sd.groupby("surface_id", sort=False)["strength_anomaly_gain_frac"].first()
            t300 = t.loc[(t["temperature_K"] - 300.0).abs() < 1.0]
            if t300.empty:
                # fall back to the fatigue temperature closest to 300 K
                Tn = t["temperature_K"].astype(float)
                Tstar = float(Tn.iloc[(Tn - 300.0).abs().argmin()])
                t300 = t.loc[(t["temperature_K"] - Tstar).abs() < 1.0]
                note = f"nearest_T_{Tstar:g}K"
            else:
                note = "T_300K"
            for cid, g in t300.groupby("fracture_context", sort=False):
                gg = g.merge(gain.rename("anomaly_gain"), on="surface_id", how="inner")
                add(cid, "DKth_nearRT_vs_strength_anomaly_gain",
                    gg["DeltaK_th_MPa_sqrtm"], gg["anomaly_gain"], note)
    return pd.DataFrame(rows)


def make_overview_plot(stats: pd.DataFrame, out: Path) -> None:
    if stats.empty:
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    piv = stats.pivot(index="context_id", columns="analysis_family", values="cramers_V")
    ax = piv.plot(kind="bar", figsize=(13, 5.5))
    ax.set_ylabel("Cramér's V")
    ax.set_title("Response-class association strength by fracture context "
                 "(strength / fracture / S-N phenotype / rate-defined DeltaK_th)")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(fontsize=7)
    plt.tight_layout()
    plt.savefig(out / "fourway_class_association_overview_v5_7.png", dpi=220)
    plt.close()


def make_threshold_link_plot(thr: pd.DataFrame, mono: pd.DataFrame, out: Path,
                             primary_criterion: float) -> None:
    if thr.empty or mono.empty:
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    t = thr[(thr["threshold_status"] == "bracketed")
            & np.isclose(thr["rate_criterion_m_per_cycle"].astype(float),
                         primary_criterion, rtol=1e-6, atol=0.0)]
    m = mono.rename(columns={"context_id": "fracture_context", "T_K": "temperature_K"})
    j = t.merge(m[["surface_id", "fracture_context", "temperature_K",
                   "Kc_first_MPa_sqrtm"]],
                on=["surface_id", "fracture_context", "temperature_K"], how="inner")
    if j.empty:
        return
    fig, ax = plt.subplots(figsize=(7.0, 6.0))
    for cid, g in j.groupby("fracture_context", sort=False):
        sc = ax.scatter(g["Kc_first_MPa_sqrtm"], g["DeltaK_th_MPa_sqrtm"],
                        s=14, alpha=0.6, label=str(cid),
                        c=g["temperature_K"], cmap="viridis")
    cb = fig.colorbar(sc, ax=ax)
    cb.set_label("T [K]")
    ax.set_xlabel(r"$K_c$ (monotonic first fire) [MPa$\sqrt{m}$]")
    ax.set_ylabel(rf"$\Delta K_{{th}}$ at $da/dN={primary_criterion:g}$ m/cycle "
                  rf"[MPa$\sqrt{{m}}$]")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.grid(alpha=0.25, which="both")
    ax.legend(fontsize=7, title="fracture context")
    ax.set_title("Rate-defined fatigue threshold vs monotonic toughness (matched T)")
    plt.tight_layout()
    plt.savefig(out / "DKth_vs_Kc_matched_T_v5_7.png", dpi=220)
    plt.close()


# ----------------------------------------------------------------------------
# surface selection
# ----------------------------------------------------------------------------

def select_fatigue_surfaces(design: pd.DataFrame, phenotype: pd.DataFrame, args
                            ) -> pd.DataFrame:
    """Surfaces receiving the (serial, per-front) fatigue-growth calculation.

    Priority: explicit --fatigue-surface-ids, then a seeded stratified
    representative draw of --n-fatigue-surfaces balanced over
    (strength_response_class x fatigue_temperature_pattern), then all selected
    monotonic surfaces.  Stage 2 is per-front and NOT vectorized, so the
    representative draw is the practical default at V5.6 atlas scale.
    """
    if args.fatigue_surface_ids:
        wanted = set()
        for token in args.fatigue_surface_ids:
            wanted.update(x.strip() for x in str(token).split(",") if x.strip())
        sel = design[design.surface_id.astype(str).isin(wanted)]
        missing = wanted - set(sel.surface_id.astype(str))
        if missing:
            raise SystemExit(f"Unknown fatigue surface ids: {sorted(missing)}")
        return sel
    if args.n_fatigue_surfaces is not None and args.n_fatigue_surfaces < len(design):
        rng = np.random.default_rng(int(args.fatigue_surface_seed))
        ph = phenotype[["surface_id", "strength_response_class",
                        "fatigue_temperature_pattern"]].copy()
        d = design.merge(ph, on="surface_id", how="left")
        d["stratum"] = (d["strength_response_class"].astype(str) + "|"
                        + d["fatigue_temperature_pattern"].astype(str))
        strata = list(d.groupby("stratum", sort=True))
        n_target = int(args.n_fatigue_surfaces)
        picked: List[pd.DataFrame] = []
        remaining = n_target
        # Round-robin one per stratum until the budget is spent.
        pools = {name: g.sample(frac=1.0, random_state=int(args.fatigue_surface_seed))
                 for name, g in strata}
        while remaining > 0 and any(len(p) for p in pools.values()):
            for name in sorted(pools):
                if remaining <= 0:
                    break
                p = pools[name]
                if len(p):
                    picked.append(p.iloc[[0]])
                    pools[name] = p.iloc[1:]
                    remaining -= 1
        sel = pd.concat(picked, ignore_index=False)
        return design.loc[design.surface_id.isin(sel.surface_id)]
    return design


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--atlas-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--context-table", default="fracture_contexts_v5_7.csv")
    ap.add_argument("--context-filter", nargs="*", default=[])
    ap.add_argument("--temperatures", nargs="+", type=float,
                    default=[100, 200, 300, 400, 500, 600, 700, 800, 900])
    ap.add_argument("--surface-start", type=int, default=0)
    ap.add_argument("--surface-stop", type=int, default=None)
    ap.add_argument("--monotonic-Kmax-MPa", type=float, default=40.0)
    ap.add_argument("--monotonic-dK-MPa", type=float, default=0.10)
    ap.add_argument("--Kdot-MPa-sqrtm-per-s", type=float, default=0.005)
    ap.add_argument("--skip-monotonic", action="store_true")
    ap.add_argument("--analysis-only", action="store_true")
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--overwrite", action="store_true")

    # --- stage 2: genuine V1 crack-growth thresholds ---
    ap.add_argument("--skip-fatigue-growth", action="store_true")
    ap.add_argument("--map-runner", default="run_v1_two_barrier_dbtt_fatigue_map_corrected.py",
                    help="Corrected two-barrier map runner supplying the established "
                         "1-D fatigue block-stepping path (_map_cycle_step).")
    ap.add_argument("--adaptive-study", default="run_adaptive_two_barrier_threshold_study.py",
                    help="Adaptive threshold study supplying the established bracketing "
                         "numerics (locate_bracket / crossing_estimate / threshold_record).")
    ap.add_argument("--fatigue-temperatures", nargs="+", type=float,
                    default=[100, 300, 500, 700, 900])
    ap.add_argument("--fatigue-surface-ids", nargs="*", default=[])
    ap.add_argument("--n-fatigue-surfaces", type=int, default=96,
                    help="Stratified representative draw over (strength class x S-N "
                         "phenotype).  Set >= n selected surfaces to run all.")
    ap.add_argument("--fatigue-surface-seed", type=int, default=42)
    ap.add_argument("--rate-criteria", nargs="+", type=float, default=[1e-10, 1e-12],
                    help="da/dN crossing criteria [m/cycle]: 1e-10 primary, 1e-12 sensitivity.")
    ap.add_argument("--primary-rate-criterion", type=float, default=1e-10)
    ap.add_argument("--DeltaK-seeds", nargs="+", type=float,
                    default=[0.05, 0.10, 0.20, 0.40, 0.80, 1.60, 3.20, 6.40, 12.80])
    ap.add_argument("--DeltaK-min", type=float, default=0.025)
    ap.add_argument("--DeltaK-max", type=float, default=20.0)
    ap.add_argument("--threshold-abs-tol", type=float, default=0.05)
    ap.add_argument("--threshold-rel-tol", type=float, default=0.03)
    ap.add_argument("--max-refine-iters", type=int, default=10)

    # Physics/controller settings copied from the adaptive two-barrier study defaults.
    ap.add_argument("--R", type=float, default=0.1)
    ap.add_argument("--frequency-Hz", type=float, default=1000.0)
    ap.add_argument("--cycles-max", type=float, default=2e14)
    ap.add_argument("--max-blocks", type=int, default=10000)
    ap.add_argument("--max-block-cycles", type=float, default=float("inf"))
    ap.add_argument("--min-block-cycles", type=float, default=1e-10)
    ap.add_argument("--target-state-fraction", type=float, default=0.01)
    ap.add_argument("--saturation-tol-fraction", type=float, default=1e-4)
    ap.add_argument("--target-dN-store-unbounded", type=float, default=5.0)
    ap.add_argument("--n-advances", type=int, default=5)
    ap.add_argument("--da-m", type=float, default=20e-6)
    ap.add_argument("--n-phase", type=int, default=96)
    ap.add_argument("--target-dB", type=float, default=0.02)
    ap.add_argument("--target-dN-store", type=float, default=0.01)
    args = ap.parse_args()

    atlas = Path(args.atlas_dir)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    design_path = find_one(atlas, ["independent_exp_floor_design_v5_6.csv", "independent_exp_floor_design_v5_*.csv"])
    phenotype_path = find_one(atlas, ["surface_phenotype_summary_v5_6.csv", "surface_phenotype_summary_v5_*.csv"])
    sn_path = find_one(atlas, ["sn_initiation_points_multiT_paired.csv"])
    summary_path = find_one(atlas, ["surface_temperature_descriptors_v5_6.csv", "surface_temperature_descriptors_v5_*.csv"])
    design = pd.read_csv(design_path).sort_values("surface_index")
    design = design[(design.surface_index >= args.surface_start)]
    if args.surface_stop is not None:
        design = design[design.surface_index < args.surface_stop]
    contexts = pd.read_csv(args.context_table)
    if args.context_filter:
        wanted = set()
        for token in args.context_filter:
            wanted.update(x.strip() for x in str(token).split(",") if x.strip())
        contexts = contexts[contexts.context_id.astype(str).isin(wanted)]
        missing = wanted - set(contexts.context_id.astype(str))
        if missing:
            raise SystemExit(f"Unknown context ids: {sorted(missing)}")
    if contexts.empty:
        raise SystemExit("No fracture contexts selected")
    phenotype = pd.read_csv(phenotype_path)

    # ---------------- stage 1: monotonic Kc(T) ----------------
    mono_path = out / "fracture_monotonic_points_v5_7.csv"
    if mono_path.exists() and not args.resume and not args.analysis_only:
        if args.overwrite:
            mono_path.unlink()
        else:
            raise SystemExit(f"{mono_path} exists. Use --resume to extend/continue, or --overwrite for a fresh extension run.")

    existing = pd.read_csv(mono_path) if mono_path.exists() else pd.DataFrame()
    done = set()
    if not existing.empty:
        done = {(str(r.surface_id), str(r.context_id), round(float(r.T_K), 9))
                for r in existing.itertuples(index=False)}
        print(f"resume database: {len(done)} completed monotonic task keys")

    if not args.analysis_only and not args.skip_monotonic:
        total = len(design) * len(contexts) * len(set(args.temperatures))
        for _, cser in contexts.iterrows():
            context = cser.to_dict()
            cid = str(context["context_id"])
            for T in sorted(set(float(x) for x in args.temperatures)):
                mask = [((str(r.surface_id), cid, round(T, 9)) not in done)
                        for r in design.itertuples(index=False)]
                missing_design = design.loc[np.asarray(mask, bool)]
                if missing_design.empty:
                    print(f"skip complete context={cid} T={T:g} K")
                    continue
                rows = run_monotonic_batch(
                    missing_design, context, T,
                    args.monotonic_Kmax_MPa,
                    args.monotonic_dK_MPa,
                    args.Kdot_MPa_sqrtm_per_s,
                )
                append_rows(mono_path, rows)
                for r in rows:
                    done.add((str(r["surface_id"]), cid, round(T, 9)))
                print(f"fracture context={cid} T={T:g} K: +{len(rows)} rows; total keys={len(done)}/{total}")

    mono = pd.read_csv(mono_path) if mono_path.exists() else pd.DataFrame()
    frac_desc = fracture_curve_descriptors(mono)
    frac_desc.to_csv(out / "fracture_curve_descriptors_v5_7.csv", index=False)

    # ---------------- stage 2: rate-defined DeltaK_th(T) ----------------
    thr = pd.DataFrame()
    thr_desc = pd.DataFrame()
    growth_points = pd.DataFrame()
    if not args.skip_fatigue_growth:
        map_mod = load_map_runner(Path(args.map_runner))
        adapt_mod = load_adaptive_study(Path(args.adaptive_study))
        fatigue_design = select_fatigue_surfaces(design, phenotype, args)
        print(f"fatigue-growth stage: {len(fatigue_design)} surfaces x "
              f"{len(contexts)} contexts x {len(set(args.fatigue_temperatures))} temperatures")
        growth_points, thr = fatigue_growth_and_thresholds(
            fatigue_design, contexts, args.fatigue_temperatures, args,
            map_mod, adapt_mod, out, analysis_only=args.analysis_only,
        )
        desc_frames = [fatigue_threshold_curve_descriptors(thr, c)
                       for c in sorted(float(x) for x in args.rate_criteria)]
        desc_frames = [d for d in desc_frames if not d.empty]
        thr_desc = pd.concat(desc_frames, ignore_index=True) if desc_frames else pd.DataFrame()
        thr_desc.to_csv(out / "fatigue_threshold_curve_descriptors_v5_7.csv", index=False)

    # ---------------- association / correlation analysis ----------------
    fourway, stats, cells = build_fourway(frac_desc, thr_desc, phenotype,
                                          args.primary_rate_criterion)
    fourway.to_csv(out / "fourway_phenotype_summary_v5_7.csv", index=False)
    stats.to_csv(out / "fourway_class_association_statistics_v5_7.csv", index=False)
    cells.to_csv(out / "fourway_class_association_cells_v5_7.csv", index=False)
    make_overview_plot(stats, out)

    surface_desc = pd.read_csv(summary_path)
    corr = scalar_threshold_correlations(thr, mono, surface_desc,
                                         args.primary_rate_criterion)
    corr.to_csv(out / "fatigue_threshold_correlations_v5_7.csv", index=False)
    make_threshold_link_plot(thr, mono, out, args.primary_rate_criterion)

    config = vars(args).copy()
    config["atlas_dir_resolved"] = str(atlas.resolve())
    config["design_csv"] = str(design_path)
    config["phenotype_csv"] = str(phenotype_path)
    config["sn_csv"] = str(sn_path)
    config["surface_descriptor_csv"] = str(summary_path)
    config["n_selected_surfaces"] = int(len(design))
    config["n_contexts"] = int(len(contexts))
    with (out / "extension_config_v5_7.json").open("w") as f:
        json.dump(config, f, indent=2, default=str)
    print(f"Wrote V5.7 extension outputs to {out}")


if __name__ == "__main__":
    main()
