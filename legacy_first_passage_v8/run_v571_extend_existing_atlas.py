"""V5.7.1 incremental three-way extension of an existing V5.6 atlas.

The driver reuses the existing strength-temperature and blunt-feature S-N atlas
without recomputing either. It adds two new projections on the exact same
barrier surfaces:

1. monotonic K-ramp fracture response Kc(T) under several fracture contexts;
2. rate-defined fatigue crack-growth thresholds DeltaK_th(T) evaluated by the
   existing one-dimensional sharp-front fatigue-growth driver.

No life-to-DeltaK or notch-equivalent conversion is used. The fatigue threshold
is defined only by the simulated da/dN(DeltaK) crossing of a specified crack-
growth-rate criterion.

The workflow is appendable and restartable. Monotonic task keys are
(surface_id, context_id, T_K). Fatigue-growth point keys are
(surface_id, context_id, T_K, DeltaK). Existing raw points are reused by the
adaptive threshold search, so interrupted refinement resumes from the current
bracket rather than restarting the state.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from types import SimpleNamespace
import importlib
import shutil
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.special import gammainc
from scipy.stats import chi2_contingency

KBEV = 8.617333262145e-5
EV_TO_J = 1.602176634e-19

# W-like elastic constants used by the existing V1 sharp-front model.
E_PA = 410.0e9
NU = 0.28
G_PA = E_PA / (2.0 * (1.0 + NU))
B_M = 2.74e-10


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


def interp_at(x, y, xq):
    x = np.asarray(x, float); y = np.asarray(y, float)
    good = np.isfinite(x) & np.isfinite(y)
    if good.sum() < 2 or xq < np.min(x[good]) or xq > np.max(x[good]):
        return np.nan
    return float(np.interp(xq, x[good], y[good]))


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
        rows.append({
            **base,
            "fracture_response_class": cls,
            "Kc_ref_MPa_sqrtm": float(Kref),
            "Kc_lowT_MPa_sqrtm": low,
            "Kc_highT_MPa_sqrtm": high,
            "Kc_high_over_low": high_low_ratio,
            "Kc_total_range_norm": total_range,
            "Kc_peak_T_K": float(Tg[ipeak]),
            "Kc_peak_prominence_norm": float(peak_prom),
            "Kc_post_peak_drop_norm": float(post_drop),
            "Kc_max_positive_slope_norm_per_100K": float(np.max(slope)) if len(slope) else np.nan,
            "Kc_max_negative_slope_norm_per_100K": float(np.min(slope)) if len(slope) else np.nan,
            "fracture_transition_T_K": float(Ttrans),
            "fracture_positive_slope_fraction": pos_frac,
            "fracture_negative_slope_fraction": neg_frac,
        })
    return pd.DataFrame(rows)


# -----------------------------------------------------------------------------
# Existing V1 fatigue-growth driver adapter
# -----------------------------------------------------------------------------

def _load_existing_v1_driver(module_name: str):
    """Import the existing corrected V1 K-controlled fatigue-growth driver.

    The adapter intentionally uses the existing FrontEngine/controller point
    evaluator instead of creating a new fatigue-growth surrogate. The only
    adaptation is replacing the plastic barrier family with the selected V5.6
    EXP-floor surface while retaining the native one-dimensional crack-growth
    state update, renewal semantics, waveform integration, and da/dN reporting.
    """
    try:
        mod = importlib.import_module(module_name)
    except Exception as exc:
        raise RuntimeError(
            f"Could not import existing V1 fatigue-growth driver '{module_name}'. "
            "Copy run_v1_two_barrier_dbtt_fatigue_map_corrected.py into pwd "
            "before running V5.7.1."
        ) from exc
    required = ["run_fatigue", "ExpFloorBarrierParams", "ScaledExpFloorBarrier"]
    missing = [x for x in required if not hasattr(mod, x)]
    if missing:
        raise RuntimeError(f"Existing V1 driver is missing required symbols: {missing}")
    return mod


class _DirectScaledBarrierAdapter:
    """Expose a ScaledExpFloorBarrier through the FrontEngine barrier API."""

    def __init__(self, scaled, ev_to_j: float):
        self.scaled = scaled
        self.ev_to_j = float(ev_to_j)
        self.mechanism = getattr(scaled, "mechanism", "scaled_exp_floor")
        self.rate_prefactor = float(getattr(scaled, "rate_prefactor", 1.0e11))

    def deltaG_eV(self, sigma_Pa, T_K: float):
        return self.scaled.deltaG_eV(sigma_Pa, T_K)

    def rate(self, sigma_Pa, T_K: float):
        return self.scaled.rate(sigma_Pa, T_K)

    def G_barrier(self, sigma, T: float = 0.0, b: float = 2.74e-10):
        return np.asarray(self.deltaG_eV(sigma, T), dtype=float) * self.ev_to_j

    def G_J(self, sigma, T: float = 0.0):
        return self.G_barrier(sigma, T)

    def entropy_over_kB_numeric(self, sigma_Pa: float, T_K: float, dT: float = 1.0):
        return self.scaled.entropy_over_kB_numeric(sigma_Pa, T_K, dT)

    def S(self, sigma, T: float = 0.0):
        s = np.asarray(sigma, dtype=float)
        vals = np.empty_like(s, dtype=float)
        flat = vals.reshape(-1)
        sf = s.reshape(-1)
        for i, sig in enumerate(sf):
            flat[i] = self.entropy_over_kB_numeric(float(sig), float(T)) * 1.380649e-23
        return vals

    def diagnostics(self, sigma_Pa, T_K: float, b: float = 2.74e-10):
        # The current FrontEngine does not require this for lambda_emit, but
        # newer trees may audit it. Reuse the compatibility helper when present.
        sigma = np.asarray(sigma_Pa, dtype=float)
        G = self.G_barrier(sigma, T_K, b)
        return {
            "G_J": G,
            "G_eV": G / self.ev_to_j,
            "DeltaG_J": G,
            "DeltaG_eV": G / self.ev_to_j,
        }

    def __getattr__(self, name):
        return getattr(self.scaled, name)


def _surface_plastic_barriers(driver, surface: pd.Series):
    """Build the V5.6 emission/Peierls/Taylor family for the V1 driver.

    Scaling matches the V5.6 S-N atlas: emission 0.75, Peierls 0.00375,
    Taylor 0.015, with all three mechanisms derived from the same sampled
    EXP-floor base surface.
    """
    base = driver.ExpFloorBarrierParams(
        name=str(surface["surface_id"]),
        G00_eV=float(surface["exp_G00_eV"]),
        gT_eV_per_K=float(surface["exp_gT_eV_per_K"]),
        sigc0_Pa=float(surface["exp_sigc0_GPa"]) * 1.0e9,
        sT_Pa_per_K=float(surface["exp_sT_MPa_per_K"]) * 1.0e6,
        Tref_K=float(surface["exp_Tref_K"]),
        a=float(surface["exp_a"]),
        n=float(surface["exp_n"]),
        Gfloor_fraction=float(surface["exp_floor_frac"]),
        Gfloor_min_eV=1.0e-4,
        Gfloor_max_fraction=0.95,
    )
    emit = driver.ScaledExpFloorBarrier(
        base=base,
        mechanism="crack_tip_dislocation_emission",
        energy_scale=0.75,
        entropy_scale=0.75,
        stress_scale=1.0,
        rate_prefactor=1.0e11,
    )
    peierls = driver.ScaledExpFloorBarrier(
        base=base,
        mechanism="peierls_glide_escape",
        energy_scale=0.00375,
        entropy_scale=0.00375,
        stress_scale=1.0,
        rate_prefactor=1.0e11,
    )
    taylor = driver.ScaledExpFloorBarrier(
        base=base,
        mechanism="taylor_junction_depinning_escape",
        energy_scale=0.015,
        entropy_scale=0.015,
        stress_scale=1.0,
        rate_prefactor=1.0e11,
    )
    ev_to_j = float(getattr(driver, "EV_TO_J", 1.602176634e-19))
    return (
        _DirectScaledBarrierAdapter(emit, ev_to_j),
        _DirectScaledBarrierAdapter(peierls, ev_to_j),
        _DirectScaledBarrierAdapter(taylor, ev_to_j),
    )


def _context_case(context: dict) -> dict:
    """Convert a V5.7 fracture-context row into the native V1 case dictionary."""
    return {
        "case_label": str(context["context_id"]),
        "response_regime": str(context.get("source_response_regime", "")),
        "cleave_G00_eV": float(context["cleave_G00_eV"]),
        "cleave_sigc0_GPa": float(context["cleave_sigc0_GPa"]),
        "cleave_exp_a": float(context["cleave_exp_a"]),
        "cleave_exp_n": float(context["cleave_exp_n"]),
        "cleave_floor_frac": float(context["cleave_floor_frac"]),
        "chi_shield": float(context["chi_shield"]),
        "N_sat": context["N_sat"],
        "emb_sat_frac": float(context["emb_sat_frac"]),
        "c_blunt": float(context["c_blunt"]),
        "v_emb_b3": float(context["v_emb_b3"]),
        "wake_retain": float(context["wake_retain"]),
        "recover_k": float(context["recover_k"]),
        "map_storage_model": str(context.get("map_storage_model", "all_retained")),
    }


def _threshold_driver_args(args) -> SimpleNamespace:
    """Namespace expected by the existing corrected V1 run_fatigue function."""
    return SimpleNamespace(
        T_anchor_K=300.0,
        R=float(args.R),
        frequency_Hz=float(args.frequency_Hz),
        cycles_max=float(args.threshold_cycles_max),
        max_blocks=int(args.threshold_max_blocks),
        max_block_cycles=float("inf"),
        min_block_cycles=float(args.threshold_min_block_cycles),
        target_state_fraction=float(args.threshold_target_state_fraction),
        saturation_tol_fraction=float(args.threshold_saturation_tol_fraction),
        target_dN_store_unbounded=float(args.threshold_target_dN_store_unbounded),
        n_advances=int(args.threshold_n_advances),
        da_m=float(args.threshold_da_m),
        n_phase=int(args.threshold_n_phase),
        target_dB=float(args.threshold_target_dB),
        target_dN_store=0.01,
    )


def run_existing_v1_growth_point(driver, surface: pd.Series, context: dict,
                                 T_K: float, DeltaK_MPa_sqrtm: float,
                                 args) -> dict:
    """Evaluate one da/dN point through the existing V1 fatigue-growth driver."""
    case = _context_case(context)
    native_args = _threshold_driver_args(args)
    Kmax = float(DeltaK_MPa_sqrtm) / max(1.0 - float(args.R), 1.0e-12)

    barriers = _surface_plastic_barriers(driver, surface)
    original_builder = driver.make_plastic_barriers

    def _surface_builder(_case, _S_e_kB, _T_anchor):
        return barriers

    # The existing run_fatigue point evaluator is sequential. Temporarily bind
    # the selected surface family to its native barrier-builder hook.
    driver.make_plastic_barriers = _surface_builder
    try:
        rec = driver.run_fatigue(
            case,
            S_e_kB=0.0,
            S_c_kB=float(context.get("cleave_S_kB", 0.0)),
            T=float(T_K),
            Kmax_MPa=float(Kmax),
            args=native_args,
        )
    finally:
        driver.make_plastic_barriers = original_builder

    rec.update({
        "surface_id": str(surface["surface_id"]),
        "surface_index": int(surface["surface_index"]),
        "context_id": str(context["context_id"]),
        "source_case_label": str(context["source_case_label"]),
        "source_response_regime": str(context["source_response_regime"]),
        "T_K": float(T_K),
        "DeltaK_MPa_sqrtm": float(DeltaK_MPa_sqrtm),
        "Kmax_MPa_sqrtm": float(Kmax),
        "threshold_point_engine": "existing_V1_run_fatigue",
    })
    return rec


# -----------------------------------------------------------------------------
# Rate-defined threshold bracketing and restartable refinement
# -----------------------------------------------------------------------------

def _point_rate_state(row: pd.Series, criterion: float) -> Tuple[str, float]:
    """Classify a raw da/dN point relative to a rate criterion.

    Censored no-growth points are usable only when their upper bound is below
    the criterion. They are never treated as zero growth.
    """
    status = str(row.get("status", ""))
    rate = row.get("da_dN_m_per_cycle", np.nan)
    ub = row.get("da_dN_upper_bound_m_per_cycle", np.nan)
    if pd.notna(rate) and np.isfinite(float(rate)) and float(rate) > 0.0:
        val = float(rate)
        return ("above" if val >= criterion else "below", val)
    if status == "censored_cycle_horizon" and pd.notna(ub) and np.isfinite(float(ub)):
        val = float(ub)
        if val < criterion:
            return "below", val
    return "unknown", np.nan


def _threshold_bracket(points: pd.DataFrame, criterion: float,
                       grid_min: float, grid_max: float) -> dict:
    if points.empty or "DeltaK_MPa_sqrtm" not in points.columns:
        return {
            "status": "unresolved",
            "lo_K": np.nan,
            "lo_rate": np.nan,
            "hi_K": np.nan,
            "hi_rate": np.nan,
        }
    g = points.sort_values("DeltaK_MPa_sqrtm").drop_duplicates("DeltaK_MPa_sqrtm", keep="last")
    below, above = [], []
    for _, r in g.iterrows():
        state, rate = _point_rate_state(r, criterion)
        k = float(r["DeltaK_MPa_sqrtm"])
        if state == "below":
            below.append((k, rate))
        elif state == "above":
            above.append((k, rate))

    lo = max(below, key=lambda x: x[0]) if below else None
    hi_candidates = [x for x in above if lo is None or x[0] > lo[0]]
    hi = min(hi_candidates, key=lambda x: x[0]) if hi_candidates else None

    min_row = g.iloc[0] if len(g) else None
    max_row = g.iloc[-1] if len(g) else None
    min_state = _point_rate_state(min_row, criterion)[0] if min_row is not None else "unknown"
    max_state = _point_rate_state(max_row, criterion)[0] if max_row is not None else "unknown"

    if lo and hi and lo[0] < hi[0]:
        status = "bracketed"
    elif min_state == "above" and float(min_row["DeltaK_MPa_sqrtm"]) <= grid_min * (1 + 1e-9):
        status = "below_grid"
    elif max_state == "below" and float(max_row["DeltaK_MPa_sqrtm"]) >= grid_max * (1 - 1e-9):
        status = "above_grid_or_no_growth"
    else:
        status = "unresolved"
    return {
        "status": status,
        "lo_K": lo[0] if lo else np.nan,
        "lo_rate": lo[1] if lo else np.nan,
        "hi_K": hi[0] if hi else np.nan,
        "hi_rate": hi[1] if hi else np.nan,
    }


def _next_refinement_K(bracket: dict, criterion: float,
                       abs_tol: float, rel_tol: float) -> Tuple[float | None, bool]:
    if bracket["status"] != "bracketed":
        return None, False
    lo = float(bracket["lo_K"]); hi = float(bracket["hi_K"])
    width = hi - lo
    mid = 0.5 * (lo + hi)
    converged = width <= max(float(abs_tol), float(rel_tol) * max(mid, 1e-12))
    if converged:
        return None, True
    rlo = max(float(bracket["lo_rate"]), 1e-300)
    rhi = max(float(bracket["hi_rate"]), 1e-300)
    y0 = math.log10(rlo); y1 = math.log10(rhi); yt = math.log10(float(criterion))
    if np.isfinite(y0) and np.isfinite(y1) and abs(y1 - y0) > 1e-12:
        k = lo + (yt - y0) * (hi - lo) / (y1 - y0)
    else:
        k = mid
    # Keep the next point away from a numerically duplicated endpoint.
    margin = 0.10 * width
    k = min(max(k, lo + margin), hi - margin)
    return float(k), False


def _has_K(points: pd.DataFrame, K: float, tol: float = 1e-9) -> bool:
    if points.empty:
        return False
    vals = points["DeltaK_MPa_sqrtm"].to_numpy(float)
    return bool(np.any(np.isclose(vals, float(K), rtol=0.0, atol=tol)))


def _choose_seed_K(points: pd.DataFrame, seeds: Sequence[float],
                   criteria: Sequence[float], grid_min: float, grid_max: float) -> float | None:
    """Choose the next missing coarse seed needed to establish threshold brackets."""
    ordered = sorted(set(float(x) for x in seeds if grid_min <= float(x) <= grid_max))
    if grid_min not in ordered:
        ordered.insert(0, grid_min)
    if grid_max not in ordered:
        ordered.append(grid_max)
    # Endpoints first, then interior points from coarse to fine.
    priority = [ordered[0], ordered[-1]]
    interior = ordered[1:-1]
    if interior:
        center = 0.5 * (grid_min + grid_max)
        interior = sorted(interior, key=lambda x: abs(math.log(max(x,1e-12)) - math.log(max(center,1e-12))))
        priority.extend(interior)
    for K in priority:
        if not _has_K(points, K):
            return float(K)
    # All seeds exist; no additional seed required.
    return None


def summarize_rate_thresholds(points: pd.DataFrame, criteria: Sequence[float],
                              grid_min: float, grid_max: float,
                              abs_tol: float, rel_tol: float) -> pd.DataFrame:
    rows = []
    if points.empty:
        return pd.DataFrame()
    keys = ["surface_id", "context_id", "T_K"]
    for key, g in points.groupby(keys, sort=False):
        sid, cid, T = key
        for crit in criteria:
            b = _threshold_bracket(g, float(crit), grid_min, grid_max)
            converged = False
            estimate = np.nan
            if b["status"] == "bracketed":
                _, converged = _next_refinement_K(b, crit, abs_tol, rel_tol)
                lo, hi = float(b["lo_K"]), float(b["hi_K"])
                rlo = max(float(b["lo_rate"]), 1e-300)
                rhi = max(float(b["hi_rate"]), 1e-300)
                y0, y1, yt = math.log10(rlo), math.log10(rhi), math.log10(float(crit))
                if abs(y1-y0) > 1e-12:
                    estimate = lo + (yt-y0)*(hi-lo)/(y1-y0)
                    estimate = float(np.clip(estimate, lo, hi))
                else:
                    estimate = 0.5*(lo+hi)
                status = "bracketed_converged" if converged else "bracketed_not_converged"
            else:
                status = b["status"]
            rows.append({
                "surface_id": str(sid),
                "context_id": str(cid),
                "T_K": float(T),
                "rate_criterion_m_per_cycle": float(crit),
                "DeltaK_th_MPa_sqrtm": estimate,
                "DeltaK_lower_MPa_sqrtm": b["lo_K"],
                "DeltaK_upper_MPa_sqrtm": b["hi_K"],
                "threshold_status": status,
                "threshold_converged": bool(converged),
                "n_DeltaK_evaluations": int(g["DeltaK_MPa_sqrtm"].nunique()),
                "DeltaK_grid_min_MPa_sqrtm": float(grid_min),
                "DeltaK_grid_max_MPa_sqrtm": float(grid_max),
                "threshold_definition": "existing_V1_da_dN_rate_crossing",
            })
    return pd.DataFrame(rows)


def threshold_temperature_descriptors(thresholds: pd.DataFrame) -> pd.DataFrame:
    rows = []
    if thresholds.empty:
        return pd.DataFrame()
    for (sid, cid, crit), g in thresholds.groupby(
            ["surface_id", "context_id", "rate_criterion_m_per_cycle"], sort=False):
        g = g.sort_values("T_K")
        good = g["DeltaK_th_MPa_sqrtm"].notna()
        Tg = g.loc[good, "T_K"].to_numpy(float)
        Kg = g.loc[good, "DeltaK_th_MPa_sqrtm"].to_numpy(float)
        base = {
            "surface_id": str(sid),
            "context_id": str(cid),
            "rate_criterion_m_per_cycle": float(crit),
            "n_temperature_points": int(len(g)),
            "n_resolved_thresholds": int(len(Tg)),
        }
        if len(Tg) < 3:
            rows.append({**base, "threshold_temperature_class": "insufficient_or_unresolved"})
            continue
        order = np.argsort(Tg); Tg, Kg = Tg[order], Kg[order]
        ref = interp_at(Tg, Kg, 300.0)
        if not np.isfinite(ref) or ref <= 0:
            ref = float(np.median(Kg))
        low = float(np.mean(Kg[:min(2,len(Kg))]))
        high = float(np.mean(Kg[-min(2,len(Kg)):]))
        ipeak = int(np.argmax(Kg))
        interior = 0 < ipeak < len(Kg)-1
        peak_prom = (Kg[ipeak] - max(Kg[0], Kg[-1])) / max(ref,1e-12)
        total_range = (np.max(Kg)-np.min(Kg)) / max(ref,1e-12)
        slope = np.diff(Kg)/np.diff(Tg)/max(ref,1e-12)*100.0
        pos_frac = float(np.mean(slope > 0.01)) if len(slope) else np.nan
        neg_frac = float(np.mean(slope < -0.01)) if len(slope) else np.nan
        ratio = high/max(low,1e-12)
        if interior and peak_prom >= 0.12:
            cls = "peak_shaped"
        elif ratio >= 1.35 and pos_frac >= 0.35:
            cls = "increasing_threshold"
        elif total_range <= 0.20:
            cls = "weak_temperature"
        elif ratio <= 0.80 and neg_frac >= 0.50:
            cls = "decreasing_threshold"
        else:
            cls = "mixed_temperature"
        rows.append({
            **base,
            "threshold_temperature_class": cls,
            "DeltaK_th_ref_MPa_sqrtm": float(ref),
            "DeltaK_th_lowT_MPa_sqrtm": low,
            "DeltaK_th_highT_MPa_sqrtm": high,
            "DeltaK_th_high_over_low": ratio,
            "DeltaK_th_total_range_norm": float(total_range),
            "DeltaK_th_peak_T_K": float(Tg[ipeak]),
            "DeltaK_th_peak_prominence_norm": float(peak_prom),
            "DeltaK_th_positive_slope_fraction": pos_frac,
            "DeltaK_th_negative_slope_fraction": neg_frac,
        })
    return pd.DataFrame(rows)


def run_adaptive_threshold_state(driver, surface: pd.Series, context: dict, T_K: float,
                                 raw_path: Path, existing_points: pd.DataFrame,
                                 args) -> pd.DataFrame:
    """Resume/refine one (surface, context, T) da/dN threshold state."""
    points = existing_points.copy()
    criteria = [float(x) for x in args.rate_criteria]
    grid_min = float(args.threshold_DeltaK_min)
    grid_max = float(args.threshold_DeltaK_max)

    def evaluate(K: float):
        nonlocal points
        rec = run_existing_v1_growth_point(driver, surface, context, T_K, K, args)
        append_rows(raw_path, [rec])
        points = pd.concat([points, pd.DataFrame([rec])], ignore_index=True, sort=False)
        rate = rec.get("da_dN_m_per_cycle", np.nan)
        ub = rec.get("da_dN_upper_bound_m_per_cycle", np.nan)
        print(
            f"  growth {surface['surface_id']} {context['context_id']} T={T_K:g} "
            f"DeltaK={K:.6g} status={rec.get('status')} rate={rate} ub={ub}"
        )

    # Coarse bracket establishment. Endpoints are always included.
    max_seed_evals = max(len(args.threshold_seed_DeltaK) + 2, 2)
    for _ in range(max_seed_evals):
        brackets = [_threshold_bracket(points, c, grid_min, grid_max) for c in criteria]
        if all(b["status"] in {"bracketed", "below_grid", "above_grid_or_no_growth"} for b in brackets):
            break
        K = _choose_seed_K(points, args.threshold_seed_DeltaK, criteria, grid_min, grid_max)
        if K is None:
            break
        evaluate(K)

    # Adaptive refinement. Candidate points from both criteria share the same
    # raw database, so one point can tighten both thresholds.
    for _ in range(int(args.threshold_max_refine)):
        candidates = []
        all_done = True
        for crit in criteria:
            b = _threshold_bracket(points, crit, grid_min, grid_max)
            if b["status"] == "bracketed":
                K, conv = _next_refinement_K(
                    b, crit, args.threshold_abs_tol_MPa_sqrtm,
                    args.threshold_rel_tol,
                )
                if not conv and K is not None and not _has_K(points, K, tol=1e-8):
                    candidates.append(float(K))
                    all_done = False
            elif b["status"] == "unresolved":
                all_done = False
        if all_done or not candidates:
            break
        # Evaluate the most central unique proposal first; recompute brackets on
        # the next iteration rather than over-evaluating stale proposals.
        K = sorted(set(round(x, 10) for x in candidates))[0]
        evaluate(float(K))
    return points


# -----------------------------------------------------------------------------
# Cross-phenomenon class associations
# -----------------------------------------------------------------------------

def cramers_v_table(df: pd.DataFrame, row_col: str, col_col: str,
                    context_id: str, family: str) -> Tuple[dict, List[dict]]:
    tab = pd.crosstab(df[row_col], df[col_col])
    if tab.shape[0] < 2 or tab.shape[1] < 2:
        return ({"analysis_family": family, "context_id": context_id,
                 "cramers_V": np.nan, "chi2_p": np.nan,
                 "n": int(tab.values.sum())}, [])
    chi2, p, _, exp = chi2_contingency(tab, correction=False)
    n = tab.values.sum(); k = min(tab.shape)
    V = math.sqrt(chi2 / max(n * (k - 1), 1.0))
    obs = tab.to_numpy(float)
    z = (obs - exp) / np.sqrt(np.maximum(exp, 1e-30))
    cells = []
    for i, rv in enumerate(tab.index):
        for j, cv in enumerate(tab.columns):
            cells.append({
                "analysis_family": family,
                "context_id": context_id,
                "row_class": str(rv),
                "column_class": str(cv),
                "observed": int(obs[i, j]),
                "expected": float(exp[i, j]),
                "standardized_residual": float(z[i, j]),
            })
    return ({
        "analysis_family": family,
        "context_id": context_id,
        "cramers_V": float(V),
        "chi2_p": float(p),
        "n": int(n),
        "n_row_classes": tab.shape[0],
        "n_column_classes": tab.shape[1],
    }, cells)


def build_threeway(frac_desc: pd.DataFrame, phenotype: pd.DataFrame,
                   threshold_desc: pd.DataFrame,
                   primary_rate_criterion: float) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if frac_desc.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    merged = frac_desc.merge(phenotype, on="surface_id", how="left")
    if not threshold_desc.empty:
        td = threshold_desc[np.isclose(
            threshold_desc["rate_criterion_m_per_cycle"].to_numpy(float),
            float(primary_rate_criterion), rtol=0.0, atol=max(1e-18, 1e-6*float(primary_rate_criterion))
        )].copy()
        merged = merged.merge(td, on=["surface_id", "context_id"], how="left")

    stats, cells = [], []
    pairs = [
        ("strength_response_class", "fracture_response_class", "strength_vs_fracture_class"),
        ("fracture_response_class", "fatigue_temperature_pattern", "fracture_vs_SN_temperature_phenotype"),
        ("strength_response_class", "fatigue_temperature_pattern", "strength_vs_SN_temperature_phenotype"),
        ("fracture_response_class", "threshold_temperature_class", "fracture_vs_DeltaKth_temperature_class"),
        ("strength_response_class", "threshold_temperature_class", "strength_vs_DeltaKth_temperature_class"),
        ("fatigue_temperature_pattern", "threshold_temperature_class", "SN_vs_DeltaKth_temperature_class"),
    ]
    for cid, g in merged.groupby("context_id", sort=False):
        for rcol, ccol, fam in pairs:
            if rcol not in g.columns or ccol not in g.columns:
                continue
            gg = g.dropna(subset=[rcol, ccol])
            st, ce = cramers_v_table(gg, rcol, ccol, cid, fam)
            stats.append(st); cells.extend(ce)
    return merged, pd.DataFrame(stats), pd.DataFrame(cells)


def make_overview_plot(stats: pd.DataFrame, out: Path) -> None:
    if stats.empty:
        return
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    piv = stats.pivot(index="context_id", columns="analysis_family", values="cramers_V")
    ax = piv.plot(kind="bar", figsize=(14, 6.2))
    ax.set_ylabel("Cramér's V")
    ax.set_title("Cross-phenomenon response-class association strength by fracture context")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(fontsize=7, ncol=2)
    plt.tight_layout()
    plt.savefig(out / "threeway_class_association_overview_v5_7_1.png", dpi=220)
    plt.close()


# -----------------------------------------------------------------------------
# Main incremental workflow
# -----------------------------------------------------------------------------

def _parse_filter(tokens: Sequence[str]) -> set[str]:
    out = set()
    for token in tokens:
        out.update(x.strip() for x in str(token).split(",") if x.strip())
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--atlas-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--context-table", default="fracture_contexts_v5_7.csv")
    ap.add_argument("--context-filter", nargs="*", default=[])
    ap.add_argument("--temperatures", nargs="+", type=float,
                    default=[100,200,300,400,500,600,700,800,900])
    ap.add_argument("--surface-start", type=int, default=0)
    ap.add_argument("--surface-stop", type=int, default=None)
    ap.add_argument("--monotonic-Kmax-MPa", type=float, default=40.0)
    ap.add_argument("--monotonic-dK-MPa", type=float, default=0.10)
    ap.add_argument("--Kdot-MPa-sqrtm-per-s", type=float, default=0.005)

    # Existing V1 rate-defined threshold settings.
    ap.add_argument("--v1-driver-module", default="run_v1_two_barrier_dbtt_fatigue_map_corrected")
    ap.add_argument("--threshold-temperatures", nargs="+", type=float,
                    default=[100,200,300,400,500,600,700,800,900])
    ap.add_argument("--threshold-context-filter", nargs="*", default=[])
    ap.add_argument("--threshold-surface-start", type=int, default=None)
    ap.add_argument("--threshold-surface-stop", type=int, default=None)
    ap.add_argument("--rate-criteria", nargs="+", type=float, default=[1e-10,1e-12])
    ap.add_argument("--threshold-DeltaK-min", type=float, default=0.025)
    ap.add_argument("--threshold-DeltaK-max", type=float, default=20.0)
    ap.add_argument("--threshold-seed-DeltaK", nargs="+", type=float,
                    default=[0.025,0.05,0.1,0.2,0.4,0.8,1.5,2.5,4,6,8,10,14,20])
    ap.add_argument("--threshold-abs-tol-MPa-sqrtm", type=float, default=0.05)
    ap.add_argument("--threshold-rel-tol", type=float, default=0.03)
    ap.add_argument("--threshold-max-refine", type=int, default=12)
    ap.add_argument("--threshold-cycles-max", type=float, default=2e14)
    ap.add_argument("--threshold-max-blocks", type=int, default=10000)
    ap.add_argument("--threshold-min-block-cycles", type=float, default=1e-10)
    ap.add_argument("--threshold-target-state-fraction", type=float, default=0.01)
    ap.add_argument("--threshold-saturation-tol-fraction", type=float, default=1e-4)
    ap.add_argument("--threshold-target-dN-store-unbounded", type=float, default=5.0)
    ap.add_argument("--threshold-n-advances", type=int, default=5)
    ap.add_argument("--threshold-da-m", type=float, default=20e-6)
    ap.add_argument("--threshold-n-phase", type=int, default=96)
    ap.add_argument("--threshold-target-dB", type=float, default=0.02)
    ap.add_argument("--R", type=float, default=0.1)
    ap.add_argument("--frequency-Hz", type=float, default=1000.0)

    ap.add_argument("--skip-monotonic", action="store_true")
    ap.add_argument("--skip-thresholds", action="store_true")
    ap.add_argument("--analysis-only", action="store_true")
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    atlas = Path(args.atlas_dir)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    design_path = find_one(atlas, ["independent_exp_floor_design_v5_6.csv", "independent_exp_floor_design_v5_*.csv"])
    phenotype_path = find_one(atlas, ["surface_phenotype_summary_v5_6.csv", "surface_phenotype_summary_v5_*.csv"])
    design_all = pd.read_csv(design_path).sort_values("surface_index")
    design = design_all[design_all.surface_index >= args.surface_start].copy()
    if args.surface_stop is not None:
        design = design[design.surface_index < args.surface_stop].copy()

    contexts_all = pd.read_csv(args.context_table)
    wanted = _parse_filter(args.context_filter)
    contexts = contexts_all.copy()
    if wanted:
        contexts = contexts[contexts.context_id.astype(str).isin(wanted)].copy()
        missing = wanted - set(contexts.context_id.astype(str))
        if missing:
            raise SystemExit(f"Unknown context ids: {sorted(missing)}")
    if contexts.empty:
        raise SystemExit("No fracture contexts selected")

    # ------------------------------------------------------------------
    # Monotonic Kc(T) extension
    # ------------------------------------------------------------------
    mono_path = out / "fracture_monotonic_points_v5_7_1.csv"
    legacy_mono_path = out / "fracture_monotonic_points_v5_7.csv"
    if (not mono_path.exists()) and legacy_mono_path.exists():
        shutil.copy2(legacy_mono_path, mono_path)
        print(f"reused legacy V5.7 monotonic database: {legacy_mono_path} -> {mono_path}")
    legacy_notch_equiv = out / "notch_equivalent_thresholds_v5_7.csv"
    if legacy_notch_equiv.exists():
        print(f"NOTE: ignoring obsolete life-derived notch-equivalent file: {legacy_notch_equiv}")
    if mono_path.exists() and not args.resume and not args.analysis_only:
        if args.overwrite:
            mono_path.unlink()
        else:
            raise SystemExit(f"{mono_path} exists. Use --resume or --overwrite.")
    mono_existing = pd.read_csv(mono_path) if mono_path.exists() else pd.DataFrame()
    done_m = set()
    if not mono_existing.empty:
        done_m = {(str(r.surface_id), str(r.context_id), round(float(r.T_K),9))
                  for r in mono_existing.itertuples(index=False)}
        print(f"resume monotonic database: {len(done_m)} completed keys")

    if not args.analysis_only and not args.skip_monotonic:
        total = len(design) * len(contexts) * len(set(args.temperatures))
        for _, cser in contexts.iterrows():
            context = cser.to_dict(); cid = str(context["context_id"])
            for T in sorted(set(float(x) for x in args.temperatures)):
                mask = [((str(r.surface_id), cid, round(T,9)) not in done_m)
                        for r in design.itertuples(index=False)]
                missing_design = design.loc[np.asarray(mask, bool)]
                if missing_design.empty:
                    print(f"skip complete monotonic context={cid} T={T:g} K")
                    continue
                rows = run_monotonic_batch(
                    missing_design, context, T,
                    args.monotonic_Kmax_MPa,
                    args.monotonic_dK_MPa,
                    args.Kdot_MPa_sqrtm_per_s,
                )
                append_rows(mono_path, rows)
                for r in rows:
                    done_m.add((str(r["surface_id"]), cid, round(T,9)))
                print(f"monotonic context={cid} T={T:g}: +{len(rows)} rows; keys={len(done_m)}/{total}")

    mono = pd.read_csv(mono_path) if mono_path.exists() else pd.DataFrame()
    frac_desc = fracture_curve_descriptors(mono)
    frac_desc.to_csv(out / "fracture_curve_descriptors_v5_7_1.csv", index=False)

    # ------------------------------------------------------------------
    # Existing V1 da/dN -> rate-defined DeltaK_th extension
    # ------------------------------------------------------------------
    raw_growth_path = out / "fatigue_growth_points_v5_7_1.csv"
    if raw_growth_path.exists() and not args.resume and not args.analysis_only and not args.skip_thresholds:
        if args.overwrite:
            raw_growth_path.unlink()
        else:
            raise SystemExit(f"{raw_growth_path} exists. Use --resume or --overwrite.")

    if not args.analysis_only and not args.skip_thresholds:
        driver = _load_existing_v1_driver(args.v1_driver_module)
        tstart = args.threshold_surface_start if args.threshold_surface_start is not None else args.surface_start
        tstop = args.threshold_surface_stop if args.threshold_surface_stop is not None else args.surface_stop
        threshold_design = design_all[design_all.surface_index >= tstart].copy()
        if tstop is not None:
            threshold_design = threshold_design[threshold_design.surface_index < tstop].copy()

        twanted = _parse_filter(args.threshold_context_filter)
        threshold_contexts = contexts_all.copy()
        if twanted:
            threshold_contexts = threshold_contexts[threshold_contexts.context_id.astype(str).isin(twanted)].copy()
            missing = twanted - set(threshold_contexts.context_id.astype(str))
            if missing:
                raise SystemExit(f"Unknown threshold context ids: {sorted(missing)}")
        if threshold_contexts.empty:
            raise SystemExit("No threshold contexts selected")

        raw_all = pd.read_csv(raw_growth_path) if raw_growth_path.exists() else pd.DataFrame()
        total_states = len(threshold_design) * len(threshold_contexts) * len(set(args.threshold_temperatures))
        state_i = 0
        for _, cser in threshold_contexts.iterrows():
            context = cser.to_dict(); cid = str(context["context_id"])
            for T in sorted(set(float(x) for x in args.threshold_temperatures)):
                for _, surface in threshold_design.iterrows():
                    state_i += 1
                    sid = str(surface["surface_id"])
                    if raw_all.empty:
                        state_points = pd.DataFrame()
                    else:
                        mask = (
                            raw_all.surface_id.astype(str).eq(sid)
                            & raw_all.context_id.astype(str).eq(cid)
                            & np.isclose(raw_all.T_K.to_numpy(float), T, rtol=0.0, atol=1e-9)
                        )
                        state_points = raw_all.loc[mask].copy()
                    print(f"threshold state {state_i}/{total_states}: {sid} {cid} T={T:g} K; existing points={len(state_points)}")
                    updated = run_adaptive_threshold_state(
                        driver, surface, context, T, raw_growth_path,
                        state_points, args,
                    )
                    # Keep in-memory database synchronized so later states and
                    # analysis see points appended during this process.
                    if len(updated) > len(state_points):
                        new_rows = updated.iloc[len(state_points):].copy()
                        raw_all = pd.concat([raw_all, new_rows], ignore_index=True, sort=False)

    raw_growth = pd.read_csv(raw_growth_path) if raw_growth_path.exists() else pd.DataFrame()
    thresholds = summarize_rate_thresholds(
        raw_growth,
        args.rate_criteria,
        args.threshold_DeltaK_min,
        args.threshold_DeltaK_max,
        args.threshold_abs_tol_MPa_sqrtm,
        args.threshold_rel_tol,
    )
    thresholds.to_csv(out / "fatigue_rate_thresholds_v5_7_1.csv", index=False)
    threshold_desc = threshold_temperature_descriptors(thresholds)
    threshold_desc.to_csv(out / "fatigue_threshold_temperature_descriptors_v5_7_1.csv", index=False)

    # ------------------------------------------------------------------
    # Cross-phenomenon analysis
    # ------------------------------------------------------------------
    phenotype = pd.read_csv(phenotype_path)
    primary_crit = float(args.rate_criteria[0])
    threeway, stats, cells = build_threeway(frac_desc, phenotype, threshold_desc, primary_crit)
    threeway.to_csv(out / "threeway_phenotype_summary_v5_7_1.csv", index=False)
    stats.to_csv(out / "threeway_class_association_statistics_v5_7_1.csv", index=False)
    cells.to_csv(out / "threeway_class_association_cells_v5_7_1.csv", index=False)
    make_overview_plot(stats, out)

    config = vars(args).copy()
    config["atlas_dir_resolved"] = str(atlas.resolve())
    config["design_csv"] = str(design_path)
    config["phenotype_csv"] = str(phenotype_path)
    config["n_selected_monotonic_surfaces"] = int(len(design))
    config["n_selected_monotonic_contexts"] = int(len(contexts))
    config["threshold_definition"] = "existing_V1_da_dN_rate_crossing"
    with (out / "extension_config_v5_7_1.json").open("w") as f:
        json.dump(config, f, indent=2)
    print(f"Wrote V5.7.1 incremental extension outputs to {out}")


if __name__ == "__main__":
    main()
