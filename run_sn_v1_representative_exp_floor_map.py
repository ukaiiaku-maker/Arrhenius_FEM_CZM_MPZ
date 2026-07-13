#!/usr/bin/env python3
"""V5.6 physically stratified expandable barrier-phenomena atlas for strength and S-N initiation.

Purpose
-------
Build the paper's cross-phenomenon section around *independent complete EXP-floor
free-energy surfaces*.  Each sampled surface is projected under:

1. fixed-rate plastic flow -> one strength-versus-temperature curve;
2. cyclic loading of a blunt stress concentrator -> S-N initiation response;
3. no-shield and shielded crack-opening couplings -> a paired mechanistic control.

The same emission free-energy surface is scaled to the Peierls and Taylor
barriers using the production fatigue-model ratios.  The crack-opening barrier
is held fixed in this map so that correlations isolate the plastic free-energy
landscape and its cyclic state evolution.

V5.6 preserves the V5.6 expandable architecture and adds physically constrained anomaly-amplitude stratification, denser fatigue-temperature sampling, class-association robustness versus anomaly-amplitude cap, and explicit global phenotype association outputs.

V5.6 remains intentionally faster than the earlier multirate maps:
* no matched thermal variants with duplicate 300 K fatigue curves;
* one strength rate;
* two fatigue temperatures by default;
* vectorized phase quadrature and stress-grid integration;
* adaptive stress refinement only around target-life brackets;
* checkpointing once per completed surface rather than rewriting the entire CSV
  after every stress point.
"""
from __future__ import annotations

import argparse
import csv
from collections import defaultdict
import json
import math
from pathlib import Path
from types import SimpleNamespace
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.optimize import brentq
from scipy.stats import chi2_contingency, pearsonr, qmc, spearmanr

try:
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.inspection import permutation_importance
    from sklearn.model_selection import KFold, cross_val_score
    HAVE_SKLEARN = True
    SKLEARN_IMPORT_ERROR = ""
except Exception as exc:
    HAVE_SKLEARN = False
    SKLEARN_IMPORT_ERROR = f"{type(exc).__name__}: {exc}"

from arrhenius_fracture.sn_arrhenius_chain import build_chain_from_namespace
from arrhenius_fracture.sn_v1_arrhenius import SNCase, build_parser as sn_parser
from arrhenius_fracture.sn_v1_arrhenius_batch import run_stress_grid

KBEV = 8.617333262145e-5
FIXED_LIFE_EXPONENTS = list(range(4, 13))
DESIGN_VERSION = "v5_6"

# Complete independent surface design.  Thermal coordinates are sampled in a
# dimensionless form that naturally spans softening, compensation, and anomaly
# regimes, but final selection is made in emergent strength-curve space.
SURFACE_RANGES = {
    "exp_G00_eV": (0.8, 2.5),
    "exp_sigc0_GPa": (0.5, 4.0),
    "exp_a": (0.05, 0.80),
    "exp_n": (0.50, 1.50),
    "exp_floor_frac": (0.005, 0.050),
    "eta_G_Tref_over_G00": (-0.50, 3.00),
    "eta_sigc_Tref_over_sigc0": (-0.80, 0.50),
}

# Per-batch target composition for the paper-facing atlas.  The primary anomaly
# strata are capped at realistic strength gains; a small sensitivity stratum
# retains moderately stronger responses without allowing the extreme >2x tails
# present in the broad V5.6 mathematical atlas.
PHYSICAL_STRATUM_WEIGHTS = [
    ("strong_softening", 0.125),
    ("ordinary_softening", 0.125),
    ("plateau", 0.1875),
    ("lowT_anomaly_realistic", 0.15625),
    ("midT_anomaly_realistic", 0.15625),
    ("highT_anomaly_realistic", 0.15625),
    ("anomaly_sensitivity", 0.09375),
]


def write_csv(path: Path | str, rows: Sequence[dict]) -> None:
    rows = list(rows)
    if not rows:
        return
    keys: List[str] = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def apply_barrier(ns, r: dict):
    ns.exp_system = "W[100]"  # fallback container; all surface fields overridden
    for k in [
        "exp_G00_eV", "exp_gT_eV_per_K", "exp_sigc0_GPa",
        "exp_sT_MPa_per_K", "exp_a", "exp_n", "exp_floor_frac",
        "exp_Tref_K",
    ]:
        setattr(ns, k, float(r[k]))
    # Established scaling to the surface-dislocation-nucleation barrier family.
    ns.emit_energy_scale = 0.75
    ns.emit_entropy_scale = 0.75
    ns.peierls_energy_scale = 0.00375
    ns.peierls_entropy_scale = 0.00375
    ns.taylor_energy_scale = 0.015
    ns.taylor_entropy_scale = 0.015
    return ns


def barrier_override_audit(ns, r: dict) -> dict:
    chain = build_chain_from_namespace(ns, ns.b_m)
    b = chain.emit.base
    expected = {
        "G00_eV": float(r["exp_G00_eV"]),
        "gT_eV_per_K": float(r["exp_gT_eV_per_K"]),
        "sigc0_GPa": float(r["exp_sigc0_GPa"]),
        "sT_MPa_per_K": float(r["exp_sT_MPa_per_K"]),
        "Tref_K": float(r["exp_Tref_K"]),
        "a": float(r["exp_a"]),
        "n": float(r["exp_n"]),
        "floor_frac": float(r["exp_floor_frac"]),
    }
    actual = {
        "G00_eV": float(b.G00_eV),
        "gT_eV_per_K": float(b.gT_eV_per_K),
        "sigc0_GPa": float(b.sigc0_Pa / 1e9),
        "sT_MPa_per_K": float(b.sT_Pa_per_K / 1e6),
        "Tref_K": float(b.Tref_K),
        "a": float(b.a),
        "n": float(b.n),
        "floor_frac": float(b.Gfloor_fraction),
    }
    errs = {f"abs_error_{k}": abs(actual[k] - expected[k]) for k in actual}
    ok = max(errs.values()) < 1e-10
    return {**{f"actual_{k}": v for k, v in actual.items()}, **errs, "override_audit_pass": bool(ok)}


def strength_at_T(ns, T: float, target_rate: float, max_sigma_GPa: float = 100.0):
    chain = build_chain_from_namespace(ns, ns.b_m)
    rho = np.array([ns.rho0], float)

    def rate(sig):
        return float(chain.rates(np.array([sig]), rho, T)["dot_ep"][0])

    if rate(0.0) >= target_rate:
        return 0.0, "zero_stress_flow"
    hi = 0.25e9
    while hi < max_sigma_GPa * 1e9 and rate(hi) < target_rate:
        hi *= 2.0
    if rate(hi) < target_rate:
        return float("nan"), "above_search_limit"
    root = brentq(
        lambda s: math.log(max(rate(s), 1e-300)) - math.log(target_rate),
        0.0, hi, maxiter=160,
    )
    return root / 1e6, "resolved"


def interp_at(x, y, xq):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    good = np.isfinite(x) & np.isfinite(y)
    if good.sum() < 2 or xq < x[good].min() or xq > x[good].max():
        return float("nan")
    return float(np.interp(xq, x[good], y[good]))


def contiguous_segments(mask: np.ndarray) -> List[Tuple[int, int]]:
    mask = np.asarray(mask, bool)
    out = []
    start = None
    for i, ok in enumerate(mask):
        if ok and start is None:
            start = i
        if start is not None and ((not ok) or i == len(mask) - 1):
            end = i if ok and i == len(mask) - 1 else i - 1
            out.append((start, end))
            start = None
    return out


def longest_width(T: np.ndarray, segmask: np.ndarray) -> float:
    best = 0.0
    for i, j in contiguous_segments(segmask):
        best = max(best, float(T[j + 1] - T[i]))
    return best


def curve_shape_descriptors(T, s, reference_T: float = 300.0) -> dict:
    T = np.asarray(T, float)
    s = np.asarray(s, float)
    good = np.isfinite(T) & np.isfinite(s)
    T, s = T[good], s[good]
    if len(T) < 5:
        return {"valid": False}
    order = np.argsort(T)
    T, s = T[order], s[order]
    sref = interp_at(T, s, reference_T)
    if not np.isfinite(sref) or sref <= 1e-9:
        return {"valid": False}

    ds = np.diff(s)
    dT = np.diff(T)
    midT = 0.5 * (T[:-1] + T[1:])
    slope100 = ds / np.maximum(dT, 1e-30) * 100.0 / sref
    active = (s[:-1] >= 0.05 * sref) & (s[1:] >= 0.05 * sref)
    pos = (slope100 > 0.0) & active
    plateau = (np.abs(slope100) <= 0.02) & active
    nonsoft = (slope100 >= -0.02) & active

    pos_segments = contiguous_segments(pos)
    gains = []
    for i, j in pos_segments:
        gains.append(((s[j + 1] - s[i]) / sref, T[j + 1] - T[i], i, j))
    if gains:
        gain, width, i, j = max(gains, key=lambda z: (z[0], z[1]))
        onset = float(T[i])
        peak = float(T[j + 1])
        end = float(T[j + 1])
        gain = float(max(gain, 0.0))
        width = float(width)
    else:
        gain = 0.0
        onset = peak = end = float("nan")
        width = 0.0

    pos_weights = np.maximum(ds, 0.0)
    pos_centroid = float(np.sum(midT * pos_weights) / np.sum(pos_weights)) if np.sum(pos_weights) > 0 else float("nan")
    global_peak_i = int(np.argmax(s))
    out = {
        "valid": True,
        "sigma_ref_MPa": float(sref),
        "anomaly_gain_frac": gain,
        "anomaly_onset_K": onset,
        "anomaly_peak_T_K": peak,
        "anomaly_end_K": end,
        "anomaly_width_K": width,
        "positive_slope_area_frac": float(np.sum(np.maximum(ds, 0.0)) / sref),
        "positive_slope_width_total_K": float(np.sum(dT[pos])) if len(dT) else 0.0,
        "positive_slope_width_longest_K": longest_width(T, pos),
        "positive_slope_T_centroid_K": pos_centroid,
        "nonsoftening_width_K": longest_width(T, nonsoft),
        "plateau_width_K": longest_width(T, plateau),
        "max_norm_slope_per_100K": float(np.max(slope100)),
        "min_norm_slope_per_100K": float(np.min(slope100)),
        "global_peak_T_K": float(T[global_peak_i]),
        "global_peak_amp_vs_ref_frac": float((s[global_peak_i] - sref) / sref),
    }
    for tq in [50, 100, 150, 200, 300, 500, 700, 900, 1100]:
        sq = interp_at(T, s, tq)
        out[f"retention_{tq}"] = sq / sref if np.isfinite(sq) else float("nan")
    return out


def relative_temperature_descriptors(T, s, Tf: float, global_desc: dict) -> dict:
    T = np.asarray(T, float)
    s = np.asarray(s, float)
    order = np.argsort(T)
    T, s = T[order], s[order]
    sTf = interp_at(T, s, Tf)
    if not np.isfinite(sTf) or sTf <= 0:
        return {"valid": False}
    ds = np.diff(s)
    dT = np.diff(T)
    slope100 = ds / np.maximum(dT, 1e-30) * 100.0 / sTf
    midT = 0.5 * (T[:-1] + T[1:])
    above = midT >= Tf
    active = (s[:-1] >= 0.05 * sTf) & (s[1:] >= 0.05 * sTf)
    pos = (slope100 > 0) & active & above
    nonsoft = (slope100 >= -0.02) & active & above
    plateau = (np.abs(slope100) <= 0.02) & active & above
    peakT = global_desc.get("anomaly_peak_T_K", np.nan)
    onset = global_desc.get("anomaly_onset_K", np.nan)
    return {
        "valid": True,
        "sigma_Tf_MPa": float(sTf),
        "peak_minus_Tf_K": float(peakT - Tf) if np.isfinite(peakT) else np.nan,
        "onset_minus_Tf_K": float(onset - Tf) if np.isfinite(onset) else np.nan,
        "positive_slope_area_above_Tf_frac": float(np.sum(np.maximum(ds[above], 0.0)) / sTf),
        "positive_slope_width_above_Tf_K": float(np.sum(dT[pos])) if len(dT) else 0.0,
        "nonsoftening_width_above_Tf_K": longest_width(T, nonsoft),
        "plateau_width_above_Tf_K": longest_width(T, plateau),
        "Tf_below_anomaly_peak": bool(np.isfinite(peakT) and Tf < peakT),
        "Tf_below_anomaly_onset": bool(np.isfinite(onset) and Tf < onset),
    }


def response_class(desc: dict) -> str:
    gain = float(desc.get("anomaly_gain_frac", 0.0))
    peak = desc.get("anomaly_peak_T_K", np.nan)
    plat = float(desc.get("plateau_width_K", 0.0))
    nons = float(desc.get("nonsoftening_width_K", 0.0))
    ret900 = float(desc.get("retention_900", np.nan))
    if gain >= 0.02 and np.isfinite(peak):
        if peak < 300:
            return "lowT_anomaly"
        if peak <= 700:
            return "midT_anomaly"
        return "highT_anomaly"
    if plat >= 150:
        return "plateau"
    if nons >= 200:
        return "compensated"
    if np.isfinite(ret900) and ret900 < 0.4:
        return "strong_softening"
    return "ordinary_softening"


def sample_candidate_surfaces(n: int, seed: int, Tref: float = 300.0) -> List[dict]:
    """Generate a mixed candidate bank enriched around useful response regions.

    The final design is still selected from emergent strength-curve topology.
    These banks only improve candidate efficiency so realistic low-, mid-, and
    high-temperature anomalies are adequately represented without requiring an
    enormous rejection pool.
    """
    keys = list(SURFACE_RANGES)
    sampler = qmc.LatinHypercube(d=len(keys), seed=seed)
    u = sampler.random(n)
    lo = np.array([SURFACE_RANGES[k][0] for k in keys])
    hi = np.array([SURFACE_RANGES[k][1] for k in keys])
    vals = qmc.scale(u, lo, hi)
    rows = []
    rng = np.random.default_rng(seed + 7919)
    for i, row in enumerate(vals):
        r = {k: float(v) for k, v in zip(keys, row)}
        q = i / max(n, 1)
        if q < 0.25:
            r["candidate_bank"] = "broad"
        elif q < 0.45:
            r["eta_G_Tref_over_G00"] = float(rng.uniform(0.70, 1.25))
            r["eta_sigc_Tref_over_sigc0"] = float(rng.uniform(-0.20, 0.20))
            r["candidate_bank"] = "compensation_enriched"
        elif q < 0.60:
            r["eta_G_Tref_over_G00"] = float(rng.uniform(0.95, 1.40))
            r["eta_sigc_Tref_over_sigc0"] = float(rng.uniform(-0.70, -0.05))
            r["candidate_bank"] = "lowT_anomaly_enriched"
        elif q < 0.80:
            r["eta_G_Tref_over_G00"] = float(rng.uniform(1.05, 2.30))
            r["eta_sigc_Tref_over_sigc0"] = float(rng.uniform(-0.65, -0.02))
            r["candidate_bank"] = "midT_anomaly_enriched"
        else:
            r["eta_G_Tref_over_G00"] = float(rng.uniform(0.72, 1.12))
            r["eta_sigc_Tref_over_sigc0"] = float(rng.uniform(0.00, 0.45))
            r["candidate_bank"] = "highT_anomaly_enriched"
        r["exp_Tref_K"] = float(Tref)
        r["exp_gT_eV_per_K"] = r["eta_G_Tref_over_G00"] * r["exp_G00_eV"] / Tref
        r["exp_sT_MPa_per_K"] = r["eta_sigc_Tref_over_sigc0"] * (1000.0 * r["exp_sigc0_GPa"]) / Tref
        r["implied_S0_kB"] = -r["exp_gT_eV_per_K"] / KBEV
        r["candidate_index"] = i
        rows.append(r)
    return rows


def physical_design_stratum(cls: str, desc: dict, gain_min: float, main_gain_max: float, sensitivity_gain_max: float) -> str | None:
    gain = float(desc.get("anomaly_gain_frac", 0.0) or 0.0)
    if cls in {"strong_softening", "ordinary_softening", "plateau"}:
        return cls
    if cls in {"lowT_anomaly", "midT_anomaly", "highT_anomaly"}:
        if gain_min <= gain <= main_gain_max:
            return f"{cls}_realistic"
        if main_gain_max < gain <= sensitivity_gain_max:
            return "anomaly_sensitivity"
    return None

def shape_vector(T: Sequence[float], s: Sequence[float], d: dict) -> np.ndarray:
    T = np.asarray(T, float)
    s = np.asarray(s, float)
    sref = d.get("sigma_ref_MPa", np.nan)
    if not np.isfinite(sref) or sref <= 0:
        return np.full(len(T) + 7, np.nan)
    y = np.clip(np.nan_to_num(s / sref, nan=5.0, posinf=5.0, neginf=0.0), 0.0, 5.0)
    extra = np.array([
        d.get("anomaly_gain_frac", 0.0),
        d.get("anomaly_peak_T_K", 0.0) / 1200.0 if np.isfinite(d.get("anomaly_peak_T_K", np.nan)) else 0.0,
        d.get("positive_slope_area_frac", 0.0),
        d.get("positive_slope_width_longest_K", 0.0) / 1200.0,
        d.get("nonsoftening_width_K", 0.0) / 1200.0,
        d.get("plateau_width_K", 0.0) / 1200.0,
        d.get("retention_900", 0.0),
    ], float)
    return np.concatenate([y, np.nan_to_num(extra)])


def farthest_select(Z: np.ndarray, candidates: List[int], k: int, already: List[int] | None = None) -> List[int]:
    candidates = list(candidates)
    already = list(already or [])
    if len(candidates) <= k:
        return candidates
    selected = []
    if not already:
        selected.append(candidates[0])
    while len(selected) < k:
        refs = already + selected
        if not refs:
            nxt = candidates[len(selected)]
        else:
            D = np.full(len(candidates), np.inf)
            for j in refs:
                D = np.minimum(D, np.linalg.norm(Z[candidates] - Z[j], axis=1))
            for s in selected:
                D[candidates.index(s)] = -np.inf
            nxt = candidates[int(np.argmax(D))]
        if nxt in selected:
            break
        selected.append(nxt)
    return selected[:k]


def build_independent_design(base_args, n_surfaces: int, candidate_pool: int, seed: int,
                             prescreen_T: Sequence[float], strength_rate: float,
                             surface_offset: int = 0, design_batch: int = 0,
                             anomaly_gain_min: float = 0.05,
                             anomaly_gain_main_max: float = 0.50,
                             anomaly_gain_sensitivity_max: float = 1.00) -> List[dict]:
    candidates = sample_candidate_surfaces(candidate_pool, seed)
    valid_rows, vectors, descs, classes, strata = [], [], [], [], []
    for i, r in enumerate(candidates, 1):
        if i % 100 == 0 or i == 1:
            print(f"DESIGN prescreen {i}/{len(candidates)}")
        ns = apply_barrier(SimpleNamespace(**vars(base_args)), r)
        svals = [strength_at_T(ns, float(T), strength_rate)[0] for T in prescreen_T]
        if np.mean(np.isfinite(np.asarray(svals, float))) < 0.85:
            continue
        d = curve_shape_descriptors(prescreen_T, svals)
        if not d.get("valid", False):
            continue
        v = shape_vector(prescreen_T, svals, d)
        if not np.all(np.isfinite(v)):
            continue
        cls = response_class(d)
        stratum = physical_design_stratum(cls, d, anomaly_gain_min, anomaly_gain_main_max, anomaly_gain_sensitivity_max)
        if stratum is None:
            continue
        valid_rows.append(r)
        vectors.append(v)
        descs.append(d)
        classes.append(cls)
        strata.append(stratum)

    if len(valid_rows) < n_surfaces:
        raise RuntimeError(f"Only {len(valid_rows)} physically admissible candidate surfaces; need {n_surfaces}. Increase CANDIDATE_POOL.")
    X = np.vstack(vectors)
    mu, sd = X.mean(axis=0), X.std(axis=0)
    sd[sd < 1e-12] = 1.0
    Z = (X - mu) / sd

    # Convert target weights to integer quotas while preserving the requested total.
    raw = np.array([w * n_surfaces for _, w in PHYSICAL_STRATUM_WEIGHTS], float)
    quotas = np.floor(raw).astype(int)
    remainder = int(n_surfaces - quotas.sum())
    if remainder > 0:
        order = np.argsort(-(raw - quotas))
        for k in order[:remainder]:
            quotas[k] += 1

    selected: List[int] = []
    shortfall = 0
    for (stratum_name, _), quota in zip(PHYSICAL_STRATUM_WEIGHTS, quotas):
        inds = [i for i, s in enumerate(strata) if s == stratum_name and i not in selected]
        take = min(int(quota), len(inds))
        if take:
            selected.extend(farthest_select(Z, inds, take, selected))
        shortfall += int(quota) - take
        print(f"DESIGN stratum {stratum_name}: available={len(inds)} selected={take} target={int(quota)}")

    if shortfall > 0:
        pool = [i for i in range(len(valid_rows)) if i not in selected]
        if len(pool) < shortfall:
            raise RuntimeError(f"Physical design shortfall {shortfall} cannot be filled from {len(pool)} remaining admissible candidates. Increase CANDIDATE_POOL.")
        selected.extend(farthest_select(Z, pool, shortfall, selected))
    selected = selected[:n_surfaces]

    out = []
    for j, idx in enumerate(selected):
        r = dict(valid_rows[idx])
        global_index = int(surface_offset + j)
        r["surface_id"] = f"EXP_{global_index:05d}"
        r["surface_index"] = global_index
        r["design_batch"] = int(design_batch)
        r["design_batch_seed"] = int(seed)
        r["design_version"] = DESIGN_VERSION
        r["design_profile"] = "physical_anomaly_capped"
        r["design_stratum"] = strata[idx]
        r["anomaly_gain_min_design"] = float(anomaly_gain_min)
        r["anomaly_gain_main_max_design"] = float(anomaly_gain_main_max)
        r["anomaly_gain_sensitivity_max_design"] = float(anomaly_gain_sensitivity_max)
        r["prescreen_response_class"] = classes[idx]
        d = descs[idx]
        for k in ["anomaly_gain_frac", "anomaly_peak_T_K", "positive_slope_area_frac",
                  "positive_slope_width_longest_K", "nonsoftening_width_K", "plateau_width_K",
                  "retention_900"]:
            r[f"prescreen_{k}"] = d.get(k, np.nan)
        out.append(r)
    return out



def extend_independent_design(base_args, existing: Sequence[dict], target_surfaces: int,
                              batch_size: int, candidate_pool_per_batch: int, seed: int,
                              prescreen_T: Sequence[float], strength_rate: float,
                              anomaly_gain_min: float, anomaly_gain_main_max: float,
                              anomaly_gain_sensitivity_max: float) -> List[dict]:
    """Extend a persistent design in deterministic independent batches.

    Existing rows are never regenerated or reselected.  Increasing the requested
    target therefore appends new barrier surfaces without invalidating previously
    completed strength or fatigue calculations.
    """
    design = [dict(r) for r in existing]
    if target_surfaces <= len(design):
        return design
    batch_size = max(int(batch_size), 1)
    candidate_pool_per_batch = max(int(candidate_pool_per_batch), batch_size)
    if design:
        old_batches = [int(float(r.get("design_batch", -1))) for r in design if pd.notna(r.get("design_batch", np.nan))]
        batch_index = (max(old_batches) + 1) if old_batches else int(math.ceil(len(design) / batch_size))
    else:
        batch_index = 0
    while len(design) < target_surfaces:
        n_add = min(batch_size, target_surfaces - len(design))
        batch_seed = int(seed + 100003 * batch_index)
        print(f"DESIGN batch {batch_index}: selecting {n_add} surfaces from {candidate_pool_per_batch} candidates (seed={batch_seed})")
        new_rows = build_independent_design(
            base_args, n_add, candidate_pool_per_batch, batch_seed,
            prescreen_T, strength_rate, surface_offset=len(design), design_batch=batch_index,
            anomaly_gain_min=anomaly_gain_min,
            anomaly_gain_main_max=anomaly_gain_main_max,
            anomaly_gain_sensitivity_max=anomaly_gain_sensitivity_max
        )
        design.extend(new_rows)
        batch_index += 1
    return design

def numerical_barrier_shape_descriptors(ns, Tf: float) -> dict:
    chain = build_chain_from_namespace(ns, ns.b_m)
    b = chain.emit
    sigc = float(b.sigc_Pa(Tf))
    G0 = float(b.G0_eV(Tf))
    Gf = float(b.Gfloor_eV(Tf))
    x = np.linspace(0.0, 3.0, 301)
    sig = x * sigc
    G = np.asarray(b.deltaG_eV(sig, Tf), float)
    gnorm = G / max(G0, 1e-30)
    dGdx = np.gradient(gnorm, x)
    d2Gdx2 = np.gradient(dGdx, x)
    dGdsig = np.gradient(G, sig)  # eV / Pa
    hz_sens_GPa = -dGdsig / max(KBEV * Tf, 1e-30) * 1e9

    def at(arr, xq):
        return float(np.interp(xq, x, arr))

    abs_slope = np.abs(dGdx)
    maxs = max(float(np.max(abs_slope)), 1e-30)
    flat05 = abs_slope <= 0.05 * maxs
    flat10 = abs_slope <= 0.10 * maxs
    # longest width on node mask
    def width_mask(mask):
        best = 0.0
        start = None
        for i, ok in enumerate(mask):
            if ok and start is None:
                start = i
            if start is not None and ((not ok) or i == len(mask)-1):
                end = i if ok and i == len(mask)-1 else i-1
                best = max(best, x[end] - x[start])
                start = None
        return float(best)

    # Thermal derivative of reduced barrier G/(kBT) at fixed physical sigma=sigc(Tf).
    dtemp = max(1.0, 0.005 * Tf)
    T1, T2 = max(5.0, Tf - dtemp), Tf + dtemp
    sig_fixed = sigc
    red1 = float(b.deltaG_eV(sig_fixed, T1)) / max(KBEV * T1, 1e-30)
    red2 = float(b.deltaG_eV(sig_fixed, T2)) / max(KBEV * T2, 1e-30)
    dred_dT = (red2 - red1) / (T2 - T1)

    return {
        "G0_eV": G0,
        "sigc_GPa": sigc / 1e9,
        "floor_ratio": Gf / max(G0, 1e-30),
        "drop_fraction_x0_to3": float((G[0] - G[-1]) / max(G0, 1e-30)),
        "stress_slope_abs_min_norm_per_x": float(np.min(abs_slope)),
        "stress_slope_abs_median_norm_per_x": float(np.median(abs_slope)),
        "stress_slope_abs_max_norm_per_x": float(np.max(abs_slope)),
        "stress_flat_width_x_5pct_of_maxslope": width_mask(flat05),
        "stress_flat_width_x_10pct_of_maxslope": width_mask(flat10),
        "curvature_abs_integral_norm": float(np.trapezoid(np.abs(d2Gdx2), x)),
        "curvature_abs_max_norm": float(np.max(np.abs(d2Gdx2))),
        "hazard_sensitivity_median_per_GPa": float(np.median(hz_sens_GPa)),
        "hazard_sensitivity_min_per_GPa": float(np.min(hz_sens_GPa)),
        "hazard_sensitivity_max_per_GPa": float(np.max(hz_sens_GPa)),
        "hazard_sensitivity_x0p5_per_GPa": at(hz_sens_GPa, 0.5),
        "hazard_sensitivity_x1_per_GPa": at(hz_sens_GPa, 1.0),
        "hazard_sensitivity_x1p5_per_GPa": at(hz_sens_GPa, 1.5),
        "hazard_sensitivity_x2_per_GPa": at(hz_sens_GPa, 2.0),
        "reduced_barrier_thermal_derivative_at_fixed_sigc_per_K": float(dred_dT),
    }


def barrier_descriptors(r: dict, ns, fatigue_temperatures: Sequence[float]) -> dict:
    out = {
        "barrier_G00_eV": float(r["exp_G00_eV"]),
        "barrier_gT_eV_per_K": float(r["exp_gT_eV_per_K"]),
        "barrier_sigc0_GPa": float(r["exp_sigc0_GPa"]),
        "barrier_sT_MPa_per_K": float(r["exp_sT_MPa_per_K"]),
        "barrier_a": float(r["exp_a"]),
        "barrier_n": float(r["exp_n"]),
        "barrier_floor_frac": float(r["exp_floor_frac"]),
        "barrier_eta_G": float(r["eta_G_Tref_over_G00"]),
        "barrier_eta_sigc": float(r["eta_sigc_Tref_over_sigc0"]),
        "barrier_compensation_margin": float(r["eta_G_Tref_over_G00"] - 1.0),
        "barrier_implied_S0_kB": float(r["implied_S0_kB"]),
    }
    for Tf in fatigue_temperatures:
        tag = f"T{int(round(Tf))}"
        d = numerical_barrier_shape_descriptors(ns, Tf)
        for k, v in d.items():
            out[f"barrier_{tag}_{k}"] = v
    return out


def life_class_at_N(row: pd.Series, N: float) -> str:
    status = str(row.get("status", ""))
    life = row.get("cycles_to_nucleation", np.nan)
    total = float(row.get("cycles_total", 0.0))
    if status == "failed" and pd.notna(life):
        return "fail" if float(life) <= N else "survive"
    if total >= N and status == "right_censored":
        return "survive"
    return "unknown"


def stress_bracket_at_life(points: Sequence[dict], N: float) -> dict:
    df = pd.DataFrame(points).sort_values("sigma_a_MPa")
    safe, fail = [], []
    for _, row in df.iterrows():
        c = life_class_at_N(row, N)
        if c == "survive":
            safe.append(float(row.sigma_a_MPa))
        elif c == "fail":
            fail.append(float(row.sigma_a_MPa))
    low = max(safe) if safe else np.nan
    highc = [s for s in fail if not np.isfinite(low) or s > low]
    high = min(highc) if highc else np.nan
    mid = 0.5 * (low + high) if np.isfinite(low) and np.isfinite(high) else np.nan
    width = high - low if np.isfinite(low) and np.isfinite(high) else np.nan
    return {"low": low, "high": high, "mid": mid, "width": width}


def refinement_stresses(points: Sequence[dict], target_lives: Sequence[float], min_ratio: float = 1.15) -> List[float]:
    new = []
    for N in target_lives:
        b = stress_bracket_at_life(points, N)
        lo, hi = b["low"], b["high"]
        if np.isfinite(lo) and np.isfinite(hi) and lo > 0 and hi / lo >= min_ratio:
            new.append(float(math.sqrt(lo * hi)))
    existing = np.array([float(p["sigma_a_MPa"]) for p in points], float)
    out = []
    for s in sorted(set(new)):
        if len(existing) == 0 or np.min(np.abs(existing - s)) > 1e-8 * max(s, 1.0):
            out.append(s)
    return out


def sn_case_descriptors(points: Sequence[dict], sigma_ref: float, cycles_max: float) -> dict:
    df = pd.DataFrame(points)
    out: Dict[str, object] = {}
    mids = []
    for e in FIXED_LIFE_EXPONENTS:
        b = stress_bracket_at_life(points, 10.0 ** e)
        for field in ["low", "high", "mid", "width"]:
            v = b[field]
            out[f"sn_sigma_N{e}_{field}_norm"] = v / sigma_ref if np.isfinite(v) and sigma_ref > 0 else np.nan
        mids.append(out[f"sn_sigma_N{e}_mid_norm"])
    exps = np.array(FIXED_LIFE_EXPONENTS, float)
    mids = np.asarray(mids, float)
    good = np.isfinite(mids) & (mids > 0)
    high = good & (exps >= 8)
    low = good & (exps <= 7)
    out["sn_fixedlife_coverage_count"] = int(good.sum())
    out["sn_highcycle_slope_abs_norm_per_decade"] = float(abs(np.polyfit(exps[high], mids[high], 1)[0])) if high.sum() >= 2 else np.nan
    out["sn_highcycle_frac_drop_per_decade"] = float(abs(np.polyfit(exps[high], np.log(mids[high]), 1)[0])) if high.sum() >= 2 else np.nan
    out["sn_lowcycle_slope_abs_norm_per_decade"] = float(abs(np.polyfit(exps[low], mids[low], 1)[0])) if low.sum() >= 2 else np.nan
    if np.isfinite(out["sn_highcycle_slope_abs_norm_per_decade"]) and np.isfinite(out["sn_lowcycle_slope_abs_norm_per_decade"]):
        out["sn_flattening_ratio_high_over_low"] = out["sn_highcycle_slope_abs_norm_per_decade"] / max(out["sn_lowcycle_slope_abs_norm_per_decade"], 1e-12)
    else:
        out["sn_flattening_ratio_high_over_low"] = np.nan

    plateau_width = 0.0
    cur = 0.0
    for i in range(len(exps)-1):
        if exps[i] < 6 or not (good[i] and good[i+1]):
            cur = 0.0
            continue
        frac = abs(mids[i+1] - mids[i]) / max(0.5*(mids[i+1]+mids[i]), 1e-12)
        if frac <= 0.05:
            cur += 1.0
            plateau_width = max(plateau_width, cur)
        else:
            cur = 0.0
    out["sn_plateau_width_decades"] = float(plateau_width)
    out["sn_endurance_like_5pct"] = bool(plateau_width >= 2.0)

    usable = df[df.status.isin(["failed", "right_censored", "block_limited"])].copy()
    fail = usable[(usable.status == "failed") & usable.cycles_to_nucleation.notna()]
    cens = usable[usable.status == "right_censored"]
    block = usable[usable.status == "block_limited"]
    den = max(len(usable), 1)
    out["sn_failed_fraction"] = len(fail) / den
    out["sn_censored_fraction"] = len(cens) / den
    out["sn_block_limited_fraction"] = len(block) / den
    out["sn_max_censored_stress_norm"] = float(cens.sigma_a_MPa.max()/sigma_ref) if len(cens) and sigma_ref > 0 else np.nan
    out["sn_min_failed_stress_norm"] = float(fail.sigma_a_MPa.min()/sigma_ref) if len(fail) and sigma_ref > 0 else np.nan
    horizon = stress_bracket_at_life(points, cycles_max)
    out["sn_horizon_transition_mid_norm"] = horizon["mid"]/sigma_ref if np.isfinite(horizon["mid"]) else np.nan
    out["sn_horizon_transition_width_norm"] = horizon["width"]/sigma_ref if np.isfinite(horizon["width"]) else np.nan
    return out


def pair_descriptors(no_pts, sh_pts, no_desc, sh_desc) -> dict:
    out = {}
    deltas, exps = [], []
    for e in FIXED_LIFE_EXPONENTS:
        n = no_desc.get(f"sn_sigma_N{e}_mid_norm", np.nan)
        s = sh_desc.get(f"sn_sigma_N{e}_mid_norm", np.nan)
        if np.isfinite(n) and np.isfinite(s):
            out[f"pair_delta_sigma_N{e}_norm"] = float(s-n)
            out[f"pair_ratio_sigma_N{e}"] = float(s/max(n,1e-12))
            deltas.append(s-n); exps.append(e)
        else:
            out[f"pair_delta_sigma_N{e}_norm"] = np.nan
            out[f"pair_ratio_sigma_N{e}"] = np.nan
    out["pair_area_delta_sigma_fixedlife_norm_decades"] = float(np.trapezoid(deltas, exps)) if len(deltas) >= 2 else np.nan
    out["pair_mean_delta_sigma_fixedlife_norm"] = float(np.mean(deltas)) if deltas else np.nan
    out["pair_max_delta_sigma_fixedlife_norm"] = float(np.max(deltas)) if deltas else np.nan
    nsl = no_desc.get("sn_highcycle_slope_abs_norm_per_decade", np.nan)
    ssl = sh_desc.get("sn_highcycle_slope_abs_norm_per_decade", np.nan)
    out["pair_highcycle_slope_reduction"] = float(nsl-ssl) if np.isfinite(nsl) and np.isfinite(ssl) else np.nan
    out["pair_plateau_width_gain_decades"] = float(sh_desc.get("sn_plateau_width_decades",0)-no_desc.get("sn_plateau_width_decades",0))

    no = pd.DataFrame(no_pts); sh = pd.DataFrame(sh_pts)
    m = no.merge(sh, on="sigma_a_MPa", suffixes=("_no","_sh"))
    exact, lower = [], []
    for _, r in m.iterrows():
        ln, ls = r.get("cycles_to_nucleation_no", np.nan), r.get("cycles_to_nucleation_sh", np.nan)
        if pd.notna(ln) and float(ln)>0:
            if pd.notna(ls) and float(ls)>0:
                v=math.log10(float(ls)/float(ln)); exact.append(v); lower.append(v)
            else:
                ts=float(r.get("cycles_total_sh",0.0))
                if ts>0: lower.append(math.log10(ts/float(ln)))
    out["pair_median_loglife_benefit_exact"] = float(np.median(exact)) if exact else np.nan
    out["pair_median_loglife_benefit_lower_bound"] = float(np.median(lower)) if lower else np.nan
    out["pair_max_loglife_benefit_lower_bound"] = float(np.max(lower)) if lower else np.nan
    out["pair_fraction_stresses_benefit_gt_1decade"] = float(np.mean(np.asarray(lower)>=1.0)) if lower else np.nan
    return out


def bh_qvalues(pvals):
    p=np.asarray(pvals,float); q=np.full_like(p,np.nan); good=np.isfinite(p)
    pv=p[good]
    if len(pv)==0: return q
    order=np.argsort(pv); ranked=pv[order]; m=len(ranked)
    raw=ranked*m/np.arange(1,m+1)
    adj=np.minimum.accumulate(raw[::-1])[::-1]; adj=np.minimum(adj,1.0)
    tmp=np.empty_like(adj); tmp[order]=adj; q[good]=tmp
    return q


def partial_spearman(x, y, controls):
    x=np.asarray(x,float); y=np.asarray(y,float); C=np.asarray(controls,float)
    good=np.isfinite(x)&np.isfinite(y)&np.all(np.isfinite(C),axis=1)
    x,y,C=x[good],y[good],C[good]
    if len(x)<20 or len(np.unique(x))<4 or len(np.unique(y))<4: return np.nan,np.nan,len(x)
    rx=pd.Series(x).rank().to_numpy(float); ry=pd.Series(y).rank().to_numpy(float)
    RC=np.column_stack([pd.Series(C[:,j]).rank().to_numpy(float) for j in range(C.shape[1])])
    A=np.column_stack([np.ones(len(x)),RC])
    bx=np.linalg.lstsq(A,rx,rcond=None)[0]; by=np.linalg.lstsq(A,ry,rcond=None)[0]
    ex=rx-A@bx; ey=ry-A@by
    r,p=pearsonr(ex,ey)
    return float(r),float(p),len(x)


def correlation_family(df: pd.DataFrame, predictors: Sequence[str], responses: Sequence[str], Tf: float, family: str) -> List[dict]:
    rows=[]
    base_controls=["barrier_G00_eV","barrier_sigc0_GPa"]
    for x in predictors:
        for y in responses:
            controls=[c for c in base_controls if c not in {x,y}]
            cols=list(dict.fromkeys([x,y]+controls))
            d=df[cols].replace([np.inf,-np.inf],np.nan).dropna()
            if len(d)<24 or d[x].nunique()<4 or d[y].nunique()<4: continue
            sr,sp=spearmanr(d[x],d[y]); pr,pp=pearsonr(d[x],d[y])
            if controls:
                rr,rp,rn=partial_spearman(d[x].to_numpy(),d[y].to_numpy(),d[controls].to_numpy())
            else:
                rr,rp,rn=float(sr),float(sp),len(d)
            rows.append({"analysis_family":family,"fatigue_T_K":Tf,"predictor_metric":x,"response_metric":y,
                         "n":len(d),"spearman_rho":float(sr),"spearman_p":float(sp),
                         "pearson_r":float(pr),"pearson_p":float(pp),
                         "partial_spearman_r_controlling_G00_sigc0":rr,
                         "partial_spearman_p":rp,"partial_spearman_n":rn,
                         "abs_spearman":abs(float(sr)),"abs_partial_spearman":abs(rr) if np.isfinite(rr) else np.nan})
    if rows:
        q=bh_qvalues([r["spearman_p"] for r in rows]); qp=bh_qvalues([r["partial_spearman_p"] for r in rows])
        for r,a,b in zip(rows,q,qp):
            r["spearman_q_BH"]=float(a) if np.isfinite(a) else np.nan
            r["partial_spearman_q_BH"]=float(b) if np.isfinite(b) else np.nan
    return rows


def analyze_correlations(summary_rows: Sequence[dict], fatigue_temperatures: Sequence[float]) -> List[dict]:
    """Targeted correlation families for the paper's common-landscape test.

    Scale variables that would create denominator/shared-reference tautologies
    are excluded from the strength-shape predictor set.  A separate exhaustive
    table can be built later from the raw descriptor CSV if needed.
    """
    df = pd.DataFrame(summary_rows)
    all_rows = []
    for Tf in fatigue_temperatures:
        g = df[np.isclose(df.fatigue_T_K.astype(float), float(Tf))].copy()
        tag = f"T{int(round(Tf))}"

        barrier_shape = [c for c in [
            "barrier_eta_G", "barrier_eta_sigc", "barrier_compensation_margin",
            "barrier_a", "barrier_n", "barrier_floor_frac",
            f"barrier_{tag}_floor_ratio",
            f"barrier_{tag}_drop_fraction_x0_to3",
            f"barrier_{tag}_stress_slope_abs_min_norm_per_x",
            f"barrier_{tag}_stress_slope_abs_median_norm_per_x",
            f"barrier_{tag}_stress_slope_abs_max_norm_per_x",
            f"barrier_{tag}_stress_flat_width_x_5pct_of_maxslope",
            f"barrier_{tag}_stress_flat_width_x_10pct_of_maxslope",
            f"barrier_{tag}_curvature_abs_integral_norm",
            f"barrier_{tag}_curvature_abs_max_norm",
            f"barrier_{tag}_hazard_sensitivity_median_per_GPa",
            f"barrier_{tag}_hazard_sensitivity_min_per_GPa",
            f"barrier_{tag}_hazard_sensitivity_max_per_GPa",
            f"barrier_{tag}_hazard_sensitivity_x0p5_per_GPa",
            f"barrier_{tag}_hazard_sensitivity_x1_per_GPa",
            f"barrier_{tag}_hazard_sensitivity_x1p5_per_GPa",
            f"barrier_{tag}_hazard_sensitivity_x2_per_GPa",
            f"barrier_{tag}_reduced_barrier_thermal_derivative_at_fixed_sigc_per_K",
        ] if c in g.columns]

        strength_shape = [c for c in [
            "strength_anomaly_gain_frac", "strength_anomaly_peak_T_K",
            "strength_anomaly_width_K", "strength_positive_slope_area_frac",
            "strength_positive_slope_width_total_K",
            "strength_positive_slope_width_longest_K",
            "strength_positive_slope_T_centroid_K",
            "strength_nonsoftening_width_K", "strength_plateau_width_K",
            "strength_max_norm_slope_per_100K", "strength_min_norm_slope_per_100K",
            "strength_global_peak_T_K", "strength_global_peak_amp_vs_ref_frac",
            "strength_retention_100", "strength_retention_150",
            "strength_retention_200", "strength_retention_500",
            "strength_retention_700", "strength_retention_900",
            "strength_rel_peak_minus_Tf_K", "strength_rel_onset_minus_Tf_K",
            "strength_rel_positive_slope_area_above_Tf_frac",
            "strength_rel_positive_slope_width_above_Tf_K",
            "strength_rel_nonsoftening_width_above_Tf_K",
            "strength_rel_plateau_width_above_Tf_K",
        ] if c in g.columns]

        fatigue_core = [c for c in [
            "no_shield_sn_highcycle_slope_abs_norm_per_decade",
            "no_shield_sn_highcycle_frac_drop_per_decade",
            "no_shield_sn_flattening_ratio_high_over_low",
            "no_shield_sn_plateau_width_decades",
            "no_shield_sn_horizon_transition_mid_norm",
            "shielded_sn_highcycle_slope_abs_norm_per_decade",
            "shielded_sn_highcycle_frac_drop_per_decade",
            "shielded_sn_flattening_ratio_high_over_low",
            "shielded_sn_plateau_width_decades",
            "shielded_sn_horizon_transition_mid_norm",
            "pair_highcycle_slope_reduction",
            "pair_plateau_width_gain_decades",
            "pair_area_delta_sigma_fixedlife_norm_decades",
            "pair_mean_delta_sigma_fixedlife_norm",
            "pair_max_delta_sigma_fixedlife_norm",
            "pair_median_loglife_benefit_exact",
            "pair_median_loglife_benefit_lower_bound",
            "pair_max_loglife_benefit_lower_bound",
            "pair_fraction_stresses_benefit_gt_1decade",
            "pair_delta_sigma_N8_norm", "pair_delta_sigma_N10_norm",
            "pair_delta_sigma_N12_norm",
        ] if c in g.columns]

        all_rows += correlation_family(g, barrier_shape, strength_shape, Tf, "barrier_shape_to_strength_shape")
        all_rows += correlation_family(g, barrier_shape, fatigue_core, Tf, "barrier_shape_to_fatigue")
        all_rows += correlation_family(g, strength_shape, fatigue_core, Tf, "strength_shape_to_fatigue")
    return all_rows



def fatigue_response_class(desc: dict) -> str:
    """Transparent three-class S-N topology classifier."""
    plateau = float(desc.get("sn_plateau_width_decades", 0.0) or 0.0)
    flatten = desc.get("sn_flattening_ratio_high_over_low", np.nan)
    hc = desc.get("sn_highcycle_slope_abs_norm_per_decade", np.nan)
    if plateau >= 2.0:
        return "endurance_like"
    if np.isfinite(flatten) and np.isfinite(hc) and flatten <= 0.60:
        return "knee_threshold_like"
    return "continuous_SN"


def relative_strength_regime(global_desc: dict, Tf: float) -> str:
    gain = float(global_desc.get("anomaly_gain_frac", 0.0) or 0.0)
    onset = global_desc.get("anomaly_onset_K", np.nan)
    peak = global_desc.get("anomaly_peak_T_K", np.nan)
    if gain >= 0.02 and np.isfinite(onset) and np.isfinite(peak):
        if Tf < onset:
            return "below_anomaly"
        if Tf <= peak:
            return "inside_anomaly_rise"
        return "above_anomaly_peak"
    if float(global_desc.get("plateau_width_K", 0.0) or 0.0) >= 150.0 or float(global_desc.get("nonsoftening_width_K", 0.0) or 0.0) >= 200.0:
        return "nonsoftening_family"
    return "softening_family"


def rebuild_summary_rows(design: Sequence[dict], strength_rows: Sequence[dict], sn_rows: Sequence[dict],
                         fatigue_temperatures: Sequence[float], args, base) -> List[dict]:
    """Recompute all derived descriptors from raw checkpoints.

    This makes strength-temperature and fatigue-temperature grids additive: new
    temperatures can be added later without rerunning old S-N points and without
    leaving stale descriptors in an existing summary table.
    """
    sr = pd.DataFrame(strength_rows)
    nr = pd.DataFrame(sn_rows)
    summaries: List[dict] = []
    if sr.empty or nr.empty:
        return summaries
    for ii, r in enumerate(design, 1):
        sid = str(r["surface_id"])
        g = sr[sr.surface_id.astype(str) == sid].sort_values("T_K")
        if len(g) < 5:
            continue
        T = g.T_K.to_numpy(float)
        S = g.strength_MPa.to_numpy(float)
        sd = curve_shape_descriptors(T, S)
        cls = response_class(sd) if sd.get("valid", False) else "invalid"
        ns = apply_barrier(SimpleNamespace(**vars(base)), r)
        bd = barrier_descriptors(r, ns, fatigue_temperatures)
        for Tf in fatigue_temperatures:
            h = nr[(nr.surface_id.astype(str) == sid) & np.isclose(nr.T_K.astype(float), float(Tf))]
            no_pts = h[h.case == "no_shield"].to_dict("records")
            sh_pts = h[h.case == "shielded"].to_dict("records")
            if not no_pts or not sh_pts:
                continue
            sigma_ref_vals = pd.to_numeric(h.sigma_ref_MPa, errors="coerce").dropna() if "sigma_ref_MPa" in h else pd.Series(dtype=float)
            if len(sigma_ref_vals):
                sigma_ref = float(sigma_ref_vals.iloc[0])
            else:
                sref = interp_at(T, S, args.stress_reference_temperature)
                sigma_ref = max(sref / max(base.Kt, 1e-12), args.min_stress_MPa) if np.isfinite(sref) and sref > 0 else args.min_stress_MPa
            rel = relative_temperature_descriptors(T, S, float(Tf), sd)
            nd = sn_case_descriptors(no_pts, sigma_ref, args.cycles_max)
            shd = sn_case_descriptors(sh_pts, sigma_ref, args.cycles_max)
            pair = pair_descriptors(no_pts, sh_pts, nd, shd)
            summary = {
                **r, **bd,
                "surface_id": sid,
                "fatigue_T_K": float(Tf),
                "strength_rate_s-1": args.strength_rate,
                "stress_reference_temperature_K": float(args.stress_reference_temperature),
                "strength_response_class": cls,
                "relative_strength_regime_at_fatigue_T": relative_strength_regime(sd, float(Tf)),
                "strength_sigma_ref_MPa": float(sd.get("sigma_ref_MPa", np.nan)),
                "sigma_ref_nominal_MPa": sigma_ref,
                **{f"strength_{k}": v for k, v in sd.items() if k != "valid"},
                **{f"strength_rel_{k}": v for k, v in rel.items() if k != "valid"},
                **{f"no_shield_{k}": v for k, v in nd.items()},
                **{f"shielded_{k}": v for k, v in shd.items()},
                **pair,
            }
            summary["no_shield_fatigue_response_class"] = fatigue_response_class(nd)
            summary["shielded_fatigue_response_class"] = fatigue_response_class(shd)
            summaries.append(summary)
        if ii % 250 == 0:
            print(f"ANALYSIS descriptors {ii}/{len(design)}")
    return summaries


def build_surface_phenotypes(summary_rows: Sequence[dict], fatigue_temperatures: Sequence[float]) -> List[dict]:
    df = pd.DataFrame(summary_rows)
    out = []
    if df.empty:
        return out
    temps = sorted(float(x) for x in fatigue_temperatures)
    for sid, g in df.groupby("surface_id"):
        g = g.sort_values("fatigue_T_K")
        row = {
            "surface_id": sid,
            "strength_response_class": str(g.strength_response_class.iloc[0]),
            "n_fatigue_temperatures": len(g),
        }
        classes = []
        flags = []
        used_t = []
        for T in temps:
            h = g[np.isclose(g.fatigue_T_K.astype(float), T)]
            if h.empty:
                continue
            c = str(h.shielded_fatigue_response_class.iloc[0])
            row[f"shielded_fatigue_class_T{int(round(T))}"] = c
            classes.append(c); flags.append(c == "endurance_like"); used_t.append(T)
        if not flags:
            pattern = "unresolved"
        elif all(flags):
            pattern = "persistent_endurance"
        elif not any(flags):
            pattern = "no_endurance"
        else:
            idx = np.flatnonzero(flags)
            contiguous = len(idx) == (idx[-1] - idx[0] + 1)
            if contiguous and idx[0] == 0:
                pattern = "low_temperature_endurance"
            elif contiguous and idx[-1] == len(flags) - 1:
                pattern = "high_temperature_endurance"
            elif contiguous:
                pattern = "intermediate_temperature_window"
            else:
                pattern = "mixed_or_reentrant"
        row["fatigue_temperature_pattern"] = pattern
        row["endurance_temperature_fraction"] = float(np.mean(flags)) if flags else np.nan
        row["endurance_T_min_K"] = float(min(t for t, f in zip(used_t, flags) if f)) if any(flags) else np.nan
        row["endurance_T_max_K"] = float(max(t for t, f in zip(used_t, flags) if f)) if any(flags) else np.nan
        out.append(row)
    return out


def cramers_v_and_cells(df: pd.DataFrame, row_metric: str, col_metric: str, Tf: float, family: str):
    d = df[[row_metric, col_metric]].dropna().astype(str)
    if len(d) < 20:
        return None, []
    tab = pd.crosstab(d[row_metric], d[col_metric])
    if tab.shape[0] < 2 or tab.shape[1] < 2:
        return None, []
    chi2, p, dof, expected = chi2_contingency(tab.to_numpy(), correction=False)
    n = tab.to_numpy().sum()
    denom = max(min(tab.shape[0]-1, tab.shape[1]-1), 1)
    V = math.sqrt(max(chi2, 0.0) / max(n * denom, 1e-30))
    stat = {
        "analysis_family": family,
        "fatigue_T_K": float(Tf),
        "row_metric": row_metric,
        "column_metric": col_metric,
        "n": int(n),
        "cramers_V": float(V),
        "chi2": float(chi2),
        "chi2_p": float(p),
        "dof": int(dof),
    }
    cells = []
    obs = tab.to_numpy(float)
    expected = np.asarray(expected, float)
    z = (obs - expected) / np.sqrt(np.maximum(expected, 1e-30))
    enrich = np.log2((obs + 0.5) / (expected + 0.5))
    for i, rv in enumerate(tab.index):
        for j, cv in enumerate(tab.columns):
            cells.append({
                "analysis_family": family,
                "fatigue_T_K": float(Tf),
                "row_class": str(rv),
                "column_class": str(cv),
                "observed": int(obs[i, j]),
                "expected": float(expected[i, j]),
                "standardized_residual": float(z[i, j]),
                "log2_observed_over_expected": float(enrich[i, j]),
            })
    return stat, cells


def analyze_class_associations(summary_rows: Sequence[dict], fatigue_temperatures: Sequence[float]):
    df = pd.DataFrame(summary_rows)
    stats, cells = [], []
    if df.empty:
        return stats, cells
    for Tf in fatigue_temperatures:
        g = df[np.isclose(df.fatigue_T_K.astype(float), float(Tf))]
        for row_metric, family in [
            ("strength_response_class", "strength_class_vs_fatigue_class"),
            ("relative_strength_regime_at_fatigue_T", "relative_strength_regime_vs_fatigue_class"),
        ]:
            stat, cell = cramers_v_and_cells(g, row_metric, "shielded_fatigue_response_class", Tf, family)
            if stat is not None:
                stats.append(stat); cells.extend(cell)
    if stats:
        q = bh_qvalues([r["chi2_p"] for r in stats])
        for r, qq in zip(stats, q):
            r["chi2_q_BH"] = float(qq)
    return stats, cells


def analyze_association_vs_anomaly_cap(summary_rows: Sequence[dict], fatigue_temperatures: Sequence[float], caps: Sequence[float]) -> List[dict]:
    """Sensitivity of strength-class/fatigue-class association to anomaly-gain envelope."""
    df = pd.DataFrame(summary_rows)
    rows = []
    if df.empty:
        return rows
    for Tf in fatigue_temperatures:
        g0 = df[np.isclose(df.fatigue_T_K.astype(float), float(Tf))].copy()
        for cap in caps:
            g = g0[pd.to_numeric(g0["strength_anomaly_gain_frac"], errors="coerce").fillna(0.0) <= float(cap)].copy()
            stat, _ = cramers_v_and_cells(g, "strength_response_class", "shielded_fatigue_response_class", Tf, "strength_class_vs_fatigue_class_by_anomaly_cap")
            if stat is not None:
                stat["anomaly_gain_cap_frac"] = float(cap)
                rows.append(stat)
    if rows:
        q = bh_qvalues([r["chi2_p"] for r in rows])
        for r, qq in zip(rows, q):
            r["chi2_q_BH"] = float(qq)
    return rows


def global_phenotype_association(phenotype_rows: Sequence[dict]):
    df = pd.DataFrame(phenotype_rows)
    if df.empty:
        return [], []
    stat, cells = cramers_v_and_cells(df, "strength_response_class", "fatigue_temperature_pattern", -1.0, "strength_class_vs_global_fatigue_phenotype")
    if stat is None:
        return [], []
    stat["chi2_q_BH"] = stat["chi2_p"]
    return [stat], cells


def _pretty_metric(name: str) -> str:
    replacements = {
        "drop_fraction_x0_to3": "Barrier drop",
        "stress_slope_abs_median_norm_per_x": "Median |barrier slope|",
        "stress_flat_width_x_10pct_of_maxslope": "Flat-barrier width",
        "curvature_abs_integral_norm": "Integrated barrier curvature",
        "barrier_eta_G": "Thermal barrier coordinate",
        "barrier_eta_sigc": "Thermal stress coordinate",
        "shielded_sn_plateau_width_decades": "Shielded plateau width",
        "pair_median_loglife_benefit_lower_bound": "Median shielding life benefit",
        "pair_plateau_width_gain_decades": "Plateau-width gain",
        "pair_highcycle_slope_reduction": "High-cycle slope reduction",
        "shielded_sn_highcycle_slope_abs_norm_per_decade": "Shielded high-cycle slope",
    }
    for k, v in replacements.items():
        if name.endswith(k) or name == k:
            return v
    return name.replace("_", " ")


def make_plots(out: Path, strength_rows, sn_rows, summaries, correlations, class_stats, class_cells, fatigue_temperatures, phenotypes=None, cap_sensitivity=None):
    import matplotlib.pyplot as plt
    sdf = pd.DataFrame(summaries); sr = pd.DataFrame(strength_rows); nr = pd.DataFrame(sn_rows); cr = pd.DataFrame(correlations)

    # Linked representative projections: each row is one barrier surface.
    if not sdf.empty and not sr.empty and not nr.empty:
        target_T = 300.0 if any(np.isclose(np.asarray(fatigue_temperatures, float), 300.0)) else float(sorted(fatigue_temperatures)[0])
        reps = []
        for cls in ["strong_softening", "ordinary_softening", "plateau", "compensated", "lowT_anomaly", "midT_anomaly", "highT_anomaly"]:
            g = sdf[(sdf.strength_response_class == cls) & np.isclose(sdf.fatigue_T_K.astype(float), target_T)]
            if len(g):
                reps.append(str(g.iloc[len(g)//2].surface_id))
        reps = list(dict.fromkeys(reps))[:6]
        if reps:
            fig, axs = plt.subplots(len(reps), 2, figsize=(10.5, 2.3*len(reps)), squeeze=False)
            for i, sid in enumerate(reps):
                gs = sr[sr.surface_id.astype(str) == sid].sort_values("T_K")
                axs[i,0].plot(gs.T_K, gs.strength_MPa, lw=1.5)
                cls = str(gs.response_class.iloc[0]) if "response_class" in gs and len(gs) else ""
                axs[i,0].set_title(f"{sid}: {cls}", fontsize=9)
                axs[i,0].set_xlabel("Temperature (K)"); axs[i,0].set_ylabel("Strength (MPa)"); axs[i,0].grid(True, alpha=.25)
                gn = nr[(nr.surface_id.astype(str) == sid) & np.isclose(nr.T_K.astype(float), target_T)]
                for case in ["no_shield", "shielded"]:
                    h = gn[gn.case == case].sort_values("sigma_a_MPa")
                    if len(h):
                        N = h.cycles_to_nucleation.fillna(h.cycles_total)
                        axs[i,1].plot(N, h.sigma_a_MPa, marker="o", ms=2.5, lw=1, label=case)
                axs[i,1].set_xscale("log"); axs[i,1].set_xlabel("Cycles"); axs[i,1].set_ylabel("Stress amplitude (MPa)"); axs[i,1].grid(True, which="both", alpha=.25)
            axs[0,1].legend(fontsize=7)
            fig.suptitle(f"Linked barrier projections: strength trajectory and S-N response at {target_T:g} K")
            fig.tight_layout(); fig.savefig(out/"linked_strength_and_SN_representatives_v5_6.png", dpi=220); plt.close(fig)

    # Class-association heatmaps: counts and standardized residuals show enrichment/depletion.
    cdf = pd.DataFrame(class_cells); st = pd.DataFrame(class_stats)
    if not cdf.empty:
        fam = "strength_class_vs_fatigue_class"
        ts = sorted(cdf[cdf.analysis_family == fam].fatigue_T_K.unique())
        ncol = min(3, max(1, len(ts))); nrow = int(math.ceil(len(ts)/ncol))
        fig, axs = plt.subplots(nrow, ncol, figsize=(4.7*ncol, 4.0*nrow), squeeze=False)
        row_order = ["strong_softening", "ordinary_softening", "plateau", "compensated", "lowT_anomaly", "midT_anomaly", "highT_anomaly"]
        col_order = ["continuous_SN", "knee_threshold_like", "endurance_like"]
        im = None
        for ax, Tf in zip(axs.ravel(), ts):
            g = cdf[(cdf.analysis_family == fam) & np.isclose(cdf.fatigue_T_K.astype(float), float(Tf))]
            rows = [r for r in row_order if r in set(g.row_class)]
            cols = [c for c in col_order if c in set(g.column_class)]
            M = np.full((len(rows), len(cols)), np.nan); O = np.zeros_like(M)
            for i, rv in enumerate(rows):
                for j, cv in enumerate(cols):
                    h = g[(g.row_class == rv) & (g.column_class == cv)]
                    if len(h): M[i,j] = float(h.standardized_residual.iloc[0]); O[i,j] = float(h.observed.iloc[0])
            im = ax.imshow(M, aspect="auto", cmap="coolwarm", vmin=-4, vmax=4)
            ax.set_xticks(range(len(cols))); ax.set_xticklabels([c.replace("_", " ") for c in cols], rotation=25, ha="right", fontsize=8)
            ax.set_yticks(range(len(rows))); ax.set_yticklabels([r.replace("_", " ") for r in rows], fontsize=8)
            for i in range(len(rows)):
                for j in range(len(cols)):
                    if np.isfinite(M[i,j]): ax.text(j, i, f"n={int(O[i,j])}\nz={M[i,j]:.1f}", ha="center", va="center", fontsize=7)
            hs = st[(st.analysis_family == fam) & np.isclose(st.fatigue_T_K.astype(float), float(Tf))]
            subtitle = f"T = {Tf:g} K"
            if len(hs): subtitle += f"; Cramer's V={hs.cramers_V.iloc[0]:.2f}, q={hs.chi2_q_BH.iloc[0]:.2g}"
            ax.set_title(subtitle, fontsize=9)
        for ax in axs.ravel()[len(ts):]: ax.axis("off")
        fig.suptitle("Strength-response class versus fatigue-response class", y=.995)
        fig.subplots_adjust(top=.84 if nrow == 1 else .90, right=.88, wspace=.45, hspace=.60)
        if im is not None:
            cax = fig.add_axes([.91, .20, .018, .60])
            fig.colorbar(im, cax=cax, label="Standardized residual: enrichment (+), depletion (-)")
        fig.savefig(out/"strength_vs_fatigue_class_association_v5_6.png", dpi=220, bbox_inches="tight"); plt.close(fig)

    # Global strength-class versus fatigue-temperature phenotype association.
    gdf = pd.DataFrame(class_cells)
    gst = pd.DataFrame(class_stats)
    gfam = "strength_class_vs_global_fatigue_phenotype"
    gg = gdf[gdf.analysis_family == gfam] if not gdf.empty else pd.DataFrame()
    if not gg.empty:
        row_order = ["strong_softening", "ordinary_softening", "plateau", "lowT_anomaly", "midT_anomaly", "highT_anomaly"]
        col_order = ["no_endurance", "low_temperature_endurance", "intermediate_temperature_window", "high_temperature_endurance", "persistent_endurance", "mixed_or_reentrant"]
        rows = [r for r in row_order if r in set(gg.row_class)]
        cols = [c for c in col_order if c in set(gg.column_class)]
        M = np.full((len(rows), len(cols)), np.nan); O = np.zeros_like(M)
        for i, rv in enumerate(rows):
            for j, cv in enumerate(cols):
                h = gg[(gg.row_class == rv) & (gg.column_class == cv)]
                if len(h): M[i,j] = float(h.standardized_residual.iloc[0]); O[i,j] = float(h.observed.iloc[0])
        fig, ax = plt.subplots(figsize=(1.65*max(len(cols),3)+4.0, 0.65*max(len(rows),3)+2.5))
        im = ax.imshow(M, aspect="auto", cmap="coolwarm", vmin=-4, vmax=4)
        ax.set_xticks(range(len(cols))); ax.set_xticklabels([c.replace("_"," ") for c in cols], rotation=25, ha="right")
        ax.set_yticks(range(len(rows))); ax.set_yticklabels([r.replace("_"," ") for r in rows])
        for i in range(len(rows)):
            for j in range(len(cols)):
                if np.isfinite(M[i,j]): ax.text(j,i,f"n={int(O[i,j])}\nz={M[i,j]:.1f}",ha="center",va="center",fontsize=8)
        hs = gst[gst.analysis_family == gfam]
        title = "Strength phenotype versus fatigue-temperature phenotype"
        if len(hs): title += f"\nCramer's V={hs.cramers_V.iloc[0]:.2f}, p={hs.chi2_p.iloc[0]:.2g}"
        ax.set_title(title)
        fig.colorbar(im, ax=ax, label="Standardized residual: enrichment (+), depletion (-)")
        fig.tight_layout(); fig.savefig(out/"global_strength_vs_fatigue_phenotype_v5_6.png", dpi=220, bbox_inches="tight"); plt.close(fig)

    # Robustness of categorical association to anomaly-amplitude envelope.
    capdf = pd.DataFrame(cap_sensitivity or [])
    if not capdf.empty:
        fig, ax = plt.subplots(figsize=(7.2,4.8))
        for Tf, g in capdf.groupby("fatigue_T_K"):
            g = g.sort_values("anomaly_gain_cap_frac")
            ax.plot(100*g.anomaly_gain_cap_frac, g.cramers_V, marker="o", label=f"{Tf:g} K")
        ax.set_xlabel("Maximum anomaly gain included (%)")
        ax.set_ylabel("Cramer's V")
        ax.set_title("Robustness of strength-class/fatigue-class association")
        ax.grid(True, alpha=.25); ax.legend(title="Fatigue temperature")
        fig.tight_layout(); fig.savefig(out/"class_association_vs_anomaly_cap_v5_6.png", dpi=220); plt.close(fig)

    # Curated correlation matrices.  Spearman rho_s is shown explicitly; ns marks q >= 0.05.
    if not cr.empty:
        ts = sorted(float(x) for x in fatigue_temperatures)
        ncol = min(2, max(1, len(ts))); nrow = int(math.ceil(len(ts)/ncol))
        fig, axs = plt.subplots(nrow, ncol, figsize=(8.0*ncol, 4.8*nrow), squeeze=False)
        responses = [
            "shielded_sn_plateau_width_decades",
            "pair_median_loglife_benefit_lower_bound",
            "pair_plateau_width_gain_decades",
            "pair_highcycle_slope_reduction",
            "shielded_sn_highcycle_slope_abs_norm_per_decade",
        ]
        im = None
        for ax, Tf in zip(axs.ravel(), ts):
            tag = f"T{int(round(Tf))}"
            predictors = [
                "barrier_eta_G", "barrier_eta_sigc",
                f"barrier_{tag}_drop_fraction_x0_to3",
                f"barrier_{tag}_stress_slope_abs_median_norm_per_x",
                f"barrier_{tag}_stress_flat_width_x_10pct_of_maxslope",
                f"barrier_{tag}_curvature_abs_integral_norm",
            ]
            M = np.full((len(predictors), len(responses)), np.nan)
            Q = np.full_like(M, np.nan)
            for i, pred in enumerate(predictors):
                for j, resp in enumerate(responses):
                    h = cr[(cr.analysis_family == "barrier_shape_to_fatigue") & np.isclose(cr.fatigue_T_K.astype(float), Tf) & (cr.predictor_metric == pred) & (cr.response_metric == resp)]
                    if len(h): M[i,j] = float(h.spearman_rho.iloc[0]); Q[i,j] = float(h.spearman_q_BH.iloc[0])
            im = ax.imshow(M, aspect="auto", cmap="coolwarm", vmin=-1, vmax=1)
            ax.set_xticks(range(len(responses))); ax.set_xticklabels([_pretty_metric(x) for x in responses], rotation=30, ha="right", fontsize=8)
            ax.set_yticks(range(len(predictors))); ax.set_yticklabels([_pretty_metric(x) for x in predictors], fontsize=8)
            for i in range(len(predictors)):
                for j in range(len(responses)):
                    if np.isfinite(M[i,j]):
                        tagtxt = "*" if np.isfinite(Q[i,j]) and Q[i,j] < 0.05 else "ns"
                        ax.text(j, i, f"{M[i,j]:+.2f}\n{tagtxt}", ha="center", va="center", fontsize=7)
            ax.set_title(f"T = {Tf:g} K; Spearman $\\rho_s$ (*: BH q<0.05)")
        for ax in axs.ravel()[len(ts):]: ax.axis("off")
        fig.suptitle("Selected barrier-shape correlations with fatigue topology and shielding benefit", y=.995)
        fig.subplots_adjust(top=.84 if nrow == 1 else .91, bottom=.18, left=.22, right=.88, wspace=.40, hspace=.60)
        if im is not None:
            cax = fig.add_axes([.91, .20, .018, .60])
            fig.colorbar(im, cax=cax, label="Spearman $\\rho_s$")
        fig.savefig(out/"curated_barrier_fatigue_correlations_v5_6.png", dpi=220, bbox_inches="tight"); plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default="runs/sn_v1_barrier_phenomena_map_v5_6")
    ap.add_argument("--n-surfaces", type=int, default=3840, help="Target total surface count in the persistent design.")
    ap.add_argument("--candidate-pool", type=int, default=2048, help="Candidate count per design batch.")
    ap.add_argument("--design-batch-size", type=int, default=256)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--design-csv", default="", help="Optional starting design CSV; it may be extended to --n-surfaces.")
    ap.add_argument("--temperatures", nargs="+", type=float, default=[25,40,50,75,100,125,150,175,200,225,250,275,300,325,350,375,400,450,500,550,600,650,700,750,800,850,900,1000,1100,1200])
    ap.add_argument("--prescreen-temperatures", nargs="+", type=float, default=[25,50,75,100,150,200,250,300,400,500,600,750,900,1100,1200])
    ap.add_argument("--strength-rate", type=float, default=1e-4)
    ap.add_argument("--anomaly-gain-min", type=float, default=0.05, help="Minimum gain for the main anomaly strata.")
    ap.add_argument("--anomaly-gain-main-max", type=float, default=0.50, help="Maximum anomaly gain in the primary physically representative strata.")
    ap.add_argument("--anomaly-gain-sensitivity-max", type=float, default=1.00, help="Upper bound for the small stronger-anomaly sensitivity stratum.")
    ap.add_argument("--association-anomaly-caps", nargs="+", type=float, default=[0.25,0.50,1.00], help="Amplitude envelopes used for association robustness analysis.")
    ap.add_argument("--fatigue-temperatures", nargs="+", type=float, default=[100,200,300,400,500,600,700])
    ap.add_argument("--stress-reference-temperature", type=float, default=300.0,
                    help="Common strength temperature used only to scale the S-N stress grid at every fatigue temperature.")
    ap.add_argument("--cycles-max", type=float, default=1e12)
    ap.add_argument("--max-blocks", type=int, default=12000)
    ap.add_argument("--block-cycles", type=float, default=1e10)
    ap.add_argument("--n-phase", type=int, default=32)
    ap.add_argument("--target-dP", type=float, default=0.03)
    ap.add_argument("--target-dD", type=float, default=0.03)
    ap.add_argument("--target-rho-rel-block", type=float, default=0.15)
    ap.add_argument("--target-dB-nuc", type=float, default=0.20)
    ap.add_argument("--stress-fractions", nargs="+", type=float, default=[0.025,0.05,0.10,0.18,0.30,0.48,0.72,1.00,1.35])
    ap.add_argument("--refine-target-lives", nargs="+", type=float, default=[1e6,1e8,1e10,1e12])
    ap.add_argument("--refine-rounds", type=int, default=1)
    ap.add_argument("--min-stress-MPa", type=float, default=5.0)
    ap.add_argument("--max-stress-MPa", type=float, default=5000.0)
    ap.add_argument("--checkpoint-every", type=int, default=10)
    ap.add_argument("--resume", action="store_true", default=False)
    ap.add_argument("--analysis-only", action="store_true", default=False)
    ap.add_argument("--run-random-forest", action="store_true", default=False, help="Optional nonlinear screening; off by default for large expandable maps.")
    args = ap.parse_args()

    args.temperatures = sorted(set(float(x) for x in list(args.temperatures) + [float(args.stress_reference_temperature)]))
    args.fatigue_temperatures = sorted(set(float(x) for x in args.fatigue_temperatures))
    out = Path(args.out); out.mkdir(parents=True, exist_ok=True)
    design_path = out/"independent_exp_floor_design_v5_6.csv"

    base = sn_parser().parse_args([])
    base.cycles_max = args.cycles_max; base.max_blocks = args.max_blocks; base.block_cycles = args.block_cycles
    base.n_phase = args.n_phase; base.target_rho_rel_block = args.target_rho_rel_block; base.target_dB_nuc = args.target_dB_nuc
    base.target_dP = args.target_dP; base.target_dD = args.target_dD

    if design_path.exists():
        existing_design = pd.read_csv(design_path).to_dict("records")
    elif args.design_csv:
        existing_design = pd.read_csv(args.design_csv).to_dict("records")
    else:
        existing_design = []
    if existing_design:
        versions = {str(r.get("design_version", "")) for r in existing_design}
        if versions != {DESIGN_VERSION}:
            raise RuntimeError(f"Existing design version {versions} does not match {DESIGN_VERSION}. Use a new OUT directory for V5.6.")
        for fld, expected in [
            ("anomaly_gain_min_design", args.anomaly_gain_min),
            ("anomaly_gain_main_max_design", args.anomaly_gain_main_max),
            ("anomaly_gain_sensitivity_max_design", args.anomaly_gain_sensitivity_max),
        ]:
            vals = {round(float(r.get(fld, np.nan)), 12) for r in existing_design}
            if vals != {round(float(expected), 12)}:
                raise RuntimeError(f"Existing design field {fld}={vals} does not match requested {expected}. Use a new OUT directory or preserve the original design settings.")
    design = extend_independent_design(
        base, existing_design, args.n_surfaces, args.design_batch_size,
        args.candidate_pool, args.seed, args.prescreen_temperatures, args.strength_rate,
        args.anomaly_gain_min, args.anomaly_gain_main_max, args.anomaly_gain_sensitivity_max
    )
    write_csv(design_path, design)
    print(f"DESIGN total surfaces: {len(design)}")

    audit = []
    audit_ids = np.linspace(0, len(design)-1, min(32, len(design))).astype(int) if design else []
    for idx in audit_ids:
        r = design[int(idx)]
        ns = apply_barrier(SimpleNamespace(**vars(base)), r)
        audit.append({"surface_id": r["surface_id"], **barrier_override_audit(ns, r)})
    write_csv(out/"barrier_override_audit_v5_6.csv", audit)
    if audit and not all(bool(r["override_audit_pass"]) for r in audit):
        raise RuntimeError("Barrier override audit failed; refusing to run.")

    strength_path = out/"strength_temperature_single_rate.csv"
    sn_path = out/"sn_initiation_points_multiT_paired.csv"
    summary_path = out/"surface_temperature_descriptors_v5_6.csv"

    strength_rows = pd.read_csv(strength_path).to_dict("records") if args.resume and strength_path.exists() else []
    sn_rows = pd.read_csv(sn_path).to_dict("records") if args.resume and sn_path.exists() else []
    strength_index = {(str(r["surface_id"]), float(r["T_K"])): r for r in strength_rows}
    sn_index = {(str(r["surface_id"]), float(r["T_K"]), str(r["case"]), round(float(r["sigma_a_MPa"]), 9)): r for r in sn_rows}
    sn_groups = defaultdict(list)
    for r in sn_rows:
        sn_groups[(str(r["surface_id"]), float(r["T_K"]), str(r["case"]))].append(r)

    cases = {"no_shield": SNCase("no_shield", 0.0, 0.0), "shielded": SNCase("shielded", base.shield_chi, base.Gshield_eV)}

    if not args.analysis_only:
        for ii, r in enumerate(design, 1):
            sid = str(r["surface_id"])
            print(f"=== {ii}/{len(design)} {sid} ===")
            ns = apply_barrier(SimpleNamespace(**vars(base)), r)
            chain = build_chain_from_namespace(ns, ns.b_m)

            sr_local = []
            for T in args.temperatures:
                key = (sid, float(T))
                if key in strength_index:
                    rec = strength_index[key]
                else:
                    sval, status = strength_at_T(ns, float(T), args.strength_rate)
                    rec = {**r, "surface_id": sid, "strain_rate_s-1": args.strength_rate, "T_K": float(T), "strength_MPa": sval, "strength_status": status}
                    strength_rows.append(rec); strength_index[key] = rec
                sr_local.append(rec)
            sT = np.array([float(x["T_K"]) for x in sr_local]); ss = np.array([float(x["strength_MPa"]) for x in sr_local])
            sd = curve_shape_descriptors(sT, ss)
            cls = response_class(sd) if sd.get("valid", False) else "invalid"
            for rec in sr_local: rec["response_class"] = cls

            sigma_strength_ref = interp_at(sT, ss, args.stress_reference_temperature)
            if np.isfinite(sigma_strength_ref) and sigma_strength_ref > 0:
                sigma_ref = max(sigma_strength_ref/max(base.Kt, 1e-12), args.min_stress_MPa)
                sigma_ref_source = "strength_at_reference_temperature"
            else:
                b = chain.emit.base
                sigma_ref = max(float(b.sigc_Pa(args.stress_reference_temperature))/1e6/max(base.Kt,1e-12), args.min_stress_MPa)
                sigma_ref_source = "sigc_fallback"

            for Tf in args.fatigue_temperatures:
                coarse = sorted(set(float(np.clip(f*sigma_ref, args.min_stress_MPa, args.max_stress_MPa)) for f in args.stress_fractions))
                q = SimpleNamespace(**vars(ns)); q.T = float(Tf)
                for cname, case in cases.items():
                    group_key = (sid, float(Tf), cname)
                    pts = sn_groups[group_key]
                    have = {round(float(x["sigma_a_MPa"]), 9) for x in pts}
                    missing = [s for s in coarse if round(s,9) not in have]
                    if missing:
                        new = run_stress_grid(q, case, missing, chain=chain)
                        for rec in new:
                            rec.update({**r, "surface_id": sid, "sigma_ref_MPa": sigma_ref,
                                        "sigma_ref_source": sigma_ref_source,
                                        "stress_reference_temperature_K": args.stress_reference_temperature,
                                        "fatigue_T_K": float(Tf)})
                            sn_rows.append(rec); sn_index[(sid,float(Tf),cname,round(float(rec["sigma_a_MPa"]),9))] = rec; pts.append(rec)
                    for _ in range(args.refine_rounds):
                        add = refinement_stresses(pts, args.refine_target_lives)
                        if not add: break
                        new = run_stress_grid(q, case, add, chain=chain)
                        for rec in new:
                            rec.update({**r, "surface_id": sid, "sigma_ref_MPa": sigma_ref,
                                        "sigma_ref_source": sigma_ref_source,
                                        "stress_reference_temperature_K": args.stress_reference_temperature,
                                        "fatigue_T_K": float(Tf)})
                            sn_rows.append(rec); sn_index[(sid,float(Tf),cname,round(float(rec["sigma_a_MPa"]),9))] = rec; pts.append(rec)

            if ii % max(args.checkpoint_every, 1) == 0 or ii == len(design):
                print(f"CHECKPOINT surface {ii}: strength rows={len(strength_rows)}, S-N rows={len(sn_rows)}")
                write_csv(strength_path, strength_rows)
                write_csv(sn_path, sn_rows)

    # Analysis always rebuilds from raw checkpoints, allowing ANALYSIS_ONLY reruns.
    if not strength_path.exists() and strength_rows:
        write_csv(strength_path, strength_rows)
    if not sn_path.exists() and sn_rows:
        write_csv(sn_path, sn_rows)
    if args.analysis_only:
        if strength_path.exists(): strength_rows = pd.read_csv(strength_path).to_dict("records")
        if sn_path.exists(): sn_rows = pd.read_csv(sn_path).to_dict("records")

    summary_rows = rebuild_summary_rows(design, strength_rows, sn_rows, args.fatigue_temperatures, args, base)
    write_csv(summary_path, summary_rows)
    phenotypes = build_surface_phenotypes(summary_rows, args.fatigue_temperatures)
    write_csv(out/"surface_phenotype_summary_v5_6.csv", phenotypes)

    correlations = analyze_correlations(summary_rows, args.fatigue_temperatures)
    write_csv(out/"descriptor_correlations_v5_6.csv", correlations)
    if correlations:
        top = sorted(correlations, key=lambda r: r["abs_spearman"], reverse=True)
        write_csv(out/"top_descriptor_correlations_v5_6.csv", top[:300])

    class_stats, class_cells = analyze_class_associations(summary_rows, args.fatigue_temperatures)
    global_stats, global_cells = global_phenotype_association(phenotypes)
    class_stats = class_stats + global_stats
    class_cells = class_cells + global_cells
    write_csv(out/"class_association_statistics_v5_6.csv", class_stats)
    write_csv(out/"class_association_cells_v5_6.csv", class_cells)
    cap_sensitivity = analyze_association_vs_anomaly_cap(summary_rows, args.fatigue_temperatures, args.association_anomaly_caps)
    write_csv(out/"association_vs_anomaly_cap_v5_6.csv", cap_sensitivity)

    sdf = pd.DataFrame(summary_rows)
    coverage = []
    for c in sdf.columns:
        if not pd.api.types.is_numeric_dtype(sdf[c]): continue
        v = pd.to_numeric(sdf[c], errors="coerce").replace([np.inf,-np.inf], np.nan)
        coverage.append({"metric": c, "n_total": len(v), "n_finite": int(v.notna().sum()), "n_unique_finite": int(v.dropna().nunique()), "std_finite": float(v.dropna().std()) if v.notna().sum()>1 else np.nan})
    write_csv(out/"descriptor_coverage_v5_6.csv", coverage)

    rf_rows = []
    if args.run_random_forest and HAVE_SKLEARN and len(sdf) >= 80:
        predictors = [c for c in sdf.columns if (c.startswith("strength_") or c.startswith("barrier_")) and pd.api.types.is_numeric_dtype(sdf[c]) and c not in {"strength_sigma_ref_MPa", "strength_rel_sigma_Tf_MPa"}]
        targets = [c for c in ["pair_plateau_width_gain_decades", "pair_highcycle_slope_reduction", "pair_area_delta_sigma_fixedlife_norm_decades", "pair_median_loglife_benefit_lower_bound", "pair_max_loglife_benefit_lower_bound"] if c in sdf.columns]
        for Tf in args.fatigue_temperatures:
            g = sdf[np.isclose(sdf.fatigue_T_K.astype(float), float(Tf))]
            for target in targets:
                cols = [c for c in predictors if g[c].notna().mean() >= 0.8]
                d = g[cols+[target]].replace([np.inf,-np.inf], np.nan).dropna()
                if len(d) < 80 or d[target].nunique() < 6: continue
                X = d[cols].to_numpy(float); y = d[target].to_numpy(float)
                cv = KFold(n_splits=5, shuffle=True, random_state=42)
                model = RandomForestRegressor(n_estimators=300, min_samples_leaf=4, random_state=42, n_jobs=-1)
                scores = cross_val_score(model, X, y, cv=cv, scoring="r2"); model.fit(X,y)
                perm = permutation_importance(model, X, y, n_repeats=10, random_state=42, n_jobs=-1)
                for c, imp, std in zip(cols, perm.importances_mean, perm.importances_std):
                    rf_rows.append({"fatigue_T_K": Tf, "target": target, "predictor": c, "cv_R2_mean": float(np.mean(scores)), "cv_R2_std": float(np.std(scores)), "permutation_importance": float(imp), "importance_std": float(std), "n": len(d)})
    if rf_rows: write_csv(out/"random_forest_importance_v5_6.csv", rf_rows)
    else:
        if not args.run_random_forest:
            status, detail = "disabled_by_default", "Set RUN_RF=1 for optional nonlinear screening."
        elif not HAVE_SKLEARN:
            status, detail = "skipped_dependency", SKLEARN_IMPORT_ERROR
        else:
            status, detail = "no_eligible_target", ""
        write_csv(out/"optional_analysis_status.csv", [{"analysis":"random_forest", "status":status, "detail":detail}])

    make_plots(out, strength_rows, sn_rows, summary_rows, correlations, class_stats, class_cells, args.fatigue_temperatures, phenotypes=phenotypes, cap_sensitivity=cap_sensitivity)
    with (out/"run_config_v5_6.json").open("w") as f: json.dump(vars(args), f, indent=2)

    print("Top continuous correlations (rho_s = Spearman):")
    for r in sorted(correlations, key=lambda z: z["abs_spearman"], reverse=True)[:15]:
        print(r["analysis_family"], f"T={r['fatigue_T_K']:g}K", r["predictor_metric"], "<->", r["response_metric"], f"rho_s={r['spearman_rho']:.3f}", f"q={r['spearman_q_BH']:.3g}", f"n={r['n']}")
    print("Class associations:")
    for r in class_stats:
        print(r["analysis_family"], f"T={r['fatigue_T_K']:g}K", f"Cramers_V={r['cramers_V']:.3f}", f"q={r.get('chi2_q_BH', np.nan):.3g}", f"n={r['n']}")


if __name__ == "__main__":
    main()
