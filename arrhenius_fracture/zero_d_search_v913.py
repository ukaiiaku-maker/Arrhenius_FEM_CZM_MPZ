"""Sampling, vectorized proxy, and scoring helpers for the v9.13 zero-D search."""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd
from scipy.special import gammaincinv
from scipy.stats import qmc

from .emergent_gnd_contract_v913 import ACTIVE_CANDIDATE_PARAMETER_FIELDS
from .emergent_gnd_types_v912 import KB_EV_PER_K
from .zero_d_persistent_v913 import reduction_geometry

FIXED_ACTIVE_FIELDS = {
    "Tref_K": 481.33,
    "peierls_nu0_s": 1.0e12,
    "taylor_nu0_s": 1.0e11,
}
VARIABLE_FIELDS = tuple(
    name for name in ACTIVE_CANDIDATE_PARAMETER_FIELDS if name not in FIXED_ACTIVE_FIELDS
)


def _load_policy(path: Path) -> dict[str, Any]:
    policy = json.loads(path.read_text())
    dimensions = policy.get("search_dimensions")
    if not isinstance(dimensions, Mapping):
        raise RuntimeError("policy is missing search_dimensions")
    missing = [name for name in VARIABLE_FIELDS if name not in dimensions]
    extra = [name for name in dimensions if name not in VARIABLE_FIELDS]
    if missing or extra:
        raise RuntimeError(
            f"policy/current-active-field mismatch: missing={missing}, extra={extra}"
        )
    anchor_ids = policy.get("anchor_candidate_ids")
    if not isinstance(anchor_ids, list) or not anchor_ids:
        raise RuntimeError("policy requires nonempty anchor_candidate_ids")
    local_fraction = float(policy.get("local_anchor_fraction", 0.85))
    if not 0.0 <= local_fraction <= 1.0:
        raise ValueError("local_anchor_fraction must lie in [0,1]")
    return policy


def _dimension_transform(
    u: np.ndarray,
    anchor: np.ndarray,
    spec: Mapping[str, Any],
    local_mask: np.ndarray,
) -> np.ndarray:
    low = float(spec["low"])
    high = float(spec["high"])
    mode = str(spec["mode"])
    if not high > low:
        raise ValueError(f"invalid search range [{low}, {high}]")
    if mode == "linear_delta":
        half = float(spec["half_width"])
        local = anchor + (2.0 * u - 1.0) * half
        global_values = low + u * (high - low)
    elif mode == "log10_delta":
        if low <= 0.0 or high <= 0.0 or np.any(anchor <= 0.0):
            raise ValueError("log10_delta requires positive bounds and anchors")
        half = float(spec["half_width_decades"])
        local = anchor * np.power(10.0, (2.0 * u - 1.0) * half)
        global_values = np.power(
            10.0,
            math.log10(low) + u * (math.log10(high) - math.log10(low)),
        )
    else:
        raise ValueError(f"unsupported search mode {mode!r}")
    return np.clip(np.where(local_mask, local, global_values), low, high)


def _sample_rows(
    *,
    start: int,
    count: int,
    total_samples: int,
    seed: int,
    anchors: pd.DataFrame,
    policy: Mapping[str, Any],
) -> pd.DataFrame:
    dimensions = 2 + len(VARIABLE_FIELDS)
    sampler = qmc.Sobol(dimensions, scramble=True, seed=seed)
    if start:
        sampler.fast_forward(start)
    u = sampler.random(count)
    local_fraction = float(policy.get("local_anchor_fraction", 0.85))
    local_mask = u[:, 0] < local_fraction
    anchor_index = np.minimum(
        (u[:, 1] * len(anchors)).astype(int),
        len(anchors) - 1,
    )
    selected_anchor = anchors.iloc[anchor_index].reset_index(drop=True)
    data: dict[str, Any] = {
        "sample_index": np.arange(start, start + count, dtype=np.int64),
        "candidate_id": [
            f"v913_zeroD_sobol_{index:07d}" for index in range(start, start + count)
        ],
        "sample_mode": np.where(local_mask, "anchor_local", "global"),
        "anchor_candidate_id": selected_anchor["candidate_id"].astype(str).to_numpy(),
    }
    for name, value in FIXED_ACTIVE_FIELDS.items():
        data[name] = np.full(count, float(value), dtype=float)
    specs = policy["search_dimensions"]
    for column, name in enumerate(VARIABLE_FIELDS, start=2):
        anchor = pd.to_numeric(selected_anchor[name], errors="raise").to_numpy(float)
        data[name] = _dimension_transform(
            u[:, column], anchor, specs[name], local_mask
        )
    return pd.DataFrame(data)


def _surface_barrier(
    G00: np.ndarray,
    gT: np.ndarray,
    sigc0_GPa: np.ndarray,
    sT_GPa_per_K: np.ndarray,
    exp_a: np.ndarray,
    exp_n: np.ndarray,
    floor_fraction: np.ndarray,
    Tref: np.ndarray,
    stress_Pa: np.ndarray,
    temperature_K: float,
) -> np.ndarray:
    G0 = np.maximum(G00 + gT * (temperature_K - Tref), 1.0e-12)
    floor = np.minimum(0.95 * G0, np.maximum(1.0e-4, floor_fraction * G0))
    sigc = np.maximum(
        (sigc0_GPa + sT_GPa_per_K * (temperature_K - Tref)) * 1.0e9,
        1.0,
    )
    x = np.maximum(stress_Pa, 0.0) / sigc
    return np.maximum(
        floor
        + (G0 - floor)
        * np.exp(
            -np.maximum(exp_a, 0.0)
            * np.power(x, np.maximum(exp_n, 1.0e-9))
        ),
        0.0,
    )


def _neutral_cleavage_K(
    frame: pd.DataFrame,
    temperature_K: float,
    *,
    physics: Any,
    target_rate_s: float,
) -> np.ndarray:
    hits = float(physics.cleavage_hits)
    tau = float(physics.cleavage_correlation_time_s)
    if hits <= 1.0 + 1.0e-12:
        raw_target = np.full(len(frame), target_rate_s, dtype=float)
    else:
        probability = min(max(target_rate_s * tau, 1.0e-300), 1.0 - 1.0e-15)
        raw_target = np.full(
            len(frame),
            float(gammaincinv(hits, probability) / tau),
            dtype=float,
        )
    barrier_target = -KB_EV_PER_K * temperature_K * np.log(
        np.clip(raw_target / float(physics.cleavage_nu0_s), 1.0e-300, 1.0)
    )
    Tref = frame["Tref_K"].to_numpy(float)
    G0 = np.maximum(
        frame["cleave_G00_eV"].to_numpy(float)
        + frame["cleave_gT_eV_per_K"].to_numpy(float)
        * (temperature_K - Tref),
        1.0e-12,
    )
    floor_fraction = frame["cleave_floor_frac"].to_numpy(float)
    floor = np.minimum(0.95 * G0, np.maximum(1.0e-4, floor_fraction * G0))
    sigc = np.maximum(
        (
            frame["cleave_sigc0_GPa"].to_numpy(float)
            + frame["cleave_sT_GPa_per_K"].to_numpy(float)
            * (temperature_K - Tref)
        )
        * 1.0e9,
        1.0,
    )
    a = np.maximum(frame["cleave_exp_a"].to_numpy(float), 1.0e-12)
    n = np.maximum(frame["cleave_exp_n"].to_numpy(float), 1.0e-9)
    ratio = (barrier_target - floor) / np.maximum(G0 - floor, 1.0e-30)
    sigma = np.zeros(len(frame), dtype=float)
    finite = (ratio > 0.0) & (ratio < 1.0)
    sigma[finite] = sigc[finite] * np.power(
        -np.log(np.clip(ratio[finite], 1.0e-300, 1.0)) / a[finite],
        1.0 / n[finite],
    )
    sigma[ratio <= 0.0] = np.inf
    radius = float(physics.r0_m)
    return sigma * math.sqrt(2.0 * math.pi * radius) / 1.0e6


def _proxy_response_batch(
    frame: pd.DataFrame,
    temperatures: Sequence[float],
    *,
    physics: Any,
    loading_map: Any,
    target_rate_s: float,
    history_events: float,
    target_extension_m: float,
) -> pd.DataFrame:
    reduced = reduction_geometry(physics)
    projected = np.cumsum(np.asarray(loading_map.projected_advances_m, dtype=float))
    event_count = int(np.searchsorted(projected, target_extension_m, side="left") + 1)
    event_count = min(max(event_count, 1), loading_map.n_events)
    geometry_factors = np.asarray(
        loading_map.K_per_U_MPa_sqrt_m_per_m[:event_count], dtype=float
    )
    Kdot = max(
        float(np.median(geometry_factors)) * loading_map.displacement_rate_m_s,
        1.0e-12,
    )
    breakpoints = np.asarray(physics.emission_geometry_extension_m, dtype=float)
    factors = np.asarray(physics.emission_geometry_factors, dtype=float)
    if breakpoints.size:
        mask = breakpoints <= target_extension_m + 1.0e-15
        representative_factors = np.sqrt(np.mean(np.square(factors[mask]), axis=0))
    else:
        representative_factors = np.abs(
            np.asarray(physics.emission_schmid_factors, dtype=float)
        )
    nrows = len(frame)
    rho_per = np.asarray(reduced.density_increment_per_activation_m2, dtype=float)
    slip_per = np.asarray(reduced.slip_count_increment_per_activation, dtype=float)
    active_arc = float(reduced.active_arc_factor)
    resolved = 1.0 / math.sqrt(3.0)
    backstress_prefactor = (
        float(physics.persistent_backstress_scale)
        * float(physics.G_Pa)
        * abs(float(physics.b_m))
        / resolved
    )
    out = frame[
        [
            "sample_index",
            "candidate_id",
            "sample_mode",
            "anchor_candidate_id",
            *ACTIVE_CANDIDATE_PARAMETER_FIELDS,
        ]
    ].copy()
    Kcurves = []
    max_rho = np.zeros(nrows, dtype=float)
    max_radius = np.full(nrows, float(physics.r0_m), dtype=float)
    min_width = np.full(nrows, float(physics.reference_front_width_m), dtype=float)
    for temperature in temperatures:
        T = float(temperature)
        K0 = _neutral_cleavage_K(
            frame, T, physics=physics, target_rate_s=target_rate_s
        )
        sigma0 = K0 * 1.0e6 / math.sqrt(
            2.0 * math.pi * float(physics.r0_m)
        )
        exposure_s = np.maximum(K0 / Kdot, float(loading_map.nominal_dt_s))
        exposure_s *= max(float(history_events), 1.0)
        radius = np.full(nrows, float(physics.r0_m), dtype=float)
        width = np.full(nrows, float(physics.reference_front_width_m), dtype=float)
        rho_system = np.zeros((nrows, int(physics.n_systems)), dtype=float)
        slip_count = np.zeros((nrows, int(physics.n_systems)), dtype=float)
        for _iteration in range(4):
            multiplicity = (
                frame["rho_source0_m2"].to_numpy(float)
                * active_arc
                * radius
                * width
            )
            drive = sigma0[:, None] * representative_factors[None, :]
            barrier = _surface_barrier(
                frame["emit_G00_eV"].to_numpy(float)[:, None],
                frame["emit_gT_eV_per_K"].to_numpy(float)[:, None],
                frame["emit_sigc0_GPa"].to_numpy(float)[:, None],
                frame["emit_sT_GPa_per_K"].to_numpy(float)[:, None],
                frame["emit_exp_a"].to_numpy(float)[:, None],
                frame["emit_exp_n"].to_numpy(float)[:, None],
                frame["emit_floor_frac"].to_numpy(float)[:, None],
                frame["Tref_K"].to_numpy(float)[:, None],
                np.maximum(drive, 0.0),
                T,
            )
            emission_rate = float(physics.emission_nu0_s) * np.exp(
                np.clip(-barrier / (KB_EV_PER_K * T), -700.0, 0.0)
            )
            unconstrained = (
                multiplicity[:, None] * emission_rate * exposure_s[:, None]
            )
            rho_block = np.square(
                np.divide(
                    drive,
                    backstress_prefactor,
                    out=np.zeros_like(drive),
                    where=np.isfinite(drive),
                )
            )
            block_activations = np.maximum(
                rho_block / rho_per[None, :],
                0.0,
            )
            activations = np.minimum(unconstrained, block_activations)
            rho_system = activations * rho_per[None, :]
            slip_count = activations * slip_per[None, :]
            rho_width = np.maximum(
                float(physics.rho_forest_floor_m2)
                + np.sum(rho_system, axis=1),
                float(physics.reference_density_m2),
            )
            width = np.clip(
                float(physics.reference_front_width_m)
                * np.sqrt(float(physics.reference_density_m2) / rho_width),
                max(
                    float(physics.minimum_front_width_m),
                    abs(float(physics.b_m)),
                ),
                (
                    float(physics.maximum_front_width_m)
                    if float(physics.maximum_front_width_m) > 0.0
                    else float(physics.reference_front_width_m)
                ),
            )
            radius = np.maximum(
                float(physics.r0_m)
                + frame["c_blunt"].to_numpy(float)
                * abs(float(physics.b_m))
                * np.sum(slip_count, axis=1)
                * float(physics.blunting_slip_fraction),
                float(physics.r0_m),
            )
        Kproxy = K0 * np.sqrt(radius / float(physics.r0_m))
        Kproxy[~np.isfinite(Kproxy)] = np.nan
        Kproxy = np.minimum(Kproxy, 200.0)
        tag = f"{T:g}".replace(".", "p")
        out[f"proxy_K_T{tag}"] = Kproxy
        Kcurves.append(Kproxy)
        max_rho = np.maximum(max_rho, np.sum(rho_system, axis=1))
        max_radius = np.maximum(max_radius, radius)
        min_width = np.minimum(min_width, width)
    curves = np.column_stack(Kcurves)
    metrics = _curve_metrics_matrix(
        np.asarray(temperatures, dtype=float),
        curves,
        peak_min=850.0,
        peak_max=1100.0,
    )
    for name, values in metrics.items():
        out[f"proxy_{name}"] = values
    out["proxy_max_tip_density_m2"] = max_rho
    out["proxy_max_tip_radius_um"] = max_radius * 1.0e6
    out["proxy_min_front_width_um"] = min_width * 1.0e6
    return out


def _curve_metrics_matrix(
    temperatures: np.ndarray,
    curves: np.ndarray,
    *,
    peak_min: float,
    peak_max: float,
) -> dict[str, np.ndarray]:
    nrows = curves.shape[0]
    peak_index = np.full(nrows, -1, dtype=int)
    prominence = np.full(nrows, -np.inf, dtype=float)
    for index in range(1, curves.shape[1] - 1):
        local = (curves[:, index] > curves[:, index - 1]) & (
            curves[:, index] > curves[:, index + 1]
        )
        candidate_prominence = np.minimum(
            curves[:, index] - curves[:, index - 1],
            curves[:, index] - curves[:, index + 1],
        )
        replace = local & (candidate_prominence > prominence)
        peak_index[replace] = index
        prominence[replace] = candidate_prominence[replace]
    fallback = peak_index < 0
    if np.any(fallback):
        fallback_curves = np.where(
            np.isfinite(curves[fallback]), curves[fallback], -np.inf
        )
        peak_index[fallback] = np.argmax(fallback_curves, axis=1)
    rows = np.arange(nrows)
    peak_value = curves[rows, peak_index]
    peak_temperature = temperatures[peak_index]
    peak_internal = (peak_index > 0) & (peak_index < curves.shape[1] - 1)
    prominence[~peak_internal] = -np.inf
    post_min = np.full(nrows, np.nan)
    post_max = np.full(nrows, np.nan)
    for row, index in enumerate(peak_index):
        if index + 1 < curves.shape[1]:
            post = curves[row, index + 1 :]
            post_min[row] = np.nanmin(post)
            post_max[row] = np.nanmax(post)
    drop = peak_value - post_min
    rebound = post_max - peak_value
    desired = (
        peak_internal
        & (peak_temperature >= peak_min)
        & (peak_temperature <= peak_max)
    )
    return {
        "peak_temperature_K": peak_temperature,
        "peak_value_MPa_sqrt_m": peak_value,
        "two_sided_prominence_MPa_sqrt_m": prominence,
        "post_peak_drop_MPa_sqrt_m": drop,
        "high_temperature_rebound_MPa_sqrt_m": rebound,
        "peak_internal": peak_internal.astype(int),
        "peak_in_desired_window": desired.astype(int),
    }


def _score_frame(
    frame: pd.DataFrame,
    prefix: str,
    *,
    minimum_prominence: float,
    minimum_drop: float,
    maximum_rebound: float,
    peak_min: float,
    peak_max: float,
) -> pd.DataFrame:
    out = frame.copy()
    prominence = pd.to_numeric(
        out[f"{prefix}_two_sided_prominence_MPa_sqrt_m"], errors="coerce"
    ).to_numpy(float)
    drop = pd.to_numeric(
        out[f"{prefix}_post_peak_drop_MPa_sqrt_m"], errors="coerce"
    ).to_numpy(float)
    rebound = pd.to_numeric(
        out[f"{prefix}_high_temperature_rebound_MPa_sqrt_m"], errors="coerce"
    ).to_numpy(float)
    Tpeak = pd.to_numeric(
        out[f"{prefix}_peak_temperature_K"], errors="coerce"
    ).to_numpy(float)
    peak_value = pd.to_numeric(
        out[f"{prefix}_peak_value_MPa_sqrt_m"], errors="coerce"
    ).to_numpy(float)
    internal = out[f"{prefix}_peak_internal"].astype(bool).to_numpy()
    finite = (
        np.isfinite(prominence)
        & np.isfinite(drop)
        & np.isfinite(rebound)
        & np.isfinite(Tpeak)
        & np.isfinite(peak_value)
    )
    target_center = 0.5 * (peak_min + peak_max)
    score = (
        3.0 * np.maximum(minimum_prominence - prominence, 0.0)
        + 2.0 * np.maximum(minimum_drop - drop, 0.0)
        + 4.0 * np.maximum(rebound - maximum_rebound, 0.0)
        + np.abs(Tpeak - target_center) / 100.0
        + 0.5 * np.maximum(25.0 - peak_value, 0.0)
        + 0.2 * np.maximum(peak_value - 120.0, 0.0)
        + np.where(internal, 0.0, 100.0)
        + np.where(finite, 0.0, 1000.0)
    )
    passed = (
        finite
        & internal
        & (Tpeak >= peak_min)
        & (Tpeak <= peak_max)
        & (prominence >= minimum_prominence)
        & (drop >= minimum_drop)
        & (rebound <= maximum_rebound)
    )
    out[f"{prefix}_objective"] = score
    out[f"{prefix}_gate_pass"] = passed
    return out
