#!/usr/bin/env python3
"""Build physics-informed ML features for the v9.12 campaign."""
from __future__ import annotations

import argparse
import csv
from dataclasses import replace
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from arrhenius_fracture.emergent_gnd_dbtt_v912 import (
    CommonPhysics,
    EmergentGNDState,
    candidate_from_registry_row,
)
from arrhenius_fracture.emergent_gnd_types_v912 import KB_EV_PER_K


def args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--candidate-registry", required=True)
    p.add_argument("--ranking-csv")
    p.add_argument("--physics-json", required=True)
    p.add_argument("--bounds-json", required=True)
    p.add_argument("--temperatures", nargs="+", type=float, default=[300, 600, 800, 1000, 1200])
    p.add_argument("--K-values", nargs="+", type=float, default=[15, 25, 35])
    p.add_argument("--reference-time-s", type=float, default=8.4)
    p.add_argument("--out", required=True)
    return p.parse_args()


def read_csv(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open(newline="") as fp:
        return list(csv.DictReader(fp))


def load_physics(path: str | Path) -> CommonPhysics:
    payload = json.loads(Path(path).read_text()).get("common_physics", {})
    for key in ("emission_signs", "emission_schmid_factors", "shielding_orientation_factors"):
        if key in payload:
            payload[key] = tuple(payload[key])
    for key in ("forest_interaction_matrix", "gnd_stress_projection_matrix"):
        if key in payload:
            payload[key] = tuple(tuple(row) for row in payload[key])
    return replace(CommonPhysics(**payload), n_bins=1)


def log10_positive(value: float, floor: float = 1.0e-300) -> float:
    return math.log10(max(abs(float(value)), floor))


def features(row: dict[str, str], physics: CommonPhysics, bound_names: list[str], Ts: list[float], Ks: list[float], tref: float) -> dict[str, Any]:
    candidate = candidate_from_registry_row(row)
    state = EmergentGNDState(candidate, physics)
    out: dict[str, Any] = {"candidate_id": candidate.candidate_id}
    for name in bound_names:
        out[f"x_raw__{name}"] = float(row[name])

    forest = physics.rho_forest_floor_m2
    spacing = 1.0 / (2.0 * math.sqrt(forest))
    for T in Ts:
        kBT = KB_EV_PER_K * T
        for K in Ks:
            tag = f"T{int(round(T))}_K{K:g}".replace(".", "p")
            sigma = K * 1.0e6 / math.sqrt(2.0 * math.pi * physics.r0_m)
            tau = float(physics.emission_schmid_factors[0]) * sigma
            Ge = float(candidate.emission.barrier_eV(max(tau, 0.0), T))
            Gc = float(candidate.cleavage.barrier_eV(max(sigma, 0.0), T))
            p_surface = candidate.peierls.surface(candidate.emission)
            Gp = float(p_surface.barrier_eV(max(candidate.peierls.stress_fraction * tau, 0.0), T))
            t_surface = candidate.taylor.surface(candidate.emission)
            tau_t = candidate.taylor.stress_fraction * tau * spacing / physics.b_m
            Gt = float(t_surface.barrier_eV(max(tau_t, 0.0), T))
            rates = state.local_rates(K, T)
            emit = float(np.max(rates["emission_rate_s"]))
            velocity = float(np.max(np.abs(rates["velocity_m_s"])))
            encounter = float(np.max(rates["encounter_s"]))
            release = float(np.max(rates["taylor_completion_s"]))
            recovery = float(rates["recovery_rate_s"])

            out[f"x_phys__{tag}__dG_emit_minus_cleave_kBT"] = (Ge - Gc) / max(kBT, 1.0e-30)
            out[f"x_phys__{tag}__dG_peierls_minus_emit_kBT"] = (Gp - Ge) / max(kBT, 1.0e-30)
            out[f"x_phys__{tag}__dG_taylor_minus_peierls_kBT"] = (Gt - Gp) / max(kBT, 1.0e-30)
            out[f"x_phys__{tag}__log10_emit_rate_s"] = log10_positive(emit)
            out[f"x_phys__{tag}__log10_velocity_m_s"] = log10_positive(velocity)
            out[f"x_phys__{tag}__log10_Pi_store"] = log10_positive(encounter * tref)
            out[f"x_phys__{tag}__log10_Pi_release"] = log10_positive(release * tref)
            out[f"x_phys__{tag}__log10_Pi_recovery"] = log10_positive(recovery * tref)

    source_inventory = candidate.rho_source0_m2 * physics.mpz_length_m * physics.active_strip_width_m * physics.n_systems
    out["x_phys__log10_source_inventory_per_unit_thickness"] = log10_positive(source_inventory)
    out["x_phys__log10_refresh_over_mpz"] = log10_positive(candidate.source_refresh_length_m / physics.mpz_length_m)
    corr_length = candidate.taylor_corr_scale / (2.0 * math.sqrt(candidate.taylor_corr_rho_c_m2))
    out["x_phys__log10_corr_order_increment_at_floor"] = log10_positive(2.0 * corr_length * math.sqrt(forest))
    return out


def main() -> int:
    a = args()
    rows = read_csv(a.candidate_registry)
    ranking = {r["candidate_id"]: r for r in read_csv(a.ranking_csv)} if a.ranking_csv else {}
    physics = load_physics(a.physics_json)
    bounds = json.loads(Path(a.bounds_json).read_text()).get("search_bounds", {})
    records = []
    for row in rows:
        rec = features(row, physics, list(bounds), a.temperatures, a.K_values, a.reference_time_s)
        rec.update({f"y__{k}": v for k, v in ranking.get(rec["candidate_id"], {}).items() if k != "candidate_id"})
        records.append(rec)
    fields = list(records[0])
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)
    print(f"ML_TABLE rows={len(records)} features={sum(k.startswith('x_') for k in fields)} out={out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
