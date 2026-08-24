#!/usr/bin/env python3
"""Generate the candidate-independent sequential PF sharp-wake mechanics map."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_outputs/oneD_v2_mechanics_maps_and_baselines"
PF = Path("/private/tmp/oneD-v2-pf-driver-diagnostics")
sys.path.insert(0, str(PF))

from arrhenius_fracture.config import ElasticProperties, GeometryConfig, JIntegralConfig, MeshConfig
from arrhenius_fracture.crack_backend import SharpWakeBackend
from arrhenius_fracture.crystal import cubic_plane_strain_D
from arrhenius_fracture.fem import assemble_mechanics, elastic_energy_densities, solve_dirichlet
from arrhenius_fracture.j_integral import compute_J_integral
from arrhenius_fracture.mesh import make_boundary_data, make_tri_mesh


SOURCE_COMMIT = "ab6279331d050919f78e2e7f278bc332466e34f8"
OBSERVER_COMMIT = "c695b44"
TARGETS_UM = (0, 2, 5, 10, 25, 50, 75, 100, 150, 200, 300, 400, 500, 600, 700, 800, 900, 1000)
REFERENCE_OPENING = 1.0e-5
SECOND_OPENING = 5.0e-6
EVENT_LENGTH = 5.0e-6


def initial_notch_damage(mesh, geometry):
    x, y = np.asarray(mesh.nodes, float).T
    damage = np.zeros(mesh.nn, dtype=float)
    damage[(x <= geometry.a0) & (np.abs(y) <= geometry.notch_half_thickness)] = 1.0
    return damage


def digest_arrays(*arrays) -> str:
    h = hashlib.sha256()
    for array in arrays:
        value = np.ascontiguousarray(array)
        h.update(str(value.dtype).encode()); h.update(str(value.shape).encode()); h.update(value.tobytes())
    return h.hexdigest()


def solve_state(mesh, damage, opening, D, material, tip, segments, kill_radius):
    boundary = make_boundary_data(mesh, GeometryConfig())
    u0 = np.zeros(mesh.ndof); ep = np.zeros((3, mesh.ne)); rho = np.zeros(mesh.ne)
    K, R0, *_ = assemble_mechanics(mesh, u0, ep, rho, damage, D, material, kappa=1e-6)
    u, reaction = solve_dirichlet(K, R0, u0, boundary, .5 * opening, -.5 * opening)
    _, _, sigma, *_ = assemble_mechanics(mesh, u, ep, rho, damage, D, material, kappa=1e-6)
    stored, _ = elastic_energy_densities(mesh, u, ep, sigma, D)
    energy = float(np.sum(stored * mesh.area_e))
    _, _, info = compute_J_integral(
        mesh, u, sigma, stored, damage, tip, np.array([1.0, 0.0]), material,
        ell=80e-6, cfg=JIntegralConfig(), crack_segments=segments,
        exclude_radius=2.0 * kill_radius,
    )
    jsigned = float(info.get("J_signed", info.get("J", 0.0)) or 0.0)
    jeff = max(jsigned, 0.0)
    kj = float(np.sqrt(jeff * material.Eprime))
    domain_values = {key: value for key, value in info.items()
                     if isinstance(value, (int, float, np.integer, np.floating, bool))}
    return {
        "reaction_N_per_m": float(reaction), "compliance_m2_per_N": float(opening / reaction),
        "elastic_energy_J_per_m": energy, "native_J_J_per_m2": float(jeff),
        "native_J_signed_J_per_m2": float(jsigned), "native_KJ_MPa_sqrt_m": float(kj / 1e6),
        "domain_metadata_json": json.dumps(domain_values, sort_keys=True),
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    geometry = GeometryConfig(); mesh_cfg = MeshConfig(nx=36, ny=72, jitter=0.0, tip_h_fine=1e-6, tip_ratio=1.2)
    mesh = make_tri_mesh(geometry, mesh_cfg, seed=42, tip_center=np.array([geometry.a0, 0.0]))
    material = ElasticProperties()
    D = cubic_plane_strain_D(523e9, 203e9, 160e9, 0.0)
    damage = initial_notch_damage(mesh, geometry)
    backend = SharpWakeBackend(); boundary = make_boundary_data(mesh, geometry); displacement = np.zeros(mesh.ndof)
    kill_radius = max(float(mesh.hbar_tip), 0.5e-6)
    mesh_sha = digest_arrays(mesh.nodes, mesh.elems)
    requested_counts = sorted(set(int(np.ceil(target / 5.0 - 1e-14)) for target in TARGETS_UM))
    target_by_count = {count: [target for target in TARGETS_UM if int(np.ceil(target / 5.0 - 1e-14)) == count]
                       for count in requested_counts}
    segments = [(np.array([0.0, 0.0]), np.array([geometry.a0, 0.0]))]
    rows = []
    for count in range(max(requested_counts) + 1):
        if count > 0:
            p0 = np.array([geometry.a0 + (count - 1) * EVENT_LENGTH, 0.0])
            p1 = np.array([geometry.a0 + count * EVENT_LENGTH, 0.0])
            result = backend.advance(mesh=mesh, boundary=boundary, damage=damage,
                                     displacement=displacement, p0=p0, p1=p1,
                                     direction=np.array([1.0, 0.0]), front_id=0,
                                     kill_r=kill_radius)
            if not result.inserted:
                raise RuntimeError(f"sharp-wake advance {count} failed: {result.reason}")
            damage = result.damage
            segments.append((p0, p1))
        if count not in target_by_count:
            continue
        tip = np.array([geometry.a0 + count * EVENT_LENGTH, 0.0])
        reference = solve_state(mesh, damage, REFERENCE_OPENING, D, material, tip, segments, kill_radius)
        repeat = solve_state(mesh, damage, SECOND_OPENING, D, material, tip, segments, kill_radius)
        for target in target_by_count[count]:
            rows.append({
                "target_extension_um": target, "actual_extension_um": count * 5.0,
                "geometry_event_count": count, "reference_opening_m": REFERENCE_OPENING,
                **reference,
                "F_over_U_N_per_m2": reference["reaction_N_per_m"] / REFERENCE_OPENING,
                "J_native_over_U2": reference["native_J_J_per_m2"] / REFERENCE_OPENING**2,
                "KJ_native_over_U": reference["native_KJ_MPa_sqrt_m"] / REFERENCE_OPENING,
                "wake_width_m": 2.0 * kill_radius, "mesh_scale_m": float(mesh.hbar_tip),
                "tip_x_m": float(tip[0]), "tip_y_m": 0.0, "mesh_fingerprint": mesh_sha,
                "damage_wake_fingerprint": digest_arrays(damage),
                "opening_scaling_F_ratio_error": abs((repeat["reaction_N_per_m"] / reference["reaction_N_per_m"]) / .5 - 1),
                "opening_scaling_W_ratio_error": abs((repeat["elastic_energy_J_per_m"] / reference["elastic_energy_J_per_m"]) / .25 - 1),
                "opening_scaling_J_ratio_error": abs((repeat["native_J_J_per_m2"] / reference["native_J_J_per_m2"]) / .25 - 1),
                "opening_scaling_KJ_ratio_error": abs((repeat["native_KJ_MPa_sqrt_m"] / reference["native_KJ_MPa_sqrt_m"]) / .5 - 1),
            })
    frame = pd.DataFrame(rows)
    frame.to_csv(OUT / "oneD_v2_pf_mechanics_map.csv", index=False)
    unique = frame.drop_duplicates("actual_extension_um").sort_values("actual_extension_um")
    validations = []
    for column in ("F_over_U_N_per_m2", "J_native_over_U2", "KJ_native_over_U"):
        x = unique.actual_extension_um.to_numpy(); y = unique[column].to_numpy()
        for i in range(1, len(x) - 1):
            predicted = np.interp(x[i], np.delete(x, i), np.delete(y, i))
            validations.append({"provider": "PF", "quantity": column, "left_out_extension_um": x[i],
                                "observed": y[i], "predicted": predicted,
                                "absolute_error": abs(predicted-y[i]),
                                "relative_error": abs(predicted-y[i])/max(abs(y[i]), 1e-300)})
    validation = pd.DataFrame(validations)
    validation.to_csv(OUT / "oneD_v2_map_interpolation_validation.csv", index=False)
    scaling_max = float(frame[[c for c in frame if c.startswith("opening_scaling_")]].max().max())
    status = "PF_PRODUCTION_DISCRETE_MAP_QUALIFIED" if scaling_max < 1e-9 else "PF_PRODUCTION_DISCRETE_MAP_UNQUALIFIED"
    report = f"""# One-dimensional V2 PF mechanics-map qualification V3

Status: **{status}**.

The map applies the actual `SharpWakeBackend.advance` operation sequentially on
the production 36×72 graded plane-strain mesh. Each standardized straight event
is the source-defined 5 µm segment; nominal 2 µm therefore resolves to the first
complete 5 µm event and no segment is truncated. Binary killed-element damage is
accumulated, never independently reconstructed at each node.

The production cluster-domain J convention (`ell=80 µm`, exclusion radius twice
the kill radius, signed/effective aggregation) is retained. This is a
**model-native discrete map, not continuum G**. Maximum opening-scaling error is
`{scaling_max:.3e}`. Bounded piecewise-linear interpolation is valid only over
0–1000 µm; extrapolation fails closed. Maximum LOO relative error is
`{validation.relative_error.max():.3e}` and is retained as discrete interpolation
uncertainty rather than smoothed away.
"""
    (OUT / "ONE_D_V2_PF_MECHANICS_MAP_QUALIFICATION_V3.md").write_text(report)
    fig, axes = plt.subplots(3, 1, figsize=(8, 9), sharex=True, constrained_layout=True)
    for ax, column, label in zip(axes, ("F_over_U_N_per_m2", "J_native_over_U2", "KJ_native_over_U"),
                                 ("F/U (N m⁻²)", "native J/U²", "native KJ/U")):
        ax.plot(unique.actual_extension_um, unique[column], "o-"); ax.set_ylabel(label); ax.grid(alpha=.25)
    axes[-1].set_xlabel("Actual standardized extension (µm)")
    fig.suptitle("PF production-discrete sequential sharp-wake map coefficients")
    fig.savefig(OUT / "PF_MECHANICS_MAP_COEFFICIENTS.png", dpi=180); plt.close(fig)


if __name__ == "__main__":
    main()
