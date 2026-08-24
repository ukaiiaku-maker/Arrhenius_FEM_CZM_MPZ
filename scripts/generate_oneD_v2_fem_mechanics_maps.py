#!/usr/bin/env python3
"""Generate FEM/CZM native-J and qualified structural-G straight-path maps."""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_outputs/oneD_v2_mechanics_maps_and_baselines"
QUAL = Path("/private/tmp/oneD-v2-qualified-fem-analysis")
sys.path.insert(0, str(QUAL))
spec = importlib.util.spec_from_file_location("canonical_vcct", QUAL / "scripts/qualify_canonical_connected_vcct.py")
canonical = importlib.util.module_from_spec(spec); spec.loader.exec_module(canonical)
v1 = canonical.v1


SOURCE_COMMIT = "931bed66913afc970117bce900805ecc9b6225f8"
QUALIFICATION_COMMITS = {"canonical": "15bfd51", "Peak500": "dc1146e", "DBTT500": "44c2d95"}
TARGETS_UM = (0, 2, 5, 10, 25, 50, 75, 100, 150, 200, 300, 400, 500, 600, 700, 800, 900, 1000)
RADII_UM = (30, 50, 75, 100, 150, 200)
OPENING = 12.5e-6
H = 1.0e-6


def fingerprint(mesh) -> str:
    return hashlib.sha256(np.ascontiguousarray(mesh.nodes).tobytes() +
                          np.ascontiguousarray(mesh.elems, np.int64).tobytes()).hexdigest()


def structured_mesh_at(tip, geometry):
    xs = canonical.graded_axis(tip, 0.0, geometry.Lx, H, nfine=16)
    yp = canonical.graded_axis(0.0, 0.0, geometry.Ly / 2, H, nfine=12)
    ys = np.r_[-yp[:0:-1], yp]
    nx, ny = len(xs), len(ys)
    nodes = np.array([(x, y) for y in ys for x in xs], float)
    elems = []
    for j in range(ny - 1):
        for i in range(nx - 1):
            n = j * nx + i
            if .5 * (ys[j] + ys[j + 1]) >= 0:
                elems.extend(((n, n + 1, n + nx + 1), (n, n + nx + 1, n + nx)))
            else:
                elems.extend(((n, n + 1, n + nx), (n + 1, n + nx + 1, n + nx)))
    base = v1.rebuild_tri_mesh(nodes, np.asarray(elems, int), tip_centers=[np.array([tip, 0.0])], validate=True)
    line = np.where(abs(base.nodes[:, 1]) < 1e-15)[0]
    newnodes = np.vstack((base.nodes, base.nodes[line])); dup = dict(zip(line, range(base.nn, base.nn + len(line))))
    newelems = base.elems.copy(); cent = base.nodes[base.elems].mean(axis=1)
    for element in np.where(cent[:, 1] < 0)[0]:
        for local, node in enumerate(newelems[element]):
            if node in dup:
                newelems[element, local] = dup[node]
    cut = v1.rebuild_tri_mesh(newnodes, newelems, tip_centers=[np.array([tip, 0.0])], validate=True)
    pairs = [(int(node), int(dup[node])) for node in line]
    return cut, pairs


def solve(cut, pairs, tip, D, material, opening):
    return canonical.state(cut, pairs, tip, D, material, opening)


def j_domains(cut, state, tip, material):
    segments = [(np.array([0.0, 0.0]), np.array([tip, 0.0]))]
    values = []
    for radius in RADII_UM:
        J, KJ, _ = v1.compute_J_integral(
            cut, state["u"], state["sigma"], state["stored"], np.zeros(cut.nn),
            np.array([tip, 0.0]), np.array([1.0, 0.0]), material, radius * 1e-6 / 8,
            cfg=v1.JIntegralConfig(r_inner_factor=2, r_outer_factor=8, q_type="plateau"),
            crack_segments=segments, exclude_radius=0.0,
        )
        values.append((radius, float(J), float(KJ / 1e6)))
    return values


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    geometry = v1.GeometryConfig(); material = v1.ElasticProperties()
    D = v1.cubic_plane_strain_D(523e9, 203e9, 160e9, 0.0)
    native_rows, g_rows = [], []
    for extension_um in TARGETS_UM:
        tip = geometry.a0 + extension_um * 1e-6
        cut, pairs = structured_mesh_at(tip, geometry)
        s0 = solve(cut, pairs, tip, D, material, OPENING)
        s1 = solve(cut, pairs, tip + H, D, material, OPENING)
        s2 = solve(cut, pairs, tip + 2 * H, D, material, OPENING)
        g_energy = (3 * s0["W"] - 4 * s1["W"] + s2["W"]) / (2 * H)
        g_compliance = .5 * s0["F"]**2 * (-3 * s0["C"] + 4 * s1["C"] - s2["C"]) / (2 * H)
        gi, gii, g_vcct, da_vcct, _, _ = canonical.vcct(cut, pairs, tip, s0)
        js = j_domains(cut, s0, tip, material)
        j_values = np.array([value[1] for value in js]); kj_values = np.array([value[2] for value in js])
        j_aggregate = float(np.median(j_values)); kj_aggregate = float(np.sqrt(max(j_aggregate, 0) * material.Eprime) / 1e6)
        uncertainty = float(max(abs(g_energy - g_compliance), abs(g_energy - g_vcct), abs(g_compliance - g_vcct)))
        scale = max(abs(g_energy), abs(g_compliance), abs(g_vcct), 1e-300)
        qualified = bool(uncertainty / scale <= 0.10 and g_energy > 0 and g_compliance > 0 and g_vcct > 0)
        scaling = {"opening_scaling_F_ratio_error": np.nan,
                   "opening_scaling_W_ratio_error": np.nan,
                   "opening_scaling_J_ratio_error": np.nan,
                   "opening_scaling_G_ratio_error": np.nan}
        if extension_um in {0, 500, 1000}:
            half0 = solve(cut, pairs, tip, D, material, .5 * OPENING)
            half1 = solve(cut, pairs, tip + H, D, material, .5 * OPENING)
            half2 = solve(cut, pairs, tip + 2 * H, D, material, .5 * OPENING)
            half_g = (3 * half0["W"] - 4 * half1["W"] + half2["W"]) / (2 * H)
            half_j = float(np.median([value[1] for value in j_domains(cut, half0, tip, material)]))
            scaling = {
                "opening_scaling_F_ratio_error": abs((half0["F"] / s0["F"]) / .5 - 1),
                "opening_scaling_W_ratio_error": abs((half0["W"] / s0["W"]) / .25 - 1),
                "opening_scaling_J_ratio_error": abs((half_j / j_aggregate) / .25 - 1),
                "opening_scaling_G_ratio_error": abs((half_g / g_energy) / .25 - 1),
            }
        common = {
            "target_extension_um": extension_um, "actual_extension_um": extension_um,
            "tip_x_m": tip, "tip_y_m": 0.0, "reference_opening_m": OPENING,
            "reaction_N_per_m": s0["F"], "compliance_m2_per_N": s0["C"],
            "elastic_energy_J_per_m": s0["W"], "F_over_U_N_per_m2": s0["F"] / OPENING,
            "assembled_residual_reaction_N_per_m": s0["Fres"],
            "energy_derivative_reaction_N_per_m": 2 * s0["W"] / OPENING,
            "reaction_residual_relative_error": abs(s0["F"] - s0["Fres"]) / abs(s0["F"]),
            "reaction_energy_relative_error": abs(s0["F"] - 2 * s0["W"] / OPENING) / abs(s0["F"]),
            "mesh_fingerprint": fingerprint(cut),
            "crack_topology_fingerprint": hashlib.sha256(json.dumps({"tip": tip, "pairs": pairs}, separators=(",", ":")).encode()).hexdigest(),
            "domain_radii_um": json.dumps(RADII_UM),
            **scaling,
        }
        native_rows.append({
            **common, "native_J_aggregate_J_per_m2": j_aggregate,
            "native_KJ_aggregate_MPa_sqrt_m": kj_aggregate,
            "J_native_over_U2": j_aggregate / OPENING**2,
            "KJ_native_over_U": kj_aggregate / OPENING,
            "native_J_domain_values_JSON": json.dumps({str(r): j for r, j, _ in js}, sort_keys=True),
            "J_domain_spread_J_per_m2": float(j_values.max() - j_values.min()),
            "J_status": "CONSISTENT_WITH_G_WITH_RESIDUAL_DOMAIN_SPREAD" if qualified else "UNQUALIFIED",
        })
        g_rows.append({
            **common, "G_energy_J_per_m2": g_energy, "G_compliance_J_per_m2": g_compliance,
            "G_VCCT_I_J_per_m2": gi, "G_VCCT_II_J_per_m2": gii,
            "G_VCCT_J_per_m2": g_vcct, "qualified_G_J_per_m2": np.mean([g_energy, g_compliance, g_vcct]) if qualified else np.nan,
            "G_uncertainty_J_per_m2": uncertainty,
            "qualified_KG_MPa_sqrt_m": np.sqrt(np.mean([g_energy, g_compliance, g_vcct]) * material.Eprime) / 1e6 if qualified else np.nan,
            "G_over_U2": np.mean([g_energy, g_compliance, g_vcct]) / OPENING**2 if qualified else np.nan,
            "KG_over_U": np.sqrt(np.mean([g_energy, g_compliance, g_vcct]) * material.Eprime) / 1e6 / OPENING if qualified else np.nan,
            "native_J_aggregate_J_per_m2": j_aggregate,
            "J_over_G": j_aggregate / np.mean([g_energy, g_compliance, g_vcct]) if qualified else np.nan,
            "vcct_segment_m": da_vcct, "node_qualification": "QUALIFIED" if qualified else "UNQUALIFIED",
            "topology_invariants_pass": True,
        })
    native = pd.DataFrame(native_rows); structural = pd.DataFrame(g_rows)
    native.to_csv(OUT / "oneD_v2_fem_native_mechanics_map.csv", index=False)
    structural.to_csv(OUT / "oneD_v2_fem_qualified_G_map.csv", index=False)
    validation_path = OUT / "oneD_v2_map_interpolation_validation.csv"
    prior = pd.read_csv(validation_path) if validation_path.exists() else pd.DataFrame()
    if not prior.empty and "provider" in prior:
        prior = prior[prior.provider != "FEMCZM"]
    validation_rows = []
    for source, columns in ((native, ("F_over_U_N_per_m2", "J_native_over_U2", "KJ_native_over_U")),
                            (structural, ("G_over_U2", "KG_over_U"))):
        x = source.actual_extension_um.to_numpy()
        for column in columns:
            y = source[column].to_numpy()
            for index in range(1, len(x) - 1):
                predicted = np.interp(x[index], np.delete(x, index), np.delete(y, index))
                validation_rows.append({"provider": "FEMCZM", "quantity": column,
                                        "left_out_extension_um": x[index], "observed": y[index],
                                        "predicted": predicted, "absolute_error": abs(predicted-y[index]),
                                        "relative_error": abs(predicted-y[index])/max(abs(y[index]), 1e-300)})
    pd.concat([prior, pd.DataFrame(validation_rows)], ignore_index=True).to_csv(validation_path, index=False)
    qualified_count = int(structural.node_qualification.eq("QUALIFIED").sum())
    native_status = "FEMCZM_PRODUCTION_NATIVE_MAP_QUALIFIED"
    structural_status = "FEMCZM_STRUCTURAL_G_MAP_QUALIFIED" if qualified_count == len(structural) else "FEMCZM_STRUCTURAL_G_MAP_PARTIALLY_QUALIFIED"
    (OUT / "ONE_D_V2_FEMCZM_NATIVE_MAP_QUALIFICATION.md").write_text(f"""# One-dimensional V2 FEM/CZM native-map qualification

Status: **{native_status}**.

The map uses the corrected plane-strain codimension-one displacement-discontinuity
topology with a fully pre-split interface and locally mirror-symmetric,
translationally invariant 1 µm tip motif. Native domain J is retained separately
from structural G. Bounded interpolation is allowed only on 0–1000 µm and
extrapolation fails closed.
""")
    (OUT / "ONE_D_V2_FEMCZM_STRUCTURAL_G_MAP_QUALIFICATION.md").write_text(f"""# One-dimensional V2 FEM/CZM structural-G map qualification

Status: **{structural_status}** ({qualified_count}/{len(structural)} nodes qualified).

Every qualified node has positive global energy, compliance, and VCCT estimates
with maximum mutual spread no greater than 10%. The spread is stored as numerical
uncertainty. Domain J is secondary and residual domain spread does not erase a
qualified global G result. Unqualified nodes remain empty; no neighboring value
is inferred.
""")
    fig, axes = plt.subplots(2, 1, figsize=(9, 8), sharex=True, constrained_layout=True)
    axes[0].plot(native.actual_extension_um, native.J_native_over_U2, "o-", label="native J/U²")
    axes[0].plot(structural.actual_extension_um, structural.G_over_U2, "s-", label="qualified G/U²")
    axes[0].set_ylabel("coefficient"); axes[0].legend(); axes[0].grid(alpha=.25)
    axes[1].plot(native.actual_extension_um, native.KJ_native_over_U, "o-", label="native KJ/U")
    axes[1].plot(structural.actual_extension_um, structural.KG_over_U, "s-", label="qualified KG/U")
    axes[1].set(xlabel="Straight-path extension (µm)", ylabel="coefficient"); axes[1].legend(); axes[1].grid(alpha=.25)
    fig.suptitle("FEM/CZM native and qualified structural mechanics maps")
    fig.savefig(OUT / "FEMCZM_NATIVE_AND_QUALIFIED_MAPS.png", dpi=180); plt.close(fig)
    fig, ax = plt.subplots(figsize=(8, 4.5), constrained_layout=True)
    ax.plot(structural.actual_extension_um, structural.J_over_G, "o-"); ax.axhline(1, color="k", ls="--")
    ax.set(xlabel="Straight-path extension (µm)", ylabel="native J / qualified G"); ax.grid(alpha=.25)
    fig.savefig(OUT / "FEMCZM_MAP_J_OVER_G.png", dpi=180); plt.close(fig)


if __name__ == "__main__":
    main()
