"""Full 2-D field snapshots with mapped front-local MPZ density for v9.12."""
from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

import numpy as np


def _state_payload(front: dict[str, Any]) -> dict[str, Any]:
    state = front.get("state", {})
    if isinstance(state, dict) and "state" in state and isinstance(state["state"], dict):
        state = state["state"]
    return state if isinstance(state, dict) else {}


def map_mpz_density_to_elements(snap: dict[str, Any]) -> np.ndarray:
    """Map front-local mobile+retained line counts to element density [m^-2].

    The MPZ is a one-dimensional moving inventory with an effective process-zone
    width.  Each bin count is converted to density through ``dx * width`` and
    projected into a Gaussian strip oriented with the current crack front.  The
    mapping is for visualization/output only; it does not feed back into the
    mechanics or shielding calculation.
    """
    nodes = np.asarray(snap["nodes"], dtype=float)
    elems = np.asarray(snap["elems"], dtype=int)
    centroids = nodes[elems].mean(axis=1)
    mapped = np.zeros(len(elems), dtype=float)
    for front in snap.get("mpz_front_states", []):
        state = _state_payload(front)
        cfg = dict(state.get("config", {}))
        retained = np.asarray(state.get("retained", []), dtype=float)
        mobile = np.asarray(state.get("mobile", []), dtype=float)
        if retained.ndim != 2 or retained.size == 0:
            continue
        if mobile.shape != retained.shape:
            mobile = np.zeros_like(retained)
        n_bins = retained.shape[1]
        length = max(float(cfg.get("length_m", 0.0)), 1.0e-12)
        dx = length / max(n_bins, 1)
        width = max(float(cfg.get("blunting_length_m", dx)), dx, 1.0e-12)
        count_bin = np.sum(np.maximum(retained, 0.0) + np.maximum(mobile, 0.0), axis=0)
        rho_bin = count_bin / max(dx * width, 1.0e-30)
        xbin = (np.arange(n_bins, dtype=float) + 0.5) * dx

        tip = np.asarray(front.get("xy_m", [0.0, 0.0]), dtype=float).reshape(2)
        direction = np.asarray(front.get("direction", [1.0, 0.0]), dtype=float).reshape(2)
        norm = float(np.linalg.norm(direction))
        if norm <= 1.0e-30:
            direction = np.array([1.0, 0.0])
        else:
            direction = direction / norm
        normal = np.array([-direction[1], direction[0]])
        rel = centroids - tip[None, :]
        xi = rel @ direction
        eta = rel @ normal
        active = (xi >= 0.0) & (xi <= length) & (np.abs(eta) <= 4.0 * width)
        if not np.any(active):
            continue
        longitudinal = np.interp(
            np.clip(xi[active], xbin[0], xbin[-1]), xbin, rho_bin,
            left=float(rho_bin[0]), right=float(rho_bin[-1]),
        )
        transverse = np.exp(-0.5 * (eta[active] / width) ** 2)
        mapped[active] += longitudinal * transverse
    return mapped


def render_field_snapshots_v912(out, T, mesh, snaps, max_cols=5):
    """Render and persist damage, density, stress and plastic-strain fields."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.tri as mtri
    except ImportError:
        print("  matplotlib not available, skipping v9.12 field snapshots")
        return None
    if not snaps:
        return None

    if len(snaps) > max_cols:
        idx = np.linspace(0, len(snaps) - 1, max_cols, dtype=int)
        pick = [snaps[i] for i in idx]
    else:
        pick = list(snaps)

    mapped_fields = [map_mpz_density_to_elements(snap) for snap in pick]
    rho_total = [
        np.maximum(np.asarray(snap["rho_gp"], float), 0.0) + mapped
        for snap, mapped in zip(pick, mapped_fields)
    ]
    rows = ("damage", "rho", "sig1", "epeq")
    epeq_max = max(float(np.max([np.max(s["epeq_gp"]) for s in pick])), 1.0e-6)
    s1_abs = max(float(np.max([np.max(np.abs(s["s1_gp"])) for s in pick])) / 1.0e6, 1.0)
    rho_max_log = max(16.0, float(np.max([np.nanmax(np.log10(np.maximum(r, 1.0))) for r in rho_total])))
    nrow, ncol = len(rows), len(pick)
    fig, axes = plt.subplots(
        nrow, ncol, figsize=(4.5 * ncol, 3.9 * nrow),
        squeeze=False, constrained_layout=True,
    )
    last_im = {}
    manifest = {
        "schema": "full_field_snapshots_v912",
        "temperature_K": float(T),
        "density_definition": "rho_total = FEM_forest_rho_gp + visualization_map(mobile+retained MPZ counts)",
        "mpz_mapping_role": "output_only_no_mechanical_feedback",
        "rows": [
            "damage_and_explicit_crack_path",
            "log10_total_dislocation_density_m-2",
            "maximum_principal_FEM_stress_MPa",
            "equivalent_plastic_strain",
        ],
        "snapshots": [],
    }
    array_dir = Path(out) / f"field_snapshot_arrays_{int(T)}K"
    array_dir.mkdir(parents=True, exist_ok=True)

    for j, (snap, mapped, total) in enumerate(zip(pick, mapped_fields, rho_total)):
        nodes = np.asarray(snap["nodes"], dtype=float)
        elems = np.asarray(snap["elems"], dtype=int)
        tri = mtri.Triangulation(nodes[:, 0] * 1.0e3, nodes[:, 1] * 1.0e3, elems)
        array_name = f"snapshot_{j:02d}_step_{int(snap['step']):06d}.npz"
        np.savez_compressed(
            array_dir / array_name,
            step=np.asarray(int(snap["step"])),
            KJ_Pa_sqrt_m=np.asarray(float(snap["KJ"])),
            a_tip_m=np.asarray(float(snap["a_tip"])),
            nodes_m=nodes,
            elems=elems,
            damage_nodal=np.asarray(snap["d"], dtype=float),
            rho_fem_m2=np.asarray(snap["rho_gp"], dtype=float),
            rho_mpz_mapped_m2=mapped,
            rho_total_m2=total,
            sigma1_Pa=np.asarray(snap["s1_gp"], dtype=float),
            equivalent_plastic_strain=np.asarray(snap["epeq_gp"], dtype=float),
        )
        manifest["snapshots"].append({
            "column": j,
            "step": int(snap["step"]),
            "KJ_MPa_sqrt_m": float(snap["KJ"]) / 1.0e6,
            "a_tip_mm": float(snap["a_tip"]) * 1.0e3,
            "N_em": float(snap.get("N_em", 0.0)),
            "max_rho_fem_m2": float(np.max(snap["rho_gp"])),
            "max_rho_mpz_mapped_m2": float(np.max(mapped)) if mapped.size else 0.0,
            "max_rho_total_m2": float(np.max(total)),
            "array_file": str(array_dir / array_name),
        })

        for i, row in enumerate(rows):
            ax = axes[i, j]
            if row == "damage":
                im = ax.tripcolor(
                    tri, np.asarray(snap["d"], float), shading="gouraud",
                    cmap="inferno", vmin=0.0, vmax=1.0, rasterized=True,
                )
                title = "damage d + crack path"
            elif row == "rho":
                vals = np.log10(np.maximum(total, 1.0))
                im = ax.tripcolor(
                    tri, vals, shading="flat", cmap="viridis",
                    vmin=10.0, vmax=rho_max_log, rasterized=True,
                )
                title = "log10 rho total (m^-2)"
            elif row == "sig1":
                im = ax.tripcolor(
                    tri, np.asarray(snap["s1_gp"], float) / 1.0e6,
                    shading="flat", cmap="coolwarm",
                    vmin=-s1_abs, vmax=s1_abs, rasterized=True,
                )
                title = "sigma1 FEM (MPa)"
            else:
                im = ax.tripcolor(
                    tri, np.asarray(snap["epeq_gp"], float), shading="flat",
                    cmap="magma", vmin=0.0, vmax=epeq_max, rasterized=True,
                )
                title = "equivalent plastic strain"
            last_im[row] = im
            for item in snap.get("front_paths", []):
                try:
                    _, _, path = item
                    path = np.asarray(path, float)
                    if path.ndim == 2 and path.shape[0] >= 2:
                        ax.plot(path[:, 0] * 1.0e3, path[:, 1] * 1.0e3, "w-", lw=1.35)
                        ax.plot(path[-1, 0] * 1.0e3, path[-1, 1] * 1.0e3, "wo", ms=2.8)
                except Exception:
                    pass
            ax.set_aspect("equal")
            ax.set_title(
                f"{title}\nstep {int(snap['step'])}  KJ={float(snap['KJ'])/1e6:.2f}  "
                f"N_em={float(snap.get('N_em', 0.0)):.2f}  a={float(snap['a_tip'])*1e3:.3f} mm",
                fontsize=8,
            )
            ax.set_xlabel("x (mm)")
            if j == 0:
                ax.set_ylabel("y (mm)")
    for i, row in enumerate(rows):
        fig.colorbar(last_im[row], ax=axes[i, :], shrink=0.78, pad=0.01)
    fig.suptitle(f"Full 2-D FEM/CZM field snapshots — T = {float(T):.0f} K", fontsize=13)
    path = os.path.join(out, f"field_snapshots_{int(T)}K.png")
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    manifest["image_file"] = path
    (Path(out) / f"field_snapshot_manifest_{int(T)}K.json").write_text(
        json.dumps(manifest, indent=2)
    )
    return path


__all__ = ["map_mpz_density_to_elements", "render_field_snapshots_v912"]
