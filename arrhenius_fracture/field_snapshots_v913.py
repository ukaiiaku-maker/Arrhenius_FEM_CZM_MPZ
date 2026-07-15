"""Full 2-D and crack-tip zoom snapshots for v9.13.

The front-local MPZ inventory is coarse-grained onto the resolved FE mesh for
visualization only.  The Gaussian display kernel preserves line count per unit
crack-front length and is broadened to at least the local element resolution.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

import numpy as np


def _state_payload(front: dict[str, Any]) -> dict[str, Any]:
    state = front.get("state", {})
    if isinstance(state, dict) and isinstance(state.get("state"), dict):
        state = state["state"]
    return state if isinstance(state, dict) else {}


def _triangle_area(nodes: np.ndarray, elems: np.ndarray) -> np.ndarray:
    p = nodes[elems]
    return 0.5 * np.abs(
        (p[:, 1, 0] - p[:, 0, 0]) * (p[:, 2, 1] - p[:, 0, 1])
        - (p[:, 2, 0] - p[:, 0, 0]) * (p[:, 1, 1] - p[:, 0, 1])
    )


def _front_scalars(snap: dict[str, Any]) -> dict[str, float]:
    emitted = retained = mobile = 0.0
    for front in snap.get("mpz_front_states", []):
        state = _state_payload(front)
        try:
            emitted += float(state.get("emitted_total", 0.0) or 0.0)
        except (TypeError, ValueError):
            pass
        for key, target in (("retained", "retained"), ("mobile", "mobile")):
            arr = np.asarray(state.get(key, []), dtype=float)
            if arr.size:
                if target == "retained":
                    retained += float(np.sum(np.maximum(arr, 0.0)))
                else:
                    mobile += float(np.sum(np.maximum(arr, 0.0)))
    return {"emitted_total": emitted, "retained_count": retained, "mobile_count": mobile}


def map_mpz_density_to_elements(snap: dict[str, Any]) -> tuple[np.ndarray, list[dict[str, float]]]:
    """Coarse-grain mobile+retained line counts to FE element density [m^-2]."""
    nodes = np.asarray(snap["nodes"], dtype=float)
    elems = np.asarray(snap["elems"], dtype=int)
    centroids = nodes[elems].mean(axis=1)
    h_elem = np.sqrt(np.maximum(_triangle_area(nodes, elems), 1.0e-30))
    mapped = np.zeros(len(elems), dtype=float)
    metadata: list[dict[str, float]] = []

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
        physical_width = max(float(cfg.get("blunting_length_m", dx)), dx, 1.0e-12)
        count_bin = np.sum(np.maximum(retained, 0.0) + np.maximum(mobile, 0.0), axis=0)
        xbin = (np.arange(n_bins, dtype=float) + 0.5) * dx

        tip = np.asarray(front.get("xy_m", [0.0, 0.0]), dtype=float).reshape(2)
        direction = np.asarray(front.get("direction", [1.0, 0.0]), dtype=float).reshape(2)
        norm = float(np.linalg.norm(direction))
        direction = np.array([1.0, 0.0]) if norm <= 1.0e-30 else direction / norm
        normal = np.array([-direction[1], direction[0]])
        rel = centroids - tip[None, :]
        xi = rel @ direction
        eta = rel @ normal
        near = (xi >= -0.1 * length) & (xi <= 1.1 * length)
        local_h = float(np.median(h_elem[near])) if np.any(near) else float(np.median(h_elem))
        display_width = max(physical_width, 1.5 * local_h)
        active = (xi >= 0.0) & (xi <= length) & (np.abs(eta) <= 4.0 * display_width)
        if np.any(active):
            count_per_bin = np.interp(
                np.clip(xi[active], xbin[0], xbin[-1]), xbin, count_bin,
                left=float(count_bin[0]), right=float(count_bin[-1]),
            )
            # Gaussian normalized across eta: integral rho d(eta) = count/dx.
            rho = count_per_bin / max(dx * np.sqrt(2.0 * np.pi) * display_width, 1.0e-30)
            rho *= np.exp(-0.5 * (eta[active] / display_width) ** 2)
            mapped[active] += rho
        metadata.append({
            "length_m": length, "dx_m": dx,
            "physical_width_m": physical_width,
            "display_coarse_grain_width_m": display_width,
            "local_element_scale_m": local_h,
            "active_line_count": float(np.sum(count_bin)),
            "emitted_total": float(state.get("emitted_total", 0.0) or 0.0),
        })
    return mapped, metadata


def _plot_panel(out: str | Path, T: float, pick: list[dict[str, Any]], mapped_fields: list[np.ndarray],
                total_fields: list[np.ndarray], zoom: bool) -> str:
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri

    rows = ("damage", "rho", "sig1", "epeq")
    nrow, ncol = len(rows), len(pick)
    masks: list[np.ndarray] = []
    for snap in pick:
        nodes = np.asarray(snap["nodes"], float)
        elems = np.asarray(snap["elems"], int)
        cent = nodes[elems].mean(axis=1)
        if zoom and snap.get("mpz_front_states"):
            front = snap["mpz_front_states"][0]
            tip = np.asarray(front.get("xy_m", [snap.get("a_tip", 0.0), 0.0]), float)
            state = _state_payload(front)
            cfg = dict(state.get("config", {}))
            L = max(float(cfg.get("length_m", 1.0e-4)), 1.0e-6)
            mask = (
                (cent[:, 0] >= tip[0] - 0.15 * L) & (cent[:, 0] <= tip[0] + 1.05 * L)
                & (cent[:, 1] >= tip[1] - 0.20 * L) & (cent[:, 1] <= tip[1] + 0.20 * L)
            )
        else:
            mask = np.ones(len(elems), dtype=bool)
        masks.append(mask)

    epeq_max = max(
        max(float(np.max(np.asarray(s["epeq_gp"])[m])) if np.any(m) else 0.0 for s, m in zip(pick, masks)),
        1.0e-8,
    )
    stress_abs = max(
        max(float(np.max(np.abs(np.asarray(s["s1_gp"])[m]))) if np.any(m) else 0.0 for s, m in zip(pick, masks)) / 1.0e6,
        1.0,
    )
    rho_logs = [np.log10(np.maximum(r[m], 1.0)) for r, m in zip(total_fields, masks) if np.any(m)]
    rho_min = min(10.0, min(float(np.min(q)) for q in rho_logs)) if rho_logs else 10.0
    rho_max = max(14.0, max(float(np.max(q)) for q in rho_logs)) if rho_logs else 16.0

    fig, axes = plt.subplots(nrow, ncol, figsize=(4.6 * ncol, 3.9 * nrow), squeeze=False, constrained_layout=True)
    last_im: dict[str, Any] = {}
    for j, (snap, total) in enumerate(zip(pick, total_fields)):
        nodes = np.asarray(snap["nodes"], float)
        elems = np.asarray(snap["elems"], int)
        tri = mtri.Triangulation(nodes[:, 0] * 1.0e3, nodes[:, 1] * 1.0e3, elems)
        scalars = _front_scalars(snap)
        zoom_limits = None
        if zoom and snap.get("mpz_front_states"):
            front = snap["mpz_front_states"][0]
            tip = np.asarray(front.get("xy_m", [snap.get("a_tip", 0.0), 0.0]), float)
            cfg = dict(_state_payload(front).get("config", {}))
            L = max(float(cfg.get("length_m", 1.0e-4)), 1.0e-6)
            zoom_limits = (
                (tip[0] - 0.15 * L) * 1.0e3, (tip[0] + 1.05 * L) * 1.0e3,
                (tip[1] - 0.20 * L) * 1.0e3, (tip[1] + 0.20 * L) * 1.0e3,
            )
        for i, row in enumerate(rows):
            ax = axes[i, j]
            if row == "damage":
                im = ax.tripcolor(tri, np.asarray(snap["d"], float), shading="gouraud", cmap="inferno", vmin=0.0, vmax=1.0, rasterized=True)
                title = "damage d + crack path"
            elif row == "rho":
                im = ax.tripcolor(tri, np.log10(np.maximum(total, 1.0)), shading="flat", cmap="viridis", vmin=rho_min, vmax=rho_max, rasterized=True)
                title = "log10 coarse-grained rho total (m^-2)"
            elif row == "sig1":
                im = ax.tripcolor(tri, np.asarray(snap["s1_gp"], float) / 1.0e6, shading="flat", cmap="coolwarm", vmin=-stress_abs, vmax=stress_abs, rasterized=True)
                title = "maximum principal FEM stress (MPa)"
            else:
                im = ax.tripcolor(tri, np.asarray(snap["epeq_gp"], float), shading="flat", cmap="magma", vmin=0.0, vmax=epeq_max, rasterized=True)
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
            if zoom_limits is not None:
                ax.set_xlim(zoom_limits[0], zoom_limits[1])
                ax.set_ylim(zoom_limits[2], zoom_limits[3])
            ax.set_aspect("equal")
            ax.set_title(
                f"{title}\nstep {int(snap['step'])}  KJ={float(snap['KJ'])/1e6:.2f}\n"
                f"Nemit={scalars['emitted_total']:.2f}  Nret={scalars['retained_count']:.2f}  Nmob={scalars['mobile_count']:.2f}",
                fontsize=8,
            )
            ax.set_xlabel("x (mm)")
            if j == 0:
                ax.set_ylabel("y (mm)")
    for i, row in enumerate(rows):
        fig.colorbar(last_im[row], ax=axes[i, :], shrink=0.78, pad=0.01)
    kind = "Crack-tip zoom" if zoom else "Full 2-D specimen"
    fig.suptitle(f"{kind} FEM/CZM field snapshots — T = {float(T):.0f} K", fontsize=13)
    name = f"field_snapshots_tip_zoom_{int(T)}K.png" if zoom else f"field_snapshots_{int(T)}K.png"
    path = os.path.join(out, name)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return path


def render_field_snapshots_v913(out, T, mesh, snaps, max_cols=5):
    try:
        import matplotlib
        matplotlib.use("Agg")
    except ImportError:
        print("  matplotlib not available, skipping v9.13 field snapshots")
        return None
    if not snaps:
        return None
    if len(snaps) > max_cols:
        idx = np.linspace(0, len(snaps) - 1, max_cols, dtype=int)
        pick = [snaps[i] for i in idx]
    else:
        pick = list(snaps)

    mapped_with_meta = [map_mpz_density_to_elements(s) for s in pick]
    mapped_fields = [x[0] for x in mapped_with_meta]
    map_meta = [x[1] for x in mapped_with_meta]
    total_fields = [np.maximum(np.asarray(s["rho_gp"], float), 0.0) + m for s, m in zip(pick, mapped_fields)]
    array_dir = Path(out) / f"field_snapshot_arrays_{int(T)}K"
    array_dir.mkdir(parents=True, exist_ok=True)
    manifest = {
        "schema": "full_field_snapshots_v913",
        "temperature_K": float(T),
        "density_definition": "rho_total = FEM forest density + line-count-preserving Gaussian coarse-graining of mobile+retained MPZ inventory",
        "mpz_mapping_role": "output_only_no_mechanical_feedback",
        "snapshots": [],
    }
    for j, (snap, mapped, total, meta) in enumerate(zip(pick, mapped_fields, total_fields, map_meta)):
        scalars = _front_scalars(snap)
        array_name = f"snapshot_{j:02d}_step_{int(snap['step']):06d}.npz"
        np.savez_compressed(
            array_dir / array_name,
            step=np.asarray(int(snap["step"])), KJ_Pa_sqrt_m=np.asarray(float(snap["KJ"])),
            a_tip_m=np.asarray(float(snap["a_tip"])), nodes_m=np.asarray(snap["nodes"], float),
            elems=np.asarray(snap["elems"], int), damage_nodal=np.asarray(snap["d"], float),
            rho_fem_m2=np.asarray(snap["rho_gp"], float), rho_mpz_coarse_grained_m2=mapped,
            rho_total_m2=total, sigma1_Pa=np.asarray(snap["s1_gp"], float),
            equivalent_plastic_strain=np.asarray(snap["epeq_gp"], float),
            emitted_total=np.asarray(scalars["emitted_total"]), retained_count=np.asarray(scalars["retained_count"]),
            mobile_count=np.asarray(scalars["mobile_count"]),
        )
        manifest["snapshots"].append({
            "column": j, "step": int(snap["step"]), "KJ_MPa_sqrt_m": float(snap["KJ"]) / 1.0e6,
            "a_tip_mm": float(snap["a_tip"]) * 1.0e3, **scalars,
            "max_rho_fem_m2": float(np.max(snap["rho_gp"])),
            "max_rho_mpz_coarse_grained_m2": float(np.max(mapped)) if mapped.size else 0.0,
            "max_rho_total_m2": float(np.max(total)), "mapping_fronts": meta,
            "array_file": str(array_dir / array_name),
        })
    full_path = _plot_panel(out, T, pick, mapped_fields, total_fields, zoom=False)
    zoom_path = _plot_panel(out, T, pick, mapped_fields, total_fields, zoom=True)
    manifest["full_field_image"] = full_path
    manifest["tip_zoom_image"] = zoom_path
    (Path(out) / f"field_snapshot_manifest_{int(T)}K.json").write_text(json.dumps(manifest, indent=2))
    return full_path


__all__ = ["map_mpz_density_to_elements", "render_field_snapshots_v913"]
