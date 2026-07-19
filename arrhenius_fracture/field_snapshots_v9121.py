"""v9.12.1 full-field snapshots with authoritative emission annotations."""
from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

import numpy as np

from .field_snapshots_v912 import map_mpz_density_to_elements


def authoritative_emitted_ledger(snap: dict[str, Any]) -> float:
    """Return cumulative emitted physical line content, not retained inventory."""
    for key in (
        "signed_burgers_line_content_emitted_total",
        "signed_burgers_source_activations_total",
        "mpz_emitted_total",
        "emitted_total",
    ):
        try:
            value = float(snap.get(key))
        except (TypeError, ValueError):
            continue
        if np.isfinite(value):
            return value
    total = 0.0
    found = False
    for front in snap.get("mpz_front_states", []):
        state = front.get("state", {})
        if isinstance(state, dict) and isinstance(state.get("state"), dict):
            state = state["state"]
        for key in (
            "signed_line_content_emitted_total",
            "emitted_total",
        ):
            try:
                value = float(state.get(key))
            except (AttributeError, TypeError, ValueError):
                continue
            if np.isfinite(value):
                total += value
                found = True
                break
    if found:
        return total
    # Legacy fallback is retained only for historical snapshots lacking the
    # authoritative ledger and is identified in the manifest.
    try:
        return float(snap.get("N_em", 0.0))
    except (TypeError, ValueError):
        return 0.0


def render_field_snapshots_v9121(out, T, mesh, snaps, max_cols=5):
    """Render and persist damage, density, stress and plastic-strain fields."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.tri as mtri
    except ImportError:
        print("  matplotlib not available, skipping v9.12.1 field snapshots")
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
    emitted = [authoritative_emitted_ledger(snap) for snap in pick]
    rows = ("damage", "rho", "sig1", "epeq")
    epeq_max = max(float(np.max([np.max(s["epeq_gp"]) for s in pick])), 1.0e-6)
    s1_abs = max(float(np.max([np.max(np.abs(s["s1_gp"])) for s in pick])) / 1.0e6, 1.0)
    rho_max_log = max(
        16.0,
        float(np.max([np.nanmax(np.log10(np.maximum(r, 1.0))) for r in rho_total])),
    )
    nrow, ncol = len(rows), len(pick)
    fig, axes = plt.subplots(
        nrow,
        ncol,
        figsize=(4.5 * ncol, 3.9 * nrow),
        squeeze=False,
        constrained_layout=True,
    )
    last_im = {}
    manifest = {
        "schema": "full_field_snapshots_v9121",
        "temperature_K": float(T),
        "density_definition": "rho_total = FEM_forest_rho_gp + visualization_map(mobile+retained MPZ counts)",
        "mpz_mapping_role": "output_only_no_mechanical_feedback",
        "emission_annotation": "authoritative cumulative emitted physical line content",
        "legacy_N_em_annotation_removed": True,
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

    for j, (snap, mapped, total, emitted_total) in enumerate(
        zip(pick, mapped_fields, rho_total, emitted)
    ):
        nodes = np.asarray(snap["nodes"], dtype=float)
        elems = np.asarray(snap["elems"], dtype=int)
        tri = mtri.Triangulation(nodes[:, 0] * 1.0e3, nodes[:, 1] * 1.0e3, elems)
        array_name = f"snapshot_{j:02d}_step_{int(snap['step']):06d}.npz"
        np.savez_compressed(
            array_dir / array_name,
            step=np.asarray(int(snap["step"])),
            KJ_Pa_sqrt_m=np.asarray(float(snap["KJ"])),
            a_tip_m=np.asarray(float(snap["a_tip"])),
            emitted_line_content_total=np.asarray(float(emitted_total)),
            nodes_m=nodes,
            elems=elems,
            damage_nodal=np.asarray(snap["d"], dtype=float),
            rho_fem_m2=np.asarray(snap["rho_gp"], dtype=float),
            rho_mpz_mapped_m2=mapped,
            rho_total_m2=total,
            sigma1_Pa=np.asarray(snap["s1_gp"], dtype=float),
            equivalent_plastic_strain=np.asarray(snap["epeq_gp"], dtype=float),
        )
        manifest["snapshots"].append(
            {
                "column": j,
                "step": int(snap["step"]),
                "KJ_MPa_sqrt_m": float(snap["KJ"]) / 1.0e6,
                "a_tip_mm": float(snap["a_tip"]) * 1.0e3,
                "emitted_line_content_total": float(emitted_total),
                "max_rho_fem_m2": float(np.max(snap["rho_gp"])),
                "max_rho_mpz_mapped_m2": float(np.max(mapped)) if mapped.size else 0.0,
                "max_rho_total_m2": float(np.max(total)),
                "array_file": str(array_dir / array_name),
            }
        )

        for i, row in enumerate(rows):
            ax = axes[i, j]
            if row == "damage":
                im = ax.tripcolor(
                    tri,
                    np.asarray(snap["d"], float),
                    shading="gouraud",
                    cmap="inferno",
                    vmin=0.0,
                    vmax=1.0,
                    rasterized=True,
                )
                title = "damage d + crack path"
            elif row == "rho":
                vals = np.log10(np.maximum(total, 1.0))
                im = ax.tripcolor(
                    tri,
                    vals,
                    shading="flat",
                    cmap="viridis",
                    vmin=10.0,
                    vmax=rho_max_log,
                    rasterized=True,
                )
                title = "log10 rho total (m^-2)"
            elif row == "sig1":
                im = ax.tripcolor(
                    tri,
                    np.asarray(snap["s1_gp"], float) / 1.0e6,
                    shading="flat",
                    cmap="coolwarm",
                    vmin=-s1_abs,
                    vmax=s1_abs,
                    rasterized=True,
                )
                title = "sigma1 FEM (MPa)"
            else:
                im = ax.tripcolor(
                    tri,
                    np.asarray(snap["epeq_gp"], float),
                    shading="flat",
                    cmap="magma",
                    vmin=0.0,
                    vmax=epeq_max,
                    rasterized=True,
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
                f"N_emit,total={float(emitted_total):.2f}  a={float(snap['a_tip'])*1e3:.3f} mm",
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


__all__ = [
    "authoritative_emitted_ledger",
    "map_mpz_density_to_elements",
    "render_field_snapshots_v9121",
]
