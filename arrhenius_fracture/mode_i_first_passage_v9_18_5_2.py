"""v9.18.5.2 startup-safe compact Mode-I corridor mesh.

The v9.18.5 multi-center corridor is formed by merging several radial point
clouds before Delaunay triangulation.  Qhull is permitted to omit duplicate or
numerically redundant input points from the returned simplices.  Retaining
those unused points in ``TriMesh.nodes`` creates zero stiffness rows and
one-node connected components, so the first mechanics solve can fail before a
physical event is initialized.

This wrapper compacts the generated mesh to the nodes actually referenced by
bulk triangles, remaps the connectivity, rebuilds all geometric operators, and
validates the initial mesh before the first FEM assembly.  No constitutive,
hazard, cohesive, MPZ, wake, or shielding law is changed.
"""
from __future__ import annotations

import json
import math
import os
from pathlib import Path
import sys
from typing import Any

import numpy as np

from . import mesh as _mesh
from . import mode_i_first_passage_v9_18_5 as _v9185
from . import mode_i_first_passage_v9_18_5_1 as _v91851


_STARTUP_AUDIT: dict[str, Any] = {}


def _triangle_quality(nodes: np.ndarray, elems: np.ndarray) -> np.ndarray:
    """Return the standard equilateral-normalized triangle quality in [0, 1]."""
    tri = np.asarray(nodes, float)[np.asarray(elems, int)]
    e01 = np.sum((tri[:, 1] - tri[:, 0]) ** 2, axis=1)
    e12 = np.sum((tri[:, 2] - tri[:, 1]) ** 2, axis=1)
    e20 = np.sum((tri[:, 0] - tri[:, 2]) ** 2, axis=1)
    cross = ((tri[:, 1, 0] - tri[:, 0, 0]) *
             (tri[:, 2, 1] - tri[:, 0, 1]) -
             (tri[:, 1, 1] - tri[:, 0, 1]) *
             (tri[:, 2, 0] - tri[:, 0, 0]))
    area = 0.5 * np.abs(cross)
    denom = np.maximum(e01 + e12 + e20, 1.0e-300)
    return 4.0 * math.sqrt(3.0) * area / denom


def _compact_mesh(mesh: Any, tip_centers: np.ndarray) -> tuple[Any, dict[str, Any]]:
    elems = np.asarray(mesh.elems, dtype=int)
    used = np.unique(elems.ravel())
    unused = np.setdiff1d(np.arange(int(mesh.nn), dtype=int), used, assume_unique=True)

    if unused.size:
        remap = np.full(int(mesh.nn), -1, dtype=int)
        remap[used] = np.arange(len(used), dtype=int)
        nodes = np.asarray(mesh.nodes, float)[used]
        elems_new = remap[elems]
        compact = _mesh.rebuild_tri_mesh(
            nodes,
            elems_new,
            tip_centers=np.asarray(tip_centers, float),
            validate=True,
        )
    else:
        compact = mesh

    incidence = np.bincount(
        np.asarray(compact.elems, int).ravel(), minlength=int(compact.nn)
    )
    orphan = np.where(incidence <= 0)[0]
    quality = _triangle_quality(compact.nodes, compact.elems)
    qmin = float(np.min(quality)) if quality.size else 0.0
    qfloor = float(os.environ.get(
        "ARRHENIUS_MIN_INITIAL_TRIANGLE_QUALITY",
        os.environ.get("ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY", "0.035"),
    ))

    audit = {
        "input_node_count": int(mesh.nn),
        "triangle_count": int(mesh.ne),
        "unused_input_node_count": int(unused.size),
        "unused_input_nodes_first20": unused[:20].astype(int).tolist(),
        "compacted_node_count": int(compact.nn),
        "orphan_node_count_after_compaction": int(orphan.size),
        "orphan_nodes_after_compaction_first20": orphan[:20].astype(int).tolist(),
        "minimum_initial_triangle_quality": qmin,
        "minimum_initial_triangle_quality_required": qfloor,
        "initial_hbar_tip_m": float(compact.hbar_tip),
        "initial_hbar_global_m": float(compact.hbar),
        "initial_mesh_compaction_applied": bool(unused.size),
    }

    if orphan.size:
        raise RuntimeError(
            "v9.18.5.2 corridor compaction left orphan bulk nodes: "
            f"count={len(orphan)} first={orphan[:20].tolist()}"
        )
    if not np.all(np.isfinite(compact.area_e)) or np.any(compact.area_e <= 0.0):
        raise RuntimeError("v9.18.5.2 initial corridor contains nonpositive areas")
    if qmin < qfloor:
        raise RuntimeError(
            "v9.18.5.2 initial corridor triangle quality below production floor: "
            f"qmin={qmin:.6e} floor={qfloor:.6e}"
        )
    return compact, audit


def _compact_corridor_mesh(geom, mesh_cfg, seed=None, tip_center=None):
    original = _compact_corridor_mesh._original
    graded = float(getattr(mesh_cfg, "tip_h_fine", 0.0) or 0.0) > 0.0
    enabled = os.environ.get("ARRHENIUS_PREFINED_MODE_I_CORRIDOR", "1") not in {
        "0", "false", "False", "no", "NO"
    }

    if tip_center is not None or not graded or not enabled:
        centers = np.asarray(
            [[float(geom.a0), 0.0]] if tip_center is None else tip_center,
            dtype=float,
        )
        centers = centers.reshape(1, 2) if centers.ndim == 1 else centers[:, :2]
        raw = original(geom, mesh_cfg, seed=seed, tip_center=tip_center)
    else:
        centers = _v9185._corridor_centers(geom, mesh_cfg)
        _v9185._RUNTIME["corridor_centers"] = centers.tolist()
        raw = original(geom, mesh_cfg, seed=seed, tip_center=centers)

    compact, audit = _compact_mesh(raw, centers)
    audit.update({
        "schema": "compact_corridor_mesh_v91852_v1",
        "corridor_center_count": int(len(centers)),
        "corridor_centers_m": np.asarray(centers, float).tolist(),
        "constitutive_physics_changed": False,
    })
    _STARTUP_AUDIT.clear()
    _STARTUP_AUDIT.update(audit)
    _v9185._RUNTIME["mesh"] = compact
    return compact


def _write_startup_audit(argv: list[str], error: BaseException | None = None) -> None:
    out_value = _v9185._option_value(argv, "--out")
    if out_value is None:
        return
    out = Path(out_value)
    out.mkdir(parents=True, exist_ok=True)
    payload = dict(_STARTUP_AUDIT)
    payload.setdefault("schema", "compact_corridor_mesh_v91852_v1")
    payload["startup_completed"] = error is None
    payload["startup_error_type"] = None if error is None else type(error).__name__
    payload["startup_error"] = None if error is None else str(error)
    (out / "compact_corridor_mesh_v91852.json").write_text(
        json.dumps(payload, indent=2, default=str)
    )


def main(argv=None):
    user_args = list(sys.argv[1:] if argv is None else argv)
    original = _v9185._make_corridor_mesh
    _v9185._make_corridor_mesh = _compact_corridor_mesh
    error: BaseException | None = None
    try:
        results = _v91851.main(user_args)
    except BaseException as exc:
        error = exc
        print(
            "V9_18_5_2_STARTUP_FAILURE "
            f"{type(exc).__name__}: {exc}",
            file=sys.stderr,
            flush=True,
        )
        raise
    finally:
        _write_startup_audit(user_args, error)
        _v9185._make_corridor_mesh = original
    return results


if __name__ == "__main__":
    main()
