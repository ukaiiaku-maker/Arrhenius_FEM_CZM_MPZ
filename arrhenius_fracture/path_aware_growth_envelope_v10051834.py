"""Path-aware long-growth support envelope for v10.0.5.18.3.4.

This point release changes only the startup mesh policy.  The production crack
selector continuously competes directions and then applies the global ligament
gate ``t . fwd0 >= min_global_forward``.  A horizontal corridor is therefore not
a sufficient long-growth support contract.  This module constructs a conservative
reachable envelope from the same global-forward gate, resolves that complete
envelope at the process-zone scale, and retains the fine initial-tip grading.

No constitutive, hazard, stochastic event-length, source, shielding, cohesive,
moving-tip, or adaptive-CZM quality law is changed.
"""
from __future__ import annotations

import math
import os
from typing import Any

import numpy as np
from scipy.spatial import ConvexHull, Delaunay, cKDTree

from . import crystal as _crystal
from . import mesh as _mesh
from . import mode_i_first_passage_v9_18_5 as _v9185
from . import mode_i_first_passage_v9_18_5_2 as _v91852
from . import mode_i_first_passage_v9_18_5_3 as _v91853
from . import physical_refinement_mesh_v100510 as _physical

MODEL_ID = "path_aware_crystallographic_growth_envelope_v10_0_5_18_3_4"
CORRIDOR_SCHEMA = "v10.0.5.18.3.4_path_aware_crystallographic_growth_envelope"
PATH_AWARE_PHYSICAL_POLICY = (
    "global_forward_reachable_envelope_process_zone_support_with_initial_tip_grading"
)


def _float_env(name: str, default: float) -> float:
    try:
        value = float(os.environ.get(name, default))
    except (TypeError, ValueError):
        value = float(default)
    return value if math.isfinite(value) else float(default)


def _required_positive_env(name: str) -> float:
    value = _float_env(name, float("nan"))
    if not math.isfinite(value) or value <= 0.0:
        raise RuntimeError(
            f"v10.0.5.18.3.4 requires positive {name}; "
            f"observed {os.environ.get(name)!r}"
        )
    return value


def _is_production_physical_constructor(function: Any) -> bool:
    return (
        getattr(function, "__module__", "") == _physical.__name__
        and getattr(function, "__name__", "")
        == "make_physical_refinement_mesh_v100510"
    )


def _global_forward_directions(
    minimum_cosine: float,
    sample_count: int = 181,
) -> np.ndarray:
    """Sample the complete production global-forward admissibility cone."""
    gmin = float(minimum_cosine)
    if not math.isfinite(gmin) or gmin > 1.0:
        raise ValueError("minimum global-forward cosine must be finite and <= 1")
    if gmin <= -1.0:
        angles = np.linspace(-math.pi, math.pi, max(int(sample_count), 361))
    else:
        half_angle = math.acos(max(gmin, -1.0))
        angles = np.linspace(
            -half_angle,
            half_angle,
            max(int(sample_count), 3),
        )
    return np.column_stack((np.cos(angles), np.sin(angles)))


def _forward_oriented_crystallographic_traces(
    theta_deg: float,
    minimum_cosine: float,
) -> list[dict[str, Any]]:
    """Audit the nominal BCC {100} traces using the production global gate."""
    fwd = np.array([1.0, 0.0], dtype=float)
    rows: list[dict[str, Any]] = []
    for trace in _crystal.bcc_cleavage_traces(float(theta_deg)):
        direction = np.asarray(trace["t"], dtype=float)
        if float(direction @ fwd) < float((-direction) @ fwd):
            direction = -direction
        forward_cosine = float(direction @ fwd)
        rows.append(
            {
                "name": str(trace["name"]),
                "angle_deg": float(
                    math.degrees(math.atan2(direction[1], direction[0]))
                ),
                "direction": direction.tolist(),
                "global_forward_cosine": forward_cosine,
                "globally_admissible": bool(
                    minimum_cosine <= -1.0
                    or forward_cosine >= minimum_cosine
                ),
            }
        )
    return rows


def _reachable_envelope(
    geom: Any,
    path_length_m: float,
    minimum_cosine: float,
    direction_samples: int,
) -> tuple[np.ndarray, np.ndarray, ConvexHull]:
    """Return directions, sampled endpoints, and the convex reachable envelope.

    For arbitrary direction changes, a path of cumulative arclength no greater
    than ``path_length_m`` can terminate in the convex hull of the initial tip
    and the length-scaled admissible unit directions.
    """
    start = np.array([float(geom.a0), 0.0], dtype=float)
    directions = _global_forward_directions(minimum_cosine, direction_samples)
    endpoints = start[None, :] + float(path_length_m) * directions
    points = np.vstack((start[None, :], endpoints))
    hull = ConvexHull(points)
    vertices = points[hull.vertices]

    tol = 1.0e-12
    inside_domain = (
        np.all(vertices[:, 0] >= -tol)
        and np.all(vertices[:, 0] <= float(geom.Lx) + tol)
        and np.all(vertices[:, 1] >= -0.5 * float(geom.Ly) - tol)
        and np.all(vertices[:, 1] <= 0.5 * float(geom.Ly) + tol)
    )
    if not inside_domain:
        raise RuntimeError(
            "v10.0.5.18.3.4 reachable path envelope exceeds the specimen "
            "domain before the requested target plus guard is completed"
        )
    return directions, endpoints, hull


def _signed_support_distance(points: np.ndarray, hull: ConvexHull) -> np.ndarray:
    """Conservative signed distance to the convex hull supporting lines."""
    pts = np.asarray(points, dtype=float).reshape(-1, 2)
    normals = np.asarray(hull.equations[:, :2], dtype=float)
    offsets = np.asarray(hull.equations[:, 2], dtype=float)
    norms = np.linalg.norm(normals, axis=1)
    return np.max(
        (pts @ normals.T + offsets[None, :])
        / np.maximum(norms[None, :], 1.0e-300),
        axis=1,
    )


def _triangular_lattice(geom: Any, spacing_m: float) -> np.ndarray:
    """Generate one deterministic hexagonal/triangular lattice over the domain."""
    h = float(spacing_m)
    if h <= 0.0:
        raise ValueError("triangular lattice spacing must be positive")
    Lx = float(geom.Lx)
    Ly = float(geom.Ly)
    dy = 0.5 * math.sqrt(3.0) * h
    rows: list[np.ndarray] = []
    y_values = np.arange(-0.5 * Ly, 0.5 * Ly + 0.5 * dy, dy)
    for row_index, y in enumerate(y_values):
        offset = 0.5 * h if row_index % 2 else 0.0
        x_values = np.arange(offset, Lx + 0.5 * h, h)
        keep = (
            (x_values >= 0.0)
            & (x_values <= Lx)
            & (y >= -0.5 * Ly - 1.0e-15)
            & (y <= 0.5 * Ly + 1.0e-15)
        )
        if np.any(keep):
            rows.append(
                np.column_stack(
                    (
                        x_values[keep],
                        np.full(int(np.count_nonzero(keep)), y),
                    )
                )
            )
    return np.vstack(rows)


def _domain_boundary_nodes(geom: Any, spacing_m: float) -> np.ndarray:
    """Exact boundary nodes without near-coincident lattice boundary points."""
    Lx = float(geom.Lx)
    Ly = float(geom.Ly)
    h = float(spacing_m)
    nx = max(2, int(math.ceil(Lx / h)))
    ny = max(2, int(math.ceil(Ly / h)))
    rows: list[tuple[float, float]] = []
    for x in np.linspace(0.0, Lx, nx + 1):
        rows.extend(((float(x), -0.5 * Ly), (float(x), 0.5 * Ly)))
    for y in np.linspace(-0.5 * Ly, 0.5 * Ly, ny + 1):
        rows.extend(((0.0, float(y)), (Lx, float(y))))
    return np.asarray(rows, dtype=float)


def _initial_tip_graded_nodes(
    geom: Any,
    mesh_cfg: Any,
    support_spacing_m: float,
    physical_radius_m: float,
) -> tuple[np.ndarray, float]:
    """Retain the production radial size law until it meets support spacing."""
    h_fine = float(getattr(mesh_cfg, "tip_h_fine", 0.0) or 0.0)
    if h_fine <= 0.0:
        raise ValueError("path-aware envelope requires a graded tip mesh")
    ratio = float(getattr(mesh_cfg, "tip_ratio", 1.15) or 1.15)
    slope = max(ratio - 1.0, 0.02)
    support_h = max(float(support_spacing_m), h_fine)
    blend_radius = min(
        float(physical_radius_m),
        max((support_h - h_fine) / slope + support_h, 4.0 * h_fine),
    )

    start = np.array([float(geom.a0), 0.0], dtype=float)
    points: list[tuple[float, float]] = [(float(start[0]), float(start[1]))]
    radius = h_fine
    while radius < blend_radius:
        local_h = min(h_fine + slope * radius, support_h)
        count = max(8, int(math.ceil(2.0 * math.pi * radius / local_h)))
        for index in range(count):
            angle = 2.0 * math.pi * index / count
            x = float(start[0] + radius * math.cos(angle))
            y = float(start[1] + radius * math.sin(angle))
            if (
                0.0 <= x <= float(geom.Lx)
                and -0.5 * float(geom.Ly) <= y <= 0.5 * float(geom.Ly)
            ):
                points.append((x, y))
        radius += local_h
    return np.asarray(points, dtype=float), float(blend_radius)


def _nested_envelope_nodes(
    geom: Any,
    mesh_cfg: Any,
    hull: ConvexHull,
    process_zone_m: float,
    support_ratio: float,
    physical_radius_m: float,
) -> tuple[np.ndarray, dict[str, Any]]:
    """Build compatible nested triangular lattices around the reachable hull."""
    h_support = float(support_ratio) * float(process_zone_m)
    if h_support <= 0.0:
        raise ValueError("support spacing must be positive")
    h_levels = [h_support * (2.0**level) for level in range(4)]
    shell_buffers = [
        max(
            _float_env("ARRHENIUS_ENVELOPE_BUFFER_LPZ", 1.0)
            * process_zone_m,
            2.0 * h_support,
        ),
        max(3.0 * process_zone_m, 4.0 * h_levels[1]),
        max(float(physical_radius_m), 6.0 * process_zone_m),
        float("inf"),
    ]
    h_far = h_levels[-1]
    boundary_margin = 0.35 * h_far

    clouds: list[np.ndarray] = []
    for spacing, buffer_m in zip(h_levels, shell_buffers):
        points = _triangular_lattice(geom, spacing)
        if math.isfinite(buffer_m):
            distance = _signed_support_distance(points, hull)
            points = points[distance <= buffer_m + 1.0e-12]
        near_boundary = (
            (points[:, 0] < boundary_margin)
            | (points[:, 0] > float(geom.Lx) - boundary_margin)
            | (points[:, 1] < -0.5 * float(geom.Ly) + boundary_margin)
            | (points[:, 1] > 0.5 * float(geom.Ly) - boundary_margin)
        )
        clouds.append(points[~near_boundary])

    graded, blend_radius = _initial_tip_graded_nodes(
        geom,
        mesh_cfg,
        h_support,
        physical_radius_m,
    )
    start = np.array([float(geom.a0), 0.0], dtype=float)
    radial_distance = np.linalg.norm(clouds[0] - start[None, :], axis=1)
    clouds[0] = clouds[0][radial_distance >= 0.85 * blend_radius]

    boundary = _domain_boundary_nodes(geom, h_far)
    nodes = np.vstack([*clouds, boundary, graded])
    h_fine = float(getattr(mesh_cfg, "tip_h_fine", 0.0) or 0.0)
    key_scale = max(1.0e-12, 0.05 * min(h_fine, h_support))
    key = np.round(nodes / key_scale).astype(np.int64)
    _, keep = np.unique(key, axis=0, return_index=True)
    nodes = nodes[np.sort(keep)]
    return nodes, {
        "support_spacing_m": h_support,
        "support_spacing_over_L_pz": float(support_ratio),
        "nested_lattice_spacings_m": [float(value) for value in h_levels],
        "nested_lattice_shell_buffers_m": [
            None if not math.isfinite(value) else float(value)
            for value in shell_buffers
        ],
        "initial_tip_blend_radius_m": blend_radius,
        "initial_tip_h_fine_m": h_fine,
        "initial_tip_size_law": "h_fine_plus_tip_ratio_minus_one_times_radius",
    }


def _make_path_aware_mesh(
    geom: Any,
    mesh_cfg: Any,
    hull: ConvexHull,
    support_centers: np.ndarray,
    process_zone_m: float,
    support_ratio: float,
    seed: int | None,
):
    spec = _physical._ACTIVE_SPEC
    if spec is None:
        raise RuntimeError("production physical refinement mesh was not configured")
    physical_radius = float(spec.radius_m)
    spec.validate(geom, np.asarray(hull.points[hull.vertices], dtype=float))
    if seed is not None:
        np.random.seed(seed)

    nodes, construction = _nested_envelope_nodes(
        geom,
        mesh_cfg,
        hull,
        process_zone_m,
        support_ratio,
        physical_radius,
    )
    triangulation = Delaunay(nodes)
    elems = triangulation.simplices
    centroids = nodes[elems].mean(axis=1)
    inside = (
        (centroids[:, 0] >= 0.0)
        & (centroids[:, 0] <= float(geom.Lx))
        & (centroids[:, 1] >= -0.5 * float(geom.Ly))
        & (centroids[:, 1] <= 0.5 * float(geom.Ly))
    )
    mesh = _mesh.rebuild_tri_mesh(
        nodes,
        elems[inside],
        tip_centers=np.asarray(support_centers, dtype=float),
        validate=True,
    )
    mesh.production_refinement_radius_m = physical_radius
    mesh.production_refinement_policy = PATH_AWARE_PHYSICAL_POLICY
    mesh.production_refinement_centers_m = np.asarray(
        support_centers,
        dtype=float,
    ).tolist()
    mesh.production_refinement_path_envelope_vertices_m = np.asarray(
        hull.points[hull.vertices],
        dtype=float,
    ).tolist()
    mesh.production_refinement_path_aware_envelope = True
    mesh.production_refinement_swept_capsule = False
    mesh.production_refinement_support_spacing_m = float(
        construction["support_spacing_m"]
    )
    mesh.production_refinement_construction = construction
    return mesh, construction


def _resolution_sample_points(
    geom: Any,
    hull: ConvexHull,
    sample_spacing_m: float,
) -> np.ndarray:
    points = _triangular_lattice(geom, sample_spacing_m)
    points = points[_signed_support_distance(points, hull) <= 1.0e-12]
    vertices = np.asarray(hull.points[hull.vertices], dtype=float)
    points = np.vstack((points, vertices))
    key_scale = max(1.0e-12, 0.05 * sample_spacing_m)
    key = np.round(points / key_scale).astype(np.int64)
    _, keep = np.unique(key, axis=0, return_index=True)
    return points[np.sort(keep)]


def _path_envelope_resolution(
    mesh: Any,
    hull: ConvexHull,
    da_m: float,
    process_zone_m: float,
) -> dict[str, Any]:
    """Certify local mean edge length over the full reachable envelope."""
    sample_spacing = max(
        float(da_m),
        _float_env("ARRHENIUS_ENVELOPE_AUDIT_SAMPLE_LPZ", 0.5)
        * float(process_zone_m),
    )
    sample_points = _resolution_sample_points(
        type(
            "_GeometryView",
            (),
            {
                "Lx": float(np.max(mesh.nodes[:, 0])),
                "Ly": float(
                    np.max(mesh.nodes[:, 1]) - np.min(mesh.nodes[:, 1])
                ),
            },
        )(),
        hull,
        sample_spacing,
    )
    centroids = np.asarray(mesh.nodes, float)[
        np.asarray(mesh.elems, int)
    ].mean(axis=1)
    nearest_count = max(
        4,
        int(_float_env("ARRHENIUS_ENVELOPE_AUDIT_NEAREST_ELEMENTS", 16.0)),
    )
    nearest_count = min(nearest_count, len(centroids))
    tree = cKDTree(centroids)
    _, indices = tree.query(sample_points, k=nearest_count)
    if nearest_count == 1:
        indices = np.asarray(indices, dtype=int)[:, None]

    local_mean: list[float] = []
    local_maximum: list[float] = []
    elems = np.asarray(mesh.elems, int)
    nodes = np.asarray(mesh.nodes, float)
    for row in np.asarray(indices, dtype=int):
        selected = elems[np.asarray(row, dtype=int)]
        edges = np.vstack(
            (
                selected[:, [0, 1]],
                selected[:, [1, 2]],
                selected[:, [2, 0]],
            )
        )
        lengths = np.linalg.norm(
            nodes[edges[:, 1]] - nodes[edges[:, 0]],
            axis=1,
        )
        local_mean.append(float(np.mean(lengths)))
        local_maximum.append(float(np.max(lengths)))

    local_mean_array = np.asarray(local_mean, dtype=float)
    local_maximum_array = np.asarray(local_maximum, dtype=float)
    worst_index = int(np.argmax(local_mean_array))
    return {
        "resolution_metric": (
            "mean_edge_length_over_nearest_fixed_element_patch"
        ),
        "resolution_nearest_element_count": int(nearest_count),
        "resolution_sample_spacing_m": float(sample_spacing),
        "resolution_sample_count": int(len(sample_points)),
        "maximum_sampled_hbar_tip_m": float(local_mean_array[worst_index]),
        "maximum_sampled_hbar_tip_over_da": float(
            local_mean_array[worst_index] / max(float(da_m), 1.0e-300)
        ),
        "maximum_sampled_hbar_tip_over_L_pz": float(
            local_mean_array[worst_index]
            / max(float(process_zone_m), 1.0e-300)
        ),
        "maximum_sampled_hbar_tip_location_m": sample_points[
            worst_index
        ].tolist(),
        "maximum_sampled_local_edge_m": float(np.max(local_maximum_array)),
        "maximum_sampled_local_edge_over_L_pz": float(
            np.max(local_maximum_array)
            / max(float(process_zone_m), 1.0e-300)
        ),
    }


def _support_centers(
    hull: ConvexHull,
    directions: np.ndarray,
    start: np.ndarray,
    path_length_m: float,
) -> np.ndarray:
    """Small representative set used only for mesh rebuild diagnostics."""
    fractions = np.linspace(0.0, 1.0, 5)
    chosen = directions[
        np.linspace(0, len(directions) - 1, 9, dtype=int)
    ]
    points = [
        start + fraction * path_length_m * direction
        for direction in chosen
        for fraction in fractions
    ]
    points.extend(np.asarray(hull.points[hull.vertices], dtype=float))
    values = np.asarray(points, dtype=float)
    key = np.round(values / 1.0e-12).astype(np.int64)
    _, keep = np.unique(key, axis=0, return_index=True)
    return values[np.sort(keep)]


def path_aware_long_growth_envelope_mesh(
    geom,
    mesh_cfg,
    seed=None,
    tip_center=None,
):
    """Build a process-zone-resolved mesh for every globally admissible path."""
    original = path_aware_long_growth_envelope_mesh._original
    graded = float(getattr(mesh_cfg, "tip_h_fine", 0.0) or 0.0) > 0.0
    enabled = os.environ.get("ARRHENIUS_PREFINED_MODE_I_CORRIDOR", "1") not in {
        "0",
        "false",
        "False",
        "no",
        "NO",
    }
    if tip_center is not None or not graded or not enabled:
        return original(geom, mesh_cfg, seed=seed, tip_center=tip_center)

    target_um = _required_positive_env("ARRHENIUS_COMMITTED_TARGET_EXTENSION_UM")
    guard_um = max(_float_env("ARRHENIUS_CORRIDOR_GUARD_UM", 10.0), 0.0)
    da_um = _required_positive_env("ARRHENIUS_PHYSICAL_DA_UM")
    process_zone_um = _required_positive_env(
        "ARRHENIUS_CORRIDOR_PROCESS_ZONE_UM"
    )
    theta_deg = _float_env("ARRHENIUS_CRYSTAL_THETA_DEG", 45.0)
    min_global_forward = _float_env(
        "ARRHENIUS_MIN_GLOBAL_FORWARD",
        0.05,
    )
    direction_samples = max(
        21,
        int(_float_env("ARRHENIUS_ENVELOPE_DIRECTION_SAMPLES", 181.0)),
    )
    max_h_over_lpz = max(
        _float_env("ARRHENIUS_MAX_CORRIDOR_H_OVER_LPZ", 0.25),
        1.0e-6,
    )
    audit_h_over_da = max(
        _float_env("ARRHENIUS_MAX_TIP_H_OVER_DA", 0.75),
        1.0e-6,
    )
    qfloor = _float_env(
        "ARRHENIUS_MIN_INITIAL_TRIANGLE_QUALITY",
        _float_env("ARRHENIUS_MIN_ACCEPTED_TRIANGLE_QUALITY", 0.035),
    )
    initial_support_ratio = min(
        _float_env("ARRHENIUS_ENVELOPE_SUPPORT_H_OVER_LPZ", 0.245),
        0.98 * max_h_over_lpz,
    )
    refinement_factor = min(
        max(
            _float_env(
                "ARRHENIUS_ENVELOPE_SUPPORT_REFINEMENT_FACTOR",
                0.9,
            ),
            0.5,
        ),
        0.99,
    )
    candidate_count = max(
        1,
        int(_float_env("ARRHENIUS_ENVELOPE_SUPPORT_CANDIDATES", 4.0)),
    )

    path_length_m = (target_um + guard_um) * 1.0e-6
    da_m = da_um * 1.0e-6
    process_zone_m = process_zone_um * 1.0e-6
    directions, endpoints, hull = _reachable_envelope(
        geom,
        path_length_m,
        min_global_forward,
        direction_samples,
    )
    start = np.array([float(geom.a0), 0.0], dtype=float)
    centers = _support_centers(
        hull,
        directions,
        start,
        path_length_m,
    )
    physical_provider = _is_production_physical_constructor(original)
    if not physical_provider:
        raise RuntimeError(
            "v10.0.5.18.3.4 requires the validated production physical "
            "refinement provider"
        )

    candidates: list[dict[str, Any]] = []
    selected = None
    selected_audit: dict[str, Any] | None = None
    selected_ratio: float | None = None
    selected_construction: dict[str, Any] | None = None

    ratios = [
        initial_support_ratio * (refinement_factor**index)
        for index in range(candidate_count)
    ]
    for support_ratio in ratios:
        try:
            raw, construction = _make_path_aware_mesh(
                geom,
                mesh_cfg,
                hull,
                centers,
                process_zone_m,
                support_ratio,
                seed,
            )
            compact, audit = _v91853._compact_without_quality_abort(
                raw,
                centers,
            )
            resolution = _path_envelope_resolution(
                compact,
                hull,
                da_m,
                process_zone_m,
            )
            qmin = float(audit["minimum_initial_triangle_quality"])
            h_over_lpz = float(
                resolution["maximum_sampled_hbar_tip_over_L_pz"]
            )
            quality_ok = bool(qmin >= qfloor)
            process_zone_ok = bool(h_over_lpz <= max_h_over_lpz + 1.0e-12)
            ok = bool(quality_ok and process_zone_ok)
            row = {
                "support_spacing_over_L_pz": float(support_ratio),
                "support_spacing_m": float(
                    construction["support_spacing_m"]
                ),
                "node_count": int(compact.nn),
                "triangle_count": int(compact.ne),
                "minimum_initial_triangle_quality": qmin,
                "triangle_quality_required": qfloor,
                **resolution,
                "maximum_hbar_tip_over_L_pz_required": max_h_over_lpz,
                "tip_h_over_da_audit_threshold": audit_h_over_da,
                "tip_h_over_da_enforced_as_veto": False,
                "tip_h_over_da_resolution_warning": bool(
                    resolution["maximum_sampled_hbar_tip_over_da"]
                    > audit_h_over_da
                ),
                "mesh_provider": PATH_AWARE_PHYSICAL_POLICY,
                "production_refinement_radius_m": getattr(
                    compact,
                    "production_refinement_radius_m",
                    None,
                ),
                "accepted": ok,
                "error": None,
            }
            candidates.append(row)
            if ok:
                selected = compact
                selected_audit = {**audit, **resolution, **row}
                selected_ratio = float(support_ratio)
                selected_construction = construction
                break
        except Exception as exc:
            candidates.append(
                {
                    "support_spacing_over_L_pz": float(support_ratio),
                    "mesh_provider": PATH_AWARE_PHYSICAL_POLICY,
                    "accepted": False,
                    "error": f"{type(exc).__name__}: {exc}",
                }
            )

    constructor_module = getattr(original, "__module__", None)
    constructor_qualname = getattr(original, "__qualname__", None)
    envelope_vertices = np.asarray(
        hull.points[hull.vertices],
        dtype=float,
    )
    trace_audit = _forward_oriented_crystallographic_traces(
        theta_deg,
        min_global_forward,
    )
    base_payload = {
        "schema": CORRIDOR_SCHEMA,
        "model": MODEL_ID,
        "candidate_envelopes": candidates,
        "candidate_support_h_over_L_pz": [float(value) for value in ratios],
        "candidate_search_stopped_after_first_admissible": True,
        "corridor_target_extension_um": target_um,
        "corridor_guard_um": guard_um,
        "reachable_path_arclength_um": target_um + guard_um,
        "process_zone_length_um": process_zone_um,
        "maximum_hbar_tip_over_L_pz_required": max_h_over_lpz,
        "minimum_initial_triangle_quality_required": qfloor,
        "physical_da_um": da_um,
        "tip_h_over_da_role": "audit_warning_only",
        "target_propagated_before_mesh_construction": True,
        "full_requested_corridor_covered": True,
        "full_requested_reachable_envelope_covered": True,
        "path_aware_growth_envelope": True,
        "directional_support_basis": (
            "complete_global_forward_gate_conservative_over_"
            "continuous_crack_direction_competition"
        ),
        "global_forward_reference_direction": [1.0, 0.0],
        "minimum_global_forward_cosine": min_global_forward,
        "direction_sample_count": int(len(directions)),
        "reachable_envelope_vertices_m": envelope_vertices.tolist(),
        "reachable_envelope_x_min_m": float(np.min(envelope_vertices[:, 0])),
        "reachable_envelope_x_max_m": float(np.max(envelope_vertices[:, 0])),
        "reachable_envelope_y_min_m": float(np.min(envelope_vertices[:, 1])),
        "reachable_envelope_y_max_m": float(np.max(envelope_vertices[:, 1])),
        "reachable_extreme_endpoints_m": [
            endpoints[0].tolist(),
            endpoints[len(endpoints) // 2].tolist(),
            endpoints[-1].tolist(),
        ],
        "crystal_theta_deg": theta_deg,
        "nominal_BCC_100_trace_audit": trace_audit,
        "mesh_constructor_module": constructor_module,
        "mesh_constructor_qualname": constructor_qualname,
        "production_physical_provider_detected": physical_provider,
        "swept_physical_refinement_active": False,
        "path_aware_physical_refinement_active": physical_provider,
        "path_aware_physical_refinement_policy": (
            PATH_AWARE_PHYSICAL_POLICY
        ),
        "constitutive_physics_changed": False,
    }

    if (
        selected is None
        or selected_audit is None
        or selected_ratio is None
        or selected_construction is None
    ):
        _v91852._STARTUP_AUDIT.clear()
        _v91852._STARTUP_AUDIT.update(base_payload)
        raise RuntimeError(
            "v10.0.5.18.3.4 found no path-aware reachable-envelope mesh "
            f"satisfying initial triangle quality >= {qfloor:.6g} and "
            f"h_tip/L_pz <= {max_h_over_lpz:.6g}; support ratios tried: "
            + ", ".join(f"{value:.6g}" for value in ratios)
        )

    payload = {
        **selected_audit,
        **base_payload,
        "selected_support_h_over_L_pz": selected_ratio,
        "selected_support_spacing_m": float(
            selected_construction["support_spacing_m"]
        ),
        "selected_mesh_construction": selected_construction,
        "selected_support_center_count": int(len(centers)),
        "selected_support_centers_m": centers.tolist(),
    }
    _v91852._STARTUP_AUDIT.clear()
    _v91852._STARTUP_AUDIT.update(payload)
    _v9185._RUNTIME["corridor_centers"] = centers.tolist()
    _v9185._RUNTIME["mesh"] = selected
    return selected


__all__ = [
    "CORRIDOR_SCHEMA",
    "MODEL_ID",
    "PATH_AWARE_PHYSICAL_POLICY",
    "_global_forward_directions",
    "_path_envelope_resolution",
    "_reachable_envelope",
    "path_aware_long_growth_envelope_mesh",
]
