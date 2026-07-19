from __future__ import annotations

import numpy as np
import pytest

from arrhenius_fracture.config import GeometryConfig, MeshConfig
from arrhenius_fracture.fixed_grip_audit_mesh_v10058 import (
    clear_fixed_grip_audit_mesh,
    configure_fixed_grip_audit_mesh,
    make_fixed_radius_audit_mesh,
)


def test_fixed_radius_audit_mesh_contains_every_perturbed_tip_node():
    geom = GeometryConfig(
        Lx=1.0e-3,
        Ly=2.0e-3,
        a0=0.25e-3,
        notch_half_thickness=40.0e-6,
    )
    spec = configure_fixed_grip_audit_mesh(
        crack_center_m=geom.a0,
        crack_increments_m=(40.0e-6, 20.0e-6),
        refinement_radius_m=120.0e-6,
    )
    try:
        meshes = [
            make_fixed_radius_audit_mesh(
                geom,
                MeshConfig(nx=16, ny=32, jitter=0.0, tip_h_fine=h, tip_ratio=1.20),
                seed=42,
            )
            for h in (20.0e-6, 10.0e-6)
        ]
    finally:
        clear_fixed_grip_audit_mesh()

    expected = np.asarray(
        [[spec.crack_center_m + offset, 0.0] for offset in spec.crack_offsets_m]
    )
    assert len(expected) == 5
    for mesh in meshes:
        for point in expected:
            distances = np.linalg.norm(mesh.nodes - point[None, :], axis=1)
            assert float(np.min(distances)) <= 1.0e-12


def test_fixed_radius_audit_mesh_rejects_boundary_intersection():
    with pytest.raises(ValueError):
        configure_fixed_grip_audit_mesh(
            crack_center_m=0.25e-3,
            crack_increments_m=(20.0e-6,),
            refinement_radius_m=0.30e-3,
        ).validate(width_m=1.0e-3, height_m=2.0e-3)
    clear_fixed_grip_audit_mesh()
