from __future__ import annotations

import numpy as np

from arrhenius_fracture.config import GeometryConfig, MeshConfig
from arrhenius_fracture import long_growth_corridor_v10051833 as corridor
from arrhenius_fracture import mode_i_first_passage_v9_18_5_2 as v91852
from arrhenius_fracture import mode_i_first_passage_v10_0_5_18_3_3_four_class_long_growth_corridor as entry
from arrhenius_fracture import physical_refinement_mesh_v100510 as physical


def test_entrypoint_detects_active_physical_spec_through_wrapped_provider(
    monkeypatch,
    tmp_path,
):
    observed = {}

    def wrapped_physical_provider(geom, mesh_cfg, seed=None, tip_center=None):
        return physical.make_physical_refinement_mesh_v100510(
            geom,
            mesh_cfg,
            seed=seed,
            tip_center=tip_center,
        )

    def fake_base(_argv):
        geom = GeometryConfig()
        mesh_cfg = MeshConfig(
            nx=36,
            ny=72,
            tip_h_fine=2.5e-6,
            tip_ratio=1.15,
        )
        physical.configure_physical_refinement_v100510(330.0e-6)
        had_original = hasattr(corridor.target_aware_long_growth_corridor_mesh, "_original")
        saved_original = getattr(
            corridor.target_aware_long_growth_corridor_mesh,
            "_original",
            None,
        )
        corridor.target_aware_long_growth_corridor_mesh._original = (
            wrapped_physical_provider
        )
        try:
            mesh = corridor.target_aware_long_growth_corridor_mesh(
                geom,
                mesh_cfg,
                seed=42,
            )
            audit = dict(v91852._STARTUP_AUDIT)
        finally:
            physical.clear_physical_refinement_v100510()
            if had_original:
                corridor.target_aware_long_growth_corridor_mesh._original = (
                    saved_original
                )
            else:
                delattr(
                    corridor.target_aware_long_growth_corridor_mesh,
                    "_original",
                )

        observed.update({
            "mesh": mesh,
            "audit": audit,
        })
        return "ok"

    monkeypatch.setattr(entry._base, "main", fake_base)
    original_detector = corridor._is_production_physical_constructor

    result = entry.main([
        "--target-crack-extension-um", "400",
        "--da-phys", "5e-6",
        "--mpz-length-um", "50",
        "--out", str(tmp_path),
    ])

    assert result == "ok"
    assert corridor._is_production_physical_constructor is original_detector

    mesh = observed["mesh"]
    audit = observed["audit"]
    assert audit["production_physical_provider_detected"] is True
    assert audit["swept_physical_refinement_active"] is True
    assert audit["swept_physical_refinement_policy"] == corridor.SWEPT_PHYSICAL_POLICY
    assert audit["minimum_initial_triangle_quality"] >= 0.035
    assert audit["maximum_sampled_hbar_tip_over_L_pz"] <= 0.25
    assert audit["production_refinement_radius_m"] == 330.0e-6
    assert mesh.production_refinement_swept_capsule is True
    assert mesh.production_refinement_policy == corridor.SWEPT_PHYSICAL_POLICY
    assert np.all(np.isfinite(mesh.area_e))
    assert np.all(mesh.area_e > 0.0)
