from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from arrhenius_fracture import mode_i_first_passage_v9_18_5_3 as corridor


def test_compaction_preserves_only_physical_refinement_metadata(monkeypatch):
    raw = SimpleNamespace(
        production_refinement_radius_m=330.0e-6,
        production_refinement_policy="fixed_physical_radius",
        production_refinement_centers_m=[[0.5e-3, 0.0]],
        unrelated_marker="do_not_copy",
    )
    compact = SimpleNamespace()
    centers = np.array([[0.5e-3, 0.0]])

    monkeypatch.setattr(
        corridor._v91852,
        "_compact_mesh",
        lambda source, selected_centers: (
            compact,
            {"minimum_initial_triangle_quality": 0.05},
        ),
    )

    observed, audit = corridor._compact_without_quality_abort(raw, centers)

    assert observed is compact
    assert observed.production_refinement_radius_m == 330.0e-6
    assert observed.production_refinement_policy == "fixed_physical_radius"
    assert observed.production_refinement_centers_m == [[0.5e-3, 0.0]]
    assert not hasattr(observed, "unrelated_marker")
    assert audit["physical_refinement_metadata_preserved"] is True
    assert audit["physical_refinement_metadata_preservation_is_diagnostic_only"] is True
    assert "production_refinement_radius_m" in audit[
        "physical_refinement_metadata_fields"
    ]
