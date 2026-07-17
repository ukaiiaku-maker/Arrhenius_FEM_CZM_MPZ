from __future__ import annotations

from arrhenius_fracture import sharp_front
from arrhenius_fracture.kinetic_progressive_2d_v10 import (
    build_progressive_run_2d,
    progressive_runtime_payload,
    reset_progressive_runtime,
)


def test_progressive_transform_compiles_against_actual_run_2d():
    reset_progressive_runtime()
    transformed = build_progressive_run_2d(sharp_front.run_2d)
    assert transformed is not sharp_front.run_2d
    assert transformed.__name__ == sharp_front.run_2d.__name__
    assert transformed._v10_progressive_source_transform is True
    payload = progressive_runtime_payload()
    assert payload["anchor_counts"] == {
        "backend_construction": 1,
        "adaptive_target": 1,
        "single_front_step": 1,
    }
    assert payload["active"] is False
