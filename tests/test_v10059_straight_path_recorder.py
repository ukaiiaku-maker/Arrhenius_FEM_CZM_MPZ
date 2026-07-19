from __future__ import annotations

import inspect

from arrhenius_fracture import sharp_front
from arrhenius_fracture.mode_i_first_passage_v10_0_5_9_production_j_probe import (
    patch_run_2d_source_v10059,
    validate_source_transform_v10059,
)


def test_straight_progressive_recorder_does_not_read_uninitialized_kill_r():
    patched = patch_run_2d_source_v10059(inspect.getsource(sharp_front.run_2d))
    assert "kill_radius_v10059 = float(kill_r)" in patched
    assert "kill_radius_v10059 = 0.0" in patched
    assert "production_exclude_v10059 = 0.0" in patched
    assert "kill_r_m=kill_radius_v10059" in patched
    assert "kill_r_m=kill_r," not in patched
    assert "straight_progressive_cluster_no_exclusion" in patched


def test_source_transform_reports_both_production_path_semantics():
    audit = validate_source_transform_v10059()
    assert audit["root_front_production_exclusion"] is True
    assert audit["straight_progressive_no_exclusion"] is True
    assert audit["no_unconditional_kill_r_read"] is True
    assert audit["constitutive_physics_changed"] is False
