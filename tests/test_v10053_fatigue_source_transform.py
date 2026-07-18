from __future__ import annotations

import inspect
import textwrap

from arrhenius_fracture import sharp_front
from arrhenius_fracture.mode_i_first_passage_v10_0_5_3_fatigue import (
    _build_progressive_run_2d_v10053,
    _patch_run_2d_source,
)


def test_single_front_fatigue_dispatch_is_routed_to_progressive_lifecycle():
    source = textwrap.dedent(inspect.getsource(sharp_front.run_2d))
    patched = _patch_run_2d_source(source)
    assert "if fatigue_mode and not kinetic_progressive:" in patched
    assert patched.count("if fatigue_mode and not kinetic_progressive:") == 1


def test_progressive_fatigue_transform_compiles_against_current_run_2d():
    transformed = _build_progressive_run_2d_v10053(sharp_front.run_2d)
    assert transformed._v10053_progressive_fatigue is True
    assert transformed._v10053_legacy_fatigue_commit_bypassed is True
    assert transformed._v10053_constitutive_physics_changed is False
