from __future__ import annotations

import inspect
import textwrap

from arrhenius_fracture import sharp_front
from arrhenius_fracture.kinetic_progressive_2d_v1002_anisotropic_straight import (
    build_progressive_run_2d_v1002_anisotropic_straight,
    progressive_runtime_payload_v1002_anisotropic_straight,
    reset_progressive_runtime_v1002_anisotropic_straight,
)


def test_anisotropic_elasticity_is_built_before_legacy_deflect_switch():
    source = textwrap.dedent(inspect.getsource(sharp_front.run_2d))
    anisotropic_D = "if getattr(args, 'crystal_aniso', False):"
    deflect = "deflect = bool(getattr(args, 'crystal_aniso', False))"
    assert source.count(anisotropic_D) == 1
    assert source.count(deflect) == 1
    assert source.index(anisotropic_D) < source.index(deflect)


def test_v1002_adapter_compiles_and_declares_straight_anisotropic_contract():
    reset_progressive_runtime_v1002_anisotropic_straight()
    transformed = build_progressive_run_2d_v1002_anisotropic_straight(
        sharp_front.run_2d
    )
    assert transformed._v1002_event_lifecycle is True
    assert transformed._v1002_anisotropic_elasticity_preserved is True
    assert transformed._v1002_path_deflection_forced_off is True
    assert (
        transformed._v1002_crystal_compete_used_for_v911_validation_only
        is True
    )

    payload = progressive_runtime_payload_v1002_anisotropic_straight()
    assert payload["anisotropic_elasticity_preserved"] is True
    assert payload["anisotropic_J_preserved"] is True
    assert payload["path_deflection_forced_off"] is True
    assert payload["anisotropic_path_selection_active"] is False
