from arrhenius_fracture import fem, j_integral, mesh, plasticity, sharp_front
from arrhenius_fracture.kinetic_progressive_2d_v1003 import (
    _restore_live_aliases,
    _sync_live_aliases,
)
from arrhenius_fracture.kinetic_progressive_2d_v1003_source import (
    build_progressive_run_2d_v1003_source,
)


def _sentinel(*args, **kwargs):
    raise AssertionError("sentinel should not execute during transform construction")


def test_transformed_namespace_captures_current_live_bindings(monkeypatch):
    original = sharp_front.run_2d
    monkeypatch.setattr(mesh, "make_tri_mesh", _sentinel)
    monkeypatch.setattr(fem, "assemble_mechanics", _sentinel)
    monkeypatch.setattr(fem, "solve_dirichlet", _sentinel)
    monkeypatch.setattr(j_integral, "compute_J_integral", _sentinel)
    monkeypatch.setattr(plasticity, "update_plasticity", _sentinel)
    monkeypatch.setattr(sharp_front, "build_engine", _sentinel)

    saved, _ = _sync_live_aliases(original)
    try:
        transformed = build_progressive_run_2d_v1003_source(original)
        for name in (
            "make_tri_mesh",
            "assemble_mechanics",
            "solve_dirichlet",
            "compute_J_integral",
            "update_plasticity",
            "build_engine",
        ):
            assert transformed.__globals__[name] is _sentinel
    finally:
        _restore_live_aliases(original, saved)


def test_v1003_source_adapter_fixes_accounting_and_state_compatibility():
    transformed = build_progressive_run_2d_v1003_source(sharp_front.run_2d)
    assert transformed._v1003_source_adapter is True
    assert transformed._v1003_path_deflection_forced_off is True
    assert transformed._v1003_campaign_state_compatibility is True
    assert transformed._v1003_nondeflect_summary_accounting is True
