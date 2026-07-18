from arrhenius_fracture.mode_i_first_passage_v10_0_5_3_fatigue_audited import (
    patch_run_2d_source_v10053,
    validate_source_transform_v10053,
)
from arrhenius_fracture import sharp_front
import inspect


def test_v10053_patch_anchors_match_current_run_2d():
    source = inspect.getsource(sharp_front.run_2d)
    patched = patch_run_2d_source_v10053(source)
    assert "progressive_fatigue_v10053 = bool(fatigue_mode)" in patched
    assert "if fatigue_mode and not kinetic_progressive:" in patched
    assert "cycles_requested_v10053" in patched
    assert "checkpoint_now_v10053" in patched
    assert "def _diag_with_remaining_v10053" in patched
    assert "diag_single_trial = _diag_with_remaining_v10053(" in patched
    assert "diag_single_trial = _diag_with_remaining(" not in patched


def test_v10053_exact_wrapper_chain_compiles():
    payload = validate_source_transform_v10053()
    assert payload["source_transform_preflight_passed"] is True
    assert payload["v1003_source_adapter"] is True
    assert payload["v1002_event_lifecycle"] is True
    assert payload["v10053_audited_adapter"] is True
    assert payload["legacy_fatigue_commit_bypassed"] is True
    assert payload["consumed_cycle_accounting"] is True
    assert payload["single_front_cycle_scope_fixed"] is True
    assert payload["constitutive_physics_changed"] is False
