"""Regression test for the single-front geometry-veto rollback fix.

Before this fix, `sharp_front.run_2d`'s single-front (non-deflect) accepted-
event path called `eng.step(...)` (which mutates the front engine's hazard/
MPZ state -- B, N_em, a_adv, n_adv -- for every fire) and then attempted the
actual geometry commit via `crack_backend.advance(...)`. If that commit was
vetoed (`not rr.inserted`), the code only printed a message and `break`-ed
out of the event loop: the engine state mutated by the already-completed
`eng.step()` call was never rolled back. This left the crack-tip
microstructure (MPZ/hazard state) permanently desynchronized from the actual
(unadvanced) geometric crack tip -- exactly the "reinitializes/loses state
following a veto" failure mode FEM_CZM_HANDOFF.md section 7.3 requires be
prevented ("the complete physical event must fail atomically and restore
the prior state").

The multi-front/deflect path already had this rollback (`moved_now <= 0.0`
branch, calling `eng_f.restore_geometry_veto` or the manual B/N_em/a_adv/
n_adv fallback restore); this test exercises the same mirrored logic now
present in the single-front path.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from arrhenius_fracture import sharp_front
from arrhenius_fracture.config import GeometryConfig
from arrhenius_fracture.crack_backend import AdaptiveCZMBackend, CrackAdvanceResult

MESH_ARGS = ["--nx", "6", "--ny", "10"]


def test_single_front_geometry_veto_restores_engine_state(tmp_path: Path, monkeypatch):
    captured: dict = {}
    original_step = sharp_front.FrontEngine.step

    def _forced_fire_step(self, K, T, dt):
        first_call = not captured.get("forced")
        if first_call:
            captured["forced"] = True
            # Pre-load B so the REAL, unmodified step() logic deterministically
            # fires exactly once (n_fire=1) on this call -- every resulting
            # info-dict field and state mutation is therefore genuine engine
            # behavior, not a fabricated stand-in.
            self.B = 1.5
            captured["engine"] = self
            captured["B_before_fire"] = float(self.B)
            captured["a_adv_before_fire"] = float(self.a_adv)
            captured["n_adv_before_fire"] = int(self.n_adv)
        info = original_step(self, K, T, dt)
        if first_call:
            # N_em evolves continuously via emission kinetics *within* this
            # same step() call, independent of whether cleavage also fired --
            # that emission physics is real and must NOT be rolled back. Only
            # the renewal-triggered wake-shed (self.N_em = N_retained) that
            # happens after this point is what a geometry veto must undo, so
            # the correct restore target is N_em_pre_renewal (the engine's
            # own record of "N_em right before the renewal reset"), not
            # whatever N_em was before this step() call started.
            captured["N_em_before_fire"] = float(info["N_em_pre_renewal"])
        return info

    monkeypatch.setattr(sharp_front.FrontEngine, "step", _forced_fire_step)

    def _veto_advance(self, *, mesh, boundary, damage, displacement, **_kwargs):
        return CrackAdvanceResult(
            mesh=mesh,
            boundary=boundary,
            damage=damage,
            displacement=displacement,
            moved=0.0,
            inserted=False,
            reason="forced_test_veto",
        )

    monkeypatch.setattr(AdaptiveCZMBackend, "advance", _veto_advance)

    out_dir = tmp_path / "out"
    argv = [
        "--mode", "2d",
        "--steps", "1",
        *MESH_ARGS,
        "--temperatures", "1000",
        "--out", str(out_dir),
        "--print-every", "1",
        "--no-plots",
        "--dt", "8.4",
        "--crack-backend", "adaptive_czm",
    ]
    sharp_front.main(argv)

    assert captured.get("forced") is True, "the forced fire never occurred; test setup is invalid"
    eng = captured["engine"]

    # The vetoed fire's B/N_em/a_adv/n_adv mutations must be rolled back to
    # (within benign floating-point roundoff of) their pre-fire values, not
    # left half-applied -- a real rollback bug leaves O(1) differences here
    # (e.g. B off by a full n_fire, N_em off by the whole emitted increment),
    # far larger than the tolerance below.
    assert eng.B == pytest.approx(captured["B_before_fire"], rel=1.0e-6)
    assert eng.N_em == pytest.approx(captured["N_em_before_fire"], rel=1.0e-6)
    assert eng.a_adv == pytest.approx(captured["a_adv_before_fire"], abs=1.0e-15)
    assert eng.n_adv == captured["n_adv_before_fire"]

    # The recorded step row must reflect the rolled-back state, not the
    # stale pre-rollback info dict. steps_*.csv is comma-delimited with a
    # header row (see sharp_front.py's np.savetxt(..., delimiter=',',
    # header=...) call) -- read it by column name, not positional index.
    steps_csv = out_dir / "steps_1000K.csv"
    assert steps_csv.is_file()
    last_row = pd.read_csv(steps_csv).iloc[-1]
    assert last_row["B"] == pytest.approx(captured["B_before_fire"], rel=1.0e-6)
    assert last_row["N_em"] == pytest.approx(captured["N_em_before_fire"], rel=1.0e-6)
    assert last_row["crack_extension_m"] == pytest.approx(0.0, abs=1.0e-15)


def test_single_front_geometry_veto_does_not_advance_crack_tip(tmp_path: Path, monkeypatch):
    """The vetoed geometry attempt must not leave any partial tip advance --
    unchanged from before this fix, but confirmed here alongside the state
    rollback so a future regression in either half is caught together."""
    original_step = sharp_front.FrontEngine.step
    state = {"forced": False}

    def _forced_fire_step(self, K, T, dt):
        if not state["forced"]:
            state["forced"] = True
            self.B = 1.5
        return original_step(self, K, T, dt)

    monkeypatch.setattr(sharp_front.FrontEngine, "step", _forced_fire_step)

    def _veto_advance(self, *, mesh, boundary, damage, displacement, **_kwargs):
        return CrackAdvanceResult(
            mesh=mesh, boundary=boundary, damage=damage, displacement=displacement,
            moved=0.0, inserted=False, reason="forced_test_veto",
        )

    monkeypatch.setattr(AdaptiveCZMBackend, "advance", _veto_advance)

    out_dir = tmp_path / "out"
    argv = [
        "--mode", "2d", "--steps", "1", *MESH_ARGS,
        "--temperatures", "1000", "--out", str(out_dir),
        "--print-every", "1", "--no-plots", "--dt", "8.4",
        "--crack-backend", "adaptive_czm",
    ]
    sharp_front.main(argv)

    assert state["forced"] is True
    steps_csv = out_dir / "steps_1000K.csv"
    last_row = pd.read_csv(steps_csv).iloc[-1]
    assert last_row["a_tip_m"] == pytest.approx(GeometryConfig().a0, rel=1.0e-6)  # unchanged
