"""Real-engine integration tests for the PF KJ-target controller.

These tests do not require the (currently unavailable, see
CLAUDE_PROGRESS.md's "EXTERNAL BLOCKER" and PF_REFERENCE_REGENERATION_CONTRACT.md)
real PF reference artifacts. They exercise the actual `sharp_front.run_2d`
solver -- real mesh, real FEM assembly/solve, real plasticity update -- via
a small, fast, synthetic target trajectory, proving the controller wiring
(FEM_CZM_HANDOFF.md section 6, Gate 1 of FEM_PF_PARITY_SCORECARD.md) is
transactional and numerically sane against the real engine, not just
against the synthetic callbacks in test_pf_theta0_j_controlled_loading_v10051840.py.

This is Gate 1 evidence only (numerical controller correctness). It says
nothing about physical FEM/PF parity, which requires the real PF reference
data (Gate 2 onward).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from arrhenius_fracture import sharp_front

MESH_ARGS = ["--nx", "6", "--ny", "10"]


def _write_synthetic_target(path: Path, *, n: int = 10, dt_s: float = 8.4) -> None:
    step = np.arange(1, n + 1, dtype=float)
    KJ = np.linspace(1.0e5, 5.0e6, n)
    J = np.linspace(10.0, 200.0, n)
    frame = pd.DataFrame(
        {
            "step": step,
            "dt_cur_s": np.full(n, dt_s),
            "J_effective_direct_J_per_m2": J,
            "J_signed_direct_J_per_m2": J,
            "KJ_Pa_sqrtm": KJ,
            "B": np.zeros(n),
            "N_em": np.zeros(n),
            "crack_extension_m": np.zeros(n),
        }
    )
    frame.to_csv(path, index=False)


def _run(tmp_path: Path, *, steps: int, extra_args: list[str] | None = None) -> Path:
    tmp_path.mkdir(parents=True, exist_ok=True)
    csv_path = tmp_path / "steps_1000K.csv"
    _write_synthetic_target(csv_path)
    out_dir = tmp_path / "out"
    argv = [
        "--mode", "2d",
        "--steps", str(steps),
        *MESH_ARGS,
        "--temperatures", "1000",
        "--out", str(out_dir),
        "--print-every", "1",
        "--no-plots",
        "--dt", "8.4",
    ]
    if extra_args is None:
        extra_args = ["--pf-kj-target-csv", str(csv_path)]
    sharp_front.main(argv + extra_args)
    return out_dir


def _load_audit(out_dir: Path) -> dict:
    audit_path = out_dir / "pf_kj_target_controller_audit_1000K.json"
    assert audit_path.is_file(), f"expected controller audit file at {audit_path}"
    return json.loads(audit_path.read_text())


def test_real_engine_controller_tracks_synthetic_target_within_tolerance(tmp_path: Path):
    out_dir = _run(tmp_path, steps=3)
    audit = _load_audit(out_dir)

    assert audit["schema"] == "pf_kj_target_controller_audit_v10_0_5_18_4_0"
    records = audit["records"]
    assert len(records) == 3

    for i, record in enumerate(records, start=1):
        assert record["step"] == i
        assert record["converged"] is True
        assert 0 <= record["iterations"] <= 25  # default max_iterations
        target = record["KJ_target_Pa_sqrtm"]
        achieved = record["KJ_achieved_Pa_sqrtm"]
        assert achieved == pytest.approx(target, rel=2.0e-3)
        # physical time is inherited from the PF interval (8.4 s/row), not
        # from solver step count -- both fields must agree exactly with the
        # cumulative dt, matching FEM_CZM_HANDOFF.md section 6's requirement
        # that physical time intervals are preserved.
        assert record["physical_time_target_s"] == pytest.approx(8.4 * i, rel=1.0e-12)
        assert record["physical_time_accepted_s"] == pytest.approx(8.4 * i, rel=1.0e-12)


def test_real_engine_rejected_trials_are_recorded_and_precede_the_accepted_value(tmp_path: Path):
    out_dir = _run(tmp_path, steps=1)
    audit = _load_audit(out_dir)
    record = audit["records"][0]

    # A synthetic target starting from Uapp=0 requires more than one trial
    # (the elastic predictor alone under/overshoots off a cold start), so
    # this configuration is expected to exercise the safeguarded-secant
    # rejection path for real, not just converge in one shot.
    assert record["iterations"] >= 2
    assert len(record["rejected_trials"]) == record["iterations"] - 1
    for trial in record["rejected_trials"]:
        assert trial["KJ_achieved_Pa_sqrtm"] != pytest.approx(
            record["KJ_target_Pa_sqrtm"], rel=2.0e-3
        )


def test_real_engine_run_is_deterministic_across_repeats(tmp_path: Path):
    """No cross-trial or cross-run state leakage: an identical invocation
    must reproduce the identical accepted controller trajectory bit-for-bit
    (this front-engine/mesh configuration has no RNG in its path, so this
    is a legitimate exact-equality check, not merely a tolerance check)."""
    out_a = _run(tmp_path / "a", steps=3)
    out_b = _run(tmp_path / "b", steps=3)
    audit_a = _load_audit(out_a)
    audit_b = _load_audit(out_b)

    for record_a, record_b in zip(audit_a["records"], audit_b["records"]):
        assert record_a["Uapp_accepted_m"] == record_b["Uapp_accepted_m"]
        assert record_a["KJ_achieved_Pa_sqrtm"] == record_b["KJ_achieved_Pa_sqrtm"]
        assert record_a["iterations"] == record_b["iterations"]
        assert record_a["rejected_trials"] == record_b["rejected_trials"]


def test_flag_off_path_writes_no_controller_audit_file(tmp_path: Path):
    out_dir = _run(tmp_path, steps=3, extra_args=[])
    audit_path = out_dir / "pf_kj_target_controller_audit_1000K.json"
    assert not audit_path.exists()


def test_controller_driven_accepted_state_matches_direct_single_step_at_same_Uapp(
    tmp_path: Path,
):
    """The strongest available transactionality proof without exposing
    run_2d's internal u/ep_gp/rho_gp arrays directly: the FEM state after
    the controller accepts a step (having tried and discarded several
    candidate Uapp values along the way) must be numerically
    indistinguishable from a single direct step solved once at exactly
    that accepted Uapp. If a rejected trial left any residue in the
    accepted plastic/mechanical state, this would not hold.
    """
    controller_out = _run(tmp_path / "controller", steps=1)
    audit = _load_audit(controller_out)
    record = audit["records"][0]
    assert record["iterations"] >= 2  # exercised the rejection path for real
    accepted_Uapp = record["Uapp_accepted_m"]

    controller_steps = pd.read_csv(controller_out / "steps_1000K.csv")
    assert len(controller_steps) == 1

    # Reproduce the exact same physical step directly: a single fixed-ramp
    # step whose --dU equals the controller's accepted Uapp, with no PF
    # target controller involved at all.
    direct_out = tmp_path / "direct" / "out"
    direct_argv = [
        "--mode", "2d",
        "--steps", "1",
        *MESH_ARGS,
        "--temperatures", "1000",
        "--out", str(direct_out),
        "--print-every", "1",
        "--no-plots",
        "--dt", "8.4",
        "--dU", f"{accepted_Uapp:.17g}",
    ]
    sharp_front.main(direct_argv)
    direct_steps = pd.read_csv(direct_out / "steps_1000K.csv")
    assert len(direct_steps) == 1

    for column in ("Uapp_m", "Ftop_N", "KJ_Pa_sqrtm", "sigma_tip_Pa", "sigma_back_Pa"):
        controller_value = float(controller_steps.iloc[0][column])
        direct_value = float(direct_steps.iloc[0][column])
        assert controller_value == pytest.approx(direct_value, rel=1.0e-9), column


def test_stochastic_identity_is_preserved_across_rejected_trials(tmp_path: Path):
    """The mechanical/plastic proof above says nothing about the
    stochastic Arrhenius first-passage state (RNG streams, cleavage
    threshold, event-length factor, B, N_em, MPZ populations), which lives
    on the front-engine object, not in the mesh arrays. This test forces
    multiple rejected controller trials per step (already established:
    >=2 iterations for a cold start, see the rejected-trials test above)
    and inspects the stochastic_identity block sharp_front.run_2d now
    writes into the controller audit JSON.

    The `legacy_scalar` front engine used here has no RNG, threshold, or
    event-length-factor state at all (that only exists on the persistent-
    site production engine chain, unreachable without the real PF kernel
    -- see EXTERNAL BLOCKER), so those specific fields are expected to
    report "not_materialized" consistently. B and N_em DO exist on this
    engine and are the real, exercised part of this proof: they must be
    bit-identical across every candidate trial within a step.
    """
    out_dir = _run(tmp_path, steps=3)
    audit = _load_audit(out_dir)

    for record in audit["records"]:
        identity = record["stochastic_identity"]
        assert identity["all_pretrial_fingerprints_equal"] is True

        interval_start = identity["interval_start"]
        before_commit = identity["before_hazard_commit"]
        per_candidate = identity["per_candidate"]
        # Normally one fingerprint per solve_to_target-internal trial
        # (== iterations), plus one extra for the cold-start seed solve
        # at step 1 (Uapp_saved<=0), which happens outside solve_to_target
        # and so isn't counted in its own `iterations`.
        assert len(per_candidate) in (record["iterations"], record["iterations"] + 1)

        all_fingerprints = [interval_start] + per_candidate + [before_commit]
        for fingerprint in all_fingerprints:
            # This engine has none of the persistent-site stochastic
            # attributes -- confirm the helper honestly reports absence
            # rather than fabricating a value.
            assert fingerprint["hazard_rng_digest"] == "not_materialized"
            assert fingerprint["cleave_threshold_action"] == "not_materialized"
            assert fingerprint["stochastic_event_length_factor"] == "not_materialized"
            assert fingerprint["mpz"]["digest"] == "not_materialized"

        # B and N_em DO exist on this engine and must be identical across
        # every trial within the step -- the real, exercised assertion.
        b_values = {fp["B"] for fp in all_fingerprints}
        n_em_values = {fp["N_em"] for fp in all_fingerprints}
        assert len(b_values) == 1, f"B changed during trial search: {b_values}"
        assert len(n_em_values) == 1, f"N_em changed during trial search: {n_em_values}"


def test_stochastic_identity_violation_raises_immediately(monkeypatch, tmp_path: Path):
    """Fail-closed check: if a candidate trial's fingerprint were to
    differ from the interval-start fingerprint, run_2d must raise rather
    than silently continue. Simulated by monkeypatching the fingerprint
    capture to return a mutated B on the second call, mimicking a trial
    that illegitimately advanced hazard state.
    """
    from arrhenius_fracture import sharp_front as sf
    from arrhenius_fracture.pf_theta0_stochastic_fingerprint_v10051840 import (
        StochasticFingerprint,
        MPZFingerprint,
        capture_stochastic_fingerprint,
        NOT_MATERIALIZED,
    )

    calls = {"count": 0}
    real_capture = capture_stochastic_fingerprint

    def _tampering_capture(eng):
        calls["count"] += 1
        fingerprint = real_capture(eng)
        if calls["count"] == 2:
            return StochasticFingerprint(
                hazard_rng_digest=fingerprint.hazard_rng_digest,
                emission_rng_digest=fingerprint.emission_rng_digest,
                cleave_threshold_action=fingerprint.cleave_threshold_action,
                stochastic_event_length_factor=fingerprint.stochastic_event_length_factor,
                B=999.0,  # illegitimately different
                N_em=fingerprint.N_em,
                mpz=fingerprint.mpz,
            )
        return fingerprint

    monkeypatch.setattr(
        "arrhenius_fracture.pf_theta0_stochastic_fingerprint_v10051840.capture_stochastic_fingerprint",
        _tampering_capture,
    )
    # sharp_front imports the function by name into its own local scope
    # inside run_2d (a fresh `from ... import ...` each call), so patching
    # the source module's attribute above is sufficient -- no separate
    # sharp_front-level patch target exists to monkeypatch.

    with pytest.raises(RuntimeError, match="stochastic identity"):
        _run(tmp_path, steps=1)
