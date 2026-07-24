from __future__ import annotations

from types import MethodType, SimpleNamespace

import numpy as np
import pytest

from arrhenius_fracture.persistent_site_moving_tip_v100515 import (
    PersistentSitePFMovingTipFrontEngineV100515,
)
from arrhenius_fracture.persistent_site_registry_v100514 import (
    select_persistent_site_row,
)
from arrhenius_fracture.persistent_site_signed_mpz_v100514 import (
    SignedShieldingKernelV100514,
)
from arrhenius_fracture.sharp_front import (
    FrontConfig,
    default_cleavage_barrier,
    default_emission_barrier,
)


def make_engine():
    candidate = select_persistent_site_row("v912_peak_0118_persistent_sites")
    n = candidate.n_bins_recommended
    kernel = SignedShieldingKernelV100514(
        active_kernel_Pa_sqrt_m_per_signed_line=np.zeros((2, n)),
        wake_kernel_Pa_sqrt_m_per_signed_line=np.zeros((2, n)),
        activation_to_line_content_by_system=np.ones(2),
        metadata={
            "candidate_independent": True,
            "counts_are_signed_burgers_lines": True,
            "normalization_is_mechanically_derived": True,
        },
        source_path="synthetic_moving_tip_kernel",
    )
    PersistentSitePFMovingTipFrontEngineV100515.configure(candidate, kernel)
    f = FrontConfig()
    f.r0 = 1.0e-6
    f.L_pz = 50.0e-6
    f.da = 5.0e-6
    f.sigma_cap = 0.0
    mpz_cfg = SimpleNamespace(
        blunting_length_m=0.5e-6,
        max_transport_cfl=0.35,
        max_transport_substeps=2000,
    )
    engine = PersistentSitePFMovingTipFrontEngineV100515(
        f,
        default_cleavage_barrier(),
        default_emission_barrier(2.74e-10),
        160.0e9,
        0.28,
        2.74e-10,
        mpz_cfg,
    )
    engine._mm = SimpleNamespace(
        latest={
            "two_channel_drive_reliable": True,
            "two_channel_drive_factors": [0.2, 0.1],
            "two_channel_tau_signed_Pa": [1.0e8, -1.0e8],
            "two_channel_names": ["positive", "negative"],
        }
    )
    engine.max_action_substep = 1.0
    engine.max_translation_substep_m = 1.0
    return engine


def install_no_plastic_constant_cleavage(engine, rate_s):
    def no_plastic(self, **kwargs):
        return {
            "dN_emit": 0.0,
            "dN_trapped": 0.0,
            "dN_released": 0.0,
            "dN_escaped": 0.0,
            "dN_recovered": 0.0,
        }, 0.0

    def cleavage(self, sigma, T):
        return float(rate_s), float(rate_s), 1.0e-19

    engine._plastic_half_step = MethodType(no_plastic, engine)
    engine.lambda_cleave = MethodType(cleavage, engine)


def coupled(engine, dt=1.0):
    return engine._integrate_coupled(
        K_cleave=20.0e6,
        K_emit=20.0e6,
        T_K=700.0,
        dt_s=dt,
        drive_factors=np.array([0.2, 0.1]),
        tau_signed_Pa=np.array([1.0e8, -1.0e8]),
    )


def test_fractional_cleavage_progress_translates_mpz_before_checkpoint():
    engine = make_engine()
    install_no_plastic_constant_cleavage(engine, 0.25)
    result = coupled(engine, dt=1.0)
    assert result["fired"] is False
    assert result["dB"] == pytest.approx(0.25)
    assert result["da"] == pytest.approx(0.25 * engine.f.da)
    assert engine.B == pytest.approx(0.25)
    assert engine.mpz_state.advance_total_m == pytest.approx(0.25 * engine.f.da)
    assert engine.a_adv == 0.0


def test_checkpoint_does_not_repeat_already_accumulated_mpz_advance():
    engine = make_engine()
    install_no_plastic_constant_cleavage(engine, 0.25)
    results = [coupled(engine, dt=1.0) for _ in range(4)]
    assert results[-1]["fired"] is True
    assert engine.n_adv == 1
    assert engine.a_adv == pytest.approx(engine.f.da)
    assert engine.checkpoint_advance_total_m == pytest.approx(engine.f.da)
    assert engine.micro_advance_total_m == pytest.approx(engine.f.da)
    assert engine.mpz_state.advance_total_m == pytest.approx(engine.f.da)
    assert engine.B == pytest.approx(0.0, abs=1e-12)


def test_fast_microstructure_development_reduces_same_step_cleavage_progress():
    engine = make_engine()

    def develop_retained(self, **kwargs):
        self.mpz_state.retained_positive[0, 0] += 1.0
        return {
            "dN_emit": 1.0,
            "dN_trapped": 0.0,
            "dN_released": 0.0,
            "dN_escaped": 0.0,
            "dN_recovered": 0.0,
        }, 0.0

    def state_dependent_cleavage(self, sigma, T):
        rate = 1.0 if self.mpz_state.retained_count == 0.0 else 0.1
        return rate, rate, 1.0e-19

    engine._plastic_half_step = MethodType(develop_retained, engine)
    engine.lambda_cleave = MethodType(state_dependent_cleavage, engine)
    result = coupled(engine, dt=1.0)
    assert result["dB"] == pytest.approx(0.1)
    assert result["dB"] < 1.0
    assert engine.mpz_state.advance_total_m == pytest.approx(0.1 * engine.f.da)


def test_coupled_predictor_is_nonmutating_and_includes_fractional_motion():
    engine = make_engine()
    install_no_plastic_constant_cleavage(engine, 0.25)
    before_B = engine.B
    before_advance = engine.mpz_state.advance_total_m
    prediction = engine.predict_clock_increment_drives(
        20.0e6, 20.0e6, 700.0, 1.0
    )
    assert prediction == pytest.approx(0.25)
    assert engine.B == before_B
    assert engine.mpz_state.advance_total_m == before_advance
    assert engine.kinetic_prediction_calls == 1


def test_geometry_veto_restores_renewal_origin_and_fails_closed():
    engine = make_engine()
    install_no_plastic_constant_cleavage(engine, 1.0)
    result = engine.step_drives(20.0e6, 20.0e6, 700.0, 1.0)
    assert result["fired"] is True
    assert engine.mpz_state.advance_total_m == pytest.approx(engine.f.da)
    with pytest.raises(RuntimeError, match="continuous moving-tip checkpoint"):
        engine.restore_geometry_veto(1)
    assert engine.mpz_state.advance_total_m == pytest.approx(0.0)
    assert engine.a_adv == pytest.approx(0.0)
    assert engine.n_adv == 0


def test_failed_constitutive_trial_rolls_back_complete_signed_state():
    engine = make_engine()
    before = engine.mpz_state.state_dict()

    def fail_after_mutation(self, **kwargs):
        self.mpz_state.mobile_positive[0, 0] += 5.0
        raise RuntimeError("synthetic constitutive failure")

    engine._plastic_half_step = MethodType(fail_after_mutation, engine)
    engine.lambda_cleave = MethodType(
        lambda self, sigma, T: (0.1, 0.1, 1.0e-19), engine
    )
    with pytest.raises(RuntimeError, match="synthetic constitutive failure"):
        engine.step_drives(20.0e6, 20.0e6, 700.0, 1.0)
    after = engine.mpz_state.state_dict()
    assert after["mobile_positive"] == before["mobile_positive"]
    assert engine.B == 0.0
    assert engine.mpz_state.advance_total_m == 0.0
