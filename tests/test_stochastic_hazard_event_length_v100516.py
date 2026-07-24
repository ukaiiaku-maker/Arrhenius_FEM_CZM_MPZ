from __future__ import annotations

import copy
from types import MethodType, SimpleNamespace

import numpy as np
import pytest

from arrhenius_fracture.persistent_site_registry_v100514 import (
    select_persistent_site_row,
)
from arrhenius_fracture.persistent_site_signed_mpz_v100514 import (
    SignedShieldingKernelV100514,
)
from arrhenius_fracture.persistent_site_stochastic_tip_v100516 import (
    PersistentSitePFStochasticMovingTipFrontEngineV100516,
    clipped_exponential_mean,
    draw_hazard_threshold,
    threshold_event_length_factor,
)
from arrhenius_fracture.sharp_front import (
    FrontConfig,
    default_cleavage_barrier,
    default_emission_barrier,
)
from arrhenius_fracture import stochastic_geometry_v100516 as geometry


def make_engine(seed=1720):
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
        source_path="synthetic_stochastic_kernel",
    )
    cls = PersistentSitePFStochasticMovingTipFrontEngineV100516
    cls.configure(candidate, kernel)
    cls.configure_stochastic(
        hazard_mode="exponential",
        hazard_seed=seed,
        event_length_mode="threshold_scaled",
        event_minimum_factor=0.5,
        event_maximum_factor=4.0,
    )
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
    engine = cls(
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


def coupled(engine, dt):
    return engine._integrate_coupled(
        K_cleave=20.0e6,
        K_emit=20.0e6,
        T_K=700.0,
        dt_s=dt,
        drive_factors=np.array([0.2, 0.1]),
        tau_signed_Pa=np.array([1.0e8, -1.0e8]),
    )


def test_exponential_threshold_draw_is_reproducible_and_unit_mean():
    rng1 = np.random.default_rng(31415)
    rng2 = np.random.default_rng(31415)
    first = [draw_hazard_threshold("exponential", rng1) for _ in range(100)]
    second = [draw_hazard_threshold("exponential", rng2) for _ in range(100)]
    assert first == second
    rng = np.random.default_rng(2718)
    sample = np.array([draw_hazard_threshold("exponential", rng) for _ in range(200000)])
    assert float(np.mean(sample)) == pytest.approx(1.0, abs=0.01)
    assert draw_hazard_threshold("deterministic", rng) == 1.0


def test_threshold_correlated_event_length_is_bounded_and_mean_preserving():
    mean_clip = clipped_exponential_mean(0.5, 4.0)
    assert mean_clip == pytest.approx(0.5 + np.exp(-0.5) - np.exp(-4.0))
    assert threshold_event_length_factor(
        0.0, minimum_factor=0.5, maximum_factor=4.0
    ) == pytest.approx(0.5 / mean_clip)
    assert threshold_event_length_factor(
        10.0, minimum_factor=0.5, maximum_factor=4.0
    ) == pytest.approx(4.0 / mean_clip)
    rng = np.random.default_rng(97531)
    xi = rng.exponential(1.0, size=300000)
    factors = np.array(
        [
            threshold_event_length_factor(
                value, minimum_factor=0.5, maximum_factor=4.0
            )
            for value in xi
        ]
    )
    assert float(np.mean(factors)) == pytest.approx(1.0, abs=0.005)


def test_same_threshold_controls_waiting_time_and_continuous_event_length():
    engine = make_engine()
    install_no_plastic_constant_cleavage(engine, 1.0)
    engine.hazard_threshold_action = 2.0
    engine._set_current_event_length()
    expected_length = engine.f.da * threshold_event_length_factor(2.0)

    first = coupled(engine, 1.0)
    assert first["fired"] is False
    assert first["dB"] == pytest.approx(0.5)
    assert first["da"] == pytest.approx(0.5 * expected_length)
    assert engine.mpz_state.advance_total_m == pytest.approx(0.5 * expected_length)

    second = coupled(engine, 1.0)
    assert second["fired"] is True
    assert second["hazard_threshold_completed_action"] == pytest.approx(2.0)
    assert second["geometry_event_advance_m"] == pytest.approx(expected_length)
    assert engine.mpz_state.advance_total_m == pytest.approx(expected_length)
    assert engine.checkpoint_advance_total_m == pytest.approx(expected_length)
    assert engine.n_adv == 1


def test_stochastic_predictor_preserves_rng_threshold_and_state():
    engine = make_engine(seed=44)
    install_no_plastic_constant_cleavage(engine, 0.2)
    before_rng = copy.deepcopy(engine._hazard_rng.bit_generator.state)
    before_threshold = engine.hazard_threshold_action
    before_B = engine.B
    before_advance = engine.mpz_state.advance_total_m
    predicted = engine.predict_clock_increment_drives(
        20.0e6, 20.0e6, 700.0, 1.0
    )
    assert predicted >= 0.0
    assert engine._hazard_rng.bit_generator.state == before_rng
    assert engine.hazard_threshold_action == before_threshold
    assert engine.B == before_B
    assert engine.mpz_state.advance_total_m == before_advance


def test_geometry_adapter_replaces_only_requested_event_length():
    class Result:
        inserted = True
        moved = 7.5e-6
        reason = "ok"
        angle_error_deg = 0.0

    class Delegate:
        name = "adaptive_czm"

        def __init__(self):
            self.advance_log = []
            self.last = None

        def advance(self, **kwargs):
            self.last = kwargs
            p1 = np.asarray(kwargs["p1"], dtype=float)
            self.advance_log.append({"x1": float(p1[0]), "y1": float(p1[1])})
            return Result()

    geometry.clear_stochastic_geometry_state()
    geometry.enqueue_stochastic_geometry_event(
        {
            "event_advance_m": 7.5e-6,
            "event_length_factor": 1.5,
            "threshold_action": 1.2,
        }
    )
    delegate = Delegate()
    adapter = geometry._StochasticLengthBackend(delegate)
    adapter.advance(
        p0=np.array([0.5e-3, 0.0]),
        p1=np.array([0.505e-3, 0.0]),
        direction=np.array([1.0, 0.0]),
        front_id=0,
    )
    assert delegate.last is not None
    assert delegate.last["p1"][0] == pytest.approx(0.5075e-3)
    records = geometry.realized_stochastic_geometry_events()
    assert len(records) == 1
    assert records[0]["requested_fixed_length_m"] == pytest.approx(5.0e-6)
    assert records[0]["requested_stochastic_length_m"] == pytest.approx(7.5e-6)
