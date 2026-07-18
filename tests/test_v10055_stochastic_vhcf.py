from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from arrhenius_fracture.mode_i_first_passage_v10_0_5_5_stochastic_vhcf import (
    patch_run_2d_source_v10055,
)
from arrhenius_fracture.mode_i_first_passage_v10_0_5_5_stochastic_vhcf_audited import (
    patch_run_2d_source_v10055_audited,
    validate_source_transform_v10055,
)
from arrhenius_fracture.stochastic_campaign_v10055 import (
    HybridSchedulerConfigV10055,
    StochasticCampaignKineticMPZStateV10055,
    StochasticPredictorMeanFieldMixinV10055,
    hybrid_choose_block_factory_v10055,
)
from arrhenius_fracture import sharp_front


def _pred(mu_emit, escape=0.0):
    return SimpleNamespace(mu_emit=mu_emit, escape_per_cycle=escape)


def _base_choice(controller, pred, user_block_cycles=None):
    return {
        "cycles": 100.0,
        "limiter": "cleavage_clock",
        "unlimited_cycles": 100.0,
        "candidate_limits": {"cleavage_clock": 100.0},
    }


def test_hybrid_scheduler_quiet_rare_and_tau_modes():
    cfg = HybridSchedulerConfigV10055(
        rare_event_target=0.25,
        tau_leap_target=3.0,
        tau_switch_expected_events=10.0,
    )
    audit = {}
    choose = hybrid_choose_block_factory_v10055(_base_choice, cfg, audit)
    controller = SimpleNamespace(cfg=SimpleNamespace(min_block_cycles=1.0e-6))

    quiet = choose(controller, _pred(1.0e-3))
    assert quiet["cycles"] == pytest.approx(100.0)
    assert quiet["v10055_stochastic_mode"] == "quiet"

    rare = choose(controller, _pred(1.0e-2))
    assert rare["cycles"] == pytest.approx(25.0)
    assert rare["limiter"] == "stochastic_rare_event"
    assert rare["v10055_expected_state_events"] == pytest.approx(0.25)

    tau = choose(controller, _pred(1.0))
    assert tau["cycles"] == pytest.approx(3.0)
    assert tau["limiter"] == "stochastic_tau_leap"
    assert tau["v10055_expected_state_events"] == pytest.approx(3.0)
    assert audit["quiet_calls"] == 1
    assert audit["rare_event_calls"] == 1
    assert audit["tau_leap_calls"] == 1


def _minimal_stochastic_state(seed):
    state = object.__new__(StochasticCampaignKineticMPZStateV10055)
    state.event_statistics = "stochastic"
    state.stochastic_emission = True
    state._v10055_force_mean_field = False
    state.stochastic_seed = seed
    state.stochastic_stream = 17011
    state._emission_rng = np.random.default_rng(seed)
    state.stochastic_emission_events = 0
    state.n_systems = 2
    state.n_bins = 4
    state.available_sites = np.array([10.0, 10.0])
    state.site_capacity = state.available_sites.copy()
    state.mobile = np.zeros((2, 4))
    state.accumulated_slip = np.zeros((2, 4))
    state.cfg = SimpleNamespace(source_bin_count=1)
    state.manifest = SimpleNamespace(
        emission=SimpleNamespace(
            rate=lambda sigma, temperature: np.full_like(
                np.asarray(sigma, dtype=float), 0.4
            )
        )
    )
    state.emitted_total = 0.0
    state.cumulative_emitted = 0.0
    state.last_emitted_per_system = np.zeros(2)
    state.last_emission_rate_per_system_s = np.zeros(2)
    state.last_sigma_emit_per_system_Pa = np.zeros(2)
    state.last_local_density_per_system_m2 = np.zeros(2)
    state.taylor_backstress_Pa = lambda: np.zeros(2)
    return state


def test_stochastic_source_commit_is_seed_reproducible_and_bounded():
    a = _minimal_stochastic_state(17)
    b = _minimal_stochastic_state(17)
    out_a = a.emit_exact(1.0, 1.0e9, 700.0, np.ones(2))
    out_b = b.emit_exact(1.0, 1.0e9, 700.0, np.ones(2))

    assert out_a["dN_emit_per_system"] == out_b["dN_emit_per_system"]
    assert 0.0 <= out_a["dN_emit"] <= 20.0
    assert np.all(a.available_sites >= 0.0)
    assert np.all(a.available_sites <= a.site_capacity)
    assert out_a["stochastic_source_commit"] == 1.0


def test_mean_field_predictor_flag_is_transactional():
    class Base:
        def __init__(self):
            self.mpz_state = SimpleNamespace(_v10055_force_mean_field=False)

        def predict_fatigue_cycle(self, *args, **kwargs):
            assert self.mpz_state._v10055_force_mean_field is True
            return {"dN_emit_per_cycle": 0.1}

    class Engine(StochasticPredictorMeanFieldMixinV10055, Base):
        pass

    engine = Engine()
    result = engine.predict_fatigue_cycle(object(), 700.0, 96)
    assert result["v10055_predictor_mean_field"] is True
    assert engine.mpz_state._v10055_force_mean_field is False
    assert engine._v10055_predictor_mean_field_calls == 1


def test_v10055_source_transform_compiles_current_run_2d():
    original = __import__("inspect").getsource(sharp_front.run_2d)
    source = patch_run_2d_source_v10055(original)
    compile(source, "<v10055_run_2d>", "exec")
    assert "vhcf_fem_cache_v10_0_5_5.json" in source
    assert "reuse_mechanics_v10055" in source

    audited = patch_run_2d_source_v10055_audited(original)
    compile(audited, "<v10055_audited_run_2d>", "exec")
    assert "cohesive_elements" in audited
    assert "clock_sum" in audited

    result = validate_source_transform_v10055()
    assert result["v10055_source_transform_preflight_passed"] is True
    assert result["stochastic_vhcf_adapter"] is True
    assert result["fem_cache_adapter"] is True
    assert result["cohesive_element_cache_signature"] is True


def test_shell_uses_flexible_horizon_and_opt_in_cache():
    text = Path("run_v10_0_5_5_stochastic_vhcf_delta_sigma.sh").read_text()
    assert 'CYCLES_MAX="${CYCLES_MAX:-1e12}"' in text
    assert "1e14 remains supported" in text
    assert 'VHCF_FEM_CACHE="${VHCF_FEM_CACHE:-0}"' in text
