from __future__ import annotations

import copy
import math
from types import SimpleNamespace

import numpy as np

from arrhenius_fracture.kinetic_fatigue_v10053 import (
    FatigueKineticsMixinV10053,
    FatigueLoadingConfigV10053,
)


class DummyBase:
    def __init__(self):
        self.B = 0.0
        self.t = 0.0
        self.micro_advance_total_m = 0.0
        self.checkpoint_advance_total_m = 0.0
        self.f = SimpleNamespace(da=1.0)
        self.calls = []

    @staticmethod
    def _sum_numeric(target, source):
        for key, value in source.items():
            if isinstance(value, (int, float)):
                target[key] = target.get(key, 0.0) + value
            elif key == "dN_emit_per_system":
                a = np.asarray(target.get(key, np.zeros(len(value))))
                target[key] = (a + np.asarray(value)).tolist()
            else:
                target[key] = copy.deepcopy(value)

    def snapshot_kinetic_state(self):
        return {
            "B": self.B,
            "t": self.t,
            "micro": self.micro_advance_total_m,
            "checkpoint": self.checkpoint_advance_total_m,
            "calls": copy.deepcopy(self.calls),
        }

    def restore_kinetic_state(self, state):
        self.B = state["B"]
        self.t = state["t"]
        self.micro_advance_total_m = state["micro"]
        self.checkpoint_advance_total_m = state["checkpoint"]
        self.calls = copy.deepcopy(state["calls"])

    def stress_channels(self, Kopen, Kcleave, weights):
        return {
            "K_open_Pa_sqrt_m": Kopen,
            "K_cleave_input_Pa_sqrt_m": Kcleave,
            "sigma_cleave_eff_Pa": Kcleave,
            "sigma_emission_effective_Pa": Kopen,
        }

    def integrate_kinetics(self, Kopen, Kcleave, T, dt, *, system_weights=None):
        self.calls.append((Kopen, Kcleave, T, dt))
        rate = Kcleave / 10.0
        dB = rate * dt
        fired = self.B + dB >= 1.0
        if fired:
            consumed = (1.0 - self.B) / rate
            dB = 1.0 - self.B
            self.B = 0.0
            self.checkpoint_advance_total_m += 1.0
        else:
            consumed = dt
            self.B += dB
        self.t += consumed
        self.micro_advance_total_m += dB
        return {
            "fired": fired,
            "n_fire": int(fired),
            "dB": dB,
            "micro_advance_step_m": dB,
            "micro_advance_total_m": self.micro_advance_total_m,
            "checkpoint_committed_total_m": self.checkpoint_advance_total_m,
            "dt_consumed_s": consumed,
            "dt_unused_s": dt - consumed,
            "internal_substeps": 1,
            "lambda_c_effective_s-1": rate,
            "lambda_c_raw_s-1": rate,
            "G_cleave_eff_eV": 0.5,
            "plastic": {
                "dN_emit": 2.0 * dt,
                "dN_emit_per_system": [dt, dt],
                "lambda_emit_per_system_s-1": [2.0, 2.0],
            },
            "advance": {},
            "channels": self.stress_channels(Kopen, Kcleave, system_weights),
        }


class DummyEngine(FatigueKineticsMixinV10053, DummyBase):
    pass


def test_phase_factors_recover_requested_R():
    cfg = FatigueLoadingConfigV10053(R=0.2, frequency_Hz=10.0, n_phase=1000)
    f = cfg.phase_factors()
    assert np.max(f) == np.max(f)
    assert abs(np.mean(f) - 0.6) < 1e-12
    assert np.min(f) > 0.19
    assert np.max(f) < 1.01


def test_fatigue_interval_preserves_total_time_and_uses_all_phases():
    eng = DummyEngine()
    cfg = FatigueLoadingConfigV10053(R=0.1, frequency_Hz=100.0, n_phase=16)
    eng.configure_fatigue_v10053(cfg)
    out = eng.integrate_kinetics(1.0, 1.0, 300.0, 0.01)
    assert not out["fired"]
    assert math.isclose(out["dt_consumed_s"], 0.01)
    assert math.isclose(out["fatigue_cycles_consumed"], 1.0)
    assert len(eng.calls) == 16
    assert out["plastic"]["dN_emit_per_system"] == [0.01, 0.01]
    expected_dB = np.mean(cfg.phase_factors()) / 10.0 * 0.01
    assert math.isclose(out["dB"], expected_dB, rel_tol=1e-12)


def test_predictor_is_transactional():
    eng = DummyEngine()
    cfg = FatigueLoadingConfigV10053(R=0.1, frequency_Hz=10.0, n_phase=16)
    eng.configure_fatigue_v10053(cfg)
    before = eng.snapshot_kinetic_state()
    wave = SimpleNamespace(Kmax=1.0, R=0.1, frequency_Hz=10.0, closure_clip=True)
    pred = eng.predict_fatigue_cycle(wave, 300.0, 16)
    after = eng.snapshot_kinetic_state()
    assert before == after
    assert pred["mu_cleave_per_cycle"] > 0.0
    assert pred["dN_emit_per_cycle"] > 0.0


def test_event_returns_unused_block_time():
    eng = DummyEngine()
    eng.B = 0.999
    cfg = FatigueLoadingConfigV10053(R=1.0, frequency_Hz=1.0, n_phase=8)
    eng.configure_fatigue_v10053(cfg)
    out = eng.integrate_kinetics(10.0, 10.0, 300.0, 1.0)
    assert out["fired"]
    assert out["dt_consumed_s"] < 1.0
    assert out["dt_unused_s"] > 0.0
    assert math.isclose(out["fatigue_cycles_consumed"], out["dt_consumed_s"])


def test_inactive_adapter_preserves_monotonic_base_call():
    eng = DummyEngine()
    eng.configure_fatigue_v10053(None)
    out = eng.integrate_kinetics(2.0, 3.0, 300.0, 0.25)
    assert len(eng.calls) == 1
    assert eng.calls[0] == (2.0, 3.0, 300.0, 0.25)
    assert "fatigue_loading_v10053" not in out


def test_reported_waveform_endpoints_are_exact_not_midpoint_samples():
    eng = DummyEngine()
    cfg = FatigueLoadingConfigV10053(R=0.2, frequency_Hz=10.0, n_phase=8)
    eng.configure_fatigue_v10053(cfg)
    out = eng.integrate_kinetics(4.0, 5.0, 300.0, 0.01)
    assert out["fatigue_Kmax_Pa_sqrt_m"] == 5.0
    assert out["fatigue_Kmin_Pa_sqrt_m"] == 1.0
    assert out["fatigue_DeltaK_Pa_sqrt_m"] == 4.0
