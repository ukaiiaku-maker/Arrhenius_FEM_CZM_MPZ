from pathlib import Path
from types import SimpleNamespace
import math

import pytest

import run_v10_0_5_4_vhcf_delta_sigma as campaign
from arrhenius_fracture import (
    mode_i_first_passage_v10_0_5_3_fatigue as v10053,
)
from arrhenius_fracture.mode_i_first_passage_v10_0_5_4_vhcf import (
    classify_termination_v10054,
)


def _row(cycles, extension=0.0, n_fire=0):
    return {
        "fatigue_cycles": str(cycles),
        "crack_extension_m": str(extension),
        "n_fire": str(n_fire),
        "B": "0.2",
        "mpz_available_site_fraction": "0.5",
    }


def test_v10054_classifies_max_block_exhaustion_as_right_censored():
    rows = [_row(10.0) for _ in range(5)]
    result = classify_termination_v10054(
        rows,
        cycles_max=1.0e6,
        max_blocks=5,
        target_extension_um=5.0,
    )
    assert result["status"] == "right_censored"
    assert result["termination"] == "right_censored_max_blocks"
    assert result["reached_cycle_horizon"] is False


def test_v10054_classifies_physical_cycle_horizon_as_complete():
    rows = [_row(4.0e13), _row(6.0e13)]
    result = classify_termination_v10054(
        rows,
        cycles_max=1.0e14,
        max_blocks=500,
        target_extension_um=5.0,
    )
    assert result["status"] == "complete"
    assert result["termination"] == "cycle_horizon"
    assert result["cycle_horizon_fraction"] == pytest.approx(1.0)


def test_v10054_campaign_defaults_are_vhcf_safe():
    args = campaign.build_parser().parse_args(
        ["--out", "runs/test", "--delta-sigma-MPa", "300"]
    )
    assert args.cycles_max == pytest.approx(1.0e14)
    assert math.isinf(args.max_block_cycles)
    assert args.resolve_cyclic_mechanics is False

    command = campaign._base_command(args, Path("runs/case"), 700.0, 1.0e-7)
    assert campaign.ENTRY_MODULE in command
    assert "--target-da-per-block-um" not in command
    assert "--crystal-compete" in command
    assert "--no-cyclic-mechanics" in command
    assert command[command.index("--max-block-cycles") + 1] == "inf"


def test_v10053_authoritative_predictor_uses_tensor_weights(monkeypatch):
    seen = {}

    class Front:
        def predict_fatigue_cycle(
            self, waveform, temperature, n_phase, *, system_weights
        ):
            seen["weights"] = list(system_weights)
            seen["n_phase"] = n_phase
            return {
                "dN_emit_per_cycle": 0.2,
                "dN_store_per_cycle": 0.03,
                "dN_mobile_per_cycle": 0.17,
                "dN_escape_per_cycle": 0.04,
                "dN_peierls_per_cycle": 0.05,
                "dN_taylor_per_cycle": 0.06,
                "mu_cleave_per_cycle": 1.0e-8,
                "avg_sigma_tip": 1.0,
                "max_sigma_tip": 2.0,
                "avg_sigma_emit_eff": 3.0,
            }

    def forbidden_fallback(*args, **kwargs):
        raise AssertionError("legacy scalar predictor was called")

    monkeypatch.setattr(
        v10053,
        "latest_tensor_drive",
        lambda: {"slip_system_drive_factors": [0.25, 1.0]},
    )
    wrapped = v10053._fatigue_predictor_dispatch(forbidden_fallback)
    controller = SimpleNamespace(cfg=SimpleNamespace(n_phase=96))
    result = wrapped(controller, Front(), object(), 700.0)

    assert seen["weights"] == [0.25, 1.0]
    assert seen["n_phase"] == 96
    assert result.mu_emit == pytest.approx(0.2)
    assert result.store_per_cycle == pytest.approx(0.03)
    assert result.mobile_per_cycle == pytest.approx(0.17)
    assert result.mu_cleave == pytest.approx(1.0e-8)
