from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from arrhenius_fracture.kinetic_campaign_czm_v10052 import (
    ChannelDiagnosticParallelOpeningEmissionCZMFrontEngine,
)
from arrhenius_fracture.mode_i_first_passage_v10_0_5_2_parallel import (
    make_progressive_formatter_v10052,
)


def test_strang_aggregation_sums_counts_and_keeps_final_rates():
    target = {}
    aggregate = ChannelDiagnosticParallelOpeningEmissionCZMFrontEngine._sum_numeric

    aggregate(target, {
        "dN_emit": 0.3,
        "dN_emit_per_system": [0.1, 0.2],
        "lambda_emit_per_system_s-1": [1.0, 2.0],
        "sigma_emission_effective_per_system_Pa": [10.0, 20.0],
    })
    aggregate(target, {
        "dN_emit": 0.7,
        "dN_emit_per_system": [0.4, 0.3],
        "lambda_emit_per_system_s-1": [3.0, 4.0],
        "sigma_emission_effective_per_system_Pa": [30.0, 40.0],
    })

    assert target["dN_emit"] == pytest.approx(1.0)
    assert np.allclose(target["dN_emit_per_system"], [0.5, 0.5])
    assert np.allclose(target["lambda_emit_per_system_s-1"], [3.0, 4.0])
    assert np.allclose(
        target["sigma_emission_effective_per_system_Pa"], [30.0, 40.0]
    )
    assert sum(target["dN_emit_per_system"]) == pytest.approx(target["dN_emit"])


def _step_result(dN_per_system=(0.1, 0.2)):
    return SimpleNamespace(
        kinetics={
            "channels": {
                "slip_system_drive_factors": [0.4, 0.1],
                "sigma_emission_effective_per_system_Pa": [4.0e8, 1.0e8],
                "sigma_emission_backstress_per_system_Pa": [2.0e7, 2.0e7],
            },
            "plastic": {
                "lambda_emit_per_system_s-1": [3.0e6, 0.25],
                "dN_emit_per_system": list(dN_per_system),
            },
        },
        mechanics_predictor={
            "tensor_resolved_drive_active": True,
            "slip_system_tau_signed_Pa": [8.0e9, -1.0e8],
            "opening_probe_sigma1_Pa": 1.0e10,
            "opening_probe_sigma_nn_Pa": 9.0e9,
            "opening_shape_factor": 0.9,
        },
        mechanics_corrector=None,
    )


def test_formatter_persists_complete_channel_partition():
    def original_formatter(_engine, _step_result, _KJ, _N_em_pre):
        return {"dN_emit_block": 0.3}

    formatter = make_progressive_formatter_v10052(original_formatter)
    out = formatter(None, _step_result(), 20.0e6, 0.0)

    assert out["lambda_emit_s-1_0"] == pytest.approx(3.0e6)
    assert out["lambda_emit_s-1_1"] == pytest.approx(0.25)
    assert out["dN_emit_0"] == pytest.approx(0.1)
    assert out["dN_emit_1"] == pytest.approx(0.2)
    assert out["per_channel_emission_partition_verified"] is True
    assert out["per_channel_emission_partition_residual"] == pytest.approx(0.0)
    assert out["per_channel_strang_diagnostics_complete"] is True


def test_formatter_fails_closed_on_partition_mismatch():
    def original_formatter(_engine, _step_result, _KJ, _N_em_pre):
        return {"dN_emit_block": 0.5}

    formatter = make_progressive_formatter_v10052(original_formatter)
    with pytest.raises(RuntimeError, match="partition mismatch"):
        formatter(None, _step_result(), 20.0e6, 0.0)
