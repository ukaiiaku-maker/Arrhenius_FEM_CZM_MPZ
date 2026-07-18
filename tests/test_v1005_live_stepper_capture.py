from types import SimpleNamespace

from arrhenius_fracture import kinetic_progressive_2d_v1002 as v1002
from arrhenius_fracture import sharp_front
from arrhenius_fracture.kinetic_progressive_2d_v1003_source import (
    build_progressive_run_2d_v1003_source,
)
from arrhenius_fracture.mode_i_first_passage_v10_0_5_parallel import (
    make_progressive_formatter_v1005,
)
from arrhenius_fracture.tensor_resolved_coupling_v1005 import (
    TensorResolvedKineticCohesiveStepper,
)


def _base_formatter(*_args, **_kwargs):
    return {"base": True}


def test_generated_progressive_function_captures_v1005_stepper_and_formatter(
    monkeypatch,
):
    formatter = make_progressive_formatter_v1005(_base_formatter)
    monkeypatch.setattr(
        v1002,
        "KineticCohesiveStepper",
        TensorResolvedKineticCohesiveStepper,
    )
    monkeypatch.setattr(v1002, "_v10_format_progressive_info", formatter)
    transformed = build_progressive_run_2d_v1003_source(sharp_front.run_2d)
    assert (
        transformed.__globals__["KineticCohesiveStepper"]
        is TensorResolvedKineticCohesiveStepper
    )
    assert transformed.__globals__["_v10_format_progressive_info"] is formatter


def test_v1005_formatter_persists_per_system_parallel_channels():
    formatter = make_progressive_formatter_v1005(_base_formatter)
    mechanics = {
        "tensor_resolved_drive_active": True,
        "opening_probe_sigma1_Pa": 3.0,
        "opening_probe_sigma_nn_Pa": 2.0,
        "opening_shape_factor": 2.0 / 3.0,
        "slip_system_tau_signed_Pa": [1.0, -2.0],
        "slip_system_drive_factors": [0.4, 0.8],
    }
    kinetics = {
        "channels": {
            "slip_system_drive_factors": [0.4, 0.8],
            "sigma_emission_effective_per_system_Pa": [10.0, 20.0],
            "sigma_emission_backstress_per_system_Pa": [1.0, 2.0],
        },
        "plastic": {
            "lambda_emit_per_system_s-1": [4.0, 5.0],
            "dN_emit_per_system": [0.1, 0.2],
        },
    }
    result = SimpleNamespace(
        kinetics=kinetics,
        mechanics_corrector=mechanics,
        mechanics_predictor={},
    )
    out = formatter(None, result, 0.0, 0.0)
    assert out["base"] is True
    assert out["tensor_resolved_drive_active"] is True
    assert out["slip_drive_factor_0"] == 0.4
    assert out["slip_drive_factor_1"] == 0.8
    assert out["slip_tau_signed_Pa_1"] == -2.0
    assert out["sigma_emission_effective_Pa_1"] == 20.0
    assert out["lambda_emit_s-1_0"] == 4.0
    assert out["dN_emit_1"] == 0.2
    assert out["directional_multiplier_applied_after_hazard"] is False
    assert out["fit_derived_shielding_cap_active"] is False
