from arrhenius_fracture import kinetic_progressive_2d_v1002 as v1002
from arrhenius_fracture import sharp_front
from arrhenius_fracture.kinetic_progressive_2d_v1003_source import (
    build_progressive_run_2d_v1003_source,
)
from arrhenius_fracture.tensor_resolved_coupling_v1005 import (
    TensorResolvedKineticCohesiveStepper,
)


def test_generated_progressive_function_captures_v1005_stepper(monkeypatch):
    monkeypatch.setattr(
        v1002,
        "KineticCohesiveStepper",
        TensorResolvedKineticCohesiveStepper,
    )
    transformed = build_progressive_run_2d_v1003_source(sharp_front.run_2d)
    assert (
        transformed.__globals__["KineticCohesiveStepper"]
        is TensorResolvedKineticCohesiveStepper
    )
