from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np

from arrhenius_fracture.four_class_parameter_bridge_v100518 import (
    load_four_class_parameter_option,
)
from arrhenius_fracture.persistent_site_state_coupled_stochastic_emission_v1005182 import (
    PersistentSiteStateCoupledStochasticEmissionFrontEngineV1005182,
)
from arrhenius_fracture.sharp_front import (
    FrontConfig,
    default_cleavage_barrier,
    default_emission_barrier,
)
from arrhenius_fracture.signed_kernel_family_v1005141 import (
    load_signed_shielding_artifact_v1005141,
)

ROOT = Path(__file__).resolve().parents[1]
PARAMETER_ROOT = ROOT / "runtime_inputs" / "v10_0_5_18_four_class"
FAMILY = (
    ROOT
    / "runtime_inputs"
    / "v10_0_5_17_frozen_pf_inputs"
    / "signed_kernel"
    / "v10_2_14_active_only_campaign_family.json"
)


def test_trial_fraction_subdivision_reaches_exact_predictor(monkeypatch):
    monkeypatch.setenv("EMISSION_HAZARD_SEED", "444111166")
    monkeypatch.setenv("EMISSION_ADAPT_REFERENCE_TARGET", "0.05")
    monkeypatch.setenv("EMISSION_ADAPT_MAX_LOG_HAZARD_CHANGE", "0.25")
    monkeypatch.setenv("EMISSION_ADAPT_MAX_BACKSTRESS_FRACTION", "0.05")
    monkeypatch.setenv("EMISSION_ADAPT_MAX_GEOMETRY_FRACTION", "0.05")

    candidate, _ = load_four_class_parameter_option(
        PARAMETER_ROOT,
        "v913_paper_weakT01_0129902_persistent_sites",
    )
    family = load_signed_shielding_artifact_v1005141(FAMILY)
    cls = PersistentSiteStateCoupledStochasticEmissionFrontEngineV1005182
    cls.configure(candidate, family)
    cls.configure_stochastic(
        hazard_mode="exponential",
        hazard_seed=444111166,
        hazard_minimum_threshold=1.0e-12,
        event_length_mode="threshold_scaled",
        event_minimum_factor=0.5,
        event_maximum_factor=4.0,
    )
    front = FrontConfig()
    front.r0 = 1.0e-6
    front.L_pz = 50.0e-6
    front.da = 5.0e-6
    front.sigma_cap = 30.0e9
    engine = cls(
        front,
        default_cleavage_barrier(),
        default_emission_barrier(2.74e-10),
        160.0e9,
        0.28,
        2.74e-10,
        SimpleNamespace(
            blunting_length_m=0.5e-6,
            max_transport_cfl=0.35,
            max_transport_substeps=2000,
        ),
    )
    engine._mm = SimpleNamespace(
        latest={
            "two_channel_drive_reliable": True,
            "two_channel_drive_factors": [1.0, 0.5],
            "two_channel_tau_signed_Pa": [1.0e8, -1.0e8],
            "two_channel_names": ["positive", "negative"],
        }
    )

    before = engine.mpz_state.state_dict()
    fraction = 1.0
    accepted_prediction = None
    accepted_audit = None
    for _ in range(40):
        predicted = engine.predict_clock_increment_drives(
            80.0e6 * fraction,
            80.0e6 * fraction,
            300.0,
            840.0 * fraction,
        )
        assert np.isfinite(predicted)
        if predicted <= 0.05:
            accepted_prediction = predicted
            accepted_audit = dict(engine.emission_last_adaptive_prediction)
            break
        fraction *= 0.5

    assert accepted_prediction is not None
    assert fraction > 1.0e-8
    assert accepted_audit is not None
    assert accepted_audit["adaptive_increment"] <= 0.05
    assert accepted_audit["backstress_fraction"] <= 0.05 * (1.0 + 1.0e-10)
    assert accepted_audit["geometry_fraction"] <= 0.05 * (1.0 + 1.0e-10)
    assert accepted_audit["load_log_hazard_change"] <= 0.25 * (1.0 + 1.0e-10)
    assert engine.mpz_state.state_dict() == before
    assert engine.emission_event_count_total == 0
