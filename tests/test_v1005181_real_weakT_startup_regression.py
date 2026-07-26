from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np

from arrhenius_fracture.four_class_parameter_bridge_v100518 import (
    load_four_class_parameter_option,
)
from arrhenius_fracture.persistent_site_stochastic_onset_burst_v1005181 import (
    PersistentSiteStochasticOnsetBurstMovingTipFrontEngineV1005181,
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


def test_real_weakT_300K_startup_predictor_handles_stiff_emission(monkeypatch):
    monkeypatch.setenv("EMISSION_HAZARD_SEED", "444111166")
    candidate, _ = load_four_class_parameter_option(
        PARAMETER_ROOT,
        "v913_paper_weakT01_0129902_persistent_sites",
    )
    family = load_signed_shielding_artifact_v1005141(FAMILY)
    cls = PersistentSiteStochasticOnsetBurstMovingTipFrontEngineV1005181
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
    front.sigma_cap = 0.0
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
            "two_channel_drive_factors": [0.2, 0.1],
            "two_channel_tau_signed_Pa": [1.0e8, -1.0e8],
            "two_channel_names": ["positive", "negative"],
        }
    )

    before = engine.mpz_state.state_dict()
    predicted = engine.predict_clock_increment_drives(
        20.0e6,
        20.0e6,
        300.0,
        840.0,
    )

    assert np.isfinite(predicted)
    assert engine.mpz_state.state_dict() == before
    assert engine.emission_event_count_total == 0
