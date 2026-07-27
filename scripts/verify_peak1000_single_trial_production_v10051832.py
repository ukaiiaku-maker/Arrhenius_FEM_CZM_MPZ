from __future__ import annotations

import json
import os
from pathlib import Path
from types import SimpleNamespace

import numpy as np

from arrhenius_fracture.four_class_parameter_bridge_v100518 import (
    load_four_class_parameter_option,
)
from arrhenius_fracture.persistent_site_joint_K_single_trial_v10051832 import (
    PersistentSiteJointKSingleTrialFrontEngineV10051832,
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
OPTION = "v913_paper_peak01_0242980_persistent_sites"
SEED = 3100913530

os.environ["EMISSION_HAZARD_SEED"] = str(SEED)
os.environ["EMISSION_EVENT_HORIZON_FACTOR"] = "1.25"
os.environ["EMISSION_MAX_EVENTS_PER_HALF_STEP"] = "10000"

candidate, _ = load_four_class_parameter_option(PARAMETER_ROOT, OPTION)
family = load_signed_shielding_artifact_v1005141(FAMILY)
cls = PersistentSiteJointKSingleTrialFrontEngineV10051832
cls.configure(candidate, family)
cls.configure_stochastic(
    hazard_mode="exponential",
    hazard_seed=SEED,
    hazard_minimum_threshold=1.0e-12,
    event_length_mode="threshold_scaled",
    event_minimum_factor=0.5,
    event_maximum_factor=4.0,
)


def build_engine() -> PersistentSiteJointKSingleTrialFrontEngineV10051832:
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
            "cleavage_factor": 1.0,
            "emission_factor": 4.0,
        }
    )
    return engine


engine = build_engine()
before = engine._capture_state()
predicted = engine.predict_clock_increment_drives(
    16.5e6,
    16.5e6,
    1000.0,
    840.0,
)
after = engine._capture_state()

errors = []
if not np.isfinite(predicted):
    errors.append("predicted clock increment is not finite")
if not 0.0 <= predicted <= 1.0 + 1.0e-12:
    errors.append("predicted clock increment is outside [0, 1]")
if after["mpz_state"].state_dict() != before["mpz_state"].state_dict():
    errors.append("predictor mutated MPZ state")
if after["emission_event_count_total"] != before["emission_event_count_total"]:
    errors.append("predictor mutated accepted event count")
if after["hazard_rng_state"] != before["hazard_rng_state"]:
    errors.append("predictor mutated RNG state")

trial = build_engine()
trial._integrate_coupled(
    K_cleave=16.5e6,
    K_emit=16.5e6,
    T_K=1000.0,
    dt_s=840.0,
    drive_factors=np.asarray([1.0, 0.5]),
    tau_signed_Pa=np.asarray([1.0e8, -1.0e8]),
)

events = int(trial.emission_event_count_total)
if events <= 0:
    errors.append("trial produced no emission events")
if events > trial.emission_max_events_per_half_step:
    errors.append("trial exceeded the exact event guard")

payload = {
    "passed": not errors,
    "option": OPTION,
    "temperature_K": 1000,
    "K_end_MPa_sqrt_m": 16.5,
    "dt_s": 840.0,
    "seed": SEED,
    "predicted_clock_increment": float(predicted),
    "trial_exact_emission_events": events,
    "event_guard": int(trial.emission_max_events_per_half_step),
    "event_horizon_factor": float(trial.inner_event_horizon_factor),
    "state_dependent_log_rate_rejection_active": bool(
        trial.state_dependent_log_rate_rejection_active
    ),
    "PF_transport_operator_preserved": bool(
        trial.PF_transport_operator_preserved
    ),
    "errors": errors,
}
print(json.dumps(payload, indent=2, sort_keys=True))
if errors:
    raise SystemExit(1)
