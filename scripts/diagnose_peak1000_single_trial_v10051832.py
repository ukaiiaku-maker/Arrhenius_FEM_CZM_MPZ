from __future__ import annotations

import json
import os
from pathlib import Path
from types import SimpleNamespace

import arrhenius_fracture.persistent_site_joint_K_ramp_emission_v10051832 as joint32
from arrhenius_fracture.four_class_parameter_bridge_v100518 import (
    load_four_class_parameter_option,
)
from arrhenius_fracture.persistent_site_joint_K_ramp_emission_v10051832 import (
    PersistentSiteJointKRampEventDrivenFrontEngineV10051832,
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
LIMIT = int(os.environ.get("PEAK1000_DIAGNOSTIC_HAZARD_LIMIT", "50000"))

os.environ["EMISSION_HAZARD_SEED"] = str(SEED)
os.environ["EMISSION_RAMP_MAX_LOG_RATE_CHANGE"] = "0.25"
os.environ["EMISSION_RAMP_ACTION_FLOOR"] = "1e-6"
os.environ["EMISSION_INNER_MAX_LOG_HAZARD_CHANGE"] = os.environ.get(
    "EMISSION_INNER_MAX_LOG_HAZARD_CHANGE", "1e300"
)
os.environ["EMISSION_INNER_ACTION_FLOOR"] = "1e-6"
os.environ["EMISSION_EVENT_HORIZON_FACTOR"] = os.environ.get(
    "EMISSION_EVENT_HORIZON_FACTOR", "1000"
)
os.environ["EMISSION_INNER_MAX_REFINEMENTS"] = "64"
os.environ["EMISSION_MAX_EVENTS_PER_HALF_STEP"] = os.environ.get(
    "EMISSION_MAX_EVENTS_PER_HALF_STEP", "10000"
)

candidate, _ = load_four_class_parameter_option(PARAMETER_ROOT, OPTION)
family = load_signed_shielding_artifact_v1005141(FAMILY)
cls = PersistentSiteJointKRampEventDrivenFrontEngineV10051832
cls.configure(candidate, family)
cls.configure_stochastic(
    hazard_mode="exponential",
    hazard_seed=SEED,
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
        "cleavage_factor": 1.0,
        "emission_factor": 4.0,
    }
)

counts = {"hazard_calls": 0, "transport_calls": 0, "events": 0}
original_hazards = joint32._emission_hazards
original_transport = joint32.advance_pf_transport_only_v10222
original_deposit = joint32._deposit_one_emission_event


class DiagnosticStop(RuntimeError):
    pass


def counted_hazards(*args, **kwargs):
    counts["hazard_calls"] += 1
    if counts["hazard_calls"] > LIMIT:
        raise DiagnosticStop("peak-1000 diagnostic hazard-call limit reached")
    return original_hazards(*args, **kwargs)


def counted_transport(*args, **kwargs):
    counts["transport_calls"] += 1
    return original_transport(*args, **kwargs)


def counted_deposit(*args, **kwargs):
    result = original_deposit(*args, **kwargs)
    counts["events"] += 1
    return result


joint32._emission_hazards = counted_hazards
joint32.advance_pf_transport_only_v10222 = counted_transport
joint32._deposit_one_emission_event = counted_deposit

completed = False
predicted = None
error = None
try:
    predicted = engine.predict_clock_increment_drives(16.5e6, 16.5e6, 1000.0, 840.0)
    completed = True
except (DiagnosticStop, RuntimeError) as exc:
    error = f"{type(exc).__name__}: {exc}"
finally:
    joint32._emission_hazards = original_hazards
    joint32.advance_pf_transport_only_v10222 = original_transport
    joint32._deposit_one_emission_event = original_deposit

payload = {
    **counts,
    "completed": completed,
    "predicted_clock_increment": predicted,
    "hazard_limit": LIMIT,
    "error": error,
    "event_horizon_factor": engine.inner_event_horizon_factor,
    "max_log_hazard_change": engine.inner_max_log_hazard_change,
    "max_events_per_half_step": engine.emission_max_events_per_half_step,
}
print("PEAK1000_SINGLE_TRIAL_DIAGNOSTIC " + json.dumps(payload, sort_keys=True), flush=True)
if not completed and error and "DiagnosticStop" not in error:
    raise SystemExit(1)
