from __future__ import annotations

import json
import os
from collections import Counter, defaultdict
from pathlib import Path
from types import SimpleNamespace

import arrhenius_fracture.persistent_site_joint_K_action_error_v10051832 as action32
from arrhenius_fracture.four_class_parameter_bridge_v100518 import (
    load_four_class_parameter_option,
)
from arrhenius_fracture.persistent_site_joint_K_action_error_v10051832 import (
    PersistentSiteJointKActionErrorFrontEngineV10051832,
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
LIMIT = int(os.environ.get("PEAK1000_DIAGNOSTIC_HAZARD_LIMIT", "5000"))
HORIZON_FACTOR = float(os.environ.get("EMISSION_EVENT_HORIZON_FACTOR", "1000"))
ACTION_ERROR = float(os.environ.get("EMISSION_INNER_MAX_ACTION_ERROR", "0.25"))
MAX_REFINEMENTS = int(os.environ.get("EMISSION_INNER_MAX_REFINEMENTS", "64"))

os.environ["EMISSION_HAZARD_SEED"] = str(SEED)
os.environ["EMISSION_RAMP_MAX_LOG_RATE_CHANGE"] = "0.25"
os.environ["EMISSION_RAMP_ACTION_FLOOR"] = "1e-6"
os.environ["EMISSION_INNER_MAX_LOG_HAZARD_CHANGE"] = "0.25"
os.environ["EMISSION_INNER_ACTION_FLOOR"] = "1e-6"
os.environ["EMISSION_EVENT_HORIZON_FACTOR"] = str(HORIZON_FACTOR)
os.environ["EMISSION_INNER_MAX_ACTION_ERROR"] = str(ACTION_ERROR)
os.environ["EMISSION_INNER_MAX_REFINEMENTS"] = str(MAX_REFINEMENTS)

candidate, _ = load_four_class_parameter_option(PARAMETER_ROOT, OPTION)
family = load_signed_shielding_artifact_v1005141(FAMILY)
cls = PersistentSiteJointKActionErrorFrontEngineV10051832
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
start_times = defaultdict(Counter)
durations = defaultdict(list)
original_hazards = action32._emission_hazards
original_transport = action32.advance_pf_transport_only_v10222
original_deposit = action32._deposit_one_emission_event


class DiagnosticStop(RuntimeError):
    pass


def counted_hazards(*args, **kwargs):
    counts["hazard_calls"] += 1
    if counts["hazard_calls"] > LIMIT:
        raise DiagnosticStop("peak-1000 diagnostic hazard-call limit reached")
    return original_hazards(*args, **kwargs)


def counted_transport(state, *, dt_s, **kwargs):
    event_key = str(counts["events"])
    counts["transport_calls"] += 1
    start_key = f"{float(state.time_s):.12e}"
    start_times[event_key][start_key] += 1
    durations[event_key].append(float(dt_s))
    return original_transport(state, dt_s=dt_s, **kwargs)


def counted_deposit(*args, **kwargs):
    result = original_deposit(*args, **kwargs)
    counts["events"] += 1
    return result


action32._emission_hazards = counted_hazards
action32.advance_pf_transport_only_v10222 = counted_transport
action32._deposit_one_emission_event = counted_deposit

completed = False
predicted = None
error = None
try:
    predicted = engine.predict_clock_increment_drives(16.5e6, 16.5e6, 1000.0, 840.0)
    completed = True
except DiagnosticStop as exc:
    error = str(exc)
finally:
    action32._emission_hazards = original_hazards
    action32.advance_pf_transport_only_v10222 = original_transport
    action32._deposit_one_emission_event = original_deposit

stage_summary = {}
for key in sorted(set(start_times) | set(durations), key=int):
    counter = start_times[key]
    values = durations[key]
    stage_summary[key] = {
        "transport_calls": len(values),
        "distinct_transport_start_times": len(counter),
        "maximum_calls_at_one_start_time": max(counter.values(), default=0),
        "minimum_dt_s": min(values, default=None),
        "maximum_dt_s": max(values, default=None),
        "median_dt_s": sorted(values)[len(values) // 2] if values else None,
    }

payload = {
    **counts,
    "completed_before_limit": completed,
    "predicted_clock_increment": predicted,
    "limit": LIMIT,
    "error": error,
    "action_error_tolerance": engine.inner_max_action_error,
    "event_horizon_factor": engine.inner_event_horizon_factor,
    "maximum_refinements": engine.inner_max_refinements,
    "stages": stage_summary,
}
print("PEAK1000_HORIZON_FACTOR_DIAGNOSTIC " + json.dumps(payload, sort_keys=True), flush=True)
