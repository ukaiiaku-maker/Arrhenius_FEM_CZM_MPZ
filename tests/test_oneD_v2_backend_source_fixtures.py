from pathlib import Path
import json
import pandas as pd
from reduced_fracture_v2.adapters import FEMCZMProcessZoneStateAdapter,PFProcessZoneStateAdapter
from reduced_fracture_v2.kinetics import SharedBarrierHazardCore
from reduced_fracture_v2.lifecycle import FEMCZMEventLifecyclePolicy,PFEventLifecyclePolicy
from reduced_fracture_v2.mechanics import FEMCZMMechanicsProvider,PFMechanicsProvider

OUT=Path(__file__).parents[1]/"analysis_outputs/oneD_v2_backend_qualification"
def fixture(name):return json.loads((OUT/name).read_text())
def test_required_architecture_names_are_distinct():
 assert PFProcessZoneStateAdapter is not FEMCZMProcessZoneStateAdapter
 assert PFEventLifecyclePolicy is not FEMCZMEventLifecyclePolicy
 assert PFMechanicsProvider is not FEMCZMMechanicsProvider
 assert SharedBarrierHazardCore.__name__=="SharedBarrierHazardCore"
def test_common_lifecycle_removed():
 import reduced_fracture_v2.lifecycle as x;assert not hasattr(x,"CommonEventLifecyclePolicy")
def test_pf_direct_source_fixture():
 x=fixture("oneD_v2_pf_source_fixture.json");assert x["rollback_pass"] and x["cleavage_crossing_fired"] and x["post_event_renewal_pass"]
def test_fem_direct_source_fixture():
 x=fixture("oneD_v2_femczm_source_fixture.json");assert x["rollback_pass"] and x["post_event_renewal_pass"] and x["cleavage_crossing_status"].startswith("DRIVER_OWNED")
def test_source_sufficiency_is_fail_closed():
 x=pd.read_csv(OUT/"oneD_v2_source_state_sufficiency_matrix.csv");assert (x.exact_comparison_status=="REQUIRES_INSTRUMENTED_DRIVER_DIAGNOSTIC").any()
def test_mechanics_maps_not_fabricated():
 assert pd.read_csv(OUT/"oneD_v2_pf_mechanics_map.csv").availability_status.iloc[0].startswith("BLOCKED")
 assert pd.read_csv(OUT/"oneD_v2_femczm_mechanics_map.csv").availability_status.iloc[0].startswith("BLOCKED")
def test_campaigns_remain_unlaunched():
 x=fixture("oneD_v2_backend_gate_decision.json");assert not x["campaigns_launched"] and not x["canonical_parameters_changed"]
