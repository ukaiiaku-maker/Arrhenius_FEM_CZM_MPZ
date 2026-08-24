import csv,json
from pathlib import Path
import pytest
from reduced_fracture_v2.mechanics import FEMCZMMechanicsProvider,PFMechanicsProvider

ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/"analysis_outputs/oneD_v2_mechanics_maps_and_baselines"

def test_typed_pf_provider_has_native_quantities_but_no_continuum_g():
 p=PFMechanicsProvider.from_csv(OUT/"oneD_v2_pf_mechanics_map.csv","test")
 s=p.evaluate(100e-6,1e-5)
 assert s.pf_native_J_J_m2 and s.pf_native_KJ_Pa_sqrt_m
 assert s.structural_G_J_m2 is None and s.qualified_KG_Pa_sqrt_m is None

def test_typed_fem_provider_keeps_native_and_structural_quantities_distinct():
 p=FEMCZMMechanicsProvider.from_csv(OUT/"oneD_v2_fem_native_mechanics_map.csv",OUT/"oneD_v2_fem_qualified_G_map.csv","test")
 s=p.evaluate(100e-6,1e-5)
 assert s.fem_native_J_J_m2 != s.structural_G_J_m2
 assert s.fem_native_KJ_Pa_sqrt_m != s.qualified_KG_Pa_sqrt_m

def test_both_providers_fail_closed_outside_map():
 p=PFMechanicsProvider.from_csv(OUT/"oneD_v2_pf_mechanics_map.csv","test")
 with pytest.raises(ValueError): p.evaluate(1001e-6,1e-5)

def test_baselines_not_launched_without_local_stress_closure():
 d=json.loads((OUT/"oneD_v2_baseline_decision.json").read_text())
 assert d["prefixes_launched"]==d["long_cases_launched"]==0
 assert d["forward_launch_gate"].startswith("BLOCKED_")
 assert d["canonical_parameters_changed"] is False and d["refit_justified"] is False
 for name in ("oneD_v2_baseline_event_transactions.csv","oneD_v2_baseline_physical_avalanches.csv","oneD_v2_baseline_onset_candidates.csv","oneD_v2_baseline_summary.csv"):
  with (OUT/name).open() as f: assert list(csv.DictReader(f))==[]
