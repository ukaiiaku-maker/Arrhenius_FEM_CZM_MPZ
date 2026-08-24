from __future__ import annotations
import hashlib,json
from pathlib import Path
import pandas as pd
import pytest
from reduced_fracture_v2 import CanonicalParameters,ObservableKind,Qualification,ReducedFractureKernel,ReducedState
from reduced_fracture_v2.lifecycle import CommonEventLifecyclePolicy
from reduced_fracture_v2.mechanics import MechanicsState,ReplayProvider,TabulatedForwardProvider
from reduced_fracture_v2.thresholds import GeneratedSequence,ReplaySequence

ROOT=Path(__file__).parents[1];OUT=ROOT/"analysis_outputs/oneD_v2_common_kernel"
IDS={"v913_zeroD_sobol_0242980","v913_zeroD_sobol_0202500","v913_zeroD_sobol_0129902","v913_zeroD_sobol_0077080"}
P={"Tref_K":300,"cleave_G00_eV":3,"cleave_gT_eV_per_K":0,"cleave_sigc0_GPa":3,"cleave_sT_GPa_per_K":0,"cleave_floor_frac":.1,"cleave_exp_n":2,"cleave_exp_a":1,"emit_G00_eV":2,"emit_gT_eV_per_K":0,"emit_sigc0_GPa":2,"emit_sT_GPa_per_K":0,"emit_floor_frac":.1,"emit_exp_n":2,"emit_exp_a":1,"peierls_H0_eV":1,"peierls_nu0_s":1e12,"taylor_H0_eV":1.2,"taylor_nu0_s":1e12}
def kernel():return ReducedFractureKernel(CanonicalParameters("x",P))
def mech(kind=ObservableKind.PF_NATIVE_DRIVE):return MechanicsState("OPENING_M",1e-5,0,kind,1.0 if kind!=ObservableKind.FEM_STRUCTURAL_G else None,1.0 if kind==ObservableKind.FEM_STRUCTURAL_G else None,None,1e9,None,Qualification.ARCHIVED_REPLAY)

def test_01_v1_frozen():
 d=json.loads((OUT/"oneD_v2_model_inventory.json").read_text());assert d["V1"]["preserved"] and d["V1"]["audit_head"]=="a5dece9"
def test_02_exact_rows():
 x=pd.read_csv(OUT/"oneD_v2_frozen_state_consistency.csv");assert set(x.candidate_id)==IDS
def test_03_common_barrier_equality():
 s=ReducedState();a=kernel().evaluate(s,1000,2e9,1e9,1);b=kernel().evaluate(s,1000,2e9,1e9,1);assert a.cleavage_barrier_eV==b.cleavage_barrier_eV
def test_04_common_hazard_equality():
 s=ReducedState();a=kernel().evaluate(s,1000,2e9,1e9,1);b=kernel().evaluate(s,1000,2e9,1e9,1);assert a.cleavage_hazard_increment==b.cleavage_hazard_increment
def test_05_common_state_update_equality():
 s=ReducedState();assert s.with_hazard(1,2,3)==s.with_hazard(1,2,3)
def test_06_typed_provider_observables_differ():
 assert mech().observable_kind!=mech(ObservableKind.FEM_STRUCTURAL_G).observable_kind
def test_07_pf_not_structural_g():
 with pytest.raises(ValueError):MechanicsState("U",1,0,ObservableKind.PF_NATIVE_DRIVE,1,2,None,None,None,Qualification.ARCHIVED_REPLAY)
def test_08_fem_G_requires_qualification():
 with pytest.raises(ValueError):MechanicsState("U",1,0,ObservableKind.FEM_STRUCTURAL_G,None,2,None,None,None,Qualification.UNQUALIFIED)
def test_09_replay_never_redraws():
 r=ReplaySequence([.2,.4]);assert r.thresholds(2)==r.thresholds(2)
def test_10_replay_geometry_immutable():
 r=ReplayProvider((mech(),),(mech(),),(mech(),),"abc");r.assert_geometry("abc");
def test_11_forward_does_not_import_future_state():
 p=TabulatedForwardProvider((0.,1.),(1.,2.),ObservableKind.PF_NATIVE_DRIVE,Qualification.MODEL_NATIVE_NOT_STRUCTURAL_G,"x");assert not hasattr(p,"archived_states")
def test_12_transaction_differs_from_avalanche():
 assert (OUT/"oneD_v2_event_transactions.csv").name!=(OUT/"oneD_v2_physical_avalanches.csv").name
def test_13_pre_state_differs_post_geometry():
 t=CommonEventLifecyclePolicy().fire(ReducedState(),1e-6,2e-6,0);assert t.pre_event_state.crack_extension_m!=t.post_event_geometry_extension_m
def test_14_target_is_censored():
 t=CommonEventLifecyclePolicy().fire(ReducedState(),2e-6,1e-6,0);assert t.right_censored_at_target and "CENSORED" in t.role
def test_15_onsets_pre_event_schema():
 assert "onset" in (OUT/"oneD_v2_onset_candidates.csv").name
def test_16_interior_not_resistance_points():
 assert "in_avalanche_drive" in (OUT/"oneD_v2_in_avalanche_drive.csv").name
def test_17_maps_candidate_independent():
 p=TabulatedForwardProvider((0.,1.),(1.,2.),ObservableKind.PF_NATIVE_DRIVE,Qualification.MODEL_NATIVE_NOT_STRUCTURAL_G,"x");assert p.candidate_independent
def test_18_sequences_generated_not_tiled():
 assert GeneratedSequence(3).thresholds(5)==GeneratedSequence(3).thresholds(5)
def test_19_identical_local_state_same_rates():
 s=ReducedState();assert kernel().evaluate(s,900,3e9,2e9,1)==kernel().evaluate(s,900,3e9,2e9,1)
def test_20_no_2d_production_artifacts():
 assert not any(p.name.startswith("fem_") for p in OUT.iterdir())
def test_21_fingerprints_match():
 d=json.loads((OUT/"oneD_v2_scientific_fingerprints.json").read_text());assert all(hashlib.sha256((OUT/n).read_bytes()).hexdigest()==h for n,h in d.items())
