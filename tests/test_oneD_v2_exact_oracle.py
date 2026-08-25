from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd
import pytest

from reduced_fracture_v2.exact_oracle import (ElasticFieldCacheKey,ExactOracleCache,
    NormalizedElasticFieldSnapshot,OptionalRegimeAwareSurrogate,ProbeQueryKey)
from reduced_fracture_v2.local_drive import DriveQualificationError
from reduced_fracture_v2.production_oracles import (FEMCZMExactElasticFieldOracle,
    PFExactElasticFieldOracle)

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'analysis_outputs/oneD_v2_exact_oracle'


def snapshot(**changes):
    values=dict(backend='TEST',source_repository='repo',source_commit='source',factory_wrapper_identity='factory',
        mechanics_factory_hash='mechanics',bulk_constraint='PLANE_STRAIN',crack_representation='CUT',failed_region_state='FREE',
        geometry_fingerprint='geometry',mesh_fingerprint='mesh',topology_wake_fingerprint='topology',
        crack_tip_xy_m=(0.,0.),crack_tangent=(1.,0.),crack_normal=(0.,1.),
        nodes_m=np.array([[0.,0.],[1.,0.],[0.,1.]]),elements=np.array([[0,1,2]]),
        element_centroids_m=np.array([[1/3,1/3]]),element_areas_m2=np.array([.5]),
        stress_per_opening_Pa_per_m=np.array([[[2.,-1.],[-1.,3.]]]),
        strain_per_opening_per_m=np.array([[[.2,-.1],[-.1,.3]]]),damage_or_interface_state=np.zeros(3),
        reaction_per_opening_N_per_m2=4.,elastic_energy_per_opening2_J_per_m3=2.,native_J_per_opening2_J_per_m4=5.,
        native_KJ_per_opening_Pa_sqrt_m_per_m=6.,
        native_J_domain_metadata={'r':1},validity_flags=('LINEAR_ELASTIC','CANDIDATE_INDEPENDENT'),reference_opening_m=1.)
    values.update(changes);return NormalizedElasticFieldSnapshot(**values)


def test_01_normalized_field_is_candidate_independent_and_read_only():
    s=snapshot();assert 'candidate' not in s.__dict__ and not s.stress_per_opening_Pa_per_m.flags.writeable
    with pytest.raises(DriveQualificationError,match='forbidden'):
        snapshot(validity_flags=('candidate_specific_state_embedded',))


def test_02_opening_scaling_is_exact():
    s=snapshot();assert np.array_equal(s.stress_at_opening(3.),3*s.stress_at_opening(1.))
    assert s.reaction_per_opening_N_per_m2*3==12 and s.native_J_per_opening2_J_per_m4*3**2==45


def test_03_field_cache_key_includes_full_geometry_topology_and_basis():
    a=snapshot();b=replace(a,geometry_fingerprint='different');c=replace(a,topology_wake_fingerprint='different')
    assert len({a.cache_key.fingerprint(),b.cache_key.fingerprint(),c.cache_key.fingerprint()})==3


def test_04_process_zone_state_does_not_trigger_field_solve_at_fixed_geometry():
    cache=ExactOracleCache();s=snapshot();count={'n':0}
    def solve():count['n']+=1;return s
    assert cache.field(s.cache_key,solve) is cache.field(s.cache_key,solve)
    assert count['n']==1 and cache.field_hits==1 and cache.field_misses==1
    for radius in (1.,2.):
        key=ProbeQueryKey(s.snapshot_hash,radius,3.,4.,'g','p','slip','probe')
        cache.probes[key.fingerprint()]=None
    assert count['n']==1


def test_05_probe_key_uses_current_radius_front_width_and_source_geometry():
    s=snapshot();base=ProbeQueryKey(s.snapshot_hash,1.,2.,3.,'g','p','slip','probe')
    assert len({base.fingerprint(),replace(base,tip_radius_m=2.).fingerprint(),replace(base,front_width_m=3.).fingerprint(),replace(base,source_geometry_fingerprint='h').fingerprint()})==4


def test_06_backstress_is_not_embedded_in_elastic_field():
    assert 'backstress' not in snapshot().__dict__
    with pytest.raises(DriveQualificationError,match='forbidden'):snapshot(validity_flags=('backstress_embedded',))


def test_07_pf_probe_reproduces_all_existing_observer_tensors():
    data=pd.read_csv(OUT/'oneD_v2_pf_exact_oracle_validation.csv')
    assert len(data)==463 and data.status.str.startswith('PASS').all()
    assert data.probe_support_exact.all() and data.maximum_tensor_relative_error.max()<=5e-12


def test_08_pf_full_dynamic_shadow_remains_fail_closed_when_kinetics_are_absent():
    data=pd.read_csv(OUT/'oneD_v2_pf_dynamic_state_shadow.csv')
    assert len(data)==8 and data.status.eq('PARTIAL_TENSOR_DOMAIN_ONLY_FULL_KINETIC_STATE_UNAVAILABLE').all()
    assert data.mobile_state_status.eq('NOT_SERIALIZED').all()


def test_09_fem_factory_uses_authoritative_d41_kernel():
    data=pd.read_csv(OUT/'oneD_v2_fem_exact_oracle_validation.csv')
    assert data.kernel_sha256.nunique()==1 and data.kernel_sha256.iloc[0]==FEMCZMExactElasticFieldOracle.kernel_sha256


def test_10_fem_exact_probe_reproduces_corrected_factory_resolved_drives():
    data=pd.read_csv(OUT/'oneD_v2_fem_exact_oracle_validation.csv')
    assert len(data)==6 and (data.maximum_factory_resolved_drive_absolute_error==0).all() and data.reliability_exact.all()


def test_11_scalar_kj_cannot_replace_tensor_input_contract_is_preserved():
    schema=json.loads((OUT/'oneD_v2_normalized_field_schema.json').read_text())
    assert 'stress_per_opening_Pa_per_m' in schema['required']


def test_12_common_tensor_parity_remains_exact():
    data=pd.read_csv(OUT/'oneD_v2_common_tensor_parity_v2.csv')
    exact=[c for c in data if c.endswith('_exact')]
    assert len(data)==4 and data[exact].all().all() and data.status.eq('PASS_EXACT_ARCHITECTURAL_COMMON_MODE').all()


def test_13_level3_microtrajectories_cannot_claim_surrogate_or_archived_mechanics():
    data=pd.read_csv(OUT/'oneD_v2_level3_exact_microtrajectories.csv')
    assert len(data)==4 and data.accepted_events.eq(0).all() and data.status.str.startswith('NOT_LAUNCHED').all()


def test_14_surrogate_is_locked_before_exact_level3():
    with pytest.raises(DriveQualificationError,match='forbidden'):
        OptionalRegimeAwareSurrogate('NATIVE_CLOSURE_PARTIALLY_QUALIFIED').evaluate()


def test_15_out_of_domain_queries_fail_closed():
    with pytest.raises(DriveQualificationError,match='straight interface'):
        FEMCZMExactElasticFieldOracle().solve({'events':[{'p0_m':[0.,0.],'p1_m':[1e-6,1e-6]}]})


def test_16_no_future_archived_state_enters_forward_stepping_contract():
    required=set(json.loads((OUT/'oneD_v2_normalized_field_schema.json').read_text())['required'])
    assert not any('future' in field or 'archive' in field for field in required)


def test_17_canonical_parameters_remain_unchanged():
    decision=json.loads((OUT/'oneD_v2_level3_decision_v4.json').read_text())
    assert not decision['canonical_parameters_changed'] and decision['parameter_status']=='INSUFFICIENT_EVIDENCE_FOR_V2_PARAMETER_DECISION'


def test_18_no_new_stochastic_run_was_launched():
    manifest=json.loads((OUT/'oneD_v2_exact_oracle_cache_manifest.json').read_text())
    assert manifest['new_stochastic_PF_runs']==0 and manifest['new_stochastic_FEMCZM_runs']==0
    assert manifest['PF_field_solves']==7 and manifest['FEM_field_solves']==2


def test_exact_oracle_bundle_fingerprints_reproduce():
    fingerprints=json.loads((OUT/'oneD_v2_exact_oracle_fingerprints.json').read_text())
    assert all(hashlib.sha256((OUT/name).read_bytes()).hexdigest()==value for name,value in fingerprints.items())
