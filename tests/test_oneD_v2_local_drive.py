import numpy as np
import pytest
from pathlib import Path
import json
import pandas as pd

from reduced_fracture_v2.local_drive import (
    DriveQualificationError, EmissionDrive, GlobalStructuralDrive,
    NativeDriveBundle, interpolate_symmetric_tensor,
    PFTensorDriveProvider,
)


def _bundle(emission=None, radius=1e-6):
    return NativeDriveBundle(
        GlobalStructuralDrive(1e-6, 1.0, 2.0, 3.0), (0.0, 0.0),
        (1.0, 0.0), (0.0, 1.0), (), emission, radius, 2e-6,
        "topology", "source", "map", (0.0, 1e-3), (0.0, 2e-6, -1e-6, 1e-6),
    )


def test_scalar_kj_cannot_satisfy_pf_signed_emission():
    with pytest.raises(DriveQualificationError, match="scalar KJ is insufficient"):
        _bundle().require_pf_signed_emission()


def test_global_kj_and_local_tensor_are_separate_fields():
    assert _bundle().global_drive.native_KJ_Pa_sqrt_m == 3.0
    assert _bundle().emission is None


def test_tensor_interpolation_preserves_sign_and_symmetry():
    out = interpolate_symmetric_tensor(.5, 0., 1., ((2., -4.), (-4., 6.)), ((4., -2.), (-2., 8.)))
    assert np.array_equal(out, out.T) and out[0, 1] == -3.0


def test_tensor_interpolation_and_radius_fail_closed_outside_domain():
    with pytest.raises(DriveQualificationError, match="outside qualified bounds"):
        interpolate_symmetric_tensor(2., 0., 1., ((1., 0.), (0., 1.)), ((2., 0.), (0., 2.)))
    with pytest.raises(DriveQualificationError, match="outside the qualified profile"):
        _bundle(radius=3e-6).require_radius_in_domain()


def test_source_projection_retains_signed_shear_and_basis():
    root2=np.sqrt(2.0)
    drive=EmissionDrive(
        ((10., 0.), (0., 20.)),
        (((10., 4.), (4., 20.)), ((10., -4.), (-4., 20.))),
        ((1/root2, 1/root2), (1/root2, -1/root2)),
        ((-1/root2, 1/root2), (1/root2, 1/root2)),
        (5., -5.), (0.5, 0.5), True, 1,
    )
    assert _bundle(drive).require_pf_signed_emission() is drive
    _bundle(drive).validate_basis()


def test_unreliable_drive_fails_closed():
    root2=np.sqrt(2.0)
    drive=EmissionDrive(
        ((1.,0.),(0.,1.)), (((1.,0.),(0.,1.)),)*2,
        ((1/root2,1/root2),(1/root2,-1/root2)),
        ((-1/root2,1/root2),(1/root2,1/root2)),
        (0.,0.), (0.,0.), False, 0,
    )
    with pytest.raises(DriveQualificationError, match="reliable 2-D tensor drive"):
        _bundle(drive).require_pf_signed_emission()


def test_unqualified_pf_tensor_map_fails_closed(tmp_path):
    csv=tmp_path/'nodes.csv';manifest=tmp_path/'manifest.json'
    pd.DataFrame([{'actual_extension_um':0}]).to_csv(csv,index=False)
    manifest.write_text(json.dumps({'status':'PF_RAW_TENSOR_PROFILE_MAP_QUALIFIED_AT_SOURCE_NODES'}))
    with pytest.raises(DriveQualificationError,match='interpolation is not qualified'):
        PFTensorDriveProvider.from_artifacts(csv,manifest)


ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'analysis_outputs/oneD_v2_local_drive'


def test_raw_pf_tensor_nodes_scale_linearly_and_retain_source_signs():
    manifest=json.loads((OUT/'oneD_v2_pf_tensor_drive_map_manifest.json').read_text())
    validation=pd.read_csv(OUT/'oneD_v2_pf_tensor_interpolation_validation.csv')
    assert manifest['all_reliable'] and manifest['maximum_linear_scaling_error']==0.0
    scaling=validation[validation.quantity.isin(['opening_tensor_Pa','channel_tensors_Pa','tau_signed_Pa'])]
    assert (scaling.max_linear_opening_error==0.0).all()
    assert (scaling.sign_mismatch_count==0).all()


def test_sparse_pf_tensor_interpolation_fails_sign_sensitive_gate():
    manifest=json.loads((OUT/'oneD_v2_pf_tensor_drive_map_manifest.json').read_text())
    assert manifest['status']=='PF_RAW_TENSOR_PROFILE_MAP_QUALIFIED_AT_SOURCE_NODES'
    assert manifest['LOO_sign_mismatch_count']==23
    with pytest.raises(DriveQualificationError,match='interpolation is not qualified'):
        PFTensorDriveProvider.from_artifacts(
            OUT/'oneD_v2_pf_tensor_drive_node_summary.csv',
            OUT/'oneD_v2_pf_tensor_drive_map_manifest.json')


def test_bounded_pf_observer_is_trajectory_neutral_and_event_bounded():
    audit=json.loads((OUT/'oneD_v2_pf_tensor_diagnostic_neutrality.json').read_text())
    assert audit['qualification']=='PASS_OFF_ON_PHYSICAL_BYTE_IDENTITY'
    assert {k:v['accepted_geometry_events'] for k,v in audit['cases'].items()}=={'Peak':5,'DBTT':3}
    assert all(v['all_physical_artifacts_identical'] for v in audit['cases'].values())


def test_level3_stays_closed_and_parameters_unchanged():
    d=json.loads((OUT/'oneD_v2_level3_decision_v2.json').read_text())
    assert d['Level3']=='NATIVE_CLOSURE_UNQUALIFIED'
    assert not d['forward_baselines_authorized'] and not d['canonical_parameters_changed']
    assert d['new_2d_PF_diagnostics_launched']==2 and d['new_2d_FEMCZM_runs_launched']==0


def test_local_drive_fingerprints_reproduce():
    import hashlib
    fingerprints=json.loads((OUT/'oneD_v2_local_drive_fingerprints.json').read_text())
    assert all(hashlib.sha256((OUT/name).read_bytes()).hexdigest()==value
               for name,value in fingerprints.items())
