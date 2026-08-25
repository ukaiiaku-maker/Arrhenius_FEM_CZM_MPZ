#!/usr/bin/env python3
"""Generate the raw candidate-independent PF local tensor/profile map."""
from __future__ import annotations
import hashlib,json,os,sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'analysis_outputs/oneD_v2_local_drive'
PF=Path('/private/tmp/pf-v2-tensor-diagnostic')
sys.path.insert(0,str(PF))
os.environ.setdefault('PF_TENSOR_DRIVE_DIAGNOSTIC_JSONL',str(OUT/'_map_probe_enabled_sentinel.jsonl'))
from arrhenius_fracture.config import ElasticProperties,GeometryConfig,MeshConfig
from arrhenius_fracture.crack_backend import SharpWakeBackend
from arrhenius_fracture.crystal import cubic_plane_strain_D
from arrhenius_fracture.fem import assemble_mechanics,solve_dirichlet
from arrhenius_fracture.mesh import make_boundary_data,make_tri_mesh
from arrhenius_fracture.anisotropic_front_direction_fix_v10227 import install_front_direction_fix
from arrhenius_fracture.anisotropic_emission_v10174 import AnisotropicEmissionConfig
from arrhenius_fracture import anisotropic_emission_v10174 as ae

TARGETS=(0,2,5,10,25,50,75,100,150,200,300,400,500,600,700,800,900,1000)
U0=1e-5;U1=5e-6;DA=5e-6

def digest(*arrays):
 h=hashlib.sha256()
 for a in arrays:
  a=np.ascontiguousarray(a);h.update(str(a.dtype).encode());h.update(str(a.shape).encode());h.update(a.tobytes())
 return h.hexdigest()

def notch(mesh,g):
 x,y=np.asarray(mesh.nodes).T;d=np.zeros(mesh.nn);d[(x<=g.a0)&(np.abs(y)<=g.notch_half_thickness)]=1.;return d

def state(mesh,damage,U,D,mat,tip,cfg):
 bd=make_boundary_data(mesh,GeometryConfig());u0=np.zeros(mesh.ndof);ep=np.zeros((3,mesh.ne));rho=np.zeros(mesh.ne)
 K,R,*_=assemble_mechanics(mesh,u0,ep,rho,damage,D,mat,kappa=1e-6)
 u,reaction=solve_dirichlet(K,R,u0,bd,.5*U,-.5*U)
 _,_,sigma,*_=assemble_mechanics(mesh,u,ep,rho,damage,D,mat,kappa=1e-6)
 drive=ae.build_front_drive(mesh,sigma,damage,tip,cfg)
 return float(reaction),drive

def main():
 OUT.mkdir(parents=True,exist_ok=True);install_front_direction_fix()
 g=GeometryConfig();mc=MeshConfig(nx=36,ny=72,jitter=0.,tip_h_fine=1e-6,tip_ratio=1.2)
 mesh=make_tri_mesh(g,mc,seed=42,tip_center=np.array([g.a0,0.]));mat=ElasticProperties();D=cubic_plane_strain_D(523e9,203e9,160e9,0.)
 damage=notch(mesh,g);backend=SharpWakeBackend();bd=make_boundary_data(mesh,g);disp=np.zeros(mesh.ndof);kill=max(float(mesh.hbar_tip),.5e-6)
 cfg=AnisotropicEmissionConfig(crystal_theta_deg=0.,probe_radius_m=10e-6,sector_half_angle_deg=25.,damage_cutoff=.85,min_elements=3,schmid_reference=.5)
 counts=sorted(set(int(np.ceil(t/5.-1e-14)) for t in TARGETS));target_by={c:[t for t in TARGETS if int(np.ceil(t/5.-1e-14))==c] for c in counts}
 nodes=[];profile=[];scaling=[]
 for count in range(max(counts)+1):
  if count:
   p0=np.array([g.a0+(count-1)*DA,0.]);p1=np.array([g.a0+count*DA,0.])
   r=backend.advance(mesh=mesh,boundary=bd,damage=damage,displacement=disp,p0=p0,p1=p1,direction=np.array([1.,0.]),front_id=0,kill_r=kill)
   if not r.inserted:raise RuntimeError(r.reason)
   damage=r.damage
  if count not in target_by:continue
  tip=np.array([g.a0+count*DA,0.]);F0,d0=state(mesh,damage,U0,D,mat,tip,cfg);F1,d1=state(mesh,damage,U1,D,mat,tip,cfg)
  for target in target_by[count]:
   row={'target_extension_um':target,'actual_extension_um':count*5.,'reference_opening_m':U0,'reaction_N_per_m':F0,
    'tip_x_m':tip[0],'tip_y_m':tip[1],'front_tangent_x':d0['front_direction'][0],'front_tangent_y':d0['front_direction'][1],
    'front_normal_x':d0['front_normal'][0],'front_normal_y':d0['front_normal'][1],
    'opening_sigma_xx_Pa':d0['opening_tensor_Pa'][0][0],'opening_sigma_yy_Pa':d0['opening_tensor_Pa'][1][1],'opening_sigma_xy_Pa':d0['opening_tensor_Pa'][0][1],
    'channel0_sigma_xx_Pa':d0['channel_tensors_Pa'][0][0][0],'channel0_sigma_yy_Pa':d0['channel_tensors_Pa'][0][1][1],'channel0_sigma_xy_Pa':d0['channel_tensors_Pa'][0][0][1],
    'channel1_sigma_xx_Pa':d0['channel_tensors_Pa'][1][0][0],'channel1_sigma_yy_Pa':d0['channel_tensors_Pa'][1][1][1],'channel1_sigma_xy_Pa':d0['channel_tensors_Pa'][1][0][1],
    'tau0_signed_Pa':d0['tau_signed_Pa'][0],'tau1_signed_Pa':d0['tau_signed_Pa'][1],'factor0':d0['drive_factors'][0],'factor1':d0['drive_factors'][1],
    'opening_probe_elements':d0['opening_probe_elements'],'channel0_probe_elements':d0['channel_probe_elements'][0],'channel1_probe_elements':d0['channel_probe_elements'][1],
    'mesh_fingerprint':digest(mesh.nodes,mesh.elems),'wake_fingerprint':digest(damage),'reliable':d0['reliable']}
   nodes.append(row)
   for probe_name,probe,tensor in [('opening',d0['opening_probe'],d0['opening_tensor_Pa']),('channel0',d0['channel_probes'][0],d0['channel_tensors_Pa'][0]),('channel1',d0['channel_probes'][1],d0['channel_tensors_Pa'][1])]:
    for idx,xy,w,sig,dam in zip(probe['element_indices'],probe['element_centroids_m'],probe['element_weights'],probe['element_stress_components_Pa'],probe['element_damage']):
     profile.append({'target_extension_um':target,'actual_extension_um':count*5.,'probe':probe_name,'element_index':idx,'x_m':xy[0],'y_m':xy[1],'s_local_m':xy[0]-tip[0],'n_local_m':xy[1]-tip[1],'weight':w,'sigma_xx_Pa':sig[0],'sigma_yy_Pa':sig[1],'sigma_xy_Pa':sig[2],'damage':dam})
  tensor_keys=['opening_tensor_Pa','channel_tensors_Pa','tau_signed_Pa']
  for key in tensor_keys:
   a=np.asarray(d0[key],float);b=np.asarray(d1[key],float)
   scaling.append({'actual_extension_um':count*5.,'quantity':key,'max_linear_opening_error':float(np.max(np.abs(b-.5*a))/max(np.max(np.abs(.5*a)),1.)),'sign_mismatch_count':int(np.count_nonzero(np.sign(b)!=np.sign(a)))})
  scaling.append({'actual_extension_um':count*5.,'quantity':'reaction','max_linear_opening_error':abs(F1/.5/F0-1.),'sign_mismatch_count':0})
 nodes=pd.DataFrame(nodes);prof=pd.DataFrame(profile);validation_rows=list(scaling)
 unique=nodes.drop_duplicates('actual_extension_um').sort_values('actual_extension_um')
 tensor_columns=[
  'opening_sigma_xx_Pa','opening_sigma_yy_Pa','opening_sigma_xy_Pa',
  'channel0_sigma_xx_Pa','channel0_sigma_yy_Pa','channel0_sigma_xy_Pa',
  'channel1_sigma_xx_Pa','channel1_sigma_yy_Pa','channel1_sigma_xy_Pa',
  'tau0_signed_Pa','tau1_signed_Pa','factor0','factor1']
 x=unique.actual_extension_um.to_numpy(float)
 for column in tensor_columns:
  y=unique[column].to_numpy(float)
  for i in range(1,len(x)-1):
   prediction=float(np.interp(x[i],np.delete(x,i),np.delete(y,i)))
   validation_rows.append({'actual_extension_um':x[i],'quantity':column,
    'validation_kind':'LEAVE_ONE_NODE_OUT','observed':y[i],'predicted':prediction,
    'absolute_error':abs(prediction-y[i]),'relative_error':abs(prediction-y[i])/max(abs(y[i]),1.),
    'max_linear_opening_error':np.nan,'sign_mismatch_count':int(np.sign(prediction)!=np.sign(y[i]))})
 validation=pd.DataFrame(validation_rows)
 prof.to_parquet(OUT/'oneD_v2_pf_tensor_drive_nodes.parquet',index=False);nodes.to_csv(OUT/'oneD_v2_pf_tensor_drive_node_summary.csv',index=False);validation.to_csv(OUT/'oneD_v2_pf_tensor_interpolation_validation.csv',index=False)
 manifest={'schema':'oneD_v2_pf_tensor_drive_map_manifest_v2','status':'PF_RAW_TENSOR_PROFILE_MAP_QUALIFIED_AT_SOURCE_NODES',
  'model_kind':'ONE_D_FRONT_MODEL_WITH_LOCAL_TENSOR_DRIVE','source_commit':'9e884fb0b0845da621d2612bdf1042e481b8df49','observer_default_off':True,'observer_analysis_commit':'187c5ed',
  'candidate_independent':True,'bulk_constraint':'PLANE_STRAIN','crack_representation':'FINITE_WIDTH_STIFFNESS_KILLED_WAKE','failed_region_state':'BINARY_STIFFNESS_KILL',
  'reference_opening_m':U0,'extensions_um':TARGETS,'actual_extensions_um':sorted(nodes.actual_extension_um.unique().tolist()),'probe_radius_m':cfg.probe_radius_m,
  'profile_rows':len(prof),'all_reliable':bool(nodes.reliable.all()),'maximum_linear_scaling_error':float(validation.max_linear_opening_error.max()),
  'maximum_LOO_relative_error':float(validation.relative_error.max()),'LOO_sign_mismatch_count':int(validation.loc[validation.validation_kind.eq('LEAVE_ONE_NODE_OUT'),'sign_mismatch_count'].sum()),
  'interpolation':'BOUNDED_PIECEWISE_LINEAR_SIGN_PRESERVING','extrapolation':'FAIL_CLOSED','nodes_sha256':hashlib.sha256((OUT/'oneD_v2_pf_tensor_drive_nodes.parquet').read_bytes()).hexdigest(),
  'node_summary_sha256':hashlib.sha256((OUT/'oneD_v2_pf_tensor_drive_node_summary.csv').read_bytes()).hexdigest()}
 (OUT/'oneD_v2_pf_tensor_drive_map_manifest.json').write_text(json.dumps(manifest,indent=2,sort_keys=True)+'\n')

if __name__=='__main__':main()
