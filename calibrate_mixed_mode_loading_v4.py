#!/usr/bin/env python3
"""Anisotropic process-zone traction calibration for mixed-mode FEM/CZM v4."""
from __future__ import annotations
import argparse,csv,json,math
from pathlib import Path
import numpy as np

MODEL_ID='FEM_CZM_mixed_mode_calibration_v4_anisotropic_traction_basis'
def vals(s):return [float(x) for x in str(s).replace(',',' ').split() if x]
def err(v,t):return (float(v)-float(t)+180)%360-180

def main():
 p=argparse.ArgumentParser();p.add_argument('--out',default='runs/mixed_mode_fem_czm_v4_anisotropic_calibration')
 p.add_argument('--target-psi-deg',default='-60 -45 -30 -15 0 15 30 45 60');p.add_argument('--U-cal-m',type=float,default=2e-7)
 p.add_argument('--nx',type=int,default=24);p.add_argument('--ny',type=int,default=48);p.add_argument('--tip-h-fine',type=float,default=3e-6);p.add_argument('--tip-ratio',type=float,default=1.25);p.add_argument('--mesh-seed',type=int,default=42)
 p.add_argument('--crystal-theta-deg',type=float,default=0.0);p.add_argument('--crystal-C11',type=float,default=523e9);p.add_argument('--crystal-C12',type=float,default=203e9);p.add_argument('--crystal-C44',type=float,default=160e9)
 p.add_argument('--traction-probe-radius-m',type=float,default=10e-6);p.add_argument('--traction-annulus-half-width',type=float,default=.45);p.add_argument('--traction-sector-half-angle-deg',type=float,default=40);p.add_argument('--psi-tol-deg',type=float,default=.75);p.add_argument('--basis-condition-max',type=float,default=1e8)
 a=p.parse_args()
 from arrhenius_fracture.config import make_emergent_config
 from arrhenius_fracture.mesh import make_tri_mesh,make_boundary_data
 from arrhenius_fracture.fem import assemble_mechanics,solve_dirichlet
 from arrhenius_fracture.crystal import cubic_plane_strain_D,zener_ratio
 from arrhenius_fracture.j_integral import compute_J_integral
 from arrhenius_fracture.mixed_mode_first_passage_v4 import (AnisotropicContext,_mixed_solve_factory,process_zone_traction_probe,
    shear_sign_from_basis,loading_angle_from_response_basis,traction_phase_deg,energy_matrix_from_basis,energy_phase_from_matrix)
 cfg=make_emergent_config();cfg.mesh.nx=a.nx;cfg.mesh.ny=a.ny;cfg.mesh.tip_h_fine=a.tip_h_fine;cfg.mesh.tip_ratio=a.tip_ratio
 mesh=make_tri_mesh(cfg.geometry,cfg.mesh,seed=a.mesh_seed);bnd=make_boundary_data(mesh,cfg.geometry);mat=cfg.material
 D=cubic_plane_strain_D(a.crystal_C11,a.crystal_C12,a.crystal_C44,a.crystal_theta_deg)
 u0=np.zeros(mesh.ndof);ep=np.zeros((3,mesh.ne));rho=np.zeros(mesh.ne);d=np.zeros(mesh.nn);d[bnd.notch_nodes]=1
 K,R,*_=assemble_mechanics(mesh,u0,ep,rho,d,D,mat);cache={}
 def evaluate(alpha):
  key=round(float(alpha),10)
  if key in cache:return dict(cache[key])
  ctx=AnisotropicContext(alpha,0,a.crystal_theta_deg,0.3,a.traction_probe_radius_m,a.traction_annulus_half_width,a.traction_sector_half_angle_deg)
  solve=_mixed_solve_factory(solve_dirichlet,ctx);u,F=solve(K,R,u0.copy(),bnd,.5*a.U_cal_m,-.5*a.U_cal_m)
  _,_,sig,_,_,psi=assemble_mechanics(mesh,u,ep,rho,d,D,mat)
  pr=process_zone_traction_probe(mesh,sig,d,np.array([cfg.geometry.a0,0.]),np.array([1.,0.]),a.traction_probe_radius_m,a.traction_annulus_half_width,a.traction_sector_half_angle_deg)
  J,KJ,ji=compute_J_integral(mesh,u,sig,psi,d,np.array([cfg.geometry.a0,0.]),np.array([1.,0.]),mat,ell=20e-6)
  row={'loading_angle_deg':float(alpha),'sigma_nn_raw_Pa':pr.get('sigma_nn_Pa',np.nan),'tau_tn_raw_Pa':pr.get('tau_tn_Pa',np.nan),'probe_reliable':pr.get('reliable',False),'probe_n_elements':pr.get('n_elements',0),'J_J_per_m2':J,'KJ_reference_Pa_sqrt_m':KJ,'generalized_reaction_N':F}
  cache[key]=dict(row);return row
 o=evaluate(0);s=evaluate(90)
 Mraw=np.array([[o['sigma_nn_raw_Pa'],s['sigma_nn_raw_Pa']],[o['tau_tn_raw_Pa'],s['tau_tn_raw_Pa']]])
 sign=shear_sign_from_basis(Mraw);M=np.diag([1,sign])@Mraw;cond=float(np.linalg.cond(M))
 if not np.isfinite(cond) or cond>a.basis_condition_max:raise SystemExit(f'traction basis invalid cond={cond}')
 Jeq=evaluate(45)['J_J_per_m2'];G=energy_matrix_from_basis(o['J_J_per_m2'],s['J_J_per_m2'],Jeq,a.U_cal_m)
 ew=np.linalg.eigvalsh(G); energy_ok=bool(np.all(np.isfinite(ew)) and np.min(ew)>0)
 rows=[]
 for target in vals(a.target_psi_deg):
  try:alpha=loading_angle_from_response_basis(M,target)
  except Exception as ex:
   rows.append({'target_psi_deg':target,'phase_converged':False,'reason':f'basis_error:{ex}'});continue
  r=evaluate(alpha);phase=traction_phase_deg(r['sigma_nn_raw_Pa'],r['tau_tn_raw_Pa'],sign);e=err(phase,target)
  row={**r,'target_psi_deg':target,'loading_angle_deg':alpha,'traction_shear_sign':sign,'achieved_traction_phase_deg':phase,'traction_phase_error_deg':e,'phase_converged':bool(r['probe_reliable'] and abs(e)<=a.psi_tol_deg),'basis_condition':cond,'zener_ratio':zener_ratio(a.crystal_C11,a.crystal_C12,a.crystal_C44),'crystal_theta_deg':a.crystal_theta_deg,'energy_matrix_positive_definite':energy_ok,'energy_phase_deg':energy_phase_from_matrix(G,alpha),'response_11_Pa':M[0,0],'response_12_Pa':M[0,1],'response_21_Pa':M[1,0],'response_22_Pa':M[1,1],'energy_G11':G[0,0],'energy_G12':G[0,1],'energy_G22':G[1,1]}
  rows.append(row)
 out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
 with (out/'mixed_mode_loading_calibration_v4.csv').open('w',newline='') as fp:
  cols=sorted({k for r in rows for k in r});w=csv.DictWriter(fp,fieldnames=cols);w.writeheader();w.writerows(rows)
 (out/'anisotropic_calibration_basis.json').write_text(json.dumps({'model':MODEL_ID,'raw_traction_basis_Pa':Mraw.tolist(),'normalized_traction_basis_Pa':M.tolist(),'traction_shear_sign':sign,'basis_condition':cond,'energy_matrix':G.tolist(),'energy_eigenvalues':ew.tolist(),'energy_matrix_positive_definite':energy_ok,'crystal':{'C11':a.crystal_C11,'C12':a.crystal_C12,'C44':a.crystal_C44,'theta_deg':a.crystal_theta_deg,'zener_ratio':zener_ratio(a.crystal_C11,a.crystal_C12,a.crystal_C44)},'probe_radius_m':a.traction_probe_radius_m},indent=2))
 print('raw traction basis [Pa per calibration amplitude]:\n',Mraw);print('shear sign:',sign);print('normalized traction basis:\n',M);print('condition:',cond);print('energy J matrix:',G)
 for r in rows:print({'target':r.get('target_psi_deg'),'alpha':r.get('loading_angle_deg'),'phase':r.get('achieved_traction_phase_deg'),'error':r.get('traction_phase_error_deg'),'ok':r.get('phase_converged'),'n':r.get('probe_n_elements')})
 bad=[r for r in rows if not r.get('phase_converged',False)]
 if bad:raise SystemExit('anisotropic traction calibration failed: '+', '.join(str(r.get('target_psi_deg')) for r in bad))
if __name__=='__main__':main()
