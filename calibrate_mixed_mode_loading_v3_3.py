#!/usr/bin/env python3
"""Two-basis isotropic elastic calibration for mixed-mode FEM/CZM v3.3.

The authoritative amplitude in v3.3 comes from the domain J integral.  This
calibrator therefore uses the Williams projection only for the phase ratio and
never rejects a calibration because its amplitude residual is high.
"""
from __future__ import annotations
import argparse, csv, json, math
from pathlib import Path
import numpy as np

MODEL_ID = "FEM_CZM_mixed_mode_first_passage_v3_3_J_consistent_circular_phase"

def vals(text):
    return [float(x) for x in str(text).replace(',', ' ').split() if x]

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--out',default='runs/mixed_mode_fem_czm_v3_3_circular_phase_calibration')
    p.add_argument('--target-psi-deg',default='-60 -45 -30 -15 0 15 30 45 60')
    p.add_argument('--diagnostic-alpha-grid-deg',default='-70 -60 -45 -30 -15 0 15 30 45 60 70')
    p.add_argument('--U-cal-m',type=float,default=2e-7)
    p.add_argument('--nx',type=int,default=24);p.add_argument('--ny',type=int,default=48)
    p.add_argument('--tip-h-fine',type=float,default=3e-6);p.add_argument('--tip-ratio',type=float,default=1.25)
    p.add_argument('--mesh-seed',type=int,default=42)
    p.add_argument('--mode-projection-angle-deg',type=float,default=105.0)
    p.add_argument('--mode-projection-damage-cutoff',type=float,default=0.85)
    p.add_argument('--psi-tol-deg',type=float,default=0.75)
    p.add_argument('--phase-spread-tol-deg',type=float,default=20.0)
    p.add_argument('--projection-condition-max',type=float,default=1e12)
    p.add_argument('--basis-condition-max',type=float,default=1e8)
    p.add_argument('--max-refine-iters',type=int,default=5)
    p.add_argument('--refine-dalpha-deg',type=float,default=1.0)
    a=p.parse_args()

    from arrhenius_fracture.config import make_emergent_config
    from arrhenius_fracture.mesh import make_tri_mesh, make_boundary_data
    from arrhenius_fracture.fem import plane_strain_D, assemble_mechanics, solve_dirichlet
    from arrhenius_fracture.mixed_mode_first_passage_v3_3 import (
        MixedModeContext,_mixed_solve_factory,project_near_tip_modes,
        loading_angle_from_mode_basis,phase_projection_gate,angle_error_deg,\
        mode_signs_from_basis,apply_mode_signs)

    cfg=make_emergent_config();cfg.mesh.nx=a.nx;cfg.mesh.ny=a.ny
    cfg.mesh.tip_h_fine=a.tip_h_fine;cfg.mesh.tip_ratio=a.tip_ratio
    mesh=make_tri_mesh(cfg.geometry,cfg.mesh,seed=a.mesh_seed)
    bnd=make_boundary_data(mesh,cfg.geometry);mat=cfg.material;D=plane_strain_D(mat)
    u0=np.zeros(mesh.ndof);ep=np.zeros((3,mesh.ne));rho=np.zeros(mesh.ne)
    d=np.zeros(mesh.nn);d[bnd.notch_nodes]=1.0
    K,R,*_=assemble_mechanics(mesh,u0,ep,rho,d,D,mat)
    htip=float(mesh.hbar_tip or mesh.hbar)
    cache={}

    def evaluate(alpha):
        key=round(float(alpha),10)
        if key in cache:return dict(cache[key])
        ctx=MixedModeContext(float(alpha),target_mode_phase_deg=0.0,solver_seed=1)
        solve=_mixed_solve_factory(solve_dirichlet,ctx)
        u,F=solve(K,R,u0.copy(),bnd,0.5*a.U_cal_m,-0.5*a.U_cal_m)
        _K2,_R2,sig,*_=assemble_mechanics(mesh,u,ep,rho,d,D,mat)
        md=project_near_tip_modes(mesh,sig,d,np.array([cfg.geometry.a0,0.0]),np.array([1.0,0.0]),
                                  1.25*htip,6.0*htip,a.mode_projection_angle_deg,
                                  a.mode_projection_damage_cutoff)
        KI=md['KI_Pa_sqrt_m'];KII=md['KII_Pa_sqrt_m']
        row={'loading_angle_deg':float(alpha),'achieved_psi_raw_deg':md['mode_phase_deg'],
             'KI_raw_MPa_sqrt_m':KI/1e6,'KII_raw_MPa_sqrt_m':KII/1e6,
             'KI_MPa_sqrt_m':KI/1e6,'KII_MPa_sqrt_m':KII/1e6,
             'mode_ratio_KII_over_KI':KII/KI if np.isfinite(KI) and abs(KI)>1e-30 else np.nan,
             'generalized_reaction_N':F,'projection_n':md['mode_projection_n'],
             'projection_fit_count':md['mode_projection_fit_count'],
             'projection_rel_rmse':md['mode_projection_rel_rmse'],
             'projection_condition':md['mode_projection_condition'],
             'projection_psi_spread_deg':md['mode_projection_psi_spread_deg'],
             'projection_K_spread_frac':md['mode_projection_K_spread_frac'],
             'amplitude_projection_reliable':md['mode_projection_reliable']}
        cache[key]=dict(row);return row

    # Two measured basis responses retain geometry-induced cross-coupling.
    # Normalize the Williams sign convention so imposed opening has KI > 0
    # and imposed positive sliding has KII > 0.
    opening_raw=evaluate(0.0)
    sliding_raw=evaluate(90.0)
    M_raw=np.array([[opening_raw['KI_raw_MPa_sqrt_m'],sliding_raw['KI_raw_MPa_sqrt_m']],
                    [opening_raw['KII_raw_MPa_sqrt_m'],sliding_raw['KII_raw_MPa_sqrt_m']]],float)
    mode_signs=mode_signs_from_basis(M_raw)

    def normalized(row):
        out=dict(row)
        KI,KII,psi=apply_mode_signs(out['KI_raw_MPa_sqrt_m'],out['KII_raw_MPa_sqrt_m'],mode_signs)
        out['KI_MPa_sqrt_m']=KI;out['KII_MPa_sqrt_m']=KII
        out['achieved_psi_deg']=psi
        out['mode_sign_I']=float(mode_signs[0]);out['mode_sign_II']=float(mode_signs[1])
        out['mode_ratio_KII_over_KI']=KII/KI if np.isfinite(KI) and abs(KI)>1e-30 else np.nan
        usable,reasons=phase_projection_gate(out,max_condition=a.projection_condition_max,
                                              max_phase_spread_deg=a.phase_spread_tol_deg)
        out['phase_projection_usable']=usable
        out['phase_projection_reasons']=';'.join(reasons)
        return out

    opening=normalized(opening_raw)
    sliding=normalized(sliding_raw)
    M=np.diag(mode_signs) @ M_raw
    basis_cond=float(np.linalg.cond(M))
    if not np.all(np.isfinite(M)) or not np.isfinite(basis_cond) or basis_cond>a.basis_condition_max:
        raise SystemExit(f'mode-basis calibration matrix invalid: cond={basis_cond:.6g}, M={M.tolist()}')

    grid=[normalized(evaluate(x)) for x in vals(a.diagnostic_alpha_grid_deg)]
    mappings=[]
    for target in vals(a.target_psi_deg):
        try:
            alpha=loading_angle_from_mode_basis(M,target,max_abs_alpha_deg=89.9)
        except Exception as exc:
            mappings.append({'target_psi_deg':target,'loading_angle_deg':np.nan,
                             'achieved_psi_deg':np.nan,'psi_error_deg':np.nan,
                             'phase_converged':False,'phase_projection_usable':False,
                             'phase_projection_reasons':f'basis_error:{exc}',
                             'basis_condition':basis_cond,
                             'calibration_model':'two_basis_signed_phase_ratio'})
            continue
        hist=[]
        for it in range(max(1,a.max_refine_iters)):
            r=normalized(evaluate(alpha));err=angle_error_deg(r['achieved_psi_deg'],target)
            usable=bool(r['phase_projection_usable'])
            hist.append((alpha,err,usable,r))
            if usable and abs(err)<=a.psi_tol_deg:break
            da=max(float(a.refine_dalpha_deg),0.1)
            rp=normalized(evaluate(min(alpha+da,89.9)));rm=normalized(evaluate(max(alpha-da,-89.9)))
            denom=rp['loading_angle_deg']-rm['loading_angle_deg']
            slope=angle_error_deg(rp['achieved_psi_deg'],rm['achieved_psi_deg'])/max(denom,1e-12)
            if not np.isfinite(slope) or abs(slope)<0.05:
                break
            alpha=float(np.clip(alpha-err/slope,-89.9,89.9))
        best=min(hist,key=lambda z:(not z[2],abs(z[1]),float(z[3].get('projection_psi_spread_deg',np.inf))))
        alpha,err,usable,r=best
        row=dict(r)
        row.update(target_psi_deg=target,psi_error_deg=err,
                   phase_converged=bool(usable and abs(err)<=a.psi_tol_deg),
                   calibration_iterations=len(hist),basis_condition=basis_cond,
                   basis_KI_open=opening['KI_MPa_sqrt_m'],basis_KII_open=opening['KII_MPa_sqrt_m'],
                   basis_KI_shear=sliding['KI_MPa_sqrt_m'],basis_KII_shear=sliding['KII_MPa_sqrt_m'],
                   calibration_model='two_basis_signed_phase_ratio')
        mappings.append(row)

    out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
    datasets=[('elastic_alpha_grid_v3_3.csv',grid),('mixed_mode_loading_calibration_v3_3.csv',mappings)]
    for name,rows in datasets:
        fields=[]
        for row in rows:
            for key in row:
                if key not in fields:fields.append(key)
        with (out/name).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
    payload={'model':MODEL_ID,'settings':vars(a),'basis_matrix_raw':M_raw.tolist(),'basis_matrix':M.tolist(),'mode_signs':mode_signs.tolist(),
             'basis_condition':basis_cond,'opening_basis':opening,'sliding_basis':sliding,
             'mapping':mappings}
    (out/'mixed_mode_loading_calibration_v3_3.json').write_text(json.dumps(payload,indent=2))

    print('raw mode basis matrix [MPa sqrt(m) per calibration amplitude]:')
    print(M_raw)
    print('mode sign normalization [KI, KII]:',mode_signs)
    print('normalized mode basis matrix:')
    print(M)
    print(f'basis condition: {basis_cond:.6g}')
    for r in mappings:
        print({'target':r.get('target_psi_deg'),'alpha':r.get('loading_angle_deg'),
               'psi':r.get('achieved_psi_deg'),'error':r.get('psi_error_deg'),
               'phase_ok':r.get('phase_projection_usable'),
               'amplitude_fit_ok':r.get('amplitude_projection_reliable'),
               'phase_spread':r.get('projection_psi_spread_deg'),
               'reasons':r.get('phase_projection_reasons')})
    print('wrote',out/'mixed_mode_loading_calibration_v3_3.csv')
    bad=[r for r in mappings if not bool(r.get('phase_converged'))]
    if bad:
        detail=', '.join(f"{r.get('target_psi_deg')}:err={r.get('psi_error_deg')},reason={r.get('phase_projection_reasons')}" for r in bad)
        raise SystemExit('phase-ratio calibration failed: '+detail)

if __name__=='__main__':main()
