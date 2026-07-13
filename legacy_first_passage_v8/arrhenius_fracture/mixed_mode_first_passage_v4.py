"""Anisotropic traction-controlled mixed-mode first-passage wrapper.

The active sharp-front solver is retained.  This wrapper changes only:
  * boundary loading: combined opening/sliding displacement;
  * first-passage driving stresses: direct anisotropic FEM tractions sampled at
    a fixed physical process-zone radius;
  * deterministic H=1 first-passage threshold and audit output.

No isotropic KI/KII decomposition is used for kinetics.  The domain J integral
remains the authoritative energy-release audit.  Mode mixity is defined by the
reference-plane process-zone traction phase

    psi_sigma = atan2(tau_tn, sigma_nn).

The candidate cleavage direction and surface-energy weighting are supplied by
existing crystal.py functions in the active project.
"""
from __future__ import annotations

import argparse, csv, json, math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

MODEL_ID = "FEM_CZM_mixed_mode_first_passage_v4_anisotropic_traction_controlled"


def _unit(v):
    v=np.asarray(v,float); n=float(np.linalg.norm(v))
    return np.array([1.0,0.0]) if n<1e-30 else v/n


def angle_error_deg(value,target):
    return (float(value)-float(target)+180.0)%360.0-180.0


def traction_phase_deg(sigma_nn,tau_tn,shear_sign=1.0):
    return math.degrees(math.atan2(float(shear_sign)*float(tau_tn),float(sigma_nn)))


def shear_sign_from_basis(response_raw, min_diagonal=1e-12):
    M=np.asarray(response_raw,float)
    if M.shape!=(2,2) or not np.all(np.isfinite(M)):
        raise ValueError('response_raw must be a finite 2x2 matrix')
    scale=max(float(np.max(np.abs(M))),1.0)
    if M[0,0] <= min_diagonal*scale:
        raise ValueError('opening basis does not produce positive opening traction')
    if abs(M[1,1]) <= min_diagonal*scale:
        raise ValueError('sliding basis has negligible shear response')
    return 1.0 if M[1,1]>=0 else -1.0


def loading_angle_from_response_basis(response_matrix,target_phase_deg,max_abs_alpha_deg=89.9):
    M=np.asarray(response_matrix,float)
    if M.shape!=(2,2) or not np.all(np.isfinite(M)):
        raise ValueError('response_matrix must be a finite 2x2 matrix')
    cond=float(np.linalg.cond(M))
    if not np.isfinite(cond) or cond>1e12:
        raise ValueError(f'response matrix is singular or ill-conditioned: cond={cond:.6g}')
    p=math.radians(float(target_phase_deg))
    q=np.linalg.solve(M,np.array([math.cos(p),math.sin(p)]))
    if not np.all(np.isfinite(q)) or np.linalg.norm(q)<1e-30:
        raise ValueError('invalid boundary vector')
    alpha=math.degrees(math.atan2(float(q[1]),float(q[0])))
    if abs(alpha)>float(max_abs_alpha_deg):
        raise ValueError(f'required loading angle {alpha:.6g} exceeds {max_abs_alpha_deg:.6g}')
    return float(alpha)


def energy_matrix_from_basis(J_open,J_slide,J_equal,U_cal=1.0):
    """Quadratic J matrix for q=[opening,sliding] from 0,90,45 degree solves."""
    u2=max(float(U_cal)**2,1e-300)
    g11=float(J_open)/u2; g22=float(J_slide)/u2
    # q=(1/sqrt(2),1/sqrt(2)) for the equal-amplitude solve.
    g12=float(J_equal)/u2 - 0.5*(g11+g22)
    G=np.array([[g11,g12],[g12,g22]],float)
    return G


def energy_phase_from_matrix(G,alpha_deg):
    """Energy-orthonormal generalized phase for audit only."""
    G=np.asarray(G,float)
    q=np.array([math.cos(math.radians(alpha_deg)),math.sin(math.radians(alpha_deg))])
    # Eigen square root is robust to small roundoff asymmetry.
    w,V=np.linalg.eigh(0.5*(G+G.T))
    if np.min(w)<=0: return float('nan')
    L=np.diag(np.sqrt(w)) @ V.T
    k=L@q
    return math.degrees(math.atan2(float(k[1]),float(k[0])))


def _element_centroids(mesh):
    return np.asarray(mesh.nodes,float)[np.asarray(mesh.elems,int)].mean(axis=1)


def process_zone_traction_probe(mesh,sigma_gp,d,tip,crack_direction,radius_m,
                                annulus_half_width=0.45,sector_half_angle_deg=40.0,
                                damage_cutoff=0.85,min_elements=4):
    """Area-weighted anisotropic FEM traction at a fixed physical tip radius.

    Stresses are averaged in a forward annular sector and transformed into the
    crack-local frame (t,n).  This is a finite-radius process-zone quantity,
    deliberately not an isotropic Williams extrapolation.
    """
    tip=np.asarray(tip,float); t=_unit(crack_direction); n=np.array([-t[1],t[0]])
    c=_element_centroids(mesh); dx=c-tip
    x1=dx@t; x2=dx@n; rr=np.sqrt(x1*x1+x2*x2)
    th=np.degrees(np.arctan2(x2,x1))
    r0=max(float(radius_m),1e-30); hw=float(np.clip(annulus_half_width,0.05,0.95))
    dg=np.asarray(d,float)[np.asarray(mesh.elems,int)].mean(axis=1)
    sel=(rr>=(1-hw)*r0)&(rr<=(1+hw)*r0)&(np.abs(th)<=float(sector_half_angle_deg))&(dg<float(damage_cutoff))
    # Controlled expansion, preserving forward sampling, rather than silently
    # falling back to one singular element.
    expansion=1.0
    for expansion in (1.0,1.5,2.0,3.0):
        sel=(rr>=max((1-hw)*r0/expansion,0.25*r0))&(rr<=(1+hw)*r0*expansion)&(np.abs(th)<=min(float(sector_half_angle_deg)*expansion,85.0))&(dg<float(damage_cutoff))
        if int(np.count_nonzero(sel))>=int(min_elements): break
    idx=np.flatnonzero(sel)
    if len(idx)<int(min_elements):
        return {'reliable':False,'n_elements':int(len(idx)),'probe_radius_m':r0,'expansion':float(expansion)}
    area=np.maximum(np.asarray(mesh.area_e,float)[idx],1e-30); w=area/area.sum()
    sig=np.asarray(sigma_gp,float)
    sxx=float(w@sig[0,idx]); syy=float(w@sig[1,idx]); sxy=float(w@sig[2,idx])
    S=np.array([[sxx,sxy],[sxy,syy]])
    stt=float(t@S@t); snn=float(n@S@n); ttn=float(t@S@n)
    eig=np.linalg.eigvalsh(S)
    return {'reliable':True,'n_elements':int(len(idx)),'probe_radius_m':r0,
            'expansion':float(expansion),'sigma_tt_Pa':stt,'sigma_nn_Pa':snn,
            'tau_tn_Pa':ttn,'sigma1_Pa':float(eig[-1]),'sigma2_Pa':float(eig[0]),
            'sigma_xx_Pa':sxx,'sigma_yy_Pa':syy,'sigma_xy_Pa':sxy,'stress_tensor':S}


def max_resolved_slip_shear(S,theta_deg):
    from .crystal import bcc_slip_traces
    best=None
    for p in bcc_slip_traces(theta_deg):
        tau=float(np.asarray(p['t'])@S@np.asarray(p['n']))
        row={'name':p['name'],'tau_Pa':tau,'abs_tau_Pa':abs(tau),'angle_deg':float(p['angle_deg'])}
        if best is None or row['abs_tau_Pa']>best['abs_tau_Pa']: best=row
    return best or {'name':'none','tau_Pa':0.0,'abs_tau_Pa':0.0,'angle_deg':float('nan')}


@dataclass
class AnisotropicContext:
    loading_angle_deg: float
    target_traction_phase_deg: float
    crystal_theta_deg: float=0.0
    cleavage_gamma_aniso: float=0.3
    probe_radius_m: float=10e-6
    annulus_half_width: float=0.45
    sector_half_angle_deg: float=40.0
    damage_cutoff: float=0.85
    shear_sign: float=1.0
    solver_seed: int=1
    records: list[dict[str,Any]]=field(default_factory=list)
    latest: dict[str,Any]=field(default_factory=dict)
    last_loading: dict[str,float]=field(default_factory=dict)
    def __post_init__(self): self.rng=np.random.default_rng(int(self.solver_seed))
    @property
    def alpha_rad(self): return math.radians(float(self.loading_angle_deg))


class DirectTractionEngineMixin:
    """Use direct anisotropic process-zone tractions for both Arrhenius hazards."""
    def _mm_init(self,context):
        self._mm=context; self._mm_threshold=1.0; self._prev_sigma_cleave=None
    def _drives(self):
        r=self._mm.latest
        if not r or not bool(r.get('traction_probe_reliable',False)):
            raise RuntimeError('anisotropic process-zone traction probe unavailable')
        return max(float(r['sigma_cleave_drive_Pa']),0.0),max(float(r['sigma_emit_drive_Pa']),0.0),r
    def predict_clock_increment(self,K,T,dt):
        sig2,_se,_=self._drives(); lam2,_,_=self.lambda_cleave(sig2,T)
        if self._prev_sigma_cleave is not None:
            lam1,_,_=self.lambda_cleave(max(self._prev_sigma_cleave,0.0),T)
            lo,hi=sorted((max(lam1,0.0),max(lam2,0.0)))
            leff=0.5*hi if lo<=0 else (hi if abs(hi-lo)<=1e-12*max(hi,1e-300) else (hi-lo)/math.log(hi/lo))
        else: leff=max(lam2,0.0)
        return float(max(leff*float(dt),0.0))
    def step(self,K,T,dt):
        sig_c,sig_e,mm=self._drives()
        if self.f.sigma_cap>0:
            sig_c=min(sig_c,self.f.sigma_cap); sig_e=min(sig_e,self.f.sigma_cap)
        lam_e,sig_em_eff,Ge=self.lambda_emit(sig_e,T)
        prod_raw=lam_e*dt; prod=min(prod_raw,self.f.dN_cap)
        cap=bool(np.isfinite(self.f.dN_cap) and prod_raw>self.f.dN_cap)
        sat=1.0
        if np.isfinite(self.f.N_sat) and self.f.N_sat>0:
            sat=max(1.0-self.N_em/self.f.N_sat,0.0); prod*=sat
        ann=self.f.recover_k*self.N_em*dt
        self.W_emit+=sig_em_eff*self.b*self.f.L_pz*prod
        self.N_em=max(self.N_em+prod-ann,0.0)
        lam_c,lam_raw,Gc=self.lambda_cleave(sig_c,T)
        if self.f.tau_B>0 and dt>0: self.B*=np.exp(-min(dt/self.f.tau_B,80.0))
        if self._prev_sigma_cleave is not None:
            lam1,_,_=self.lambda_cleave(max(self._prev_sigma_cleave,0.0),T)
            lo,hi=sorted((max(lam1,0.0),max(lam_c,0.0)))
            leff=0.5*hi if lo<=0 else (hi if abs(hi-lo)<=1e-12*max(hi,1e-300) else (hi-lo)/math.log(hi/lo))
        else: leff=lam_c
        self.B+=leff*dt; self._prev_sigma_cleave=sig_c; self.t+=dt
        Npre=float(self.N_em); spre=float(self.sigma_back()); rpre=float(self.r_eff()); dGpre=float(self.dG_emb())
        nfire=0
        while self.B>=self._mm_threshold and nfire<100000:
            self.B-=self._mm_threshold; nfire+=1
        fired=nfire>0; retained=Npre; shed=0.0
        if fired:
            retain=float(np.clip(self.f.wake_retain,0,1))**nfire
            retained=Npre*retain; shed=Npre-retained; self.N_em=retained
            self.a_adv+=self.f.da*nfire; self.n_adv+=nfire
        info={'fired':fired,'n_fire':nfire,'v_crack':self.f.da*nfire/dt if dt>0 else 0.0,
              'sigma_tip':sig_c,'sigma_emit_tip':sig_e,'sigma_back':self.sigma_back(),
              'lambda_e':lam_e,'lambda_c':lam_c,'lambda_c_raw':lam_raw,'B':self.B,
              'N_em':self.N_em,'r_eff':self.r_eff(),'dG_emb_eV':self.dG_emb()/1.602176634e-19,
              'G_cleave_eff_eV':Gc/1.602176634e-19,**self.cleavage_diagnostics(sig_c,T),
              'G_emit_eV':Ge/1.602176634e-19,'W_emit':self.W_emit,
              'sigma_tip_uncapped':float(mm['sigma_cleave_drive_Pa']),
              'sigma_cap_active':bool(self.f.sigma_cap>0 and mm['sigma_cleave_drive_Pa']>self.f.sigma_cap),
              'dN_emit_raw':float(prod_raw),'dN_cap_active':cap,'N_sat_factor':sat,
              'N_sat_active':bool(sat<0.999999),'N_em_pre_renewal':Npre,
              'N_em_retained':retained,'N_em_shed_to_wake':shed,
              'sigma_back_pre_renewal':spre,'r_eff_pre_renewal':rpre,
              'dG_emb_pre_renewal_eV':dGpre/1.602176634e-19,
              'anisotropic_reference_phase_deg':mm.get('reference_traction_phase_deg'),
              'anisotropic_candidate_sigma_nn_Pa':mm.get('candidate_sigma_nn_Pa'),
              'anisotropic_candidate_tau_tn_Pa':mm.get('candidate_tau_tn_Pa'),
              'anisotropic_sigma_cleave_drive_Pa':mm.get('sigma_cleave_drive_Pa'),
              'anisotropic_sigma_emit_drive_Pa':mm.get('sigma_emit_drive_Pa'),
              'anisotropic_candidate_angle_deg':mm.get('candidate_angle_deg'),
              'anisotropic_gamma_rel':mm.get('candidate_gamma_rel'),
              'anisotropic_slip_system':mm.get('slip_system_name')}
        return info


def _mixed_solve_factory(original_solve,context):
    def solve_mixed(K,Rint,u,bnd,Uy_top,Uy_bot):
        from scipy.sparse.linalg import spsolve
        alpha=context.alpha_rad; Atotal=float(Uy_top-Uy_bot)
        Un=Atotal*math.cos(alpha); Us=Atotal*math.sin(alpha)
        u_open,_=original_solve(K,Rint,u,bnd,0.5*Un,-0.5*Un)
        Kc=K.tocsr(); Ropen=Rint+Kc@(u_open-u)
        if abs(Us)<=1e-30:
            unew=u_open; Rfull=Ropen
        else:
            ndof=len(u); prescribed=np.zeros(ndof,bool); target=u_open.copy()
            tn=bnd.top_nodes; bn=bnd.bot_nodes
            prescribed[2*tn]=True; prescribed[2*tn+1]=True; prescribed[2*bn]=True; prescribed[2*bn+1]=True
            target[2*tn]=u_open[2*tn]+0.5*Us; target[2*bn]=u_open[2*bn]-0.5*Us
            target[2*tn+1]=0.5*Un; target[2*bn+1]=-0.5*Un
            free=~prescribed; dup=target[prescribed]-u_open[prescribed]
            rhs=-Ropen[free]-Kc[np.ix_(free,prescribed)]@dup
            unew=u_open.copy(); unew[free]=u_open[free]+spsolve(Kc[np.ix_(free,free)],rhs); unew[prescribed]=target[prescribed]
            Rfull=Ropen+Kc@(unew-u_open)
        Fx=float(np.sum(Rfull[2*bnd.top_nodes])); Fy=float(np.sum(Rfull[2*bnd.top_nodes+1]))
        Fgen=Fx*math.sin(alpha)+Fy*math.cos(alpha)
        context.last_loading={'U_total_m':Atotal,'U_open_m':Un,'U_shear_m':Us,
                              'loading_angle_deg':context.loading_angle_deg,
                              'generalized_reaction_N':Fgen,'reaction_x_N':Fx,'reaction_y_N':Fy}
        return unew,Fgen
    return solve_mixed


def _j_wrapper_factory(original_compute,context):
    def wrapped(mesh,u,sigma_gp,psi_e_gp,d,crack_tip,crack_direction,mat,ell,cfg=None,crack_segments=None,exclude_radius=0.0):
        from .crystal import cubic_cleavage_gamma
        J,KJ,info=original_compute(mesh,u,sigma_gp,psi_e_gp,d,crack_tip,crack_direction,mat,ell,
                                   cfg=cfg,crack_segments=crack_segments,exclude_radius=exclude_radius)
        cand=process_zone_traction_probe(mesh,sigma_gp,d,crack_tip,crack_direction,
                  context.probe_radius_m,context.annulus_half_width,context.sector_half_angle_deg,context.damage_cutoff)
        ref=process_zone_traction_probe(mesh,sigma_gp,d,crack_tip,np.array([1.0,0.0]),
                  context.probe_radius_m,context.annulus_half_width,context.sector_half_angle_deg,context.damage_cutoff)
        reliable=bool(cand.get('reliable',False) and ref.get('reliable',False))
        t=_unit(crack_direction); n=np.array([-t[1],t[0]])
        normal_angle=math.atan2(float(n[1]),float(n[0]))
        gamma=float(cubic_cleavage_gamma(normal_angle,context.crystal_theta_deg,context.cleavage_gamma_aniso))
        S=np.asarray(cand.get('stress_tensor',np.zeros((2,2))),float)
        slip=max_resolved_slip_shear(S,context.crystal_theta_deg) if reliable else {'name':'none','abs_tau_Pa':0.0,'tau_Pa':0.0,'angle_deg':np.nan}
        snn=float(cand.get('sigma_nn_Pa',np.nan)); tau=float(cand.get('tau_tn_Pa',np.nan))
        sigma_c=max(snn,0.0)/math.sqrt(max(gamma,1e-12)) if reliable else float('nan')
        sigma_e=float(slip['abs_tau_Pa']) if reliable else float('nan')
        phase=traction_phase_deg(ref.get('sigma_nn_Pa',np.nan),ref.get('tau_tn_Pa',np.nan),context.shear_sign) if reliable else float('nan')
        md={'model':MODEL_ID,'J_J_per_m2':float(J),'KJ_reference_Pa_sqrt_m':float(KJ),
            'traction_probe_reliable':reliable,'traction_probe_radius_m':context.probe_radius_m,
            'reference_sigma_nn_Pa':ref.get('sigma_nn_Pa'),'reference_tau_tn_Pa':ref.get('tau_tn_Pa'),
            'reference_traction_phase_deg':phase,'target_traction_phase_deg':context.target_traction_phase_deg,
            'reference_phase_error_deg':angle_error_deg(phase,context.target_traction_phase_deg) if np.isfinite(phase) else np.nan,
            'candidate_sigma_nn_Pa':snn,'candidate_tau_tn_Pa':tau,
            'candidate_angle_deg':math.degrees(math.atan2(t[1],t[0])),
            'candidate_gamma_rel':gamma,'sigma_cleave_drive_Pa':sigma_c,'sigma_emit_drive_Pa':sigma_e,
            'slip_system_name':slip['name'],'slip_tau_signed_Pa':slip['tau_Pa'],
            'candidate_probe_n_elements':cand.get('n_elements',0),'reference_probe_n_elements':ref.get('n_elements',0),
            **context.last_loading,'tip_x_m':float(np.asarray(crack_tip)[0]),'tip_y_m':float(np.asarray(crack_tip)[1])}
        info.update(md); context.latest=md.copy(); context.records.append(md.copy())
        return J,KJ,info
    return wrapped


def _engine_factory(original_build,context,base_class):
    class Engine(DirectTractionEngineMixin,base_class):
        def __init__(self,*a,**kw): base_class.__init__(self,*a,**kw); self._mm_init(context)
    def build(args,mat):
        base=original_build(args,mat); return Engine(base.f,base.cb,base.eb,base.G,base.nu,base.b)
    return build


def _write_records(out,context):
    if not context.records:return
    cols=sorted({k for r in context.records for k in r})
    with (out/'anisotropic_traction_calls.csv').open('w',newline='') as fp:
        w=csv.DictWriter(fp,fieldnames=cols);w.writeheader();w.writerows(context.records)


def _summary(out,T,context,base_summary):
    steps=out/f'steps_{int(T):04d}K.csv'; accepted={}
    if steps.exists():
        a=np.genfromtxt(steps,delimiter=',',names=True)
        if np.size(a):
            a=np.atleast_1d(a); fire=np.flatnonzero(a['n_fire']>0); i=int(fire[0]) if len(fire) else len(a)-1
            accepted={n:float(a[n][i]) for n in a.dtype.names}
    U=accepted.get('Uapp_m',np.nan); Kacc=accepted.get('KJ_Pa_sqrtm',np.nan)
    rec=context.records[-1] if context.records else {}
    if context.records and np.isfinite(U):
        du=np.array([abs(float(r.get('U_total_m',np.nan))-U) for r in context.records])
        if np.any(np.isfinite(du)):
            best=float(np.nanmin(du)); tol=max(1e-14,1e-10*max(abs(U),1e-12))
            near=[r for r,q in zip(context.records,du) if np.isfinite(q) and q<=best+tol]
            if near and np.isfinite(Kacc):
                rec=min(near,key=lambda r: abs(float(r.get('KJ_reference_Pa_sqrt_m',np.nan))-Kacc) if np.isfinite(float(r.get('KJ_reference_Pa_sqrt_m',np.nan))) else float('inf'))
            elif near: rec=near[-1]
    b=base_summary[0] if base_summary else {}
    payload={'model':MODEL_ID,'T_K':float(T),'loading_angle_deg':context.loading_angle_deg,
      'target_traction_phase_deg':context.target_traction_phase_deg,
      'traction_phase_first_deg':rec.get('reference_traction_phase_deg'),
      'traction_phase_error_first_deg':rec.get('reference_phase_error_deg'),
      'J_first_J_per_m2':rec.get('J_J_per_m2'),'KJ_reference_first_MPa_sqrt_m':float(rec.get('KJ_reference_Pa_sqrt_m',np.nan))/1e6,
      'reference_sigma_nn_first_GPa':float(rec.get('reference_sigma_nn_Pa',np.nan))/1e9,
      'reference_tau_tn_first_GPa':float(rec.get('reference_tau_tn_Pa',np.nan))/1e9,
      'candidate_sigma_nn_first_GPa':float(rec.get('candidate_sigma_nn_Pa',np.nan))/1e9,
      'candidate_tau_tn_first_GPa':float(rec.get('candidate_tau_tn_Pa',np.nan))/1e9,
      'sigma_cleave_drive_first_GPa':float(rec.get('sigma_cleave_drive_Pa',np.nan))/1e9,
      'sigma_emit_drive_first_GPa':float(rec.get('sigma_emit_drive_Pa',np.nan))/1e9,
      'candidate_angle_first_deg':rec.get('candidate_angle_deg'),'candidate_gamma_rel':rec.get('candidate_gamma_rel'),
      'slip_system_first':rec.get('slip_system_name'),'traction_probe_reliable':rec.get('traction_probe_reliable'),
      'control_state':'first_passage' if accepted.get('n_fire',0)>0 else 'right_censored_endpoint',
      'Kc_first_existing_MPa_sqrt_m':b.get('Kc_first_MPa_sqrt_m'),'N_em_final':b.get('N_em_final'),'mode_classification':b.get('mode'),
      'crystal_aniso':True,'crystal_theta_deg':context.crystal_theta_deg,'cleavage_gamma_aniso':context.cleavage_gamma_aniso,
      'phase_definition':'process_zone_traction_at_fixed_physical_radius',
      'event_phase_within_2deg':bool(np.isfinite(float(rec.get('reference_phase_error_deg',np.nan))) and abs(float(rec.get('reference_phase_error_deg',np.nan)))<=2.0)}
    (out/'anisotropic_mixed_mode_first_passage_summary.json').write_text(json.dumps(payload,indent=2,default=str))
    with (out/'anisotropic_mixed_mode_first_passage_summary.csv').open('w',newline='') as fp:
        w=csv.DictWriter(fp,fieldnames=list(payload));w.writeheader();w.writerow(payload)
    return payload


def parser():
    p=argparse.ArgumentParser(add_help=False,allow_abbrev=False)
    p.add_argument('--mixity-loading-angle-deg',type=float,default=0.0)
    p.add_argument('--target-traction-phase-deg',type=float,required=True)
    p.add_argument('--traction-shear-sign',type=float,default=1.0)
    p.add_argument('--traction-probe-radius-m',type=float,default=10e-6)
    p.add_argument('--traction-annulus-half-width',type=float,default=0.45)
    p.add_argument('--traction-sector-half-angle-deg',type=float,default=40.0)
    p.add_argument('--traction-damage-cutoff',type=float,default=0.85)
    p.add_argument('--solver-seed',type=int,default=1)
    return p


def main(argv=None):
    from . import sharp_front as sf
    from . import fem as femmod
    from . import j_integral as jimod
    mm,remaining=parser().parse_known_args(argv); args=sf._build_parser().parse_args(remaining)
    if args.mode!='2d': raise SystemExit('v4 anisotropic mixed mode requires --mode 2d')
    if not bool(getattr(args,'crystal_aniso',False)): raise SystemExit('v4 requires --crystal-aniso')
    if not bool(getattr(args,'crystal_compete',False)): raise SystemExit('v4 requires --crystal-compete')
    if bool(getattr(args,'crystal_branch',False)) or int(getattr(args,'max_fronts',1))!=1:
        raise SystemExit('v4 first-passage screen requires branching off and --max-fronts 1')
    context=AnisotropicContext(mm.mixity_loading_angle_deg,mm.target_traction_phase_deg,
      float(getattr(args,'crystal_theta_deg',0.0) or 0.0),
      float(0.3 if getattr(args,'cleave_gamma_aniso',None) is None else getattr(args,'cleave_gamma_aniso')),mm.traction_probe_radius_m,
      mm.traction_annulus_half_width,mm.traction_sector_half_angle_deg,mm.traction_damage_cutoff,
      mm.traction_shear_sign,mm.solver_seed)
    out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    (out/'anisotropic_mixed_mode_run_config.json').write_text(json.dumps({'model':MODEL_ID,**vars(mm),
       'crystal_theta_deg':context.crystal_theta_deg,'cleavage_gamma_aniso':context.cleavage_gamma_aniso,
       'note':'anisotropic FEM + crystal-resolved direct process-zone tractions; J is energy audit only'},indent=2))
    osolve=femmod.solve_dirichlet; oJ=jimod.compute_J_integral; obuild=sf.build_engine
    try:
        femmod.solve_dirichlet=_mixed_solve_factory(osolve,context)
        jimod.compute_J_integral=_j_wrapper_factory(oJ,context)
        sf.build_engine=_engine_factory(obuild,context,sf.FrontEngine)
        base=sf.run_2d(args)
    finally:
        femmod.solve_dirichlet=osolve; jimod.compute_J_integral=oJ; sf.build_engine=obuild
    _write_records(out,context)
    vals=[]
    for T in args.temperatures: vals.append(_summary(out,T,context,base))
    print('MIXED_MODE_V4 complete:',json.dumps(vals,indent=2,default=str));return vals

if __name__=='__main__': main()
