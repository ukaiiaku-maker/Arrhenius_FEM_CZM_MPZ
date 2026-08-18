#!/usr/bin/env python3
"""Fast 100-um frozen-topology gate using archived DBTT event geometry.

The production archive retains each developed event length and the accepted step
history retains the corresponding projected x increment. Their ratio recovers
|theta_i| exactly. The transverse sign is the only missing historical datum;
this control chooses the sign that keeps the path closest to y=0 and falls back
to the opposite sign only if the unchanged production corridor rejects it.
"""
from __future__ import annotations
import csv, copy, json, math
from pathlib import Path
import numpy as np

import audit_dbtt1000_frozen_czm_topology as full
from arrhenius_fracture.crack_backend import SharpWakeBackend
from arrhenius_fracture.atomic_path_corridor_adaptive_quality_v10051839 import AdaptiveQualityAtomicPathCorridorCZMBackendV10051839
from arrhenius_fracture.atomic_path_corridor_local_scale_v10051839 import install as install_local_scale_gate, restore as restore_local_scale_gate
from arrhenius_fracture.atomic_path_corridor_czm_v10051839 import reset_audit

OUT=Path('dbtt1000_frozen_czm_topology_short100')
TARGETS=(0.0,25.0,50.0,99.0)
# Projected x increments from the exact 17 accepted FEM/CZM DBTT event rows.
DX=np.asarray([
2.2973400956249022e-06,2.2973400956249022e-06,3.2768833392917687e-06,
2.2949469871753224e-06,2.7508744695670037e-06,2.411563032448153e-06,
3.3412766523120792e-06,2.297340095624909e-06,3.426681529279359e-06,
7.567019359322075e-06,2.2906946240722984e-06,1.7861227442407224e-05,
2.8267031959075867e-06,1.8376593229120627e-05,4.7689471261807935e-06,
3.5103425910797305e-06,1.8344690041444448e-05],float)

def add(rows,name,ext,result,failed=0):
    rows.append({'representation':name,'extension_um':ext*1e6,'failed_interfaces':int(failed),**result})

def attempt_signed_event(adaptive,mesh,bnd,d,u,p0,dx,L):
    dy_abs=math.sqrt(max(float(L*L-dx*dx),0.0))
    if dy_abs <= 1e-15:
        signs=(0.0,)
    else:
        preferred=-1.0 if p0[1] > 0 else 1.0
        signs=(preferred,-preferred)
    failures=[]
    for sign in signs:
        dy=0.0 if sign==0.0 else sign*dy_abs
        p1=p0+np.array([float(dx),float(dy)])
        direction=(p1-p0)/float(L)
        # A failed advance is transactional. For a successful trial we keep the
        # returned state, so no physics gate is weakened or bypassed.
        rr=adaptive.advance(mesh=mesh,boundary=bnd,damage=d,displacement=u,
            p0=p0,p1=p1,direction=direction,front_id=0)
        if rr.inserted:
            return rr,p1,float(sign),failures
        failures.append({'sign':float(sign),'reason':rr.reason})
    return None,None,None,failures

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    ledger=json.loads(full.EVENT_LEDGER.read_text())
    lengths=np.asarray(ledger['event_lengths_m'],float)[:17]
    if len(lengths)!=len(DX): raise RuntimeError('17-event geometry ledger mismatch')
    cosines=DX/lengths
    if np.min(cosines)<0.9926003082124 or np.max(cosines)>1.0+1e-10:
        raise RuntimeError(f'archived direction reconstruction failed: min/max={cosines.min()}/{cosines.max()}')
    cfg_a,mesh_a,bnd_a,d_a,u_a=full.initial_state()
    cfg_s,mesh_s,bnd_s,d_s,u_s=full.initial_state()
    startup={'nodes':int(mesh_a.nn),'triangles':int(mesh_a.ne),'hbar_tip_m':float(mesh_a.hbar_tip)}
    if startup['nodes']!=20321 or startup['triangles']!=40502:
        raise RuntimeError(f'production corridor mismatch: {startup}')
    reset_audit(); saved=install_local_scale_gate()
    adaptive=AdaptiveQualityAtomicPathCorridorCZMBackendV10051839(geom=cfg_a.geometry,event_damage=1.0,max_angle_error_deg=35.0)
    sharp=SharpWakeBackend(); rows=[]; path=[]
    add(rows,'adaptive_czm_failed',0,full.frozen_solve(mesh_a,bnd_a,d_a,adaptive.cohesive_network,cfg_a))
    add(rows,'adaptive_split_no_network',0,full.frozen_solve(mesh_a,bnd_a,d_a,None,cfg_a))
    add(rows,'sharp_wake',0,full.frozen_solve(mesh_s,bnd_s,d_s,None,cfg_s))
    p_a=np.array([full.A0_M,0.0]); p_s=p_a.copy(); ext=0.0; pending=list(TARGETS[1:]); failure=None
    try:
        for iev,(L,dx) in enumerate(zip(lengths,DX),1):
            prev=ext*1e6
            rr,p1,sign,trial_failures=attempt_signed_event(adaptive,mesh_a,bnd_a,d_a,u_a,p_a,float(dx),float(L))
            if rr is None:
                failure={'event':iev,'extension_before_um':prev,'requested_length_um':float(L*1e6),
                         'dx_um':float(dx*1e6),'angle_abs_deg':float(np.degrees(np.arccos(np.clip(dx/L,-1,1)))),
                         'sign_trials':trial_failures}; break
            mesh_a,bnd_a,d_a,u_a=rr.mesh,rr.boundary,rr.damage,rr.displacement
            p_a=p1; ext=float(p_a[0]-full.A0_M)
            path.append({'event':iev,'dx_um':float(dx*1e6),'ds_um':float(L*1e6),
                'angle_abs_deg':float(np.degrees(np.arccos(np.clip(dx/L,-1,1)))),
                'selected_sign':sign,'x_um':float((p_a[0]-full.A0_M)*1e6),'y_um':float(p_a[1]*1e6),
                'opposite_sign_failed_first':bool(trial_failures)})
            # Sharp-wake uses the same reconstructed physical segment.
            dy=float(p_a[1]-p_s[1]); p1s=p_s+np.array([float(dx),dy])
            direction_s=(p1s-p_s)/float(L)
            rs=sharp.advance(mesh=mesh_s,boundary=bnd_s,damage=d_s,displacement=u_s,p0=p_s,p1=p1s,
                direction=direction_s,front_id=0,kill_r=max(mesh_s.hbar_tip,0.5e-6))
            if not rs.inserted: raise RuntimeError(f'sharp wake failed event {iev}: {rs.reason}')
            mesh_s,bnd_s,d_s,u_s=rs.mesh,rs.boundary,rs.damage,rs.displacement; p_s=p1s
            for cp in [x for x in pending if prev<x<=ext*1e6+1e-9]:
                add(rows,'adaptive_czm_failed',ext,full.frozen_solve(mesh_a,bnd_a,d_a,adaptive.cohesive_network,cfg_a),adaptive.cohesive_network.failed_count())
                add(rows,'adaptive_split_no_network',ext,full.frozen_solve(mesh_a,bnd_a,d_a,None,cfg_a),adaptive.cohesive_network.failed_count())
                add(rows,'sharp_wake',ext,full.frozen_solve(mesh_s,bnd_s,d_s,None,cfg_s)); pending.remove(cp)
        # Always include attained final state if successful.
        if failure is None:
            last_ad=[r for r in rows if r['representation']=='adaptive_czm_failed'][-1]
            if abs(last_ad['extension_um']-ext*1e6)>1e-8:
                add(rows,'adaptive_czm_failed',ext,full.frozen_solve(mesh_a,bnd_a,d_a,adaptive.cohesive_network,cfg_a),adaptive.cohesive_network.failed_count())
                add(rows,'adaptive_split_no_network',ext,full.frozen_solve(mesh_a,bnd_a,d_a,None,cfg_a),adaptive.cohesive_network.failed_count())
                add(rows,'sharp_wake',ext,full.frozen_solve(mesh_s,bnd_s,d_s,None,cfg_s))
    finally:
        restore_local_scale_gate(saved)
    initial={}
    for r in rows: initial.setdefault(r['representation'],r['compliance_m2_per_N'])
    for r in rows: r['normalized_compliance']=r['compliance_m2_per_N']/initial[r['representation']]
    A=[r for r in rows if r['representation']=='adaptive_czm_failed']; N=[r for r in rows if r['representation']=='adaptive_split_no_network']
    identity=[abs(a['abs_Ftop_N_per_m']-n['abs_Ftop_N_per_m'])/max(a['abs_Ftop_N_per_m'],n['abs_Ftop_N_per_m'],1e-300) for a,n in zip(A,N)]
    with (OUT/'compliance.csv').open('w',newline='') as fp:
        w=csv.DictWriter(fp,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    (OUT/'reconstructed_path.json').write_text(json.dumps(path,indent=2))
    summary={'schema':'dbtt1000_frozen_topology_short100_v2','startup':startup,'events_requested':17,
        'projected_extension_target_um':float(DX.sum()*1e6),'developed_extension_um':float(lengths.sum()*1e6),
        'minimum_reconstructed_forward_cosine':float(cosines.min()),'failure':failure,
        'max_failed_network_identity_rel_diff':max(identity),'rows':rows,'path':path}
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))
if __name__=='__main__': main()
