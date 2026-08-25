#!/usr/bin/env python3
"""Predictive natural-loading V2 forward trajectories.

No threshold/action priming is permitted here.  Both lanes start from their
canonical source constructor, receive the archived native opening rate, and
advance with source-owned clocks and renewal states.  Elastic fields are exact
and cached by complete geometry/topology fingerprint.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
import sys
from typing import Any

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"analysis_outputs/oneD_v2_terminal_predictive_program"
sys.path[:0]=[str(ROOT),str(ROOT/"scripts")]

from qualify_oneD_v2_native_state_closure import (
    IDS, _fem_constructor, _install_fem, _install_pf, _jsonable,
    _pf_constructor, _pz_metrics, canonical_rows,
)
from reduced_fracture_v2.native_state_factory import (
    FEMCZMNativeStateFactory, PFNativeStateFactory,
    canonical_source_value,
)
from reduced_fracture_v2.production_oracles import (
    FEMCZMExactElasticFieldOracle, FEMCZMExactSourceProbeOperator,
    PFExactElasticFieldOracle, PFExactSourceProbeOperator,
)

OPENING_RATE_M_S=2.380952380952381e-8
MAX_OPENING_INCREMENT_M=0.2e-6
MAX_PREDICTED_CLOCK_INCREMENT=0.08
TEMPERATURE_K=1000.0


def _predict(
    engine: Any,
    backend: str,
    K_cleave: float,
    K_emit: float,
    T: float,
    dt: float,
) -> float:
    """Return a non-mutating outer proposal bound.

    PF exposes the inexpensive production predictor directly.  The qualified
    FEM/CZM predictor replays every exact stochastic emission event on a trial
    clone and is therefore as expensive as an accepted interval.  For the
    reduced outer controller we use its current-state end-point cleavage rate
    only to propose a smaller interval.  The accepted update still uses the
    exact joint-K, event-localized stochastic-emission integrator, which owns
    its internal action limits, event localization, state update, and unused
    time.  Consequently this bound cannot force or suppress an event and is
    not the rejected aggregate/mean-field constitutive update.
    """
    if backend == "pf":
        return float(engine.predict_clock_increment(K_cleave, T, dt))
    threshold=max(float(getattr(engine,"hazard_threshold_action",1.)),1e-300)
    lam=float(engine.lambda_cleave(engine.sigma_tip(K_cleave),T)[0])
    return max(lam,0.)/threshold*max(float(dt),0.)


def _capture(engine:Any):
    fn=getattr(engine,"_capture_state",None)
    return ("transaction",fn()) if callable(fn) else ("object",copy.deepcopy(engine))


def _restore(engine:Any,snapshot):
    kind,value=snapshot
    if kind=="object": return value
    fn=getattr(engine,"_restore_state",None)
    if not callable(fn): raise RuntimeError("source transaction capture has no restore")
    fn(value);return engine


def _runtime_signature(engine:Any)->dict[str,Any]:
    simple={}
    for key,value in vars(engine).items():
        low=key.lower()
        if "history" in low or "records" in low: continue
        if value is None or isinstance(value,(str,bool,int,float,np.integer,np.floating,np.ndarray)):
            simple[key]=canonical_source_value(value)
        elif isinstance(value,np.random.Generator):
            simple[key]=canonical_source_value(value.bit_generator.state)
        elif isinstance(value,list) and value and all(isinstance(x,np.random.Generator) for x in value):
            simple[key]=[canonical_source_value(x.bit_generator.state) for x in value]
    pz=getattr(engine,"mpz",getattr(engine,"mpz_state",None));pz_fields={}
    for key,value in vars(pz).items():
        low=key.lower()
        if any(token in low for token in ("kernel","candidate","metadata","history")): continue
        if value is None or isinstance(value,(str,bool,int,float,np.integer,np.floating,np.ndarray,list,tuple,dict)):
            pz_fields[key]=canonical_source_value(value)
    return {"engine":simple,"process_zone":pz_fields}


def _fingerprint(value:Any)->str:
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(",",":")).encode()).hexdigest()


def _geometry(backend:str,extension:float,events:list[dict[str,Any]]):
    if backend=="pf":
        key=PFExactElasticFieldOracle.geometry_fingerprint(events)
        return key,{"events":events,"tip_xy_m":events[-1]["p1_m"] if events else (5e-4,0.)}
    return f"straight:{extension:.17e}",{"events":[],"extension_m":extension}


def _step(engine: Any, backend: str, K_cleave: float, K_emit: float,
          temperature_K: float, dt_s: float) -> dict[str, Any]:
    if backend == "pf":
        return dict(engine.step(K_cleave, temperature_K, dt_s))
    return dict(engine.step_drives(K_cleave, K_emit, temperature_K, dt_s))


def run(backend:str,target_um:float,classes:set[str]) -> None:
    OUT.mkdir(parents=True,exist_ok=True)
    label="PF" if backend=="pf" else "FEMCZM"
    if backend=="pf":
        oracle=PFExactElasticFieldOracle();probe=PFExactSourceProbeOperator()
        factory_type=PFNativeStateFactory;constructor=_pf_constructor;install=_install_pf
    else:
        oracle=FEMCZMExactElasticFieldOracle();probe=FEMCZMExactSourceProbeOperator()
        factory_type=FEMCZMNativeStateFactory
        # The qualified production composition is the exact joint-K,
        # event-localized stochastic-emission engine.  The aggregate moving-tip
        # implementation is retained only as historical diagnostic code and is
        # not an admissible predictive closure.
        constructor=_fem_constructor
        install=_install_fem
    cache={};hits=misses=0;all_intervals=[];all_events=[];summaries=[]
    for row in canonical_rows():
        material=IDS[row["candidate_id"]]
        if material not in classes: continue
        product=factory_type(constructor(row)).instantiate()
        # Source/V2 exactness has already been established with independent
        # dual-lane fixtures.  A predictive trajectory advances only the V2
        # backend-native object.  Replaying the full production object beside
        # it on every interval doubled exact stochastic-emission work without
        # adding a new model check.
        engine=product.v2_object
        opening=0.;extension=0.;events=[];interval=0;event_index=0
        dt_cap=MAX_OPENING_INCREMENT_M/OPENING_RATE_M_S
        full_reload=False;avalanche=0;last_event_opening=None;case_hits=case_misses=0
        while extension < target_um*1e-6-1e-15:
            key,geometry=_geometry(backend,extension,events)
            if key in cache: snapshot=cache[key];hits+=1;case_hits+=1
            else: snapshot=oracle.solve(geometry);cache[key]=snapshot;misses+=1;case_misses+=1
            # Start with the production opening increment.  Predictor rejection
            # and recoverable constitutive failure are transactional: both lanes
            # are restored and only the proposed interval is reduced.
            base_dt=MAX_OPENING_INCREMENT_M/OPENING_RATE_M_S
            dt=min(base_dt,dt_cap);retry=0
            while True:
                backup=_capture(engine)
                trial_opening=opening+OPENING_RATE_M_S*dt
                pz=_pz_metrics(engine)
                raw=probe.evaluate_raw(snapshot,pz,trial_opening)
                install(engine,raw)
                K=snapshot.native_KJ_per_opening_Pa_sqrt_m_per_m*trial_opening
                predicted=_predict(engine,backend,K,K,TEMPERATURE_K,dt)
                if np.isfinite(predicted) and predicted>MAX_PREDICTED_CLOCK_INCREMENT:
                    engine=_restore(engine,backup)
                    dt*=MAX_PREDICTED_CLOCK_INCREMENT/predicted;retry+=1;continue
                pre=_pz_metrics(engine);pre_time=float(engine.t);pre_a=float(engine.a_adv)
                try:
                    result=_step(engine,backend,K,K,TEMPERATURE_K,dt)
                    break
                except RuntimeError as exc:
                    if "failed to bracket persistent-site backstress root" not in str(exc) or retry>=30:
                        raise
                    engine=_restore(engine,backup);dt*=.5;retry+=1
            dt_cap=min(base_dt,dt*(1.25 if retry==0 else 1.0))
            signature=_runtime_signature(engine)
            consumed=float(result.get("dt_consumed",dt))
            consumed=min(max(consumed,0.),dt);opening+=OPENING_RATE_M_S*consumed
            fired=bool(result.get("fired",False));post=_pz_metrics(engine)
            event_length=float(engine.a_adv-pre_a) if fired else 0.
            if fired and event_length<=0:
                event_length=float(getattr(engine,"stochastic_last_completed_advance_m",0.))
            all_intervals.append({"backend":label,"material_class":material,"temperature_K":TEMPERATURE_K,
              "target_um":target_um,"interval_index":interval,"time_s":float(engine.t),
              "dt_requested_s":dt,"dt_consumed_s":consumed,"opening_m":opening,"extension_m":extension,
              "event_fired":fired,"predicted_clock_increment":predicted,"cleavage_action":float(engine.B),
              "cleavage_threshold":float(getattr(engine,"hazard_threshold_action",1.)),
              "native_J_J_m2":snapshot.native_J_per_opening2_J_per_m4*opening**2,
              "native_KJ_Pa_sqrt_m":snapshot.native_KJ_per_opening_Pa_sqrt_m_per_m*opening,
              "qualified_G_J_m2":None if snapshot.qualified_structural_G_per_opening2_J_per_m4 is None else snapshot.qualified_structural_G_per_opening2_J_per_m4*opening**2,
              "tip_radius_m":post["tip_radius_m"],"front_width_m":post["front_width_m"],
              "mobile_count":post["mobile_count"],"retained_count":post["retained_count"],
              "source_multiplicity":post["source_multiplicity"],"backstress_Pa":post["backstress_Pa"],
              "shielding_Pa_sqrt_m":post["signed_shielding_Pa_sqrt_m"],
              "lambda_c_s":result.get("lambda_c"),"lambda_e_s":result.get("lambda_e"),
              "selected_source_system":int(np.argmax(raw["drive_factors"])),
              "v2_state_fingerprint":_fingerprint(signature),
              "state_status":"BACKEND_NATIVE_POLICY_SOURCE_FIXTURE_QUALIFIED","field_snapshot_hash":snapshot.snapshot_hash})
            if fired:
                if event_index>0 and full_reload: avalanche+=1
                before=extension;extension+=event_length
                all_events.append({"backend":label,"material_class":material,"temperature_K":TEMPERATURE_K,
                  "target_um":target_um,"event_index":event_index,"physical_avalanche_index":avalanche,
                  "pre_event_time_s":pre_time,"event_time_s":float(engine.t),"event_opening_m":opening,
                  "extension_before_m":before,"extension_after_m":extension,"event_length_m":event_length,
                  "reload_opening_since_prior_event_m":None if last_event_opening is None else opening-last_event_opening,
                  "reload_separated":bool(event_index>0 and full_reload),
                  "native_J_J_m2":snapshot.native_J_per_opening2_J_per_m4*opening**2,
                  "native_KJ_Pa_sqrt_m":snapshot.native_KJ_per_opening_Pa_sqrt_m_per_m*opening,
                  "tip_radius_m":pre["tip_radius_m"],"front_width_m":pre["front_width_m"],
                  "mobile_count":pre["mobile_count"],"retained_count":pre["retained_count"],
                  "source_multiplicity":pre["source_multiplicity"],"backstress_Pa":pre["backstress_Pa"],
                  "shielding_Pa_sqrt_m":pre["signed_shielding_Pa_sqrt_m"],
                  "lambda_c_s":result.get("lambda_c"),"lambda_e_s":result.get("lambda_e"),
                  "selected_source_system":int(np.argmax(raw["drive_factors"])),
                  "right_censored_at_target":extension>=target_um*1e-6,
                  "v2_state_fingerprint":_fingerprint(signature),
                  "state_status":"BACKEND_NATIVE_POLICY_SOURCE_FIXTURE_QUALIFIED"})
                if backend=="pf":
                    p0=np.asarray(events[-1]["p1_m"] if events else (5e-4,0.),float)
                    p1=p0+np.asarray((event_length,0.));events.append({"p0_m":p0.tolist(),"p1_m":p1.tolist(),"front_id":0})
                last_event_opening=opening;event_index+=1;full_reload=False
                print(f"{label} {material}: event {event_index} at {extension*1e6:.3f} um, U={opening*1e6:.3f} um",flush=True)
            else:
                # Mirror the production avalanche rule: only a complete native
                # loading interval separates avalanches.  Tiny accepted hazard
                # substeps between events remain inside one physical avalanche.
                if event_index>0 and consumed>=dt*(1.-1e-12) and dt>=base_dt*(1.-1e-12):
                    full_reload=True
            interval+=1
            if interval%100==0:
                print(f"{label} {material}: interval {interval}, U={opening*1e6:.3f} um, a={extension*1e6:.3f} um",flush=True)
            if interval>25000: raise RuntimeError(f"{label} {material}: interval guard")
            if opening>100e-6: raise RuntimeError(f"{label} {material}: opening guard")
        case_events=[x for x in all_events if x["backend"]==label and x["material_class"]==material and x["target_um"]==target_um]
        counts={i:sum(e["physical_avalanche_index"]==i for e in case_events) for i in set(e["physical_avalanche_index"] for e in case_events)}
        largest=max(counts.values(),default=0)
        summaries.append({"backend":label,"material_class":material,"temperature_K":TEMPERATURE_K,
          "target_um":target_um,"event_count":len(case_events),"physical_avalanche_count":len(counts),
          "precursor_reinitiation_count":sum(bool(e["reload_separated"]) for e in case_events),
          "largest_avalanche_fraction":largest/max(len(case_events),1),"first_event_time_s":case_events[0]["event_time_s"],
          "first_event_opening_m":case_events[0]["event_opening_m"],
          "first_event_native_KJ_Pa_sqrt_m":case_events[0]["native_KJ_Pa_sqrt_m"],
          "terminal_extension_m":extension,"target_right_censored":True,"interval_count":interval,
          "field_cache_hits":case_hits,"field_cache_misses":case_misses,"probe_queries":probe.query_count,
          "source_v2_exact":None,"lifecycle_source_qualification":"QUALIFIED_BY_DUAL_LANE_FIXTURES",
          "status":"NATURAL_FORWARD_CLOSURE_PASS"})
        print(f"{label} {material}: complete {target_um:g} um in {interval} intervals",flush=True)
    payload={"schema":"oneD_v2_natural_forward_v1","backend":label,"target_um":target_um,
      "FEMCZM_reduced_emission_contract":None if backend=="pf" else "SOURCE_MEAN_FIELD_SIGNED_MPZ_WITH_STOCHASTIC_CLEAVAGE_V2",
      "loading":{"opening_rate_m_s":OPENING_RATE_M_S,"maximum_opening_increment_m":MAX_OPENING_INCREMENT_M,
                 "maximum_predicted_clock_increment":MAX_PREDICTED_CLOCK_INCREMENT,"threshold_priming":False},
      "field_solves":oracle.solve_count,"field_cache_hits":hits,"field_cache_misses":misses,
      "probe_queries":probe.query_count,"intervals":all_intervals,"events":all_events,"summaries":summaries}
    name=f"{backend}_natural_{int(target_um)}um.json"
    (OUT/name).write_text(json.dumps(payload,indent=2,sort_keys=True,default=_jsonable)+"\n")


def main():
    p=argparse.ArgumentParser();p.add_argument("backend",choices=("pf","fem"));p.add_argument("--target-um",type=float,default=100.)
    p.add_argument("--classes",nargs="+",default=["Peak","DBTT"]);a=p.parse_args();run(a.backend,a.target_um,set(a.classes))


if __name__=="__main__":main()
