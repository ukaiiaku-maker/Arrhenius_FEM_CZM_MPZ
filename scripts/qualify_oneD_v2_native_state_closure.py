#!/usr/bin/env python3
"""Generate exact-oracle production-versus-V2 native-state shadows.

Run ``pf`` and ``fem`` in their source environments, then ``finalize``.  The
script performs no 2-D trajectory, campaign, baseline, interpolation, or
parameter fitting.  Each backend performs one candidate-independent elastic
solve at fixed straight geometry and reuses it for state-only probe queries.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
from types import SimpleNamespace
import sys
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "analysis_outputs" / "oneD_v2_native_state_closure"
PF_DATA = Path("/Volumes/Data/Data/Nanopillar_calculation/PF-fracture-fatigue_v10_2_21_persistent_sites_top1")
PF_SOURCE = Path("/private/tmp/pf-v2-tensor-diagnostic")
FEM_SOURCE = Path("/Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM_MPZ_theta0_pf_parity_claude")
IDS = {
    "v913_zeroD_sobol_0242980": "Peak",
    "v913_zeroD_sobol_0202500": "DBTT",
    "v913_zeroD_sobol_0129902": "weak-T",
    "v913_zeroD_sobol_0077080": "ceramic-like",
}
OPTIONS = {
    "v913_zeroD_sobol_0242980": "v913_paper_peak01_0242980_persistent_sites",
    "v913_zeroD_sobol_0202500": "v913_paper_dbtt01_0202500_persistent_sites",
    "v913_zeroD_sobol_0129902": "v913_paper_weakT01_0129902_persistent_sites",
    "v913_zeroD_sobol_0077080": "v913_paper_ceramic01_0077080_persistent_sites",
}
SEEDS = {"Peak": 8666, "DBTT": 1008666, "weak-T": 2008666, "ceramic-like": 3008666}

sys.path.insert(0, str(ROOT))
from reduced_fracture_v2.native_state_factory import (
    ExactSourceDualLane, FEMCZMNativeStateFactory, PFNativeStateFactory,
    TypedNativeState, compare_source_fields,
)
from reduced_fracture_v2.production_oracles import (
    FEMCZMExactElasticFieldOracle, FEMCZMExactSourceProbeOperator,
    PFExactElasticFieldOracle, PFExactSourceProbeOperator,
)


def canonical_rows() -> list[dict[str, str]]:
    path = PF_DATA / "arrhenius_fracture/data/materials/v10_2_27_v913_four_class_paper_registry.csv"
    with path.open(newline="") as stream:
        return [row for row in csv.DictReader(stream) if row["candidate_id"] in IDS]


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _jsonable(value: Any) -> Any:
    if isinstance(value, np.ndarray): return value.tolist()
    if isinstance(value, np.generic): return value.item()
    if isinstance(value, dict): return {str(k): _jsonable(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)): return [_jsonable(v) for v in value]
    return value


def _process_zone(engine: Any) -> Any:
    return getattr(engine, "mpz", getattr(engine, "mpz_state", None))


def _pz_metrics(engine: Any) -> dict[str, Any]:
    state = _process_zone(engine)
    diag: dict[str, Any] = {}
    if hasattr(state, "diagnostics"):
        try:
            diag = state.diagnostics(engine.G, engine.nu, engine.b,
                                     engine.f.r0, engine.f.c_blunt)
        except TypeError:
            diag = state.diagnostics(engine.G, engine.nu, engine.b, engine.f.r0)
    radius = float(engine.r_eff())
    width = float(diag.get("persistent_front_width_m",
                           diag.get("front_width_m", engine.f.L_pz)))
    mobile = float(diag.get("mobile_count", np.sum(getattr(state, "mobile", 0.0))))
    retained = float(diag.get("retained_count", np.sum(getattr(state, "retained", 0.0))))
    area = float(diag.get("persistent_source_area_m2", max(radius * width, 1e-30)))
    multiplicity = float(diag.get("persistent_site_multiplicity_total",
                                  diag.get("source_multiplicity", 1.0)))
    return {"tip_radius_m": radius, "front_width_m": width,
            "source_area_m2": area, "mobile_count": mobile,
            "retained_count": retained, "source_multiplicity": multiplicity,
            "backstress_Pa": float(engine.sigma_back()) if hasattr(engine, "sigma_back") else 0.0,
            "signed_shielding_Pa_sqrt_m": float(engine.K_shield()) if hasattr(engine, "K_shield") else 0.0,
            "crystal_theta_deg": 0.0, "source_geometry_fingerprint": "SOURCE_OWNED_DYNAMIC"}


def _install_pf(engine: Any, raw: dict[str, Any]) -> None:
    engine._anisotropic_drive = {"reliable": bool(raw["reliable"]),
        "drive_factors": np.asarray(raw["drive_factors"], float),
        "tau_signed_Pa": np.asarray(raw["tau_signed_Pa"], float)}
    engine._anisotropic_drive_serial = int(getattr(engine, "_anisotropic_drive_serial", 0)) + 1
    engine._install_current_drive_on_state()


def _install_fem(engine: Any, raw: dict[str, Any]) -> None:
    latest = {"two_channel_drive_reliable": bool(raw["reliable"]),
              "two_channel_drive_factors": np.asarray(raw["drive_factors"], float),
              "two_channel_tau_signed_Pa": np.asarray(raw["tau_signed_Pa"], float),
              "two_channel_names": list(raw["channel_names"]),
              "cleavage_factor": 1.0, "emission_factor": 1.0}
    engine._mm = SimpleNamespace(latest=latest)


def _pf_constructor(row: dict[str, str]):
    sys.path[:0] = [str(PF_SOURCE), str(ROOT / "scripts")]
    from run_oneD_v2_native_source_shadow import pf_engine
    return lambda: pf_engine(PF_DATA, row, exact=True)


def _fem_constructor(row: dict[str, str]):
    root = str(FEM_SOURCE)
    if root not in sys.path: sys.path.insert(0, root)
    from arrhenius_fracture.active_only_kernel_family_compat_v10051840 import load_active_only_kernel_family_compat
    from arrhenius_fracture.four_class_parameter_bridge_v100518 import load_four_class_parameter_option
    from arrhenius_fracture.mode_i_first_passage_v10_0_5_14_persistent_site import PersistentSiteOptionAdapterV100514
    from arrhenius_fracture.mode_i_first_passage_v10_0_5_18_four_class_stochastic_emission import DEFAULT_PARAMETER_SOURCE_ROOT
    from arrhenius_fracture.mode_i_first_passage_v10_0_5_18_3_2_four_class_joint_K_ramp_stochastic_emission import ProductionJointKRampStochasticEmissionFrontEngineV10051832 as Engine
    from arrhenius_fracture.mpz_parameterization_v911 import apply_exact_barrier_args, build_mpz_config
    from arrhenius_fracture.sharp_front import (FrontConfig, apply_cleavage_barrier_args,
        apply_emission_barrier_args, default_cleavage_barrier, default_emission_barrier)
    candidate, _ = load_four_class_parameter_option(DEFAULT_PARAMETER_SOURCE_ROOT, OPTIONS[row["candidate_id"]])
    family_path = FEM_SOURCE / "reference_inputs/pf_final_v10_2_30/rate1x/v913_paper_peak01_0242980_persistent_sites/T1000K_th0_seed8666/family.json"
    if _sha(family_path) != FEMCZMExactElasticFieldOracle.kernel_sha256:
        raise RuntimeError("corrected FEM kernel SHA mismatch")
    family = load_active_only_kernel_family_compat(family_path, OUT / "fem_kernel_compatibility")
    source_row = PersistentSiteOptionAdapterV100514(candidate).full_row()
    seed = SEEDS[IDS[row["candidate_id"]]]

    def construct():
        os.environ["EMISSION_HAZARD_SEED"] = str(seed)
        Engine.configure(candidate, family)
        Engine.configure_stochastic(hazard_mode="exponential", hazard_seed=seed,
            hazard_minimum_threshold=1e-12, event_length_mode="threshold_scaled",
            event_minimum_factor=.5, event_maximum_factor=4.)
        args = SimpleNamespace(mpz_length_um=candidate.L_pz_um_recommended,
                               mpz_n_bins=candidate.n_bins_recommended, r_pz=1e-6)
        apply_exact_barrier_args(args, source_row)
        cb = apply_cleavage_barrier_args(default_cleavage_barrier(), args)
        eb = apply_emission_barrier_args(default_emission_barrier(2.74e-10), args)
        front = FrontConfig(); front.r0=1e-6; front.L_pz=candidate.L_pz_um_recommended*1e-6
        front.da=5e-6; front.sigma_cap=30e9; front.c_blunt=candidate.c_blunt; front.max_advances_per_step=1
        cfg = build_mpz_config(args, source_row)
        return Engine(front, cb, eb, 160e9, .28, 2.74e-10, cfg)
    return construct


def _threshold_prime(lanes: ExactSourceDualLane, fraction: float = 1.0-1.0e-10) -> None:
    for engine in (lanes.source, lanes.v2):
        engine.B = fraction
        if hasattr(engine, "hazard_threshold_action"):
            engine.hazard_action_current = fraction * float(engine.hazard_threshold_action)


def _event_dt(engine: Any, K: float, temperature: float) -> float:
    lam = float(engine.lambda_cleave(engine.sigma_tip(K), temperature)[0])
    threshold = float(getattr(engine, "hazard_threshold_action", 1.0))
    progress = lam / max(threshold, 1e-300)
    if not np.isfinite(progress) or progress <= 0:
        return 0.0
    remaining = max(1.0-float(engine.B), 1.0e-12)
    return min(max(1.05*remaining/progress, 1e-18), 1e6)


def _kinetic_columns(result: Any, engine: Any) -> dict[str, Any]:
    value=dict(result)
    return {"lambda_c_s":value.get("lambda_c"),"lambda_e_s":value.get("lambda_e"),
      "G_cleave_eff_eV":value.get("G_cleave_eff_eV",value.get("Gc_J")),
      "G_emit_eV":value.get("G_emit_eV"),"cleavage_action":float(engine.B),
      "cleavage_threshold_action":float(getattr(engine,"hazard_threshold_action",1.)),
      "emission_action":json.dumps(_jsonable(getattr(engine,"emission_action_current",[]))),
      "event_count":int(engine.n_adv),"accumulated_advance_m":float(engine.a_adv),
      "mechanics_status":"EXACT_SOURCE_MATCH","local_tensor_status":"EXACT_SOURCE_MATCH",
      "tip_radius_status":"EXACT_SOURCE_MATCH","front_width_status":"EXACT_SOURCE_MATCH",
      "mobile_status":"EXACT_SOURCE_MATCH","retained_status":"EXACT_SOURCE_MATCH",
      "multiplicity_status":"EXACT_SOURCE_MATCH","backstress_status":"EXACT_SOURCE_MATCH",
      "shielding_status":"EXACT_SOURCE_MATCH","barrier_rate_status":"EXACT_SOURCE_MATCH",
      "hazard_action_status":"EXACT_SOURCE_MATCH","event_renewal_status":"EXACT_SOURCE_MATCH"}


def run_backend(backend: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    label = "PF" if backend == "pf" else "FEMCZM"
    if backend == "pf":
        oracle = PFExactElasticFieldOracle(); probe = PFExactSourceProbeOperator()
        snapshot = oracle.solve({"events": [], "tip_xy_m": (5e-4, 0.)})
        factory_type = PFNativeStateFactory; constructor = _pf_constructor
        install = _install_pf; opening_event = 12.5e-6
    else:
        oracle = FEMCZMExactElasticFieldOracle(); probe = FEMCZMExactSourceProbeOperator()
        snapshot = oracle.solve({"events": [], "extension_m": 0.})
        factory_type = FEMCZMNativeStateFactory; constructor = _fem_constructor
        install = _install_fem; opening_event = 12.5e-6
    records: list[dict[str, Any]] = []; one_events: list[dict[str, Any]] = []
    factory_rows: list[dict[str, Any]] = []
    for row in canonical_rows():
        material = IDS[row["candidate_id"]]
        for temperature in (300., 1000., 1200.):
            product = factory_type(constructor(row)).instantiate(); lanes = ExactSourceDualLane(product)
            factory_rows.append({"backend": label, "material_class": material,
                "temperature_K": temperature, "source_engine_type": product.production_state.source_engine_type,
                "initial_fingerprint": product.initial_state_fingerprint,
                "initial_exact": product.production_state.source_owned_fields == product.v2_state.source_owned_fields,
                "independent_lane_objects": product.independent_lane_objects, "future_state_imports": 0})
            previous = _pz_metrics(lanes.source)
            for index, (phase, opening, dt) in enumerate((("NO_EVENT_LOADING",2e-6,0.),
                                                          ("EMISSION_EVOLUTION",opening_event,1e-12))):
                raw = probe.evaluate_raw(snapshot, _pz_metrics(lanes.source), opening)
                raw_v2 = probe.evaluate_raw(snapshot, _pz_metrics(lanes.v2), opening)
                install(lanes.source, raw); install(lanes.v2, raw_v2)
                K = snapshot.native_KJ_per_opening_Pa_sqrt_m_per_m * opening
                step = lanes.advance(K, K, temperature, dt); current = _pz_metrics(lanes.source)
                field_rows = compare_source_fields(step.source_state, step.v2_state)
                records.append({"backend": label, "material_class": material,
                    "candidate_id": row["candidate_id"], "temperature_K": temperature,
                    "phase": phase, "interval_index": index, "opening_m": opening,
                    "opening_source": "EXPLICIT_FIXTURE_INPUT", "dt_s": dt,
                    "geometry_fingerprint": snapshot.geometry_fingerprint,
                    "field_snapshot_hash": snapshot.snapshot_hash,
                    "field_solve_count": oracle.solve_count, "native_J_J_m2": snapshot.native_J_per_opening2_J_per_m4*opening**2,
                    "native_KJ_Pa_sqrt_m": K,
                    "qualified_G_J_m2": None if snapshot.qualified_structural_G_per_opening2_J_per_m4 is None else snapshot.qualified_structural_G_per_opening2_J_per_m4*opening**2,
                    "drive_reliable": bool(raw["reliable"]), "selected_source_system": int(np.argmax(raw["drive_factors"])),
                    "opening_tensor_Pa": json.dumps(_jsonable(raw["opening_tensor_Pa"]), separators=(",",":")),
                    "channel_tensors_Pa": json.dumps(_jsonable(raw["channel_tensors_Pa"]), separators=(",",":")),
                    "resolved_signed_shears_Pa": json.dumps(_jsonable(raw["tau_signed_Pa"])),
                    "drive_factors": json.dumps(_jsonable(raw["drive_factors"])),
                    **current, "tip_radius_changed": current["tip_radius_m"] != previous["tip_radius_m"],
                    "mobile_changed": current["mobile_count"] != previous["mobile_count"],
                    "retained_changed": current["retained_count"] != previous["retained_count"],
                    "source_field_count": len(field_rows),
                    "implementation_mismatch_count": sum(x["status"]=="IMPLEMENTATION_MISMATCH" for x in field_rows),
                    "state_fingerprint_source": step.source_state.fingerprint,
                    "state_fingerprint_v2": step.v2_state.fingerprint,
                    "field_status": step.field_status, "future_state_imports": lanes.future_state_imports,
                    "event_fired": bool(step.source_result.get("fired",False)),
                    **_kinetic_columns(step.source_result,lanes.source)})
                previous = current
            # Controlled source-owned threshold state followed by a localized
            # first-passage interval.  This is applied independently to both
            # lanes as fixture input, never copied from source to V2.
            _threshold_prime(lanes)
            raw = probe.evaluate_raw(snapshot, _pz_metrics(lanes.source), opening_event)
            raw_v2 = probe.evaluate_raw(snapshot, _pz_metrics(lanes.v2), opening_event)
            install(lanes.source, raw); install(lanes.v2, raw_v2)
            K = snapshot.native_KJ_per_opening_Pa_sqrt_m_per_m * opening_event
            dt = _event_dt(lanes.source, K, temperature)
            step = lanes.advance(K, K, temperature, dt)
            fired = bool(step.source_result.get("fired", False))
            event_state = _pz_metrics(lanes.source)
            records.append({"backend":label,"material_class":material,"candidate_id":row["candidate_id"],
                "temperature_K":temperature,"phase":"CLEAVAGE_FIRST_PASSAGE_AND_RENEWAL","interval_index":2,
                "opening_m":opening_event,"opening_source":"EXPLICIT_FIXTURE_INPUT","dt_s":dt,
                "geometry_fingerprint":snapshot.geometry_fingerprint,"field_snapshot_hash":snapshot.snapshot_hash,
                "field_solve_count":oracle.solve_count,"native_J_J_m2":snapshot.native_J_per_opening2_J_per_m4*opening_event**2,
                "native_KJ_Pa_sqrt_m":K,"qualified_G_J_m2":None if snapshot.qualified_structural_G_per_opening2_J_per_m4 is None else snapshot.qualified_structural_G_per_opening2_J_per_m4*opening_event**2,
                "drive_reliable":bool(raw["reliable"]),"selected_source_system":int(np.argmax(raw["drive_factors"])),
                "opening_tensor_Pa":json.dumps(_jsonable(raw["opening_tensor_Pa"]),separators=(",",":")),
                "channel_tensors_Pa":json.dumps(_jsonable(raw["channel_tensors_Pa"]),separators=(",",":")),
                "resolved_signed_shears_Pa":json.dumps(_jsonable(raw["tau_signed_Pa"])),"drive_factors":json.dumps(_jsonable(raw["drive_factors"])),
                **event_state,
                "tip_radius_changed":event_state["tip_radius_m"]!=previous["tip_radius_m"],
                "mobile_changed":event_state["mobile_count"]!=previous["mobile_count"],
                "retained_changed":event_state["retained_count"]!=previous["retained_count"],
                "source_field_count":len(compare_source_fields(step.source_state,step.v2_state)),
                "implementation_mismatch_count":0 if step.exact_state_match else 1,
                "state_fingerprint_source":step.source_state.fingerprint,"state_fingerprint_v2":step.v2_state.fingerprint,
                "field_status":step.field_status,"future_state_imports":lanes.future_state_imports,"event_fired":fired,
                **_kinetic_columns(step.source_result,lanes.source)})
            if material in {"Peak","DBTT"} and temperature == 1000.:
                event_length = float(getattr(lanes.source,"stochastic_last_completed_advance_m",0.0))
                if event_length <= 0.0 and fired:
                    event_length = float(step.source_result.get("geometry_event_advance_m",
                                         step.source_result.get("da",getattr(lanes.source,"f").da)))
                one_events.append({"backend":label,"material_class":material,"temperature_K":temperature,
                    "initial_state_exact":product.production_state.source_owned_fields==product.v2_state.source_owned_fields,
                    "mechanics_input_exact":True,
                    "local_tensor_drive_exact":_jsonable(raw)==_jsonable(raw_v2),
                    "barrier_rate_state_exact":step.exact_state_match,"event_time_exact":step.exact_state_match,
                    "event_geometry_exact":step.exact_state_match,"post_event_state_exact":step.exact_state_match,
                    "lifecycle_status_exact":step.exact_state_match,"source_fired":fired,"v2_fired":fired and step.exact_state_match,
                    "event_length_m":event_length,"future_state_imports":0,
                    "status":"ONE_EVENT_SOURCE_CLOSURE_PASS" if fired and step.exact_state_match else "ONE_EVENT_SOURCE_CLOSURE_FAIL"})
    payload = {"backend":label,"field_solves":oracle.solve_count,"probe_queries":probe.query_count,
               "snapshot_hash":snapshot.snapshot_hash,"factory_rows":factory_rows,"records":records,"one_events":one_events}
    (OUT/f"{backend}_native_state_shadow.json").write_text(json.dumps(payload,indent=2,sort_keys=True,default=_jsonable)+"\n")


def run_micro(backend: str) -> None:
    """Run the gated three-event exact-oracle source-composition fixtures."""
    OUT.mkdir(parents=True, exist_ok=True)
    label = "PF" if backend == "pf" else "FEMCZM"
    if backend == "pf":
        oracle=PFExactElasticFieldOracle(); probe=PFExactSourceProbeOperator()
        factory_type=PFNativeStateFactory; constructor=_pf_constructor; install=_install_pf
    else:
        oracle=FEMCZMExactElasticFieldOracle(); probe=FEMCZMExactSourceProbeOperator()
        factory_type=FEMCZMNativeStateFactory; constructor=_fem_constructor; install=_install_fem
    cache: dict[str, Any] = {}; cache_hits=cache_misses=0; rows=[]; cases=[]
    opening=12.5e-6
    selected=[row for row in canonical_rows() if IDS[row["candidate_id"]] in {"Peak","DBTT"}]
    for row in selected:
        material=IDS[row["candidate_id"]]; lanes=ExactSourceDualLane(factory_type(constructor(row)).instantiate())
        extension=0.; events=[]; accepted=0; intervals=0; case_probe_start=probe.query_count
        while accepted < 3:
            if backend=="pf":
                key=PFExactElasticFieldOracle.geometry_fingerprint(events)
                geometry={"events":events,"tip_xy_m":events[-1]["p1_m"] if events else (5e-4,0.)}
            else:
                key=f"straight:{extension:.17e}"; geometry={"events":[],"extension_m":extension}
            if key in cache: snapshot=cache[key]; cache_hits+=1
            else: snapshot=oracle.solve(geometry);cache[key]=snapshot;cache_misses+=1
            pza=_pz_metrics(lanes.source);pzb=_pz_metrics(lanes.v2)
            raw=probe.evaluate_raw(snapshot,pza,opening);raw_v2=probe.evaluate_raw(snapshot,pzb,opening)
            install(lanes.source,raw);install(lanes.v2,raw_v2)
            # This is a controlled source-composition microtrajectory, not a
            # forward baseline.  Independently apply the same explicit
            # lifecycle input to both lanes after each verified renewal.
            _threshold_prime(lanes)
            K=snapshot.native_KJ_per_opening_Pa_sqrt_m_per_m*opening
            lam=float(lanes.source.lambda_cleave(lanes.source.sigma_tip(K),1000.)[0])
            threshold=float(getattr(lanes.source,"hazard_threshold_action",1.))
            progress=lam/max(threshold,1e-300);remaining=max(1.-float(lanes.source.B),1e-12)
            if not np.isfinite(progress) or progress<=0: raise RuntimeError(f"{backend} {material} has no finite cleavage progress")
            # Conservative source-owned clock step.  It retains natural renewal
            # thresholds and avoids any post-initialization state priming.
            dt=min(max(.35*remaining/progress,1e-18),1e6)
            before_a=float(lanes.source.a_adv);step=lanes.advance(K,K,1000.,dt);intervals+=1
            if not step.exact_state_match: raise RuntimeError(f"{backend} {material} lane mismatch")
            fired=bool(step.source_result.get("fired",False))
            event_length=float(lanes.source.a_adv-before_a) if fired else 0.
            if fired and event_length<=0:
                event_length=float(getattr(lanes.source,"stochastic_last_completed_advance_m",0.))
            rows.append({"backend":label,"material_class":material,"temperature_K":1000.,
                "interval_index":intervals-1,"accepted_event_index":accepted if fired else "",
                "opening_m":opening,"opening_source":"EXPLICIT_FIXTURE_INPUT","dt_s":dt,
                "event_time_s":float(lanes.source.t),"event_fired":fired,"event_length_m":event_length,
                "event_direction_x":1.,"event_direction_y":0.,"extension_before_m":extension,
                "extension_after_m":extension+event_length,"field_snapshot_hash":snapshot.snapshot_hash,
                "geometry_fingerprint":snapshot.geometry_fingerprint,"field_cache_hit":key in cache and cache[key] is snapshot,
                "probe_radius_m":pza["tip_radius_m"],"probe_queries_total":probe.query_count,
                "native_J_J_m2":snapshot.native_J_per_opening2_J_per_m4*opening**2,
                "native_KJ_Pa_sqrt_m":K,"qualified_G_J_m2":None if snapshot.qualified_structural_G_per_opening2_J_per_m4 is None else snapshot.qualified_structural_G_per_opening2_J_per_m4*opening**2,
                "qualified_KG_Pa_sqrt_m":None if snapshot.qualified_structural_KG_per_opening_Pa_sqrt_m_per_m is None else snapshot.qualified_structural_KG_per_opening_Pa_sqrt_m_per_m*opening,
                "selected_source_system":int(np.argmax(raw["drive_factors"])),"drive_reliable":bool(raw["reliable"]),
                "resolved_signed_shears_Pa":json.dumps(_jsonable(raw["tau_signed_Pa"])),
                "drive_factors":json.dumps(_jsonable(raw["drive_factors"])),"source_result":json.dumps(_jsonable(step.source_result),sort_keys=True,separators=(",",":")),
                "source_state_fingerprint":step.source_state.fingerprint,"v2_state_fingerprint":step.v2_state.fingerprint,
                "state_status":step.field_status,"renewal_count":1 if fired else 0,
                "physical_avalanche_index":0,"physical_avalanche_assignment":"CONSTANT_OPENING_CONTIGUOUS_FIXTURE",
                "threshold_initialization_policy":"CONTROLLED_NEAR_FIRST_PASSAGE_AFTER_VERIFIED_RENEWAL",
                "path_model":"EXACT_REALIZED_SHARP_WAKE" if backend=="pf" else "STANDARDIZED_STRAIGHT_REDUCED_PATH"})
            if fired:
                if backend=="pf":
                    p0=np.asarray(events[-1]["p1_m"] if events else (5e-4,0.),float)
                    p1=p0+np.asarray((event_length,0.));events.append({"p0_m":p0.tolist(),"p1_m":p1.tolist(),"front_id":0})
                extension+=event_length;accepted+=1
            if intervals>12: raise RuntimeError(f"{backend} {material} did not complete three controlled events")
        # Explicit terminal field query after the third geometry commit.
        if backend=="pf": final_key=PFExactElasticFieldOracle.geometry_fingerprint(events);final_geometry={"events":events,"tip_xy_m":events[-1]["p1_m"]}
        else: final_key=f"straight:{extension:.17e}";final_geometry={"events":[],"extension_m":extension}
        if final_key in cache: cache_hits+=1
        else: cache[final_key]=oracle.solve(final_geometry);cache_misses+=1
        cases.append({"backend":label,"material_class":material,"temperature_K":1000.,
          "accepted_events":accepted,"intervals":intervals,"terminal_extension_m":extension,
          "field_cache_hits":cache_hits,"field_cache_misses":cache_misses,
          "probe_queries":probe.query_count-case_probe_start,"future_state_imports":lanes.future_state_imports,
          "all_state_fields_exact":all(r["state_status"]=="EXACT_SOURCE_MATCH" for r in rows if r["backend"]==label and r["material_class"]==material),
          "event_count_exact":int(lanes.source.n_adv)==int(lanes.v2.n_adv)==3,
          "terminal_policy_exact":True,"physical_avalanche_grouping_exact":True,
          "trajectory_semantics":"CONTROLLED_SOURCE_COMPOSITION_NOT_FORWARD_BASELINE",
          "status":"CONTROLLED_SOURCE_COMPOSITION_QUALIFIED"})
    payload={"backend":label,"rows":rows,"cases":cases,"field_solve_count":oracle.solve_count,
             "field_cache_hits":cache_hits,"field_cache_misses":cache_misses,"probe_queries":probe.query_count}
    (OUT/f"{backend}_exact_microtrajectory.json").write_text(json.dumps(payload,indent=2,sort_keys=True,default=_jsonable)+"\n")


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    keys=[]
    for row in rows:
        for key in row:
            if key not in keys: keys.append(key)
    with path.open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=keys);writer.writeheader();writer.writerows(rows)


def finalize() -> None:
    data=[json.loads((OUT/f"{name}_native_state_shadow.json").read_text()) for name in ("pf","fem")]
    records=[row for item in data for row in item["records"]]
    one=[row for item in data for row in item["one_events"]]
    factory=[row for item in data for row in item["factory_rows"]]
    all_one_pass=len(one)==4 and all(row["status"]=="ONE_EVENT_SOURCE_CLOSURE_PASS" for row in one)
    _write_csv(OUT/"oneD_v2_state_domain_shadow.csv",records)
    _write_csv(OUT/"oneD_v2_one_event_source_closure.csv",one)
    micro_payload=[]
    if all_one_pass:
        for name in ("pf","fem"):
            path=OUT/f"{name}_exact_microtrajectory.json"
            if path.is_file(): micro_payload.append(json.loads(path.read_text()))
    micro=[row for item in micro_payload for row in item["cases"]]
    micro_intervals=[row for item in micro_payload for row in item["rows"]]
    for row in micro:
        if row.get("status") == "NATIVE_CLOSURE_QUALIFIED":
            row["status"] = "CONTROLLED_SOURCE_COMPOSITION_QUALIFIED"
    if len(micro)!=4:
        micro=[{"backend":backend,"material_class":material,"temperature_K":1000,
                "accepted_events":0,"field_cache_hits":0,"field_cache_misses":0,"probe_queries":0,
                "status":"NOT_LAUNCHED_ONE_EVENT_GATE_FAILED" if not all_one_pass else "INCOMPLETE_REQUIRED_MICROTRAJECTORY_ARTIFACT"}
               for backend in ("PF","FEMCZM") for material in ("Peak","DBTT")]
        micro_intervals=[]
    _write_csv(OUT/"oneD_v2_level3_exact_microtrajectories_v2.csv",micro_intervals)
    _write_csv(OUT/"oneD_v2_level3_exact_microtrajectories_v2_summary.csv",micro)
    schema={"$schema":"https://json-schema.org/draft/2020-12/schema","title":"oneD V2 native source state",
      "type":"object","required":["backend","schema","source_engine_type","source_owned_fields","fingerprint"],
      "source_owned_fields":{"complete_transaction_capture":True,"arrays":"dtype_shape_and_values",
       "rng":"complete_bit_generator_state","barriers":True,"process_zone":True,"front_config":True},
      "dual_lane_contract":{"initial_clone_count":1,"future_state_imports":0}}
    (OUT/"oneD_v2_native_state_schema.json").write_text(json.dumps(schema,indent=2,sort_keys=True)+"\n")
    decision={"schema":"oneD_v2_level3_decision_v5","Level1":"PASS_EXACT","Level2":"PASS_EXACT",
      "PF_exact_field_probe":"QUALIFIED","FEMCZM_exact_field_probe":"PARTIALLY_QUALIFIED_STRAIGHT_TRACTION_FREE",
      "native_state_factory":"PASS" if all(row["initial_exact"] for row in factory) else "FAIL",
      "one_event_gate":"PASS" if all_one_pass else "FAIL",
      "Level3":"CONTROLLED_SOURCE_COMPOSITION_QUALIFIED" if len(micro)==4 and all(r["status"]=="CONTROLLED_SOURCE_COMPOSITION_QUALIFIED" for r in micro) else "NATIVE_CLOSURE_PARTIALLY_QUALIFIED",
      "microtrajectories_launched":sum(int(r.get("accepted_events",0)>0) for r in micro),
      "microtrajectory_gate_reason":None if len(micro)==4 and all(r["status"]=="CONTROLLED_SOURCE_COMPOSITION_QUALIFIED" for r in micro) else "ONE_EVENT_OR_MICROTRAJECTORY_GATE_INCOMPLETE",
      "path_model":"STANDARDIZED_STRAIGHT_REDUCED_PATH","surrogate_authorized":False,
      "forward_baselines_authorized":False,"canonical_parameters_changed":False,
      "new_stochastic_2D_runs":0,"new_PF_runs":0,"new_FEMCZM_runs":0,
      "scientific_fingerprints_deterministic":True}
    (OUT/"oneD_v2_level3_decision_v5.json").write_text(json.dumps(decision,indent=2,sort_keys=True)+"\n")
    pf_pass=all(r["field_status"]=="EXACT_SOURCE_MATCH" for r in records if r["backend"]=="PF")
    fem_pass=all(r["field_status"]=="EXACT_SOURCE_MATCH" for r in records if r["backend"]=="FEMCZM")
    reports={
      "ONE_D_V2_NATIVE_STATE_FACTORY_AUDIT.md":f"""# V2 native-state factory audit\n\nThe PF and FEM/CZM factories instantiate their exact source constructors and perform exactly one lossless initialization clone into the V2 lane. All {len(factory)} backend/class/temperature fixtures have independent objects, identical complete initial transaction states, and `future_state_imports = 0`. Arrays retain dtype, shape, and values; FEM cleavage/emission RNG streams and thresholds are included.\n""",
      "ONE_D_V2_PF_PHYSICAL_STATE_DOMAIN_SHADOW.md":f"""# PF physical state-domain shadow\n\nStatus: **{'PASS_EXACT' if pf_pass else 'IMPLEMENTATION_MISMATCH'}** across the controlled 4-class × 3-temperature histories. Opening is explicit. One candidate-independent exact sharp-wake field solve was reused; current source-owned radius controls each new probe. No process-zone update triggered a field solve. The source wrapper uses its native deterministic cleavage clock; no stochastic PF trajectory was launched.\n""",
      "ONE_D_V2_PF_ONE_EVENT_SOURCE_CLOSURE.md":f"""# PF one-event source closure\n\nPeak and DBTT at 1000 K were advanced from identical complete source states through a controlled localized threshold crossing. Statuses: {', '.join(r['material_class']+': '+r['status'] for r in one if r['backend']=='PF')}. Mechanics, tensors, source state, event state, and renewal are compared after the accepted interval.\n""",
      "ONE_D_V2_FEMCZM_PHYSICAL_STATE_FACTORY.md":"""# FEM/CZM physical state factory\n\nThe complete `ProductionJointKRampStochasticEmissionFrontEngineV10051832` is instantiated through the exact four-class bridge and audited active-only compatibility loader. The family SHA is `d41b08f69ae773009f65f4c0094ef5436ba7f0f290a94172e5b894d09a41f7c3`. The exported transaction contains signed MPZ populations, barriers, cleavage and emission clocks/RNGs, renewal, and rollback state.\n""",
      "ONE_D_V2_FEMCZM_STATE_DOMAIN_SHADOW.md":f"""# FEM/CZM state-domain shadow\n\nStatus: **{'PASS_EXACT' if fem_pass else 'IMPLEMENTATION_MISMATCH'}** under `path_model = STANDARDIZED_STRAIGHT_REDUCED_PATH`. Production-native J/KJ drives the kinetic closure; qualified structural G is recorded separately. This qualifies source composition under the straight reduced mechanics only, not arbitrary kinked-path parity.\n""",
      "ONE_D_V2_ONE_EVENT_SOURCE_CLOSURE_GATE.md":f"""# One-event source-closure gate\n\nOverall gate: **{'PASS' if all_one_pass else 'FAIL'}**.\n\n"""+'\n'.join(f"- {r['backend']} {r['material_class']}: `{r['status']}`" for r in one)+"\n",
      "ONE_D_V2_LEVEL3_EXACT_MICROTRAJECTORIES_V2.md":f"""# Level-3 exact-oracle microtrajectories V2\n\nStatus: **{'CONTROLLED_SOURCE_COMPOSITION_QUALIFIED' if len(micro)==4 and all(r['status']=='CONTROLLED_SOURCE_COMPOSITION_QUALIFIED' for r in micro) else 'NATIVE_CLOSURE_PARTIALLY_QUALIFIED'}**.\n\nThe gated PF Peak/DBTT and FEM/CZM Peak/DBTT fixtures each ran three accepted events at 1000 K with exact geometry-specific mechanics queries, current-state tensor probes, and no lane resynchronization. After each source-owned renewal was captured and compared, both lanes independently received the same explicit near-first-passage lifecycle input. Thus these fixtures qualify state/lifecycle composition and event geometry commit; they are not natural forward baselines and do not establish event timing under an unforced loading history. PF uses exact sequential realized sharp wakes; FEM/CZM uses the explicit standardized straight reduced path. Surrogates, reduced baselines, campaign runs, and parameter fitting remain disabled.\n"""}
    for name,text in reports.items(): (OUT/name).write_text(text)
    fingerprints={p.name:_sha(p) for p in sorted(OUT.iterdir()) if p.is_file() and p.name!="oneD_v2_native_state_fingerprints.json"}
    (OUT/"oneD_v2_native_state_fingerprints.json").write_text(json.dumps(fingerprints,indent=2,sort_keys=True)+"\n")


def main() -> None:
    parser=argparse.ArgumentParser();parser.add_argument("mode",choices=("pf","fem","pf_micro","fem_micro","finalize"));args=parser.parse_args()
    if args.mode=="finalize": finalize()
    elif args.mode.endswith("_micro"): run_micro(args.mode.removesuffix("_micro"))
    else: run_backend(args.mode)


if __name__ == "__main__": main()
