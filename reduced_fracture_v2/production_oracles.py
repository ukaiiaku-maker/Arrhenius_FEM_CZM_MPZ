"""Concrete deterministic PF and FEM/CZM elastic-field oracles.

Imports are delayed so the V2 package remains usable without either production
repository on ``sys.path``.  These classes perform mechanics/probe evaluations
only; they never integrate hazard or advance stochastic state.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace
from typing import Any, Mapping, Sequence
import hashlib
import importlib.util
import json
import os
import sys

import numpy as np

from .exact_oracle import NormalizedElasticFieldSnapshot, ProbeQueryKey
from .local_drive import (DriveQualificationError, EmissionDrive,
                          GlobalStructuralDrive, LocalTensorSample,
                          NativeDriveBundle)


def _sha_arrays(*values: Any) -> str:
    h = hashlib.sha256()
    for value in values:
        array = np.ascontiguousarray(value)
        h.update(str(array.dtype).encode()); h.update(str(array.shape).encode()); h.update(array.tobytes())
    return h.hexdigest()


def _json_sha(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def _tensor_field(sigma_gp: np.ndarray) -> np.ndarray:
    sigma = np.asarray(sigma_gp, float)
    out = np.zeros((sigma.shape[1], 2, 2), float)
    out[:, 0, 0], out[:, 1, 1] = sigma[0], sigma[1]
    out[:, 0, 1] = out[:, 1, 0] = sigma[2]
    return out


def _gp_field(tensors: np.ndarray) -> np.ndarray:
    value = np.asarray(tensors, float)
    return np.vstack((value[:, 0, 0], value[:, 1, 1], value[:, 0, 1]))


def _emission(drive: Mapping[str, Any]) -> EmissionDrive:
    return EmissionDrive(
        tuple(tuple(float(x) for x in row) for row in drive["opening_tensor_Pa"]),
        tuple(tuple(tuple(float(x) for x in row) for row in tensor) for tensor in drive["channel_tensors_Pa"]),
        tuple(tuple(float(x) for x in row) for row in drive["trace_directions"]),
        tuple(tuple(float(x) for x in row) for row in drive["trace_normals"]),
        tuple(float(x) for x in drive["tau_signed_Pa"]),
        tuple(float(x) for x in drive["drive_factors"]),
        bool(drive["reliable"]), int(drive.get("mechanics_serial", 0)),
    )


class PFExactElasticFieldOracle:
    source_commit = "9e884fb0b0845da621d2612bdf1042e481b8df49"
    mechanics_factory_hash = "PF_SHARP_WAKE_NX36_NY72_SEED42_THETA0_V1"

    def __init__(self, source_repository: str | Path = "/private/tmp/pf-v2-tensor-diagnostic",
                 reference_opening_m: float = 1.0e-5):
        self.source_repository = Path(source_repository)
        self.reference_opening_m = float(reference_opening_m)
        self.solve_count = 0

    def _imports(self):
        root = str(self.source_repository)
        if root not in sys.path: sys.path.insert(0, root)
        os.environ.setdefault("PF_TENSOR_DRIVE_DIAGNOSTIC_JSONL", "/dev/null")
        from arrhenius_fracture.config import ElasticProperties, GeometryConfig, JIntegralConfig, MeshConfig
        from arrhenius_fracture.crack_backend import SharpWakeBackend
        from arrhenius_fracture.crystal import cubic_plane_strain_D
        from arrhenius_fracture.fem import assemble_mechanics, elastic_energy_densities, solve_dirichlet
        from arrhenius_fracture.j_integral import compute_J_integral
        from arrhenius_fracture.mesh import make_boundary_data, make_tri_mesh
        from arrhenius_fracture.anisotropic_front_direction_fix_v10227 import install_front_direction_fix
        return SimpleNamespace(**locals())

    @staticmethod
    def geometry_fingerprint(events: Sequence[Mapping[str, Any]]) -> str:
        canonical = [{k: list(map(float, event[k])) for k in ("p0_m", "p1_m")} for event in events]
        return _json_sha(canonical)

    def solve(self, geometry: Mapping[str, Any]) -> NormalizedElasticFieldSnapshot:
        self.solve_count += 1; q = self._imports(); q.install_front_direction_fix()
        events = list(geometry.get("events", ()))
        g = q.GeometryConfig(); mc = q.MeshConfig(nx=36, ny=72, jitter=0., tip_h_fine=1e-6, tip_ratio=1.2)
        mesh = q.make_tri_mesh(g, mc, seed=42, tip_center=np.array([g.a0, 0.]))
        mat = q.ElasticProperties(); D = q.cubic_plane_strain_D(523e9, 203e9, 160e9, 0.)
        x, y = np.asarray(mesh.nodes).T
        damage = np.zeros(mesh.nn); damage[(x <= g.a0) & (np.abs(y) <= g.notch_half_thickness)] = 1.
        backend = q.SharpWakeBackend(); boundary = q.make_boundary_data(mesh, g); displacement = np.zeros(mesh.ndof)
        kill = max(float(mesh.hbar_tip), .5e-6)
        segments = [(np.array([0., 0.]), np.array([g.a0, 0.]))]
        for index, event in enumerate(events):
            p0, p1 = np.asarray(event["p0_m"], float), np.asarray(event["p1_m"], float)
            direction = p1 - p0; length = float(np.linalg.norm(direction))
            if length <= 0: raise DriveQualificationError("PF realized event has nonpositive length")
            result = backend.advance(mesh=mesh, boundary=boundary, damage=damage,
                displacement=displacement, p0=p0, p1=p1, direction=direction/length,
                front_id=int(event.get("front_id", 0)), kill_r=kill)
            if not result.inserted: raise DriveQualificationError(f"PF wake reconstruction failed: {result.reason}")
            damage = result.damage; segments.append((p0, p1))
        tip = np.asarray(events[-1]["p1_m"] if events else geometry.get("tip_xy_m", (g.a0, 0.)), float)
        tangent = np.asarray(geometry.get("crack_tangent", (1., 0.)), float); tangent /= np.linalg.norm(tangent)
        normal = np.array([-tangent[1], tangent[0]])
        u0 = np.zeros(mesh.ndof); ep = np.zeros((3, mesh.ne)); rho = np.zeros(mesh.ne)
        K, R, *_ = q.assemble_mechanics(mesh, u0, ep, rho, damage, D, mat, kappa=1e-6)
        u, reaction = q.solve_dirichlet(K, R, u0, boundary, .5*self.reference_opening_m, -.5*self.reference_opening_m)
        _, _, sigma, *_ = q.assemble_mechanics(mesh, u, ep, rho, damage, D, mat, kappa=1e-6)
        stored, _ = q.elastic_energy_densities(mesh, u, ep, sigma, D)
        energy = float(np.sum(stored * mesh.area_e))
        _, _, info = q.compute_J_integral(mesh, u, sigma, stored, damage, tip, tangent, mat,
            ell=80e-6, cfg=q.JIntegralConfig(), crack_segments=segments, exclude_radius=2*kill)
        jsigned = float(info.get("J_signed", info.get("J", 0.)) or 0.); jvalue = max(jsigned, 0.)
        strain_gp = np.linalg.solve(D, sigma)
        centroids = np.asarray(mesh.nodes)[np.asarray(mesh.elems)].mean(axis=1)
        return NormalizedElasticFieldSnapshot(
            "PF", str(self.source_repository), self.source_commit,
            "anisotropic_emission_v10174.build_front_drive", self.mechanics_factory_hash,
            "PLANE_STRAIN", "FINITE_WIDTH_STIFFNESS_KILLED_WAKE", "BINARY_STIFFNESS_KILL",
            self.geometry_fingerprint(events), _sha_arrays(mesh.nodes, mesh.elems), _sha_arrays(damage),
            tuple(tip), tuple(tangent), tuple(normal), mesh.nodes, mesh.elems, centroids, mesh.area_e,
            _tensor_field(sigma)/self.reference_opening_m,
            _tensor_field(strain_gp)/self.reference_opening_m, damage,
            float(reaction)/self.reference_opening_m, energy/self.reference_opening_m**2,
            jvalue/self.reference_opening_m**2,
            np.sqrt(jvalue*mat.Eprime)/self.reference_opening_m,
            {"J_signed_reference": jsigned, "exclude_radius_m": 2*kill, "ell_m": 80e-6},
            validity_flags=("LINEAR_ELASTIC", "EXACT_REALIZED_WAKE", "CANDIDATE_INDEPENDENT"),
            reference_opening_m=self.reference_opening_m,
        )


class PFExactSourceProbeOperator:
    probe_convention = "anisotropic_emission_v10174"
    def __init__(self, source_repository: str | Path = "/private/tmp/pf-v2-tensor-diagnostic"):
        self.source_repository = Path(source_repository); self.query_count = 0

    def evaluate_raw(self, snapshot: NormalizedElasticFieldSnapshot,
                     process_zone_state: Mapping[str, Any], opening_m: float) -> Mapping[str, Any]:
        if snapshot.backend != "PF": raise DriveQualificationError("PF probe requires a PF field snapshot")
        root = str(self.source_repository)
        if root not in sys.path: sys.path.insert(0, root)
        from arrhenius_fracture.anisotropic_emission_v10174 import AnisotropicEmissionConfig
        from arrhenius_fracture import anisotropic_emission_v10174 as ae
        radius = float(process_zone_state["tip_radius_m"])
        if radius <= 0: raise DriveQualificationError("tip radius must be positive")
        self.query_count += 1
        mesh = SimpleNamespace(nodes=snapshot.nodes_m, elems=snapshot.elements,
            area_e=snapshot.element_areas_m2, ne=len(snapshot.elements), nn=len(snapshot.nodes_m))
        cfg = AnisotropicEmissionConfig(crystal_theta_deg=float(process_zone_state.get("crystal_theta_deg", 0.)),
            probe_radius_m=radius, sector_half_angle_deg=float(process_zone_state.get("sector_half_angle_deg", 25.)),
            damage_cutoff=float(process_zone_state.get("damage_cutoff", .85)),
            min_elements=int(process_zone_state.get("min_elements", 3)), schmid_reference=float(process_zone_state.get("schmid_reference", .5)))
        return ae.build_front_drive(mesh, _gp_field(snapshot.stress_at_opening(opening_m)),
            snapshot.damage_or_interface_state, np.asarray(snapshot.crack_tip_xy_m), cfg)

    def evaluate(self, snapshot: NormalizedElasticFieldSnapshot,
                 process_zone_state: Mapping[str, Any], opening_m: float) -> NativeDriveBundle:
        raw = self.evaluate_raw(snapshot, process_zone_state, opening_m)
        supports = {"opening": tuple(map(int, raw["opening_probe"]["element_indices"]))}
        weights = {"opening": tuple(map(float, raw["opening_probe"]["element_weights"]))}
        for index, probe in enumerate(raw["channel_probes"]):
            supports[f"channel{index}"] = tuple(map(int, probe["element_indices"]))
            weights[f"channel{index}"] = tuple(map(float, probe["element_weights"]))
        radius = float(process_zone_state["tip_radius_m"]); width = float(process_zone_state["front_width_m"])
        source_area = float(process_zone_state.get("source_area_m2", radius*width))
        query = ProbeQueryKey(snapshot.snapshot_hash, radius, width, source_area,
            str(process_zone_state.get("source_geometry_fingerprint", "UNSPECIFIED")),
            _json_sha(process_zone_state.get("source_positions_m", ())), "BCC_THETA0_TWO_CHANNELS", self.probe_convention)
        samples = tuple(LocalTensorSample(float(c[0]-snapshot.crack_tip_xy_m[0]), float(c[1]-snapshot.crack_tip_xy_m[1]),
            ((float(s[0]), float(s[2])), (float(s[2]), float(s[1]))))
            for probe in (raw["opening_probe"], *raw["channel_probes"])
            for c, s in zip(probe["element_centroids_m"], probe["element_stress_components_Pa"]))
        emission = _emission(raw); chosen = int(np.argmax(np.asarray(emission.drive_factors)))
        return NativeDriveBundle(GlobalStructuralDrive(opening_m,
            snapshot.reaction_per_opening_N_per_m2*opening_m,
            snapshot.native_J_per_opening2_J_per_m4*opening_m**2,
            snapshot.native_KJ_per_opening_Pa_sqrt_m_per_m*opening_m),
            snapshot.crack_tip_xy_m, tuple(raw["front_direction"]), tuple(raw["front_normal"]), samples,
            emission, radius, width, snapshot.topology_wake_fingerprint, snapshot.source_commit,
            snapshot.snapshot_hash, (0., float("inf")), (-float("inf"), float("inf"), -float("inf"), float("inf")),
            str(raw["channel_names"][chosen]), supports, weights, dict(process_zone_state),
            snapshot.snapshot_hash, query.fingerprint())


class FEMCZMExactElasticFieldOracle:
    kernel_sha256 = "d41b08f69ae773009f65f4c0094ef5436ba7f0f290a94172e5b894d09a41f7c3"
    mechanics_factory_hash = "FEMCZM_CONNECTED_CUT_H1UM_THETA0_D41"

    def __init__(self, source_repository: str | Path = "/Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM_MPZ_theta0_pf_parity_claude",
                 qualification_repository: str | Path = "/private/tmp/oneD-v2-qualified-fem-analysis",
                 reference_opening_m: float = 12.5e-6):
        self.source_repository = Path(source_repository); self.qualification_repository = Path(qualification_repository)
        self.reference_opening_m = float(reference_opening_m); self.solve_count = 0

    def _canonical(self):
        path = self.qualification_repository/"scripts/qualify_canonical_connected_vcct.py"
        spec = importlib.util.spec_from_file_location("oneD_v2_exact_canonical_vcct", path)
        module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module); return module

    def solve(self, geometry: Mapping[str, Any]) -> NormalizedElasticFieldSnapshot:
        events = list(geometry.get("events", ()))
        if any(abs(float(event["p1_m"][1])-float(event["p0_m"][1])) > 1e-15 for event in events):
            raise DriveQualificationError("qualified FEM/CZM exact oracle currently admits straight interface cuts only")
        self.solve_count += 1; canonical = self._canonical(); v1 = canonical.v1
        extension = float(geometry.get("extension_m", sum(float(np.linalg.norm(np.asarray(e["p1_m"])-np.asarray(e["p0_m"]))) for e in events)))
        material = v1.ElasticProperties(); g = v1.GeometryConfig(); tip = g.a0 + extension
        # Import the already-qualified deterministic connected-cut constructor.
        script = self.qualification_repository.parent/"oneD-v2-common-kernel/scripts/generate_oneD_v2_fem_mechanics_maps.py"
        spec = importlib.util.spec_from_file_location("oneD_v2_fem_map_builder", script)
        builder = importlib.util.module_from_spec(spec); spec.loader.exec_module(builder)
        cut, pairs = builder.structured_mesh_at(tip, g)
        D = v1.cubic_plane_strain_D(523e9, 203e9, 160e9, 0.); state = canonical.state(cut, pairs, tip, D, material, self.reference_opening_m)
        h = 1.0e-6
        state1 = canonical.state(cut, pairs, tip+h, D, material, self.reference_opening_m)
        state2 = canonical.state(cut, pairs, tip+2*h, D, material, self.reference_opening_m)
        g_energy = (3*state["W"]-4*state1["W"]+state2["W"])/(2*h)
        g_compliance = .5*state["F"]**2*(-3*state["C"]+4*state1["C"]-state2["C"])/(2*h)
        gi, gii, g_vcct, da_vcct, _, _ = canonical.vcct(cut, pairs, tip, state)
        g_values = np.asarray((g_energy, g_compliance, g_vcct), float)
        spread = float(np.ptp(g_values)); g_scale = max(float(np.max(np.abs(g_values))), 1e-300)
        g_qualified = bool(np.all(g_values > 0) and spread/g_scale <= .10)
        g_mean = float(np.mean(g_values)) if g_qualified else None
        js = builder.j_domains(cut, state, tip, material); jvalue = float(np.median([row[1] for row in js]))
        strain_gp = np.linalg.solve(D, state["sigma"]); centroids = cut.nodes[cut.elems].mean(axis=1)
        topology = _json_sha({"tip":tip,"pairs":pairs}); geom = _json_sha({"extension_m":extension,"events":events})
        return NormalizedElasticFieldSnapshot("FEMCZM", str(self.source_repository), self.kernel_sha256,
            "active_only_kernel_family_compat_v10051840+d41+robust_probe", self.mechanics_factory_hash,
            "PLANE_STRAIN", "CODIMENSION_ONE_DISPLACEMENT_DISCONTINUITY", "TRACTION_FREE",
            geom, _sha_arrays(cut.nodes, cut.elems), topology, (tip,0.), (1.,0.), (0.,1.),
            cut.nodes, cut.elems, centroids, cut.area_e, _tensor_field(state["sigma"])/self.reference_opening_m,
            _tensor_field(strain_gp)/self.reference_opening_m, np.zeros(cut.nn),
            float(state["F"])/self.reference_opening_m, float(state["W"])/self.reference_opening_m**2,
            jvalue/self.reference_opening_m**2, np.sqrt(max(jvalue,0.)*material.Eprime)/self.reference_opening_m,
            {"domain_values":{str(r):j for r,j,_ in js}},
            None if g_mean is None else g_mean/self.reference_opening_m**2,
            None if g_mean is None else np.sqrt(max(g_mean,0.)*material.Eprime)/self.reference_opening_m,
            {"G_energy_reference":float(g_energy),"G_compliance_reference":float(g_compliance),
             "G_VCCT_reference":float(g_vcct),"G_VCCT_I_reference":float(gi),"G_VCCT_II_reference":float(gii),
             "vcct_segment_m":float(da_vcct),"uncertainty_reference":spread,"qualified":g_qualified},
            validity_flags=("LINEAR_ELASTIC","QUALIFIED_CONNECTED_CUT","CANDIDATE_INDEPENDENT","KERNEL_D41_RESOLVED"),
            reference_opening_m=self.reference_opening_m)


class FEMCZMExactSourceProbeOperator:
    def __init__(self, source_repository: str | Path = "/Volumes/Data/Data/Nanopillar_calculation/Arrhenius_FEM_CZM_MPZ_theta0_pf_parity_claude"):
        self.source_repository = Path(source_repository); self.query_count = 0

    def evaluate_raw(self, snapshot: NormalizedElasticFieldSnapshot,
                     process_zone_state: Mapping[str, Any], opening_m: float) -> Mapping[str, Any]:
        if snapshot.backend != "FEMCZM": raise DriveQualificationError("FEM/CZM probe requires FEM/CZM field")
        root = str(self.source_repository)
        if root not in sys.path: sys.path.insert(0, root)
        from arrhenius_fracture import mixed_mode_first_passage_v8 as mm
        from arrhenius_fracture.numerical_resilience_v1005183 import make_robust_process_zone_traction_probe
        from arrhenius_fracture.crystal import bcc_slip_traces
        from arrhenius_fracture.anisotropic_two_channel_drive_v100514 import resolve_two_channel_drive_from_tensors
        probe = make_robust_process_zone_traction_probe(mm.process_zone_traction_probe)
        self.query_count += 1; radius=float(process_zone_state["tip_radius_m"]); direction=np.asarray(snapshot.crack_tangent)
        mesh=SimpleNamespace(nodes=snapshot.nodes_m,elems=snapshot.elements,area_e=snapshot.element_areas_m2)
        sigma=_gp_field(snapshot.stress_at_opening(opening_m)); damage=snapshot.damage_or_interface_state
        kwargs=dict(radius_m=radius,annulus_half_width=float(process_zone_state.get("annulus_half_width",.45)),
            sector_half_angle_deg=float(process_zone_state.get("sector_half_angle_deg",40.)),damage_cutoff=float(process_zone_state.get("damage_cutoff",.85)),min_elements=int(process_zone_state.get("min_elements",4)))
        opening=probe(mesh,sigma,damage,snapshot.crack_tip_xy_m,direction,**kwargs)
        channels=[]; probes=[]
        for trace in bcc_slip_traces(float(process_zone_state.get("crystal_theta_deg",0.))):
            ray=np.asarray(trace["t"],float)
            if ray@direction<0: ray=-ray
            item=probe(mesh,sigma,damage,snapshot.crack_tip_xy_m,ray,**kwargs); probes.append(item); channels.append(np.asarray(item["stress_tensor"]))
        reliable=bool(opening.get("reliable")) and all(bool(x.get("reliable")) for x in probes)
        if not reliable: raise DriveQualificationError("corrected FEM/CZM robust tensor probe is unreliable")
        drive=resolve_two_channel_drive_from_tensors(np.asarray(opening["stress_tensor"]),channels,direction,float(process_zone_state.get("crystal_theta_deg",0.)))
        return {"opening_tensor_Pa":np.asarray(opening["stress_tensor"]).tolist(),"channel_tensors_Pa":[x.tolist() for x in channels],
            "tau_signed_Pa":drive["two_channel_tau_signed_Pa"],"drive_factors":drive["two_channel_drive_factors"],
            "trace_directions":drive["two_channel_trace_directions"],"trace_normals":drive["two_channel_trace_normals"],
            "channel_names":drive["two_channel_names"],"reliable":True,"opening_probe":opening,"channel_probes":probes}

    def evaluate(self, snapshot: NormalizedElasticFieldSnapshot,
                 process_zone_state: Mapping[str, Any], opening_m: float) -> NativeDriveBundle:
        raw=self.evaluate_raw(snapshot,process_zone_state,opening_m); emission=_emission(raw)
        radius=float(process_zone_state["tip_radius_m"]);width=float(process_zone_state["front_width_m"])
        source_area=float(process_zone_state.get("source_area_m2",radius*width))
        query=ProbeQueryKey(snapshot.snapshot_hash,radius,width,source_area,
            str(process_zone_state.get("source_geometry_fingerprint","UNSPECIFIED")),
            _json_sha(process_zone_state.get("source_positions_m",())),"BCC_THETA0_TWO_CHANNELS",
            "FEMCZM_D41_ROBUST_PROCESS_ZONE_TRACTION_PROBE")
        j=snapshot.native_J_per_opening2_J_per_m4*opening_m**2
        g=None if snapshot.qualified_structural_G_per_opening2_J_per_m4 is None else snapshot.qualified_structural_G_per_opening2_J_per_m4*opening_m**2
        selected=int(np.argmax(np.asarray(emission.drive_factors)))
        return NativeDriveBundle(GlobalStructuralDrive(opening_m,
            snapshot.reaction_per_opening_N_per_m2*opening_m,j,snapshot.native_KJ_per_opening_Pa_sqrt_m_per_m*opening_m,
            g,None if g is None else snapshot.qualified_structural_KG_per_opening_Pa_sqrt_m_per_m*opening_m),
            snapshot.crack_tip_xy_m,snapshot.crack_tangent,snapshot.crack_normal,(),emission,radius,width,
            snapshot.topology_wake_fingerprint,snapshot.source_commit,snapshot.snapshot_hash,(0.,float("inf")),
            (-float("inf"),float("inf"),-float("inf"),float("inf")),str(raw["channel_names"][selected]),
            {},{},dict(process_zone_state),snapshot.snapshot_hash,query.fingerprint())
