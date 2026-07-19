"""
Diagnostics, history tracking, and plotting for Arrhenius fracture simulations.
"""

import numpy as np
import json
import csv
from dataclasses import dataclass, field
from typing import Dict, List, Optional
import os


@dataclass
class StepDiagnostics:
    """Per-step diagnostic quantities."""
    step: int = 0
    t: float = 0.0
    Uapp: float = 0.0
    Freact: float = 0.0
    Ftop: float = 0.0
    Fbot: float = 0.0
    Fpair_abs: float = 0.0

    # Energetics
    Wext: float = 0.0         # cumulative external work
    Uel: float = 0.0          # stored degraded elastic energy used in thermodynamic audit
    Uel_drive: float = 0.0    # undegraded positive tensile-energy diagnostic
    Uel_undegraded: float = 0.0  # total undegraded elastic energy diagnostic
    Wext_top: float = 0.0     # cumulative work using top reaction times opening
    Wext_pair: float = 0.0    # cumulative signed work from top and bottom moving boundaries
    Wext_abs: float = 0.0     # cumulative positive boundary-work magnitude
    Wp: float = 0.0           # cumulative gross plastic work
    Wp_tip: float = 0.0       # cumulative crack-tip emission dissipation (sub-grid tip plasticity)
    Dp_eff: float = 0.0       # plastic dissipation after stored toughening partition
    Etough: float = 0.0       # stored process-zone toughening energy
    Dtough: float = 0.0       # dissipated process-zone toughening work
    Epf_surf: float = 0.0     # PF surface energy

    # Incremental thermodynamic energy audit.  The residual is
    # dWext - dUel - dEpf - dWp - dEmem - dDmem.  Negative residual means the
    # accepted update spent more energy than external work plus released
    # stored energy can support.
    dWext: float = 0.0
    dWext_top: float = 0.0
    dWext_pair: float = 0.0
    dWext_abs: float = 0.0
    dUel: float = 0.0
    dUel_drive: float = 0.0
    dEpf: float = 0.0
    dWp_step: float = 0.0
    dDp_eff: float = 0.0
    dEtough: float = 0.0
    dDtough: float = 0.0
    Emem: float = 0.0
    Dmem: float = 0.0
    Dfrac: float = 0.0
    dEmem: float = 0.0
    dDmem: float = 0.0
    dDfrac: float = 0.0
    energy_residual: float = 0.0
    energy_residual_rel: float = 0.0
    energy_residual_absWext: float = 0.0
    energy_residual_topWext: float = 0.0
    energy_cumulative_residual: float = 0.0
    energy_cumulative_residual_absWext: float = 0.0
    energy_units_ratio_Uel_over_Wext_abs: float = 0.0
    energy_balance_ok: float = 1.0
    J_tearing: float = 0.0
    KJ_tearing: float = 0.0

    # Crack metrics
    crack_len: float = 0.0    # Epf/Gc
    Da_projected: float = 0.0
    Gamma_total: float = 0.0
    branch_factor: float = 1.0

    # Toughness
    J_domain: float = 0.0     # domain-integral J [J/m²]
    KJ_domain: float = 0.0    # K from domain J [Pa*sqrt(m)]
    J_global: float = 0.0     # global energy balance J
    KJ_global: float = 0.0    # K from global J
    K_force: float = 0.0      # LEFM K from reaction force and initial crack geometry [Pa*sqrt(m)]

    # Field statistics
    rho_mean: float = 0.0
    rho_p95: float = 0.0
    rho_p99: float = 0.0
    rho_max: float = 0.0
    rho_gt_1e14_frac: float = 0.0
    rho_gt_1e15_frac: float = 0.0
    rho_gt_1e16_frac: float = 0.0
    rho_cap_frac: float = 0.0
    dotep_mean: float = 0.0
    dotep_max: float = 0.0
    d_frac: float = 0.0       # fraction of nodes with d > threshold
    plast_frac: float = 0.0   # fraction with active plasticity

    # Thermodynamic plastic-work diagnostics
    dWp_requested: float = 0.0
    dWp_accepted: float = 0.0
    dep_eq_requested_max: float = 0.0
    dep_eq_accepted_max: float = 0.0
    dep_eq_uncapped_max: float = 0.0
    dep_limited_frac: float = 0.0
    thermo_scale_min: float = 1.0
    thermo_scale_mean: float = 1.0
    thermo_admissible_frac: float = 0.0
    thermo_hazard_max: float = 0.0
    thermo_substeps: float = 1.0
    thermo_dt_min: float = 0.0
    thermo_retry_count: float = 0.0

    # Memory energetics / conjugate diagnostics
    memory_energy_increment: float = 0.0
    memory_dissipation_increment: float = 0.0
    memory_A_r_mean: float = 0.0
    memory_A_z_mean: float = 0.0

    # Plastic flow threshold diagnostics
    sigma_eq_mean: float = 0.0
    sigma_eq_max: float = 0.0
    sigma_y_min: float = 0.0
    sigma_y_mean: float = 0.0
    sigma_y_max: float = 0.0
    sigma_T_min: float = 0.0
    sigma_T_mean: float = 0.0
    sigma_T_max: float = 0.0
    sigma_Peierls: float = 0.0
    sigma_eq_over_sigma_y_max: float = 0.0
    yield_frac: float = 0.0
    flow_dgamma_uncapped_max: float = 0.0
    flow_dgamma_cap: float = 0.0
    flow_cap_frac: float = 0.0
    flow_phi_mean: float = 0.0
    flow_phi_max: float = 0.0
    flow_Gtarget_eV_min: float = 0.0
    flow_Gtarget_eV_mean: float = 0.0
    flow_Gtarget_eV_max: float = 0.0
    flow_DG0_eV: float = 0.0
    flow_DGfloor_eV: float = 0.0
    flow_vstar_ref_b3: float = 0.0
    flow_status_zero_stress_frac: float = 0.0
    flow_status_solved_frac: float = 0.0
    flow_status_floor_limited_frac: float = 0.0

    # Tip memory
    rtip_mean: float = 0.0
    shield_mean: float = 0.0
    rtip_amp_mean: float = 0.0
    rtip_min: float = 0.0
    rtip_max: float = 0.0
    shield_max: float = 0.0
    M_fract_mean: float = 1.0
    M_fract_max: float = 1.0
    tip_emit_prob_mean: float = 0.0
    tip_emit_prob_max: float = 0.0
    pz_emit_prob_mean: float = 0.0
    pz_emit_prob_max: float = 0.0
    pz_mobility_prob_mean: float = 0.0
    pz_mobility_prob_max: float = 0.0
    pz_mobile_prob_max: float = 0.0
    pz_escape_prob_max: float = 0.0
    pz_store_prob_mean: float = 0.0
    pz_store_prob_max: float = 0.0
    pz_storage_fraction_mean: float = 0.0
    pz_mobility_hazard_max: float = 0.0
    pz_mobility_hazard_raw_max: float = 0.0
    pz_sigma_mobility_eff_max: float = 0.0
    pz_sigma_tip_eff_max: float = 0.0
    pz_sigma_back_max: float = 0.0
    pz_sigma_back_disl_max: float = 0.0
    pz_sigma_back_mem_max: float = 0.0
    pz_sigma_back_crack_max: float = 0.0
    pz_G_shield_max: float = 0.0
    pz_G_stored_release_max: float = 0.0
    pz_G_stored_release_p99: float = 0.0
    pz_e_stored_max: float = 0.0
    pz_Gc_net_min: float = 0.0
    pz_Gc_net_p01: float = 0.0
    pz_G_app_max: float = 0.0
    pz_G_eff_max: float = 0.0
    pz_crack_R_max: float = 0.0
    pz_crack_R_p99: float = 0.0
    pz_H_eff_drive_max: float = 0.0
    pz_H_eff_drive_p99: float = 0.0
    pz_front_mask_frac: float = 0.0
    pz_H_eff_masked_max: float = 0.0
    pz_H_eff_unmasked_max: float = 0.0
    pz_crack_barrier_min_eV: float = 0.0
    pz_crack_hazard_max: float = 0.0
    pz_crack_prob_mean: float = 0.0
    pz_crack_prob_max: float = 0.0
    pz_crack_hazard_raw_max: float = 0.0
    pz_crack_B_mean: float = 0.0
    pz_crack_B_max: float = 0.0
    pz_emit_B_max: float = 0.0
    pz_emit_rho_max: float = 0.0
    pz_emit_Gshield_max: float = 0.0
    pz_crack_sigma_tip_max: float = 0.0
    pz_emission_hazard_max: float = 0.0
    pz_emission_hazard_raw_max: float = 0.0
    pz_multihit_n_emit_mean: float = 1.0
    pz_multihit_n_emit_max: float = 1.0
    pz_multihit_n_mobility_mean: float = 1.0
    pz_multihit_n_mobility_max: float = 1.0
    pz_multihit_spacing_nm_min: float = 0.0
    pz_multihit_log_suppression_emit_min: float = 0.0
    pz_multihit_log_suppression_mobility_min: float = 0.0
    pz_drho_emit_max: float = 0.0
    pz_drho_rec_max: float = 0.0
    pz_recovery_rate_max: float = 0.0
    dwp_norm_front_mean: float = 0.0

    # Process-zone toughening state driven by accepted plastic work
    q_tough_mean_front: float = 0.0
    q_tough_mean_all: float = 0.0
    q_tough_p95: float = 0.0
    q_tough_p99: float = 0.0
    q_tough_max: float = 0.0
    dqtough_max: float = 0.0
    toughening_weight_mean: float = 0.0
    toughening_energy_increment: float = 0.0
    toughening_dissipation_increment: float = 0.0

    # Cohesive gate
    cohesive_gate_max: float = 0.0
    cohesive_ratio_max: float = 0.0

    # Emergent mode
    Gc_local_mean: float = 0.0       # spatial mean local Gc
    Gc_local_p95: float = 0.0        # 95th percentile local Gc
    Gc_local_p99: float = 0.0        # 99th percentile local Gc
    Gc_local_max: float = 0.0        # max local Gc from plastic dissipation
    Gc_local_mean_front: float = 0.0 # mean Gc at crack front
    pz_rho_source_sat_max: float = 0.0
    pz_rho_source_sat_min: float = 0.0
    pz_source_availability_mean: float = 0.0
    pz_source_availability_min: float = 0.0
    pz_storage_capacity_mean: float = 0.0
    pz_storage_capacity_min: float = 0.0


@dataclass
class SimulationHistory:
    """Complete simulation history for one temperature."""
    T: float = 0.0
    Gc_eff: float = 0.0
    Kc_input: float = 0.0
    ell: float = 0.0

    # Per-step arrays (filled during simulation)
    steps: List[StepDiagnostics] = field(default_factory=list)

    # Saved fields (optional)
    u_fields: Dict[int, np.ndarray] = field(default_factory=dict)
    d_fields: Dict[int, np.ndarray] = field(default_factory=dict)
    rho_fields: Dict[int, np.ndarray] = field(default_factory=dict)
    rtip_fields: Dict[int, np.ndarray] = field(default_factory=dict)
    shield_fields: Dict[int, np.ndarray] = field(default_factory=dict)
    wp_fields: Dict[int, np.ndarray] = field(default_factory=dict)
    Gc_fields: Dict[int, np.ndarray] = field(default_factory=dict)
    M_fields: Dict[int, np.ndarray] = field(default_factory=dict)
    # Elastic / hazard diagnostics: raw FEM max principal stress, de-smeared
    # crack-tip drive stress, accumulated first-passage action, front weight.
    sig1_fields: Dict[int, np.ndarray] = field(default_factory=dict)
    sigtip_fields: Dict[int, np.ndarray] = field(default_factory=dict)
    B_fields: Dict[int, np.ndarray] = field(default_factory=dict)
    fw_fields: Dict[int, np.ndarray] = field(default_factory=dict)

    # Summary (filled at end)
    KJ_final: float = 0.0
    J_final: float = 0.0
    n_steps_used: int = 0

    def add_step(self, diag: StepDiagnostics):
        self.steps.append(diag)

    def get_array(self, field_name: str) -> np.ndarray:
        """Extract a field as a numpy array from step history."""
        return np.array([getattr(s, field_name) for s in self.steps])

    def get_Uapp(self) -> np.ndarray:
        return self.get_array('Uapp')

    def get_Freact(self) -> np.ndarray:
        return self.get_array('Freact')


def save_history(hist: SimulationHistory, filepath: str,
                 mesh_nodes: np.ndarray = None,
                 mesh_elems: np.ndarray = None):
    """Save history to npz file.

    The file contains:
        - Scalar time series (Uapp, Freact, J_domain, KJ_domain, ...)
        - Metadata (T, Gc_eff, Kc_input, ell)
        - Mesh coordinates and connectivity (if provided)
        - Per-step field snapshots: d_step_N, u_step_N, rho_step_N

    Load with:
        data = np.load('history_0700K.npz')
        nodes = data['mesh_nodes']        # (nn, 2)
        elems = data['mesh_elems']        # (ne, 3)
        d5 = data['d_step_5']             # damage at step 5
        u5 = data['u_step_5']             # displacement at step 5
        rho5 = data['rho_step_5']         # dislocation density at step 5
    """
    arrays = {}
    for field_name in [
        'Uapp', 'Freact', 'Ftop', 'Fbot', 'Fpair_abs', 'Wext', 'Wext_top', 'Wext_pair', 'Wext_abs', 'Uel', 'Uel_drive', 'Uel_undegraded', 'Wp', 'Wp_tip', 'Dp_eff', 'Etough', 'Dtough', 'Epf_surf',
        'dWext', 'dWext_top', 'dWext_pair', 'dWext_abs', 'dUel', 'dUel_drive', 'dEpf', 'dWp_step', 'dDp_eff', 'dEtough', 'dDtough',
        'Emem', 'Dmem', 'Dfrac', 'dEmem', 'dDmem', 'dDfrac',
        'energy_residual', 'energy_residual_rel', 'energy_residual_absWext', 'energy_residual_topWext', 'energy_cumulative_residual', 'energy_cumulative_residual_absWext', 'energy_units_ratio_Uel_over_Wext_abs', 'energy_balance_ok',
        'J_tearing', 'KJ_tearing',
        'J_domain', 'KJ_domain', 'J_global', 'KJ_global', 'K_force',
        'rho_mean', 'rho_p95', 'rho_p99', 'rho_max', 'rho_gt_1e14_frac', 'rho_gt_1e15_frac', 'rho_gt_1e16_frac', 'rho_cap_frac', 'dotep_mean', 'dotep_max',
        'dWp_requested', 'dWp_accepted',
        'dep_eq_requested_max', 'dep_eq_accepted_max', 'dep_eq_uncapped_max', 'dep_limited_frac',
        'thermo_scale_min', 'thermo_scale_mean',
        'thermo_admissible_frac', 'thermo_hazard_max',
        'thermo_substeps', 'thermo_dt_min', 'thermo_retry_count',
        'memory_energy_increment', 'memory_dissipation_increment',
        'memory_A_r_mean', 'memory_A_z_mean',
        'sigma_eq_mean', 'sigma_eq_max',
        'sigma_y_min', 'sigma_y_mean', 'sigma_y_max',
        'sigma_T_min', 'sigma_T_mean', 'sigma_T_max', 'sigma_Peierls',
        'sigma_eq_over_sigma_y_max', 'yield_frac',
        'flow_dgamma_uncapped_max', 'flow_dgamma_cap', 'flow_cap_frac',
        'flow_phi_mean', 'flow_phi_max',
        'flow_Gtarget_eV_min', 'flow_Gtarget_eV_mean', 'flow_Gtarget_eV_max',
        'flow_DG0_eV', 'flow_DGfloor_eV', 'flow_vstar_ref_b3',
        'flow_status_zero_stress_frac', 'flow_status_solved_frac', 'flow_status_floor_limited_frac',
        'd_frac', 'Da_projected', 'Gamma_total', 'branch_factor',
        'rtip_mean', 'shield_mean', 'rtip_amp_mean',
        'rtip_min', 'rtip_max', 'shield_max',
        'M_fract_mean', 'M_fract_max',
        'tip_emit_prob_mean', 'tip_emit_prob_max',
        'pz_emit_prob_mean', 'pz_emit_prob_max', 'pz_mobility_prob_mean', 'pz_mobility_prob_max', 'pz_mobile_prob_max', 'pz_escape_prob_max', 'pz_store_prob_mean', 'pz_store_prob_max', 'pz_storage_fraction_mean',
        'pz_mobility_hazard_max', 'pz_mobility_hazard_raw_max', 'pz_sigma_mobility_eff_max', 'pz_sigma_tip_eff_max', 'pz_sigma_back_max', 'pz_sigma_back_disl_max', 'pz_sigma_back_mem_max', 'pz_sigma_back_crack_max',
        'pz_G_shield_max', 'pz_G_stored_release_max', 'pz_G_stored_release_p99', 'pz_e_stored_max', 'pz_Gc_net_min', 'pz_Gc_net_p01',
        'pz_G_app_max', 'pz_G_eff_max', 'pz_crack_R_max', 'pz_crack_R_p99', 'pz_H_eff_drive_max', 'pz_H_eff_drive_p99', 'pz_front_mask_frac', 'pz_H_eff_masked_max', 'pz_H_eff_unmasked_max',
        'pz_crack_barrier_min_eV', 'pz_crack_hazard_max', 'pz_crack_prob_mean', 'pz_crack_prob_max', 'pz_crack_hazard_raw_max', 'pz_crack_B_mean', 'pz_crack_B_max', 'pz_crack_sigma_tip_max', 'pz_emit_B_max', 'pz_emit_rho_max', 'pz_emit_Gshield_max',
        'pz_emission_hazard_max', 'pz_emission_hazard_raw_max', 'pz_multihit_n_emit_mean', 'pz_multihit_n_emit_max', 'pz_multihit_n_mobility_mean', 'pz_multihit_n_mobility_max', 'pz_multihit_spacing_nm_min', 'pz_multihit_log_suppression_emit_min', 'pz_multihit_log_suppression_mobility_min', 'pz_drho_emit_max', 'pz_drho_rec_max', 'pz_recovery_rate_max',
        'dwp_norm_front_mean',
        'q_tough_mean_front', 'q_tough_mean_all', 'q_tough_p95', 'q_tough_p99', 'q_tough_max', 'dqtough_max',
        'toughening_weight_mean', 'toughening_energy_increment', 'toughening_dissipation_increment',
        'cohesive_gate_max', 'cohesive_ratio_max',
        'Gc_local_mean', 'Gc_local_p95', 'Gc_local_p99', 'Gc_local_mean', 'Gc_local_p95', 'Gc_local_p99', 'Gc_local_max', 'Gc_local_mean_front',
    ]:
        try:
            arrays[field_name] = hist.get_array(field_name)
        except AttributeError:
            pass

    arrays['T'] = np.array([hist.T])
    arrays['Gc_eff'] = np.array([hist.Gc_eff])
    arrays['Kc_input'] = np.array([hist.Kc_input])
    arrays['ell'] = np.array([hist.ell])

    # Mesh (saved once per file, enables field reconstruction)
    if mesh_nodes is not None:
        arrays['mesh_nodes'] = mesh_nodes
    if mesh_elems is not None:
        arrays['mesh_elems'] = mesh_elems

    # Field snapshots
    for step, d_field in hist.d_fields.items():
        arrays[f'd_step_{step}'] = d_field
    for step, u_field in hist.u_fields.items():
        arrays[f'u_step_{step}'] = u_field
    for step, rho_field in hist.rho_fields.items():
        arrays[f'rho_step_{step}'] = rho_field
    for step, rtip_field in hist.rtip_fields.items():
        arrays[f'rtip_step_{step}'] = rtip_field
    for step, shield_field in hist.shield_fields.items():
        arrays[f'shield_step_{step}'] = shield_field
    for step, wp_field in hist.wp_fields.items():
        arrays[f'wp_step_{step}'] = wp_field
    for step, Gc_field in hist.Gc_fields.items():
        arrays[f'Gc_step_{step}'] = Gc_field
    for step, M_field in hist.M_fields.items():
        arrays[f'M_step_{step}'] = M_field
    for step, f in hist.sig1_fields.items():
        arrays[f'sig1_step_{step}'] = f
    for step, f in hist.sigtip_fields.items():
        arrays[f'sigtip_step_{step}'] = f
    for step, f in hist.B_fields.items():
        arrays[f'Bcrack_step_{step}'] = f
    for step, f in hist.fw_fields.items():
        arrays[f'frontw_step_{step}'] = f

    # List of snapshot steps (for easy discovery)
    arrays['snapshot_steps'] = np.array(sorted(hist.d_fields.keys()))

    np.savez_compressed(filepath, **arrays)
    print(f"  Saved history to {filepath}")




def _safe_ratio_percent(num, den, den_min=1e-12):
    try:
        num = float(num); den = float(den)
    except Exception:
        return float('nan')
    if not np.isfinite(num) or not np.isfinite(den) or den <= den_min:
        return float('nan')
    return float(100.0 * num / den)


def _classify_failure_mode(hist: SimulationHistory, arrays: dict) -> dict:
    """Classify the run: brittle crack, soft tearing, plastic collapse,
    invalid runaway, or no resolved failure. This prevents soft/plastic
    cases from being reported as clean K_Ic values."""
    def a(name): return arrays.get(name, np.array([]))
    F=a('Freact'); Wext=a('Wext'); Wp=a('Wp'); dfrac=a('d_frac'); Da=a('Da_projected')
    rho=a('rho_max'); Gc=a('Gc_local_max'); Kd=a('KJ_domain')/1e6 if len(a('KJ_domain')) else np.array([])
    Kg=a('KJ_global')/1e6 if len(a('KJ_global')) else np.array([]); yfrac=a('yield_frac'); pfrac=a('plast_frac'); capfrac=a('flow_cap_frac')
    eok=a('energy_balance_ok')
    n=max([len(x) for x in [F,Wext,Wp,dfrac,Da,rho,Gc,Kd,Kg,yfrac,pfrac,capfrac,eok]]+[0])
    if n==0:
        return {'failure_mode':'no_data','failure_mode_reason':'no step diagnostics were recorded','valid_KIc_like':False,'soft_tearing_candidate':False,'plastic_collapse_candidate':False,'tearing_index':0.0}
    def peak(x, default=0.0): return float(np.nanmax(x)) if len(x) else default
    Fabs=np.abs(F) if len(F) else np.zeros(n); Fmax=peak(Fabs); Ffinal=float(Fabs[-1]) if len(Fabs) else 0.0
    load_drop_ratio=float(Ffinal/max(Fmax,1e-30)) if Fmax>0 else 1.0
    rho_peak=peak(rho); Gc_peak=peak(Gc, hist.Gc_eff); K_peak=max(peak(Kd),peak(Kg))
    d_peak=peak(dfrac); Da_peak=peak(Da); ell=max(float(getattr(hist,'ell',0.0)),1e-30); Da_peak_over_ell=Da_peak/ell
    yield_peak=peak(yfrac); plast_peak=peak(pfrac); cap_peak=peak(capfrac)
    wp_pct_series=np.array([])
    if len(Wp) and len(Wext):
        m=np.isfinite(Wp)&np.isfinite(Wext)&(Wext>1e-12)
        if np.any(m): wp_pct_series=100.0*Wp[m]/np.maximum(Wext[m],1e-30)
    wp_pct_peak=float(np.nanmax(wp_pct_series)) if len(wp_pct_series) else float('nan')
    wp_pct_final=_safe_ratio_percent(Wp[-1],Wext[-1]) if len(Wp) and len(Wext) else float('nan')
    energy_bad_frac = float(np.mean(eok < 0.5)) if len(eok) else 0.0
    has_plastic=(yield_peak>1e-3) or (plast_peak>1e-4) or (np.isfinite(wp_pct_peak) and wp_pct_peak>5.0)
    strong_plastic=(yield_peak>0.02) or (plast_peak>0.005) or (np.isfinite(wp_pct_peak) and wp_pct_peak>50.0)
    crack_advanced=(d_peak>0.03) or (Da_peak_over_ell>0.5); load_dropped=load_drop_ratio<0.20
    soft_tearing=bool(strong_plastic and crack_advanced and (not np.isfinite(wp_pct_peak) or wp_pct_peak<2000.0))
    plastic_collapse=bool(strong_plastic and (not crack_advanced) and load_dropped)
    mode='undetermined'; reason=''; valid_KIc_like=False
    if energy_bad_frac > 0.25:
        mode='invalid_energy_balance'; reason=f'incremental thermodynamic audit failed in {energy_bad_frac:.0%} of recorded steps'
    elif rho_peak>5e17:
        mode='invalid_rho_runaway'; reason=f'rho_max peaked at {rho_peak:.3e} m^-2'
    elif np.isfinite(wp_pct_peak) and wp_pct_peak>=2000.0:
        if crack_advanced and strong_plastic:
            mode='invalid_soft_tearing_runaway'; reason=f'plastic damage/tearing candidate, but Wp/Wext reached {wp_pct_peak:.1f}%'
        else:
            mode='invalid_plastic_work'; reason=f'Wp/Wext reached {wp_pct_peak:.1f}%'
    elif Gc_peak>200.0*max(hist.Gc_eff,1e-30):
        mode='invalid_Gc_runaway'; reason=f'local Gc peaked at {Gc_peak:.3g} J/m^2'
    elif K_peak>100.0:
        mode='invalid_K_spike'; reason=f'K diagnostic peaked at {K_peak:.3g} MPa√m'
    elif soft_tearing:
        mode='soft_tearing'; reason='plasticity and crack/damage advance co-occur without invalid work blow-up'
    elif plastic_collapse:
        mode='plastic_collapse'; reason='load dropped with strong plasticity but little projected crack advance'
    elif crack_advanced and (not has_plastic):
        mode='brittle_crack'; reason='crack/damage advanced with negligible plastic work/yielding'; valid_KIc_like=True
    elif crack_advanced and has_plastic:
        mode='mixed_fracture_plasticity'; reason='crack/damage advanced with moderate plasticity'
    elif not crack_advanced and load_dropped:
        mode='load_drop_no_resolved_crack'; reason='load dropped but crack advance metric was small'
    else:
        mode='no_failure_or_no_crack'; reason='no clear crack advance, plastic collapse, or load-drop failure detected'
    tearing_index=0.0
    if crack_advanced:
        pterm=min(1.0,max(yield_peak,plast_peak)/0.05) if max(yield_peak,plast_peak)>0 else 0.0
        wterm=min(1.0,(wp_pct_peak if np.isfinite(wp_pct_peak) else 0.0)/100.0)
        dterm=min(1.0,d_peak/0.10); tearing_index=float(max(pterm,wterm)*dterm)
    if mode.startswith('invalid') or mode in ('soft_tearing','plastic_collapse','mixed_fracture_plasticity'):
        valid_KIc_like=False
    return {'failure_mode':mode,'failure_mode_reason':reason,'valid_KIc_like':bool(valid_KIc_like),'soft_tearing_candidate':bool(mode in ('soft_tearing','invalid_soft_tearing_runaway','mixed_fracture_plasticity')),'plastic_collapse_candidate':bool(mode=='plastic_collapse'),'force_drop_ratio_final':float(load_drop_ratio),'Wp_over_Wext_peak_percent':float(wp_pct_peak) if np.isfinite(wp_pct_peak) else float('nan'),'Wp_over_Wext_final_percent':float(wp_pct_final) if np.isfinite(wp_pct_final) else float('nan'),'yield_frac_peak':float(yield_peak),'plast_frac_peak':float(plast_peak),'flow_cap_frac_peak':float(cap_peak),'d_frac_peak':float(d_peak),'Da_projected_peak_m':float(Da_peak),'Da_projected_peak_over_ell':float(Da_peak_over_ell),'tearing_index':float(tearing_index),'energy_bad_frac':float(energy_bad_frac)}

def history_summary(hist: SimulationHistory) -> dict:
    """Return compact scalar summary for one temperature run.

    The selected toughness is event/window based.  Values after full crack
    passage or numerical plastic runaway are excluded so diagnostic runs that
    continue after failure do not report artificial enormous K values.
    """
    def arr(name):
        try:
            return hist.get_array(name)
        except Exception:
            return np.array([])

    Kd = arr('KJ_domain') / 1e6
    Kg = arr('KJ_global') / 1e6
    Kf = arr('K_force') / 1e6
    F = arr('Freact')
    Wext = arr('Wext')
    Wp = arr('Wp')
    rho = arr('rho_max')
    dfrac = arr('d_frac')
    Da = arr('Da_projected')
    Gc = arr('Gc_local_max')
    M = arr('M_fract_max')
    sigy = arr('sigma_y_mean')
    seqmax = arr('sigma_eq_max')
    yfrac = arr('yield_frac')
    gtarget = arr('flow_Gtarget_eV_mean')
    dg0 = arr('flow_DG0_eV')
    dgfloor = arr('flow_DGfloor_eV')
    Jr = arr('J_tearing')
    Kr = arr('KJ_tearing') / 1e6
    eok = arr('energy_balance_ok')

    arrays_for_mode = {name: arr(name) for name in [
        'Freact','Wext','Wp','d_frac','Da_projected','rho_max','Gc_local_max',
        'KJ_domain','KJ_global','yield_frac','plast_frac','flow_cap_frac',
        'energy_balance_ok','J_tearing','KJ_tearing']}
    mode_info = _classify_failure_mode(hist, arrays_for_mode)

    n = max(len(Kd), len(Kg), len(F), len(Da), len(dfrac), 0)
    ell = max(float(getattr(hist, 'ell', 0.0)), 1e-30)
    valid = np.ones(n, dtype=bool)

    if len(F) == n:
        Fabs = np.abs(F)
        Fmax = max(float(np.nanmax(Fabs)), 1e-30)
        valid &= Fabs > 0.01 * Fmax
        ipeak = int(np.nanargmax(Fabs)) if len(Fabs) else n-1
        post = np.arange(n) > ipeak
        valid[post & (Fabs < 0.10 * Fmax)] = False
    else:
        Fabs = np.ones(n)
        Fmax = 1.0

    if len(rho) == n:
        valid &= np.isfinite(rho) & (rho < 5e17)
    if len(Wp) == n and len(Wext) == n:
        ratio = np.zeros(n)
        ok = np.abs(Wext) > 1e-30
        ratio[ok] = 100.0 * Wp[ok] / np.maximum(np.abs(Wext[ok]), 1e-30)
        valid &= (~ok) | (ratio < 2000.0)
    if len(Gc) == n and len(Gc) > 0:
        valid &= np.isfinite(Gc) & (Gc < 200.0 * max(hist.Gc_eff, 1e-30))
    if len(eok) == n:
        valid &= (eok > 0.5)

    K_upper = 100.0
    selected = 0.0
    source = 'none'

    if len(Kd) == n:
        mask = valid & np.isfinite(Kd) & (Kd > 1e-6) & (Kd < K_upper)
        if len(Da) == n:
            mask &= (Da > 0.02 * ell)
        if len(dfrac) == n:
            mask &= (dfrac > 1e-6) & (dfrac < 0.35)
        if np.any(mask):
            selected = float(np.nanmax(Kd[mask]))
            source = 'domain_event_window'

    # Do NOT use the global Wext/Da estimate as the default selected toughness.
    # It is extremely sensitive to tiny projected crack advance and has repeatedly
    # produced artificial 10-1000+ MPa√m spikes after the crack/process zone has
    # already crossed the measurement window.  Keep it in the CSV/JSON as a
    # diagnostic only.  If the domain integral is unavailable, use a simple
    # force-based LEFM estimate as a conservative elastic calibration metric.
    if selected <= 1e-6 and len(Kf) == n:
        fmask = valid & np.isfinite(Kf) & (Kf > 1e-6) & (Kf < K_upper)
        if len(dfrac) == n:
            # pre/near-onset force peak; exclude fully failed states
            fmask &= (dfrac < 0.50)
        if np.any(fmask):
            selected = float(np.nanmax(Kf[fmask]))
            source = 'force_LEFM_sanity_window'

    if selected <= 1e-6 and len(Kd) == n:
        dmask = valid & np.isfinite(Kd) & (Kd > 1e-6) & (Kd < K_upper)
        if np.any(dmask):
            selected = float(np.nanmax(Kd[dmask]))
            source = 'domain_sanity_window'

    def nanpeak(x):
        return float(np.nanmax(x)) if len(x) else 0.0

    return {
        'T_K': float(hist.T),
        'Gc_eff_J_m2': float(hist.Gc_eff),
        'Kc_input_MPa_sqrt_m': float(hist.Kc_input / 1e6),
        'KJ_selected_MPa_sqrt_m': float(selected),
        'KJ_selected_source': source,
        'KIc_like_selected_MPa_sqrt_m': float(selected) if mode_info.get('valid_KIc_like', False) else float('nan'),
        'J_tearing_peak_J_m2': nanpeak(Jr),
        'KJ_tearing_peak_MPa_sqrt_m': nanpeak(Kr),
        'energy_balance_bad_frac': float(np.mean(eok < 0.5)) if len(eok) else 0.0,
        'KJ_domain_peak_MPa_sqrt_m': nanpeak(Kd),
        'KJ_global_peak_MPa_sqrt_m': nanpeak(Kg),
        'K_force_peak_MPa_sqrt_m': nanpeak(Kf),
        'KJ_domain_final_MPa_sqrt_m': float(hist.KJ_final / 1e6),
        'Fmax_N_per_thickness': float(np.nanmax(np.abs(F))) if len(F) else 0.0,
        'Uapp_final_m': float(arr('Uapp')[-1]) if len(arr('Uapp')) else 0.0,
        'Wext_final_J_per_thickness': float(Wext[-1]) if len(Wext) else 0.0,
        'Wp_final_J_per_thickness': float(Wp[-1]) if len(Wp) else 0.0,
        'Dp_eff_final_J_per_thickness': float(arr('Dp_eff')[-1]) if len(arr('Dp_eff')) else 0.0,
        'Etough_final_J_per_thickness': float(arr('Etough')[-1]) if len(arr('Etough')) else 0.0,
        'Dtough_final_J_per_thickness': float(arr('Dtough')[-1]) if len(arr('Dtough')) else 0.0,
        'Emem_final_J_per_thickness': float(arr('Emem')[-1]) if len(arr('Emem')) else 0.0,
        'Dmem_final_J_per_thickness': float(arr('Dmem')[-1]) if len(arr('Dmem')) else 0.0,
        'energy_residual_final_J_per_thickness': float(arr('energy_residual')[-1]) if len(arr('energy_residual')) else 0.0,
        'energy_residual_rel_peak': float(np.nanmax(np.abs(arr('energy_residual_rel')))) if len(arr('energy_residual_rel')) else 0.0,
        'Wp_over_Wext_percent': (float(100 * Wp[-1] / Wext[-1])
                                  if len(Wp) and len(Wext) and Wext[-1] > 1e-12 else float('nan')),
        'rho_max_peak_m2': float(np.nanmax(rho)) if len(rho) else 0.0,
        'rho_p99_peak_m2': float(np.nanmax(arr('rho_p99'))) if len(arr('rho_p99')) else 0.0,
        'rho_cap_frac_peak': float(np.nanmax(arr('rho_cap_frac'))) if len(arr('rho_cap_frac')) else 0.0,
        'rho_gt_1e15_frac_peak': float(np.nanmax(arr('rho_gt_1e15_frac'))) if len(arr('rho_gt_1e15_frac')) else 0.0,
        'd_frac_final': float(dfrac[-1]) if len(dfrac) else 0.0,
        'Da_projected_final_m': float(Da[-1]) if len(Da) else 0.0,
        'Gc_local_max_J_m2': float(np.nanmax(Gc)) if len(Gc) else 0.0,
        'Gc_local_p99_J_m2': float(np.nanmax(arr('Gc_local_p99'))) if len(arr('Gc_local_p99')) else 0.0,
        'Gc_local_mean_front_J_m2': float(np.nanmax(arr('Gc_local_mean_front'))) if len(arr('Gc_local_mean_front')) else 0.0,
        'q_tough_max_J_m2': float(np.nanmax(arr('q_tough_max'))) if len(arr('q_tough_max')) else 0.0,
        'q_tough_p99_J_m2': float(np.nanmax(arr('q_tough_p99'))) if len(arr('q_tough_p99')) else 0.0,
        'q_tough_mean_front_J_m2': float(np.nanmax(arr('q_tough_mean_front'))) if len(arr('q_tough_mean_front')) else 0.0,
        'pz_emit_prob_max': float(np.nanmax(arr('pz_emit_prob_max'))) if len(arr('pz_emit_prob_max')) else 0.0,
        'pz_emission_hazard_max': float(np.nanmax(arr('pz_emission_hazard_max'))) if len(arr('pz_emission_hazard_max')) else 0.0,
        'pz_sigma_tip_eff_max_GPa': float(np.nanmax(arr('pz_sigma_tip_eff_max'))/1e9) if len(arr('pz_sigma_tip_eff_max')) else 0.0,
        'pz_sigma_back_max_GPa': float(np.nanmax(arr('pz_sigma_back_max'))/1e9) if len(arr('pz_sigma_back_max')) else 0.0,
        'pz_sigma_back_disl_max_GPa': float(np.nanmax(arr('pz_sigma_back_disl_max'))/1e9) if len(arr('pz_sigma_back_disl_max')) else 0.0,
        'pz_sigma_back_crack_max_GPa': float(np.nanmax(arr('pz_sigma_back_crack_max'))/1e9) if len(arr('pz_sigma_back_crack_max')) else 0.0,
        'pz_G_shield_max_J_m2': float(np.nanmax(arr('pz_G_shield_max'))) if len(arr('pz_G_shield_max')) else 0.0,
        'pz_G_stored_release_max_J_m2': float(np.nanmax(arr('pz_G_stored_release_max'))) if len(arr('pz_G_stored_release_max')) else 0.0,
        'pz_e_stored_max_J_m3': float(np.nanmax(arr('pz_e_stored_max'))) if len(arr('pz_e_stored_max')) else 0.0,
        'pz_Gc_net_min_J_m2': float(np.nanmin(arr('pz_Gc_net_min'))) if len(arr('pz_Gc_net_min')) else 0.0,
        'pz_drho_emit_max': float(np.nanmax(arr('pz_drho_emit_max'))) if len(arr('pz_drho_emit_max')) else 0.0,
        'pz_drho_rec_max': float(np.nanmax(arr('pz_drho_rec_max'))) if len(arr('pz_drho_rec_max')) else 0.0,
        'M_fract_max': float(np.nanmax(M)) if len(M) else 1.0,
        'sigma_y_mean_peak_GPa': float(np.nanmax(sigy) / 1e9) if len(sigy) else 0.0,
        'sigma_eq_max_peak_GPa': float(np.nanmax(seqmax) / 1e9) if len(seqmax) else 0.0,
        'yield_frac_peak': float(np.nanmax(yfrac)) if len(yfrac) else 0.0,
        'flow_Gtarget_eV_mean': float(gtarget[-1]) if len(gtarget) else 0.0,
        'flow_DG0_eV': float(dg0[-1]) if len(dg0) else 0.0,
        'flow_DGfloor_eV': float(dgfloor[-1]) if len(dgfloor) else 0.0,
        **mode_info,
        'n_steps_used': int(hist.n_steps_used),
        'snapshot_steps': [int(x) for x in sorted(hist.d_fields.keys())],
    }

def save_step_table(hist: SimulationHistory, filepath: str):
    """Save all scalar per-step diagnostics to CSV."""
    fields = [
        'step', 't', 'Uapp', 'Freact', 'Ftop', 'Fbot', 'Fpair_abs', 'Wext', 'Wext_top', 'Wext_pair', 'Wext_abs', 'Uel', 'Uel_drive', 'Uel_undegraded', 'Wp', 'Wp_tip', 'Dp_eff', 'Etough', 'Dtough', 'Epf_surf',
        'dWext', 'dWext_top', 'dWext_pair', 'dWext_abs', 'dUel', 'dUel_drive', 'dEpf', 'dWp_step', 'dDp_eff', 'dEtough', 'dDtough',
        'Emem', 'Dmem', 'Dfrac', 'dEmem', 'dDmem', 'dDfrac',
        'energy_residual', 'energy_residual_rel', 'energy_residual_absWext', 'energy_residual_topWext', 'energy_cumulative_residual', 'energy_cumulative_residual_absWext', 'energy_units_ratio_Uel_over_Wext_abs', 'energy_balance_ok',
        'J_tearing', 'KJ_tearing',
        'crack_len', 'Da_projected', 'Gamma_total', 'branch_factor',
        'J_domain', 'KJ_domain', 'J_global', 'KJ_global', 'K_force',
        'rho_mean', 'rho_p95', 'rho_p99', 'rho_max', 'rho_gt_1e14_frac', 'rho_gt_1e15_frac', 'rho_gt_1e16_frac', 'rho_cap_frac', 'dotep_mean', 'dotep_max',
        'dWp_requested', 'dWp_accepted',
        'dep_eq_requested_max', 'dep_eq_accepted_max', 'dep_eq_uncapped_max', 'dep_limited_frac',
        'thermo_scale_min', 'thermo_scale_mean',
        'thermo_admissible_frac', 'thermo_hazard_max',
        'thermo_substeps', 'thermo_dt_min', 'thermo_retry_count',
        'memory_energy_increment', 'memory_dissipation_increment',
        'memory_A_r_mean', 'memory_A_z_mean',
        'sigma_eq_mean', 'sigma_eq_max', 'sigma_y_min', 'sigma_y_mean', 'sigma_y_max',
        'sigma_T_min', 'sigma_T_mean', 'sigma_T_max', 'sigma_Peierls',
        'sigma_eq_over_sigma_y_max', 'yield_frac',
        'flow_dgamma_uncapped_max', 'flow_dgamma_cap', 'flow_cap_frac',
        'flow_phi_mean', 'flow_phi_max',
        'flow_Gtarget_eV_min', 'flow_Gtarget_eV_mean', 'flow_Gtarget_eV_max',
        'flow_DG0_eV', 'flow_DGfloor_eV', 'flow_vstar_ref_b3',
        'flow_status_zero_stress_frac', 'flow_status_solved_frac', 'flow_status_floor_limited_frac',
        'd_frac', 'plast_frac', 'rtip_mean', 'shield_mean',
        'M_fract_mean', 'M_fract_max', 'tip_emit_prob_mean',
        'tip_emit_prob_max', 'pz_emit_prob_mean', 'pz_emit_prob_max',
        'pz_mobility_prob_mean', 'pz_mobility_prob_max', 'pz_mobile_prob_max', 'pz_escape_prob_max',
        'pz_store_prob_mean', 'pz_store_prob_max', 'pz_storage_fraction_mean',
        'pz_mobility_hazard_max', 'pz_mobility_hazard_raw_max', 'pz_sigma_mobility_eff_max',
        'pz_sigma_tip_eff_max', 'pz_sigma_back_max', 'pz_sigma_back_disl_max', 'pz_sigma_back_mem_max', 'pz_sigma_back_crack_max', 'pz_G_shield_max', 'pz_G_stored_release_max', 'pz_G_stored_release_p99', 'pz_e_stored_max', 'pz_Gc_net_min', 'pz_Gc_net_p01', 'pz_G_app_max', 'pz_G_eff_max', 'pz_crack_R_max', 'pz_crack_R_p99', 'pz_H_eff_drive_max', 'pz_H_eff_drive_p99', 'pz_front_mask_frac', 'pz_H_eff_masked_max', 'pz_H_eff_unmasked_max', 'pz_crack_barrier_min_eV', 'pz_crack_hazard_max', 'pz_crack_prob_mean', 'pz_crack_prob_max', 'pz_crack_hazard_raw_max', 'pz_crack_B_mean', 'pz_crack_B_max', 'pz_crack_sigma_tip_max', 'pz_emit_B_max', 'pz_emit_rho_max', 'pz_emit_Gshield_max', 'pz_emission_hazard_max', 'pz_emission_hazard_raw_max',
        'pz_drho_emit_max', 'pz_drho_rec_max', 'pz_recovery_rate_max',
        'q_tough_mean_front', 'q_tough_mean_all', 'q_tough_p95', 'q_tough_p99', 'q_tough_max', 'dqtough_max',
        'toughening_weight_mean', 'toughening_energy_increment', 'toughening_dissipation_increment',
        'Gc_local_mean', 'Gc_local_p95', 'Gc_local_p99', 'Gc_local_mean', 'Gc_local_p95', 'Gc_local_p99', 'Gc_local_max', 'Gc_local_mean_front',
    ]
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    with open(filepath, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for st in hist.steps:
            writer.writerow({name: getattr(st, name, '') for name in fields})
    print(f"  Saved step table to {filepath}")



def save_results_summary_table(results: Dict[float, SimulationHistory], filepath: str):
    """Save one CSV row per temperature, including failure-mode classification."""
    rows=[]; keys=None
    for T in sorted(results.keys()):
        summ=history_summary(results[T]); rows.append(summ)
        if keys is None: keys=list(summ.keys())
    if not rows: return
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    with open(filepath,'w',newline='') as f:
        writer=csv.DictWriter(f, fieldnames=keys, extrasaction='ignore')
        writer.writeheader(); writer.writerows(rows)
    print(f"  Saved temperature/failure-mode summary to {filepath}")

def save_summary_json(hist: SimulationHistory, filepath: str):
    """Save compact run summary to JSON."""
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    with open(filepath, 'w') as f:
        json.dump(history_summary(hist), f, indent=2)
    print(f"  Saved summary to {filepath}")


def plot_field_snapshots(hist: SimulationHistory, outdir: str,
                         mesh_nodes: np.ndarray, mesh_elems: np.ndarray,
                         max_cols: int = 4):
    """Save PNG panel of crack damage, dislocation density, local Gc/work, and M_tip."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import matplotlib.tri as mtri
    except ImportError:
        print("  matplotlib not available, skipping field snapshots")
        return

    steps = sorted(hist.d_fields.keys())
    if not steps:
        print("  No saved field snapshots available")
        return

    if len(steps) > max_cols:
        idx = np.linspace(0, len(steps) - 1, max_cols, dtype=int)
        pick = [steps[i] for i in idx]
    else:
        pick = steps

    os.makedirs(outdir, exist_ok=True)
    nodes = mesh_nodes
    elems = mesh_elems
    tri = mtri.Triangulation(nodes[:, 0] * 1e3, nodes[:, 1] * 1e3, elems)

    # Decide third row: local Gc if available; otherwise log plastic work.
    has_Gc = any(st in hist.Gc_fields for st in pick)
    has_M = any(st in hist.M_fields for st in pick)
    rows = ['damage', 'rho', 'Gc' if has_Gc else 'Wp']
    if has_M:
        rows.append('M')
    # Elastic / hazard diagnostics: raw FEM principal stress, de-smeared tip
    # drive, accumulated first-passage action.  These make it visible whether
    # the crack drive is concentrated AHEAD of the tip (correct) or whether
    # damage is following some other field (e.g. Gc degradation).
    if any(st in hist.sig1_fields for st in pick):
        rows.append('sig1')
    if any(st in hist.sigtip_fields for st in pick):
        rows.append('sigtip')
    if any(st in hist.B_fields for st in pick):
        rows.append('B')

    nrow, ncol = len(rows), len(pick)
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.4 * ncol, 3.9 * nrow), squeeze=False,
                             constrained_layout=True)

    last_im = {}
    for j, step in enumerate(pick):
        d = hist.d_fields.get(step)
        rho = hist.rho_fields.get(step)
        wp = hist.wp_fields.get(step)
        Gc = hist.Gc_fields.get(step)
        M = hist.M_fields.get(step)

        for i, row in enumerate(rows):
            ax = axes[i, j]
            if row == 'damage' and d is not None:
                im = ax.tripcolor(tri, d, shading='flat', cmap='inferno', vmin=0, vmax=1, rasterized=True)
                title = 'damage d'
            elif row == 'rho' and rho is not None:
                vals = np.log10(np.maximum(rho, 1.0))
                im = ax.tripcolor(tri, vals, shading='flat', cmap='viridis', vmin=10, vmax=max(16, float(np.nanmax(vals))), rasterized=True)
                title = r'log10 rho'
            elif row == 'Gc' and Gc is not None:
                im = ax.tripcolor(tri, Gc, shading='flat', cmap='magma', rasterized=True)
                title = 'Gc local'
            elif row == 'Wp' and wp is not None:
                vals = np.log10(1.0 + np.maximum(wp, 0.0))
                im = ax.tripcolor(tri, vals, shading='flat', cmap='magma', rasterized=True)
                title = 'log10(1+Wp)'
            elif row == 'M' and M is not None:
                im = ax.tripcolor(tri, M, shading='flat', cmap='plasma', rasterized=True)
                title = 'M_tip'
            elif row == 'sig1' and hist.sig1_fields.get(step) is not None:
                vals = hist.sig1_fields[step] / 1e6  # MPa
                im = ax.tripcolor(tri, vals, shading='flat', cmap='magma', rasterized=True)
                title = 'sigma1 FEM (MPa)'
            elif row == 'sigtip' and hist.sigtip_fields.get(step) is not None:
                vals = hist.sigtip_fields[step] / 1e9  # GPa
                im = ax.tripcolor(tri, vals, shading='flat', cmap='magma', rasterized=True)
                title = 'sigma_tip drive (GPa)'
            elif row == 'B' and hist.B_fields.get(step) is not None:
                vals = np.log10(np.maximum(hist.B_fields[step], 1e-12))
                im = ax.tripcolor(tri, vals, shading='flat', cmap='viridis',
                                  vmin=-12, vmax=1, rasterized=True)
                title = 'log10 B_crack'
            else:
                ax.axis('off')
                continue
            last_im[row] = im
            ax.set_aspect('equal')
            ax.set_title(f'{title}, step {step}', fontsize=9)
            ax.set_xlabel('x (mm)')
            if j == 0:
                ax.set_ylabel('y (mm)')

    # One colorbar per row.
    for i, row in enumerate(rows):
        if row in last_im:
            fig.colorbar(last_im[row], ax=axes[i, :], shrink=0.78, pad=0.01)

    fig.suptitle(f'Field snapshots — T = {hist.T:.0f} K', fontsize=13)
    path = os.path.join(outdir, f'field_snapshots_{hist.T:.0f}K.png')
    fig.savefig(path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved field snapshots to {path}")

def plot_diagnostics(hist: SimulationHistory, outdir: str = None):
    """Generate diagnostic plots for a simulation run."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("  matplotlib not available, skipping plots")
        return

    if outdir and not os.path.exists(outdir):
        os.makedirs(outdir)

    Uapp = hist.get_Uapp() * 1e3  # mm
    T = hist.T

    # 1. Load-displacement
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(Uapp, hist.get_Freact() / 1e3, 'b-', linewidth=2)
    ax.set_xlabel('Applied opening (mm)')
    ax.set_ylabel('Reaction force (kN)')
    ax.set_title(f'Load-displacement (T = {T:.0f} K)')
    ax.grid(True, alpha=0.3)
    _savefig(fig, outdir, f'load_displacement_{T:.0f}K.png')

    # 2. Energetics
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(Uapp, hist.get_array('Wext'), label='$W_{ext}$')
    ax.plot(Uapp, hist.get_array('Uel'), label='$U_{el}$')
    ax.plot(Uapp, hist.get_array('Wp'), label='$W_p$')
    ax.plot(Uapp, hist.get_array('Epf_surf'), label='$E_{pf}$')
    ax.set_xlabel('Applied opening (mm)')
    ax.set_ylabel('Energy (J/m)')
    ax.set_title(f'Energy balance (T = {T:.0f} K)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    _savefig(fig, outdir, f'energetics_{T:.0f}K.png')

    # 3. J-integral comparison
    fig, ax = plt.subplots(figsize=(8, 5))
    KJ_domain = hist.get_array('KJ_domain') / 1e6
    KJ_global = hist.get_array('KJ_global') / 1e6
    try:
        K_force = hist.get_array('K_force') / 1e6
    except Exception:
        K_force = np.zeros_like(KJ_domain)
    ax.plot(Uapp, KJ_domain, 'b-', linewidth=2, label='$K_J$ (domain integral)')
    ax.plot(Uapp, KJ_global, 'r--', linewidth=1.5, label='$K_J$ (global energy)')
    ax.plot(Uapp, K_force, 'k-.', linewidth=1.5, label='$K$ (force LEFM)')
    ax.axhline(hist.Kc_input / 1e6, color='k', linestyle=':', alpha=0.5,
               label=f'$K_c^{{input}}$ = {hist.Kc_input/1e6:.2f}')
    ax.set_xlabel('Applied opening (mm)')
    ax.set_ylabel('$K_J$ (MPa·√m)')
    ax.set_title(f'Toughness measurement (T = {T:.0f} K)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    _savefig(fig, outdir, f'toughness_{T:.0f}K.png')

    # 4. Dislocation-density distribution diagnostics
    fig, ax = plt.subplots(figsize=(8, 5))
    rho_mean = hist.get_array('rho_mean')
    rho_p95 = hist.get_array('rho_p95')
    rho_p99 = hist.get_array('rho_p99')
    rho_max = hist.get_array('rho_max')
    rho_cap_frac = hist.get_array('rho_cap_frac')
    if len(rho_mean) and len(rho_max) and np.allclose(rho_mean, rho_max, rtol=1e-6, atol=0.0):
        ax.semilogy(Uapp, rho_mean, 'k--', linewidth=2, label=r'$\rho$ frozen / uniform')
        ax.text(0.02, 0.95, r'$\rho_{mean}=\rho_{max}$: no resolved evolution',
                transform=ax.transAxes, va='top', ha='left', fontsize=9,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.75))
    else:
        ax.semilogy(Uapp, rho_mean, label='mean')
        if len(rho_p95): ax.semilogy(Uapp, rho_p95, label='p95')
        if len(rho_p99): ax.semilogy(Uapp, rho_p99, label='p99')
        ax.semilogy(Uapp, rho_max, label='max', alpha=0.75)
    ax.set_xlabel('Applied opening (mm)')
    ax.set_ylabel(r'$\rho$ (m$^{-2}$)')
    ax.set_title(f'Dislocation density distribution (T = {T:.0f} K)')
    ax.grid(True, which='both', alpha=0.3)
    if len(rho_cap_frac) and np.nanmax(rho_cap_frac) > 0:
        ax2 = ax.twinx()
        ax2.plot(Uapp, 100*rho_cap_frac, 'm:', linewidth=1.5, label='cap fraction')
        ax2.set_ylabel(r'$\rho$ cap fraction (%)')
        lines1, labels1 = ax.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax.legend(lines1 + lines2, labels1 + labels2, loc='best')
    else:
        ax.legend(loc='best')
    _savefig(fig, outdir, f'dislocations_{T:.0f}K.png')

    # 4b. Plasticity / thermodynamic time-cone diagnostics
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(Uapp, hist.get_array('yield_frac'), label='yield fraction')
    ax.plot(Uapp, hist.get_array('plast_frac'), label='plastic fraction')
    ax.plot(Uapp, hist.get_array('thermo_admissible_frac'), label='thermo admissible fraction')
    ax.set_xlabel('Applied opening (mm)')
    ax.set_ylabel('Fraction')
    ax.set_ylim(bottom=0)
    ax.grid(True, alpha=0.3)
    ax2 = ax.twinx()
    Wext = hist.get_array('Wext')
    Wp = hist.get_array('Wp')
    with np.errstate(divide='ignore', invalid='ignore'):
        wp_pct = np.where(Wext > 1e-30, 100.0 * Wp / Wext, np.nan)
    ax2.plot(Uapp, wp_pct, 'k--', linewidth=1.5, label=r'$W_p/W_{ext}$ (%)')
    ax2.plot(Uapp, hist.get_array('thermo_hazard_max'), 'r:', linewidth=1.5, label='hazard max')
    ax2.set_ylabel('Work ratio (%) / hazard')
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, loc='best')
    ax.set_title(f'Plasticity / thermodynamic diagnostics (T = {T:.0f} K)')
    _savefig(fig, outdir, f'plasticity_thermo_{T:.0f}K.png')

    # 4c. Local toughening state diagnostics
    if len(hist.get_array('q_tough_max')) or len(hist.get_array('Gc_local_max')):
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(Uapp, hist.get_array('Gc_local_mean'), label=r'$G_c^{local}$ mean')
        ax.plot(Uapp, hist.get_array('Gc_local_mean_front'), label=r'$G_c^{local}$ front mean')
        ax.plot(Uapp, hist.get_array('Gc_local_p99'), label=r'$G_c^{local}$ p99')
        ax.plot(Uapp, hist.get_array('Gc_local_max'), label=r'$G_c^{local}$ max', alpha=0.75)
        ax.set_xlabel('Applied opening (mm)')
        ax.set_ylabel(r'$G_c^{local}$ (J/m$^2$)')
        ax.set_title(f'Retained process-zone toughening (T = {T:.0f} K)')
        ax.legend(loc='best')
        ax.grid(True, alpha=0.3)
        _savefig(fig, outdir, f'toughening_state_{T:.0f}K.png')

    # 5. Tip memory and amplification
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(Uapp, hist.get_array('rtip_mean') / max(hist.ell, 1e-30), label=r'$r_{tip}/\ell$ mean')
    ax.plot(Uapp, hist.get_array('shield_mean'), label=r'$z_{shield}$ mean')
    ax.plot(Uapp, hist.get_array('M_fract_max'), label=r'$M_{tip}$ max')
    ax.set_xlabel('Applied opening (mm)')
    ax.set_ylabel('Dimensionless')
    ax.set_title(f'Tip memory/amplification (T = {T:.0f} K)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    _savefig(fig, outdir, f'tip_memory_{T:.0f}K.png')

    plt.close('all')


def plot_toughness_vs_T(results: Dict[float, SimulationHistory], outdir: str = None):
    """Plot toughness vs temperature from multiple runs."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        return

    if outdir and not os.path.exists(outdir):
        os.makedirs(outdir)

    temps = sorted(results.keys())
    KJ_domain = []
    KJ_global = []
    KJ_selected = []
    Kc_input = []
    modes = []

    for T in temps:
        hist = results[T]
        KJ_d = hist.get_array('KJ_domain')
        KJ_domain.append(np.max(KJ_d) / 1e6 if len(KJ_d) > 0 else 0)
        # Global estimate
        KJ_g = hist.get_array('KJ_global')
        KJ_global.append(np.max(KJ_g) / 1e6 if len(KJ_g) > 0 else 0)
        summ = history_summary(hist)
        KJ_selected.append(float(summ.get('KJ_selected_MPa_sqrt_m', 0.0)))
        modes.append(str(summ.get('failure_mode', 'unknown')))
        Kc_input.append(hist.Kc_input / 1e6)

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(temps, KJ_domain, 'bo-', linewidth=2, markersize=8,
            label='$K_J$ (domain integral)')
    ax.plot(temps, KJ_global, 'rs--', linewidth=1.5, markersize=7,
            label='$K_J$ (global energy)')
    ax.plot(temps, Kc_input, 'k^:', linewidth=1, markersize=6,
            label='$K_c^{input}$ (Arrhenius Gc)')
    ax.plot(temps, KJ_selected, 'mo-.', linewidth=1.5, markersize=5,
            label='$K_{selected}$ (windowed)')
    for x, y, m in zip(temps, KJ_selected, modes):
        if m not in ('brittle_crack', 'no_failure_or_no_crack'):
            ax.annotate(m.replace('_','\n'), (x, max(y, 0.05)),
                        textcoords='offset points', xytext=(4, 5), fontsize=7, alpha=0.75)
    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel('$K$ (MPa·√m)')
    ax.set_title('Fracture toughness vs temperature')
    ax.legend()
    ax.grid(True, alpha=0.3)
    _savefig(fig, outdir, 'toughness_vs_temperature.png')

    plt.close('all')


def _savefig(fig, outdir, filename):
    """Save figure if outdir is provided."""
    if outdir:
        fig.savefig(os.path.join(outdir, filename), dpi=200, bbox_inches='tight')
