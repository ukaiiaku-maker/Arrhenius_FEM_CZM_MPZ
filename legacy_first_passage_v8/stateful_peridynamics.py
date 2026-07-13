"""Stateful local peridynamic initiation patch for blunt-scratch fatigue.

The patch is deliberately separate from every legacy diffuse-fracture module.
It combines:

* physical candidate-site density and finite-memory multi-hit completion,
* reversible embryos with competing stabilization and healing,
* stable-defect growth,
* cohesive softening of a nonlocal bond network,
* graph-based crack connectivity.

The first implementation is a one-way local coupling: the intact global FEM
supplies boundary motion, plastic eigenstrain, residual stress and cyclic stress
history.  Bond loss redistributes deformation inside the peridynamic patch but
is not yet fed back to the global FEM stiffness.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.spatial import cKDTree
from scipy.special import expit, gammainc

from .config import EV_TO_J, KB, ElasticProperties
from .sn_geometry import BluntNotchGeometry
from .sn_intact_fem import lumped_nodal_area

KBEV = KB / EV_TO_J


@dataclass
class StatefulPDConfig:
    patch_radius_m: float = 0.45e-3
    horizon_m: float = 90e-6
    boundary_shell_m: float = 100e-6
    kernel_power: float = 2.0
    residual_bond_stiffness: float = 1e-7

    # Physical initiation support. Candidate sites are excluded from the
    # displacement-coupling shell and smoothly tapered to zero outside the
    # root process zone.
    initiation_radius_m: float = 240e-6
    initiation_taper_m: float = 60e-6
    initiation_back_extent_m: float = 60e-6

    # Candidate-site population and finite-memory completion.
    site_density_m2: float = 5.0e10
    hit_count: float = 3.0
    hit_memory_s: float = 1.0e-6
    birth_scale: float = 1.0

    # Reversible embryo transitions.
    nu_stabilize_s: float = 5.0e2
    nu_heal_s: float = 2.0e2
    stabilize_stress_Pa: float = 1.4e9
    stabilize_width_Pa: float = 2.5e8
    stabilize_plastic_gain: float = 2.0
    heal_return_fraction: float = 0.9

    # Stable-object growth and bond linkage.
    nu_grow_s: float = 2.0e-2
    grow_stress_Pa: float = 1.2e9
    grow_width_Pa: float = 2.5e8
    stable_count_scale: float = 2.0
    nu_link_s: float = 8.0e-3
    link_stress_Pa: float = 1.1e9
    link_width_Pa: float = 2.0e8
    link_orientation_power: float = 2.0
    link_orientation_floor: float = 0.05
    neighbor_link_gain: float = 1.5

    # Numerical event controls and crack definition.
    max_transition_probability: float = 0.08
    broken_damage: float = 0.95
    root_seed_radius_m: float = 120e-6
    established_extent_m: float = 240e-6
    pd_amplification_cap: float = 4.0
    amplification_damage_scale: float = 0.05
    strain_reference: float = 1e-8
    random_seed: int = 1
    softening_damage: float = 1e-3


@dataclass
class StatefulPDState:
    available: np.ndarray
    embryo: np.ndarray
    stable: np.ndarray
    inactive: np.ndarray
    candidate_sites: np.ndarray
    available_sites: np.ndarray
    embryo_sites: np.ndarray
    stable_sites: np.ndarray
    inactive_sites: np.ndarray
    born_sites_cumulative: np.ndarray
    healed_sites_cumulative: np.ndarray
    hit_memory: np.ndarray
    completion: np.ndarray
    growth: np.ndarray
    bond_damage: np.ndarray
    healed_cumulative: np.ndarray
    born_cumulative: np.ndarray
    cycles_first_embryo: Optional[float] = None
    cycles_first_stable: Optional[float] = None
    cycles_first_expected_embryo: Optional[float] = None
    cycles_first_expected_stable: Optional[float] = None
    cycles_first_softening: Optional[float] = None
    cycles_connected: Optional[float] = None


@dataclass
class PDUpdateDiagnostics:
    max_hit_memory: float
    max_completion: float
    max_embryo: float
    max_stable: float
    max_growth: float
    max_bond_damage: float
    broken_bonds: int
    connected_extent_m: float
    connected_bonds: int
    expected_embryos: float
    expected_births_cumulative: float
    expected_stable: float
    realized_embryos: int
    realized_births_cumulative: int
    realized_stable: int
    max_effective_stress_Pa: float
    max_pd_amplification: float
    max_rate_per_cycle: float
    max_hit_rate_s: float
    max_birth_rate_per_cycle: float
    expected_candidate_sites: float
    realized_candidate_sites: int


class StatefulPDPatch:
    def __init__(
        self,
        mesh,
        geom: BluntNotchGeometry,
        root_xy,
        mat: ElasticProperties,
        cfg: StatefulPDConfig,
    ):
        self.cfg = cfg
        self.geom = geom
        self.root_xy = np.asarray(root_xy, float)
        self.mat = mat

        dist = np.linalg.norm(mesh.nodes - self.root_xy[None, :], axis=1)
        self.global_nodes = np.where(dist <= cfg.patch_radius_m)[0]
        if len(self.global_nodes) < 8:
            raise RuntimeError("peridynamic patch contains too few points")
        self.xy = np.asarray(mesh.nodes[self.global_nodes], float).copy()
        self.area = lumped_nodal_area(mesh)[self.global_nodes]
        self.area = np.maximum(self.area, np.median(self.area[self.area > 0]) * 1e-3)

        radial = np.linalg.norm(self.xy - self.root_xy[None, :], axis=1)
        shell0 = max(cfg.patch_radius_m - cfg.boundary_shell_m, 0.5 * cfg.patch_radius_m)
        self.boundary = radial >= shell0
        if np.count_nonzero(self.boundary) < 4:
            order = np.argsort(radial)
            self.boundary[order[-max(4, len(order) // 8):]] = True
        self._update_initiation_weight()

        pairs = np.asarray(list(cKDTree(self.xy).query_pairs(cfg.horizon_m)), dtype=int)
        if pairs.size == 0:
            raise RuntimeError("peridynamic horizon produced no bonds")
        pairs = pairs.reshape(-1, 2)
        keep = self._filter_bonds_outside_initial_void(pairs)
        self.bonds = pairs[keep]
        if len(self.bonds) < len(self.xy):
            raise RuntimeError("peridynamic patch is under-connected; increase horizon")

        self._stiffness_scale = None
        self._update_bond_geometry(calibrate=True)
        self.incident = [[] for _ in range(len(self.xy))]
        for b, (i, j) in enumerate(self.bonds):
            self.incident[int(i)].append(b)
            self.incident[int(j)].append(b)

    def _update_initiation_weight(self, update_site_measure=True):
        """Build a root-localized physical site measure.

        The coupling shell carries Dirichlet data only and cannot host
        nucleation sites. A cosine taper avoids a mesh-sensitive hard cutoff.
        """
        cfg = self.cfg
        radial = np.linalg.norm(self.xy - self.root_xy[None, :], axis=1)
        r_outer = max(float(cfg.initiation_radius_m), 0.0)
        taper = max(float(cfg.initiation_taper_m), 0.0)
        r_inner = max(r_outer - taper, 0.0)
        w = np.zeros(len(self.xy), dtype=float)
        if r_outer <= 0.0:
            w[:] = 1.0
        elif taper <= 0.0:
            w[radial < r_outer] = 1.0
        else:
            w[radial <= r_inner] = 1.0
            mid = (radial > r_inner) & (radial < r_outer)
            q = (radial[mid] - r_inner) / max(r_outer - r_inner, 1e-30)
            w[mid] = 0.5 * (1.0 + np.cos(np.pi * q))
        x_min = self.root_xy[0] - max(float(cfg.initiation_back_extent_m), 0.0)
        w[self.xy[:, 0] < x_min] = 0.0
        w[self.boundary] = 0.0
        self.initiation_weight = np.clip(w, 0.0, 1.0)
        if update_site_measure or not hasattr(self, "mean_candidate_sites"):
            self.mean_candidate_sites = (
                max(float(cfg.site_density_m2), 0.0) * self.area * self.initiation_weight
            )

    def _filter_bonds_outside_initial_void(self, pairs: np.ndarray) -> np.ndarray:
        p0 = self.xy[pairs[:, 0]]
        p1 = self.xy[pairs[:, 1]]
        keep = np.ones(len(pairs), dtype=bool)
        a = max(float(self.geom.depth_a), 1e-30)
        b = max(float(self.geom.half_height_b), 1e-30)
        for t in (0.25, 0.5, 0.75):
            q = (1.0 - t) * p0 + t * p1
            inside = (q[:, 0] >= -1e-15) & ((q[:, 0] / a) ** 2 + (q[:, 1] / b) ** 2 < 1.0 - 1e-9)
            keep &= ~inside
        return keep

    def _update_bond_geometry(self, calibrate=False):
        i, j = self.bonds[:, 0], self.bonds[:, 1]
        dx = self.xy[j] - self.xy[i]
        L = np.linalg.norm(dx, axis=1)
        self.L = np.maximum(L, 1e-30)
        self.n = dx / self.L[:, None]
        q = np.clip(self.L / max(self.cfg.horizon_m, 1e-30), 0.0, 1.0)
        kernel = np.exp(-(q ** self.cfg.kernel_power)) * (1.0 - q) ** 2
        raw_w = kernel * self.area[i] * self.area[j] / np.maximum(self.L**2, 1e-30)
        if calibrate or self._stiffness_scale is None:
            eps = 1.0e-6
            uy = eps * (self.xy[:, 1] - np.mean(self.xy[:, 1]))
            utest = np.column_stack([np.zeros(len(self.xy)), uy])
            ext = np.einsum("bi,bi->b", utest[j] - utest[i], self.n)
            raw_energy = 0.5 * np.sum(raw_w * ext**2)
            Dyyyy = self.mat.E * (1.0 - self.mat.nu) / ((1.0 + self.mat.nu) * (1.0 - 2.0 * self.mat.nu))
            target_energy = 0.5 * Dyyyy * eps**2 * np.sum(self.area)
            self._stiffness_scale = target_energy / max(raw_energy, 1e-300)
        self.k0 = self._stiffness_scale * raw_w

    def update_geometry(self, mesh, root_xy=None):
        self.xy = np.asarray(mesh.nodes[self.global_nodes], float).copy()
        if root_xy is not None:
            self.root_xy = np.asarray(root_xy, float)
        self._update_initiation_weight(update_site_measure=False)
        self._update_bond_geometry(calibrate=False)

    def initial_state(self) -> StatefulPDState:
        npnt, nb = len(self.xy), len(self.bonds)
        # A spatial crack graph is a single specimen realization, not an
        # ensemble average.  Draw a fixed physical candidate-site population
        # once and then use binomial state transitions.  Mean-field fractions
        # are retained in parallel for diagnostics and calibration.
        self._rng = np.random.default_rng(int(self.cfg.random_seed))
        mean_sites = np.maximum(self.mean_candidate_sites, 0.0)
        candidate_sites = self._rng.poisson(mean_sites).astype(np.int64)
        if np.sum(candidate_sites) == 0 and np.sum(mean_sites) > 0.0:
            candidate_sites[int(np.argmax(mean_sites))] = 1
        return StatefulPDState(
            available=np.ones(npnt),
            embryo=np.zeros(npnt),
            stable=np.zeros(npnt),
            inactive=np.zeros(npnt),
            candidate_sites=candidate_sites.copy(),
            available_sites=candidate_sites.copy(),
            embryo_sites=np.zeros(npnt, dtype=np.int64),
            stable_sites=np.zeros(npnt, dtype=np.int64),
            inactive_sites=np.zeros(npnt, dtype=np.int64),
            born_sites_cumulative=np.zeros(npnt, dtype=np.int64),
            healed_sites_cumulative=np.zeros(npnt, dtype=np.int64),
            hit_memory=np.zeros(npnt),
            completion=np.zeros(npnt),
            growth=np.zeros(npnt),
            bond_damage=np.zeros(nb),
            healed_cumulative=np.zeros(npnt),
            born_cumulative=np.zeros(npnt),
        )

    def _assemble_spring_system(self, bond_damage, ep_node):
        npnt = len(self.xy)
        ndof = 2 * npnt
        rows, cols, vals = [], [], []
        rhs = np.zeros(ndof)
        i_all, j_all = self.bonds[:, 0], self.bonds[:, 1]
        exx, eyy, gxy = ep_node[0], ep_node[1], ep_node[2]
        for b, (i, j) in enumerate(zip(i_all, j_all)):
            n = self.n[b]
            nn = np.outer(n, n)
            k = self.k0[b] * max(1.0 - float(bond_damage[b]), self.cfg.residual_bond_stiffness)
            dofi = np.array([2 * i, 2 * i + 1])
            dofj = np.array([2 * j, 2 * j + 1])
            for aa in range(2):
                for bb in range(2):
                    v = k * nn[aa, bb]
                    rows.extend([dofi[aa], dofi[aa], dofj[aa], dofj[aa]])
                    cols.extend([dofi[bb], dofj[bb], dofi[bb], dofj[bb]])
                    vals.extend([v, -v, -v, v])
            epn_i = n[0] ** 2 * exx[i] + n[1] ** 2 * eyy[i] + n[0] * n[1] * gxy[i]
            epn_j = n[0] ** 2 * exx[j] + n[1] ** 2 * eyy[j] + n[0] * n[1] * gxy[j]
            e0 = self.L[b] * 0.5 * (epn_i + epn_j)
            f0 = k * e0 * n
            rhs[dofi] -= f0
            rhs[dofj] += f0
        K = sparse.csr_matrix((vals, (rows, cols)), shape=(ndof, ndof))
        diag_reg = max(float(np.max(self.k0)), 1.0) * 1e-12
        K = K + sparse.eye(ndof, format="csr") * diag_reg
        return K, rhs

    def solve_local_mechanics(self, state, fem_u_global, ep_node_global):
        """Equilibrate the damaged nonlocal patch under FEM shell motion."""
        up = np.asarray(fem_u_global, float).reshape(-1, 2)[self.global_nodes]
        ep = np.asarray(ep_node_global, float)[:, self.global_nodes]
        K, rhs = self._assemble_spring_system(state.bond_damage, ep)
        prescribed_nodes = np.where(self.boundary)[0]
        prescribed = np.zeros(2 * len(self.xy), dtype=bool)
        prescribed[2 * prescribed_nodes] = True
        prescribed[2 * prescribed_nodes + 1] = True
        uvec = np.zeros(2 * len(self.xy))
        uvec[2 * prescribed_nodes] = up[prescribed_nodes, 0]
        uvec[2 * prescribed_nodes + 1] = up[prescribed_nodes, 1]
        free = ~prescribed
        rhs_free = rhs[free] - K[np.ix_(free, prescribed)] @ uvec[prescribed]
        if np.any(free):
            uvec[free] = spsolve(K[np.ix_(free, free)], rhs_free)
        u = uvec.reshape(-1, 2)

        i, j = self.bonds[:, 0], self.bonds[:, 1]
        du = u[j] - u[i]
        ext_total = np.einsum("bi,bi->b", du, self.n)
        exx, eyy, gxy = ep[0], ep[1], ep[2]
        epn_i = self.n[:, 0] ** 2 * exx[i] + self.n[:, 1] ** 2 * eyy[i] + self.n[:, 0] * self.n[:, 1] * gxy[i]
        epn_j = self.n[:, 0] ** 2 * exx[j] + self.n[:, 1] ** 2 * eyy[j] + self.n[:, 0] * self.n[:, 1] * gxy[j]
        mech_strain = ext_total / self.L - 0.5 * (epn_i + epn_j)

        ufem = up
        fem_ext = np.einsum("bi,bi->b", ufem[j] - ufem[i], self.n) / self.L - 0.5 * (epn_i + epn_j)
        amp_raw = np.maximum(mech_strain, 0.0) / np.maximum(
            np.maximum(fem_ext, 0.0), self.cfg.strain_reference
        )
        amp_raw = np.clip(amp_raw, 0.25, self.cfg.pd_amplification_cap)
        # The global FEM already resolves the intact scratch concentration.
        # Activate the PD/FEM redistribution ratio only after cohesive damage
        # develops; this removes spurious intact-boundary amplification.
        gate = np.clip(
            np.asarray(state.bond_damage, float)
            / max(float(self.cfg.amplification_damage_scale), 1e-30),
            0.0,
            1.0,
        )
        amp = 1.0 + gate * (amp_raw - 1.0)
        acc = np.zeros(len(self.xy)); w = np.zeros(len(self.xy))
        for b, (ii, jj) in enumerate(self.bonds):
            acc[ii] += amp[b]; acc[jj] += amp[b]
            w[ii] += 1.0; w[jj] += 1.0
        point_amp = np.where(w > 0, acc / np.maximum(w, 1.0), 1.0)
        return u, mech_strain, fem_ext, amp, point_amp

    def _point_drivers(self, sigma_hist_global, point_amp):
        sig = np.asarray(sigma_hist_global, float)[:, :, self.global_nodes]
        sx, sy, txy = sig[:, 0], sig[:, 1], sig[:, 2]
        szz = self.mat.nu * (sx + sy)
        hydro = (sx + sy + szz) / 3.0
        savg = 0.5 * (sx + sy)
        rad = np.sqrt((0.5 * (sx - sy)) ** 2 + txy**2)
        s1 = savg + rad
        s1_eff = np.maximum(s1, 0.0) * point_amp[None, :]
        hydro_eff = np.maximum(hydro, 0.0) * point_amp[None, :]
        return s1_eff, hydro_eff, sig

    def _bond_traction(self, sigma_hist_global, bond_amp):
        sig = np.asarray(sigma_hist_global, float)[:, :, self.global_nodes]
        i, j = self.bonds[:, 0], self.bonds[:, 1]
        sb = 0.5 * (sig[:, :, i] + sig[:, :, j])
        nx, ny = self.n[:, 0], self.n[:, 1]
        tn = nx[None, :] ** 2 * sb[:, 0] + ny[None, :] ** 2 * sb[:, 1] + 2.0 * nx[None, :] * ny[None, :] * sb[:, 2]
        return np.max(np.maximum(tn, 0.0), axis=0) * bond_amp

    def preview_rates(
        self,
        state,
        crack_barrier,
        sigma_hist_global,
        T_K,
        frequency_Hz,
        state_shift_eV_global,
        sigma_back_global,
        chi,
        plastic_state_global=None,
        point_amp=None,
        bond_amp=None,
    ):
        if point_amp is None:
            point_amp = np.ones(len(self.xy))
        if bond_amp is None:
            bond_amp = np.ones(len(self.bonds))
        s1, hydro, _ = self._point_drivers(sigma_hist_global, point_amp)
        shift = np.asarray(state_shift_eV_global)[self.global_nodes]
        back = np.asarray(sigma_back_global)[self.global_nodes]
        sig_open = np.maximum(s1 - chi * back[None, :], 0.0)
        G = np.maximum(crack_barrier.deltaG_eV(sig_open, T_K) + shift[None, :], 1e-12)
        lam_hit = crack_barrier.rate_prefactor * np.exp(
            np.clip(-G / max(KBEV * T_K, 1e-30), -700.0, 0.0)
        )
        # Convert the atomic attempt count to a bounded probability of one
        # coarse-grained delivery event during a cycle, then represent that
        # cycle-level event stream as a continuous Poisson intensity. This
        # preserves the physical one-hit-per-cycle coarse graining while making
        # a sub-cycle memory time dimensionally consistent.
        mu_raw = np.mean(lam_hit, axis=0) / max(frequency_Hz, 1e-30)
        mu_hit = 1.0 - np.exp(-np.clip(mu_raw, 0.0, 700.0))
        hit_rate_s = max(frequency_Hz, 1e-30) * mu_hit
        completion = gammainc(self.cfg.hit_count, np.maximum(state.hit_memory, 0.0))
        completion_eq = gammainc(
            self.cfg.hit_count,
            np.maximum(hit_rate_s * max(self.cfg.hit_memory_s, 0.0), 0.0),
        )
        completion_bound = np.maximum(completion, completion_eq)
        site_weight = self.initiation_weight
        mu_birth = self.cfg.birth_scale * mu_hit * completion * site_weight
        mu_birth_bound = self.cfg.birth_scale * mu_hit * completion_bound * site_weight

        smax = np.max(s1, axis=0)
        zstab = (smax - self.cfg.stabilize_stress_Pa) / max(self.cfg.stabilize_width_Pa, 1e-30)
        if plastic_state_global is not None:
            zstab = zstab + self.cfg.stabilize_plastic_gain * np.asarray(plastic_state_global)[self.global_nodes]
        mu_stab = self.cfg.nu_stabilize_s / max(frequency_Hz, 1e-30) * expit(zstab)
        mu_heal = self.cfg.nu_heal_s / max(frequency_Hz, 1e-30) * expit(-zstab)
        zgrow = (smax - self.cfg.grow_stress_Pa) / max(self.cfg.grow_width_Pa, 1e-30)
        stable_activity = 1.0 - np.exp(
            -np.maximum(state.stable_sites.astype(float), 0.0)
            / max(self.cfg.stable_count_scale, 1e-30)
        )
        mu_grow = self.cfg.nu_grow_s / max(frequency_Hz, 1e-30) * expit(zgrow) * stable_activity

        tn = self._bond_traction(sigma_hist_global, bond_amp)
        zlink = (tn - self.cfg.link_stress_Pa) / max(self.cfg.link_width_Pa, 1e-30)
        i, j = self.bonds[:, 0], self.bonds[:, 1]
        stable_count = np.maximum(state.stable_sites.astype(float), 0.0)
        activity_node = (1.0 - np.exp(-stable_count / max(self.cfg.stable_count_scale, 1e-30))) * state.growth
        activity = 0.5 * (activity_node[i] + activity_node[j])
        orient = self.cfg.link_orientation_floor + (1.0 - self.cfg.link_orientation_floor) * np.abs(self.n[:, 1]) ** self.cfg.link_orientation_power
        mu_link = self.cfg.nu_link_s / max(frequency_Hz, 1e-30) * expit(zlink) * activity * orient
        effective_birth = state.available * mu_birth_bound
        effective_embryo_transition = state.embryo * (mu_stab + mu_heal)
        effective_grow = (1.0 - state.growth) * mu_grow
        effective_link = (1.0 - state.bond_damage) * mu_link
        max_rate = float(np.max(np.r_[effective_birth, effective_embryo_transition, effective_grow, effective_link]))
        return {
            "hit_rate_s": hit_rate_s,
            "mu_hit": mu_hit,
            "mu_birth": mu_birth,
            "mu_birth_bound": mu_birth_bound,
            "mu_stab": mu_stab,
            "mu_heal": mu_heal,
            "mu_grow": mu_grow,
            "mu_link": mu_link,
            "smax": smax,
            "tn": tn,
            "max_rate_per_cycle": max_rate,
        }

    def update(
        self,
        state: StatefulPDState,
        crack_barrier,
        sigma_hist_global,
        T_K,
        frequency_Hz,
        dN,
        cycles_old,
        state_shift_eV_global,
        sigma_back_global,
        chi,
        plastic_state_global,
        point_amp,
        bond_amp,
    ) -> PDUpdateDiagnostics:
        cfg = self.cfg
        expected_births_old = float(np.sum(np.maximum(state.born_cumulative, 0.0) * self.mean_candidate_sites))
        expected_stable_old = float(np.sum(np.maximum(state.stable, 0.0) * self.mean_candidate_sites))
        realized_births_old = int(np.sum(state.born_sites_cumulative))
        realized_stable_old = int(np.sum(state.stable_sites))
        max_damage_old = float(np.max(state.bond_damage))
        rates = self.preview_rates(
            state,
            crack_barrier,
            sigma_hist_global,
            T_K,
            frequency_Hz,
            state_shift_eV_global,
            sigma_back_global,
            chi,
            plastic_state_global,
            point_amp,
            bond_amp,
        )
        tau_s = max(float(cfg.hit_memory_s), 1e-30)
        dt_s = dN / max(float(frequency_Hz), 1e-30)
        decay = np.exp(-dt_s / tau_s)
        state.hit_memory = (
            state.hit_memory * decay
            + rates["hit_rate_s"] * tau_s * (1.0 - decay)
        )
        state.completion = gammainc(cfg.hit_count, np.maximum(state.hit_memory, 0.0))

        # Re-evaluate birth with the updated completion state.
        mu_birth = cfg.birth_scale * rates["mu_hit"] * state.completion
        p_birth = 1.0 - np.exp(-np.clip(mu_birth * dN, 0.0, 700.0))
        new_embryo = state.available * p_birth
        state.available -= new_embryo
        state.embryo += new_embryo
        state.born_cumulative += new_embryo

        p_birth_clip = np.clip(p_birth, 0.0, 1.0)
        born_sites = self._rng.binomial(state.available_sites, p_birth_clip).astype(np.int64)
        state.available_sites -= born_sites
        state.embryo_sites += born_sites
        state.born_sites_cumulative += born_sites

        mu_s, mu_h = rates["mu_stab"], rates["mu_heal"]
        mu_tot = mu_s + mu_h
        p_any = 1.0 - np.exp(-np.clip(mu_tot * dN, 0.0, 700.0))
        trans = state.embryo * p_any
        stabilized = trans * mu_s / np.maximum(mu_tot, 1e-300)
        healed = trans * mu_h / np.maximum(mu_tot, 1e-300)
        state.embryo -= stabilized + healed
        state.stable += stabilized
        state.available += cfg.heal_return_fraction * healed
        state.inactive += (1.0 - cfg.heal_return_fraction) * healed
        state.healed_cumulative += healed

        p_any_clip = np.clip(p_any, 0.0, 1.0)
        transitioned_sites = self._rng.binomial(state.embryo_sites, p_any_clip).astype(np.int64)
        p_stab_cond = np.divide(mu_s, np.maximum(mu_tot, 1e-300))
        stabilized_sites = self._rng.binomial(transitioned_sites, np.clip(p_stab_cond, 0.0, 1.0)).astype(np.int64)
        healed_sites = transitioned_sites - stabilized_sites
        returned_sites = self._rng.binomial(healed_sites, np.clip(cfg.heal_return_fraction, 0.0, 1.0)).astype(np.int64)
        state.embryo_sites -= transitioned_sites
        state.stable_sites += stabilized_sites
        state.available_sites += returned_sites
        state.inactive_sites += healed_sites - returned_sites
        state.healed_sites_cumulative += healed_sites

        realized_total = state.available_sites + state.embryo_sites + state.stable_sites + state.inactive_sites
        if np.any(realized_total != state.candidate_sites):
            raise RuntimeError("realized site-population conservation failure")

        state.available = np.maximum(state.available, 0.0)
        state.embryo = np.maximum(state.embryo, 0.0)
        state.stable = np.maximum(state.stable, 0.0)
        state.inactive = np.maximum(state.inactive, 0.0)
        occupied = state.available + state.embryo + state.stable + state.inactive
        over = occupied > 1.0
        if np.any(over):
            scale = 1.0 / occupied[over]
            state.available[over] *= scale
            state.embryo[over] *= scale
            state.stable[over] *= scale
            state.inactive[over] *= scale

        p_grow = 1.0 - np.exp(-np.clip(rates["mu_grow"] * dN, 0.0, 700.0))
        state.growth += (1.0 - state.growth) * p_grow
        state.growth = np.clip(state.growth, 0.0, 1.0)

        # Linkage is promoted by already-softened neighbors, which supplies a
        # spatial front/coalescence feedback without a mean-field topology scalar.
        node_support = np.zeros(len(self.xy)); node_degree = np.zeros(len(self.xy))
        for b, (i, j) in enumerate(self.bonds):
            node_support[i] += state.bond_damage[b]
            node_support[j] += state.bond_damage[b]
            node_degree[i] += 1.0; node_degree[j] += 1.0
        node_support /= np.maximum(node_degree, 1.0)
        i, j = self.bonds[:, 0], self.bonds[:, 1]
        front_support = 0.5 * (node_support[i] + node_support[j])
        mu_link = rates["mu_link"] * (1.0 + cfg.neighbor_link_gain * front_support)
        p_link = 1.0 - np.exp(-np.clip(mu_link * dN, 0.0, 700.0))
        old_damage = state.bond_damage.copy()
        state.bond_damage += (1.0 - state.bond_damage) * p_link
        state.bond_damage = np.clip(state.bond_damage, 0.0, 1.0)

        expected_births = float(np.sum(np.maximum(state.born_cumulative, 0.0) * self.mean_candidate_sites))
        expected_embryo = float(np.sum(np.maximum(state.embryo, 0.0) * self.mean_candidate_sites))
        expected_stable = float(np.sum(np.maximum(state.stable, 0.0) * self.mean_candidate_sites))
        realized_births = int(np.sum(state.born_sites_cumulative))
        realized_embryo = int(np.sum(state.embryo_sites))
        realized_stable = int(np.sum(state.stable_sites))
        cycles_new = cycles_old + dN

        if state.cycles_first_expected_embryo is None and expected_births >= 1.0:
            frac = (1.0 - expected_births_old) / max(expected_births - expected_births_old, 1e-300)
            state.cycles_first_expected_embryo = cycles_old + float(np.clip(frac, 0.0, 1.0)) * dN
        if state.cycles_first_expected_stable is None and expected_stable >= 1.0:
            frac = (1.0 - expected_stable_old) / max(expected_stable - expected_stable_old, 1e-300)
            state.cycles_first_expected_stable = cycles_old + float(np.clip(frac, 0.0, 1.0)) * dN

        # Realized event times are block-censored.  They are intentionally not
        # inferred from fractional expected populations.
        if state.cycles_first_embryo is None and realized_births_old < 1 <= realized_births:
            state.cycles_first_embryo = cycles_new
        if state.cycles_first_stable is None and realized_stable_old < 1 <= realized_stable:
            state.cycles_first_stable = cycles_new
        soft_thr = float(np.clip(cfg.softening_damage, 0.0, 1.0))
        if state.cycles_first_softening is None and np.max(state.bond_damage) >= soft_thr:
            frac = (soft_thr - max_damage_old) / max(float(np.max(state.bond_damage)) - max_damage_old, 1e-300)
            state.cycles_first_softening = cycles_old + float(np.clip(frac, 0.0, 1.0)) * dN

        extent, nconn = self.connected_crack(state)
        if state.cycles_connected is None and extent >= cfg.established_extent_m:
            state.cycles_connected = cycles_new

        self.last_rates = {
            k: np.asarray(v).copy()
            for k, v in rates.items()
            if isinstance(v, np.ndarray)
        }
        self.last_point_amp = np.asarray(point_amp, float).copy()
        self.last_bond_amp = np.asarray(bond_amp, float).copy()

        return PDUpdateDiagnostics(
            max_hit_memory=float(np.max(state.hit_memory)),
            max_completion=float(np.max(state.completion)),
            max_embryo=float(np.max(state.embryo)),
            max_stable=float(np.max(state.stable)),
            max_growth=float(np.max(state.growth)),
            max_bond_damage=float(np.max(state.bond_damage)),
            broken_bonds=int(np.count_nonzero(state.bond_damage >= cfg.broken_damage)),
            connected_extent_m=float(extent),
            connected_bonds=int(nconn),
            expected_embryos=expected_embryo,
            expected_births_cumulative=expected_births,
            expected_stable=expected_stable,
            realized_embryos=realized_embryo,
            realized_births_cumulative=realized_births,
            realized_stable=realized_stable,
            max_effective_stress_Pa=float(np.max(rates["smax"])),
            max_pd_amplification=float(max(np.max(point_amp), np.max(bond_amp))),
            max_rate_per_cycle=float(max(rates["max_rate_per_cycle"], np.max(mu_birth), np.max(mu_link))),
            max_hit_rate_s=float(np.max(rates["hit_rate_s"])),
            max_birth_rate_per_cycle=float(np.max(mu_birth)),
            expected_candidate_sites=float(np.sum(self.mean_candidate_sites)),
            realized_candidate_sites=int(np.sum(state.candidate_sites)),
        )

    def connected_crack(self, state: StatefulPDState):
        active = np.where(state.bond_damage >= self.cfg.broken_damage)[0]
        if len(active) == 0:
            return 0.0, 0
        mids = 0.5 * (self.xy[self.bonds[:, 0]] + self.xy[self.bonds[:, 1]])
        seed = active[np.linalg.norm(mids[active] - self.root_xy[None, :], axis=1) <= self.cfg.root_seed_radius_m]
        if len(seed) == 0:
            return 0.0, 0
        active_set = set(map(int, active))
        node_to_bonds = [[] for _ in range(len(self.xy))]
        for b in active:
            i, j = self.bonds[b]
            node_to_bonds[i].append(int(b)); node_to_bonds[j].append(int(b))
        seen = set(map(int, seed)); stack = list(seen)
        while stack:
            b = stack.pop()
            for node in self.bonds[b]:
                for nb in node_to_bonds[node]:
                    if nb in active_set and nb not in seen:
                        seen.add(nb); stack.append(nb)
        ids = np.fromiter(seen, dtype=int)
        extent = max(float(np.max(mids[ids, 0]) - self.root_xy[0]), 0.0)
        return extent, len(ids)

    def state_fields_global(self, state: StatefulPDState, nn_global: int):
        out = {}
        for name in ("available", "embryo", "stable", "hit_memory", "completion", "growth"):
            arr = np.full(nn_global, np.nan)
            arr[self.global_nodes] = getattr(state, name)
            out[name] = arr
        return out

    def plot_initiation_diagnostics(self, state, out_png: Path, title: str = ""):
        """Plot the fields controlling the location of the first event."""
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        rates = getattr(self, "last_rates", {})
        hit_rate = np.asarray(rates.get("hit_rate_s", np.zeros(len(self.xy))), float)
        birth_rate = np.asarray(rates.get("mu_birth", np.zeros(len(self.xy))), float)
        expected_birth_node = (
            np.maximum(state.born_cumulative, 0.0) * self.mean_candidate_sites
        )
        fields = [
            ("initiation weight", self.initiation_weight),
            ("candidate sites", state.candidate_sites.astype(float)),
            ("log10 hit rate (s^-1)", np.log10(np.maximum(hit_rate, 1e-300))),
            ("completion Q(K,Lambda)", state.completion),
            ("log10 birth intensity / cycle", np.log10(np.maximum(birth_rate, 1e-300))),
            ("expected cumulative births", expected_birth_node),
        ]
        fig, axes = plt.subplots(2, 3, figsize=(12.5, 7.4), constrained_layout=True)
        x = self.xy[:, 0] * 1e3
        y = self.xy[:, 1] * 1e3
        for ax, (label, vals) in zip(axes.flat, fields):
            sc = ax.scatter(x, y, c=vals, s=18, edgecolors="none")
            ax.scatter(
                x[self.boundary], y[self.boundary], s=7,
                facecolors="none", edgecolors="k", linewidths=0.3
            )
            ax.plot(self.root_xy[0] * 1e3, self.root_xy[1] * 1e3, "rx", ms=7)
            ax.set_aspect("equal")
            ax.set_title(label)
            ax.set_xlabel("x (mm)")
            ax.set_ylabel("y (mm)")
            fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.03)
        fig.suptitle(title or "Stateful PD initiation diagnostics")
        fig.savefig(Path(out_png), dpi=220)
        plt.close(fig)

    def plot_snapshot(self, state, out_png: Path, title: str = ""):
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.collections import LineCollection

        out_png = Path(out_png)
        out_png.parent.mkdir(parents=True, exist_ok=True)
        fig, ax = plt.subplots(figsize=(8.5, 5.5), constrained_layout=True)
        seg = np.stack([self.xy[self.bonds[:, 0]], self.xy[self.bonds[:, 1]]], axis=1) * 1e3
        lc = LineCollection(seg, array=state.bond_damage, linewidths=0.6, cmap="viridis")
        lc.set_clim(0.0, 1.0)
        ax.add_collection(lc)
        realized = state.stable_sites.astype(float)
        sc = ax.scatter(self.xy[:, 0] * 1e3, self.xy[:, 1] * 1e3, c=realized, s=12, cmap="magma", edgecolors="none")
        active = self.initiation_weight > 0.0
        ax.scatter(self.xy[active, 0] * 1e3, self.xy[active, 1] * 1e3, s=6, facecolors="none", edgecolors="tab:blue", linewidths=0.25)
        ax.scatter(self.xy[self.boundary, 0] * 1e3, self.xy[self.boundary, 1] * 1e3, s=5, facecolors="none", edgecolors="k", linewidths=0.3)
        ax.plot(self.root_xy[0] * 1e3, self.root_xy[1] * 1e3, "rx", ms=8)
        ax.set_aspect("equal")
        ax.autoscale()
        ax.set_xlabel("x (mm)"); ax.set_ylabel("y (mm)")
        ax.set_title(title or "Stateful peridynamic initiation patch")
        cb1 = fig.colorbar(lc, ax=ax, fraction=0.045, pad=0.03)
        cb1.set_label("bond cohesive damage")
        cb2 = fig.colorbar(sc, ax=ax, fraction=0.045, pad=0.10)
        cb2.set_label("realized stable defects / point")
        fig.savefig(out_png, dpi=220)
        plt.close(fig)
