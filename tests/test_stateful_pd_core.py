from __future__ import annotations

import unittest
from collections import deque

import numpy as np

from arrhenius_fracture.config import ElasticProperties
from arrhenius_fracture.sn_geometry import BluntNotchGeometry, make_blunt_edge_notch_mesh
from arrhenius_fracture.stateful_peridynamics import StatefulPDConfig, StatefulPDPatch


class StatefulPDCoreTests(unittest.TestCase):
    def make_patch(self):
        geom = BluntNotchGeometry()
        mesh, _, root = make_blunt_edge_notch_mesh(
            geom, nx=12, ny=24, jitter=0.0, root_h_fine=60e-6, seed=1
        )
        cfg = StatefulPDConfig(
            patch_radius_m=0.40e-3,
            horizon_m=190e-6,
            boundary_shell_m=140e-6,
            initiation_radius_m=260e-6,
            initiation_taper_m=60e-6,
            initiation_back_extent_m=120e-6,
            root_seed_radius_m=220e-6,
            established_extent_m=100e-6,
        )
        patch = StatefulPDPatch(mesh, geom, root, ElasticProperties(), cfg)
        return mesh, patch

    def test_state_initialization_and_site_measure(self):
        _, patch = self.make_patch()
        state = patch.initial_state()
        self.assertGreater(len(patch.xy), 10)
        self.assertGreater(len(patch.bonds), len(patch.xy))
        self.assertTrue(np.allclose(state.available, 1.0))
        self.assertTrue(np.allclose(state.embryo + state.stable + state.inactive, 0.0))
        expected_sites = np.sum(patch.mean_candidate_sites)
        self.assertGreater(expected_sites, 1.0)
        self.assertGreater(np.sum(state.candidate_sites), 0)
        self.assertTrue(np.array_equal(state.available_sites, state.candidate_sites))
        self.assertTrue(np.all(state.candidate_sites[patch.boundary] == 0))
        self.assertTrue(np.all(patch.initiation_weight[patch.boundary] == 0.0))

    def test_graph_connectivity_on_forced_bond_path(self):
        _, patch = self.make_patch()
        state = patch.initial_state()
        start = int(np.argmin(np.linalg.norm(patch.xy - patch.root_xy[None, :], axis=1)))
        target = int(np.argmax(patch.xy[:, 0]))
        adjacency = [[] for _ in range(len(patch.xy))]
        for b, (i, j) in enumerate(patch.bonds):
            adjacency[i].append((j, b))
            adjacency[j].append((i, b))
        parent = {start: (None, None)}
        q = deque([start])
        while q and target not in parent:
            i = q.popleft()
            for j, b in adjacency[i]:
                if j not in parent:
                    parent[j] = (i, b)
                    q.append(j)
        self.assertIn(target, parent)
        node = target
        while parent[node][0] is not None:
            prev, bond = parent[node]
            state.bond_damage[bond] = 1.0
            node = prev
        extent, nconn = patch.connected_crack(state)
        self.assertGreater(nconn, 0)
        self.assertGreater(extent, 0.0)


    def test_stateful_birth_stabilization_update(self):
        mesh, patch = self.make_patch()
        patch.cfg.hit_memory_s = 1e-3
        patch.cfg.nu_stabilize_s = 1e3
        patch.cfg.nu_heal_s = 0.0
        state = patch.initial_state()

        class FlatBarrier:
            rate_prefactor = 1e3
            @staticmethod
            def deltaG_eV(sigma, T):
                return np.zeros_like(np.asarray(sigma, float))

        sigma_hist = np.zeros((4, 3, mesh.nn))
        diag = patch.update(
            state, FlatBarrier(), sigma_hist, 300.0, 1000.0, 10.0, 0.0,
            np.zeros(mesh.nn), np.zeros(mesh.nn), 0.0, np.zeros(mesh.nn),
            np.ones(len(patch.xy)), np.ones(len(patch.bonds)),
        )
        self.assertGreater(diag.expected_births_cumulative, 0.0)
        self.assertGreater(diag.expected_stable, 0.0)
        self.assertGreater(diag.realized_births_cumulative, 0)
        self.assertGreater(diag.realized_stable, 0)
        total = state.available + state.embryo + state.stable + state.inactive
        self.assertTrue(np.all(total <= 1.0 + 1e-12))
        self.assertTrue(np.all(total >= -1e-12))


    def test_no_bond_softening_without_realized_stable_defect(self):
        mesh, patch = self.make_patch()
        state = patch.initial_state()
        # Deliberately place a nonzero ensemble-average stable fraction and
        # growth state, but no realized stable objects.  The crack graph must
        # remain intact.
        state.stable[:] = 0.2
        state.growth[:] = 1.0

        class InactiveBarrier:
            rate_prefactor = 0.0
            @staticmethod
            def deltaG_eV(sigma, T):
                return np.ones_like(np.asarray(sigma, float)) * 10.0

        sigma_hist = np.zeros((4, 3, mesh.nn))
        sigma_hist[:, 0, :] = 5e9
        diag = patch.update(
            state, InactiveBarrier(), sigma_hist, 300.0, 1000.0, 1e6, 0.0,
            np.zeros(mesh.nn), np.zeros(mesh.nn), 0.0, np.zeros(mesh.nn),
            np.ones(len(patch.xy)), np.ones(len(patch.bonds)),
        )
        self.assertEqual(diag.realized_stable, 0)
        self.assertEqual(diag.max_bond_damage, 0.0)

    def test_continuous_time_hit_memory(self):
        mesh, patch = self.make_patch()
        patch.cfg.hit_memory_s = 1e-3
        state = patch.initial_state()

        class FlatBarrier:
            rate_prefactor = 1e3
            @staticmethod
            def deltaG_eV(sigma, T):
                return np.zeros_like(np.asarray(sigma, float))

        sigma_hist = np.zeros((4, 3, mesh.nn))
        patch.update(
            state, FlatBarrier(), sigma_hist, 300.0, 1000.0, 10.0, 0.0,
            np.zeros(mesh.nn), np.zeros(mesh.nn), 0.0, np.zeros(mesh.nn),
            np.ones(len(patch.xy)), np.ones(len(patch.bonds)),
        )
        active = patch.initiation_weight > 0.99
        self.assertTrue(np.any(active))
        # One raw attempt per cycle maps to p=1-exp(-1), then to the
        # bounded continuous delivery rate f*p.
        expected = (1.0 - np.exp(-1.0))
        self.assertTrue(np.allclose(state.hit_memory[active], expected, rtol=1e-3, atol=1e-3))

    def test_local_mechanics_is_finite(self):
        mesh, patch = self.make_patch()
        state = patch.initial_state()
        u = np.zeros(mesh.ndof)
        u[1::2] = 1e-4 * mesh.nodes[:, 1]
        ep = np.zeros((3, mesh.nn))
        upd, strain, fem_strain, amp, point_amp = patch.solve_local_mechanics(state, u, ep)
        self.assertTrue(np.all(np.isfinite(upd)))
        self.assertTrue(np.all(np.isfinite(strain)))
        self.assertTrue(np.all(np.isfinite(amp)))
        self.assertTrue(np.all(point_amp > 0.0))
        self.assertTrue(np.allclose(amp, 1.0))
        self.assertTrue(np.allclose(point_amp, 1.0))


if __name__ == "__main__":
    unittest.main()
