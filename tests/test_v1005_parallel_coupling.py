from __future__ import annotations

from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import pytest

from arrhenius_fracture.kinetic_campaign_czm import KineticCampaignCZMConfig
from arrhenius_fracture.kinetic_campaign_czm_v1005 import (
    ParallelOpeningEmissionCZMFrontEngine,
    TensorResolvedCampaignKineticMPZState,
)
from arrhenius_fracture.moving_process_zone import MovingProcessZoneConfig
from arrhenius_fracture.pf_equivalent_material_manifest import load_material_manifest
from arrhenius_fracture.tensor_resolved_coupling_v1005 import (
    TensorResolvedDriveConfig,
    TensorResolvedKineticCohesiveStepper,
    capture_tensor_resolved_drives,
    reset_tensor_drive_runtime,
)


def make_state():
    manifest = load_material_manifest("weakT")
    cfg = MovingProcessZoneConfig(length_m=100e-6, n_bins=100, n_systems=2)
    kinetic = KineticCampaignCZMConfig(backstress_scale=0.0)
    state = TensorResolvedCampaignKineticMPZState(
        cfg,
        manifest,
        b=2.48e-10,
        G_Pa=80e9,
        kinetic_cfg=kinetic,
    )
    return state, manifest, kinetic


def test_exact_emission_depletion_has_no_second_directional_multiplier():
    state, _manifest, _ = make_state()
    factors = np.array([1.0, 0.25])
    sigma = 2.0e9
    dt = 0.2
    available = state.available_sites.copy()

    class RecordingEmission:
        def __init__(self):
            self.last_stress = None

        def rate(self, stress, _temperature):
            self.last_stress = np.asarray(stress, dtype=float).copy()
            return np.array([2.0, 4.0])

    recording = RecordingEmission()
    state.manifest = SimpleNamespace(emission=recording)
    rates = np.array([2.0, 4.0])
    probability = 1.0 - np.exp(-rates * dt)
    expected = available * probability

    result = state.emit_exact(dt, sigma, 700.0, factors)

    assert np.array_equal(recording.last_stress, factors * sigma)
    assert np.allclose(result["dN_emit_per_system"], expected, rtol=1e-12, atol=1e-14)
    assert result["directional_multiplier_applied_after_hazard"] is False
    double_weighted = expected * factors
    assert not np.allclose(
        result["dN_emit_per_system"],
        double_weighted,
        rtol=1e-8,
        atol=1e-14,
    )


def test_tensor_drive_factors_are_not_normalized_or_clipped():
    factors = TensorResolvedCampaignKineticMPZState.drive_factors(
        np.array([0.2, 1.4]),
        2,
    )
    assert np.array_equal(factors, np.array([0.2, 1.4]))


def test_active_shielding_ignores_fit_derived_manifest_cap():
    state, manifest, kinetic = make_state()
    engine = object.__new__(ParallelOpeningEmissionCZMFrontEngine)
    engine.mpz_state = state
    engine.manifest = replace(manifest, max_K_shield_MPa_sqrt_m=1.0e-12)
    engine.kinetic_config = kinetic
    engine.G = 80e9
    engine.nu = 0.28
    engine.b = 2.48e-10
    engine.f = SimpleNamespace(r0=1e-6, c_blunt=manifest.c_blunt)
    engine._last_channels = {}

    state.retained[:, :5] = 1.0e6
    raw = engine._active_shielding_raw()
    active = engine.active_K_shielding()
    fitted_cap = engine.manifest.max_K_shield_MPa_sqrt_m * 1.0e6

    assert abs(raw) > fitted_cap
    assert active == pytest.approx(raw)


def _synthetic_mesh_and_stress():
    centers = []
    for angle_deg in (0.0, 45.0, -45.0):
        angle = np.deg2rad(angle_deg)
        direction = np.array([np.cos(angle), np.sin(angle)])
        for radius in (6e-6, 10e-6, 14e-6):
            centers.append(radius * direction)
    nodes = []
    elems = []
    for center in centers:
        start = len(nodes)
        eps = 0.2e-6
        nodes.extend([
            center + np.array([-eps, -eps]),
            center + np.array([eps, -eps]),
            center + np.array([0.0, eps]),
        ])
        elems.append([start, start + 1, start + 2])
    mesh = SimpleNamespace(
        nodes=np.asarray(nodes, dtype=float),
        elems=np.asarray(elems, dtype=int),
        area_e=np.full(len(elems), 0.5 * (0.4e-6) * (0.4e-6)),
    )
    sigma = 1.0e9
    sigma_gp = np.vstack([
        np.zeros(len(elems)),
        np.full(len(elems), sigma),
        np.zeros(len(elems)),
    ])
    damage = np.zeros(len(nodes))
    return mesh, sigma_gp, damage


def test_tensor_capture_resolves_slip_systems_from_fem_stress():
    mesh, sigma_gp, damage = _synthetic_mesh_and_stress()
    reset_tensor_drive_runtime(
        TensorResolvedDriveConfig(
            crystal_theta_deg=0.0,
            probe_radius_m=10e-6,
            sector_half_angle_deg=20.0,
            min_elements=2,
        )
    )
    record = capture_tensor_resolved_drives(
        mesh=mesh,
        sigma_gp=sigma_gp,
        damage=damage,
        crack_tip=np.zeros(2),
        crack_direction=np.array([1.0, 0.0]),
        KJ_Pa_sqrt_m=20e6,
    )
    assert record["tensor_resolved_drive_active"] is True
    assert record["opening_shape_factor"] == pytest.approx(1.0)
    assert len(record["slip_system_drive_factors"]) == 2
    assert np.allclose(record["slip_system_drive_factors"], [1.0, 1.0], rtol=1e-12)


def test_cohesive_stepper_uses_cached_tensor_factors_not_unit_weights():
    mesh, sigma_gp, damage = _synthetic_mesh_and_stress()
    reset_tensor_drive_runtime(
        TensorResolvedDriveConfig(
            crystal_theta_deg=0.0,
            probe_radius_m=10e-6,
            sector_half_angle_deg=20.0,
            min_elements=2,
        )
    )
    record = capture_tensor_resolved_drives(
        mesh=mesh,
        sigma_gp=sigma_gp,
        damage=damage,
        crack_tip=np.zeros(2),
        crack_direction=np.array([1.0, 0.0]),
        KJ_Pa_sqrt_m=20e6,
    )
    mechanics = {
        "K_open_Pa_sqrt_m": 20e6,
        "K_cleave_input_Pa_sqrt_m": 20e6,
        "slip_system_weights": np.ones(2),
    }
    _Kopen, _Kcleave, factors = TensorResolvedKineticCohesiveStepper._drives(mechanics)
    assert np.array_equal(factors, np.asarray(record["slip_system_drive_factors"]))
    assert mechanics["tensor_resolved_drive_active"] is True
    assert mechanics["drive_factor_normalization_or_clipping_active"] is False
