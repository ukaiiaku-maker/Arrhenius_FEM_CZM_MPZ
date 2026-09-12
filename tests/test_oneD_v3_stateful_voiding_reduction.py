from __future__ import annotations

from dataclasses import fields, replace
import json
import math
from pathlib import Path
import struct

import pytest
from setuptools import find_packages

from reduced_fracture_v2.state import ReducedState
from reduced_fracture_v3 import (
    AlignedVoidMechanicsBackend,
    HazardClock,
    MechanicsMapDomainError,
    MechanicsRepresentation,
    MonotonicPrototypePlan,
    OneDStateV3,
    OneDVoidState,
    ResponseMapNode,
    TensorProductResponseMap,
    TopologyRegime,
    VoidPhase,
    deserialize,
    export_v5_oracle,
    fingerprint,
    held_out_geometry_family,
    no_void_downgrade,
    no_void_upgrade,
    nucleate_downstream_child,
    promote_cavity,
    run_monotonic_prototype,
    seed_cavity,
    serialize,
)


def initial_state() -> OneDStateV3:
    fracture = ReducedState(tip_radius_m=1.5e-9)
    void = OneDVoidState(
        site_id="site-1",
        x_void_m=50.0e-6,
        candidate_weight=1.0,
        birth=HazardClock(0.0, 0.25),
        stabilization=HazardClock(0.0, 0.4),
        healing=HazardClock(0.0, 0.8),
        ligament=HazardClock(0.0, 0.3),
        downstream=HazardClock(0.0, 0.6),
        available_inventory_area_m2=1.0e-8,
        rng_state={"bit_generator": "PCG64", "state": {"state": 10, "inc": 3}},
    )
    return OneDStateV3(fracture, void)


def plan() -> MonotonicPrototypePlan:
    return MonotonicPrototypePlan(
        renewed_birth_threshold=0.75,
        rng_state_after_birth_renewal={
            "bit_generator": "PCG64",
            "state": {"state": 20, "inc": 3},
        },
        seed_radius_m=2.0e-6,
        promotion_radius_m=3.0e-6,
        connection_radius_m=5.0e-6,
        root_tip_position_m=40.0e-6,
        ligament_candidate_id="aligned-ligament",
        downstream_candidate_id="far-surface",
        child_length_m=4.0e-6,
        child_continuation_m=3.0e-6,
        renewed_child_fracture_state=ReducedState(tip_radius_m=2.0e-9),
    )


def test_01_phase_enumeration_is_exact_v5_crosswalk():
    assert {phase.value for phase in VoidPhase} == {
        "AVAILABLE_SITE",
        "EMBRYO",
        "HEALED_SITE",
        "STABLE_SUBGRID_VOID",
        "RESOLVED_VOID",
        "CONNECTED_VOID",
        "DOWNSTREAM_FRONT_ACTIVE",
        "MERGED_OR_CONSUMED",
    }


def test_02_no_void_upgrade_is_bitwise_transparent_to_v2_state_updates():
    direct = ReducedState()
    wrapped = no_void_upgrade(direct)
    assert wrapped.fracture_state is direct
    for increments in ((0.125, 0.25, 1.0), (0.375, 0.5, 2.0)):
        direct = direct.with_hazard(*increments)
        wrapped = replace(wrapped, fracture_state=wrapped.fracture_state.with_hazard(*increments))
    restored = no_void_downgrade(wrapped)
    assert restored == direct
    for item in fields(ReducedState):
        left, right = getattr(restored, item.name), getattr(direct, item.name)
        if isinstance(left, float):
            assert struct.pack(">d", left) == struct.pack(">d", right)


def test_03_checkpoint_roundtrip_and_fingerprint_preserve_owned_state():
    state = initial_state()
    payload = serialize(state)
    restored = deserialize(payload)
    assert serialize(restored) == payload
    assert fingerprint(restored) == fingerprint(state)
    assert restored.void_state.rng_state == state.void_state.rng_state


def test_04_invalid_transition_rolls_back_by_immutability_and_propagates():
    state = initial_state()
    before = fingerprint(state)
    with pytest.raises(ValueError, match="inventory"):
        seed_cavity(
            replace(
                state,
                void_state=replace(
                    state.void_state,
                    phase=VoidPhase.STABLE_SUBGRID_VOID,
                    mechanics_representation=MechanicsRepresentation.SUBGRID,
                ),
            ),
            radius_m=1.0,
        )
    assert fingerprint(state) == before


def test_05_promotion_changes_only_phase_representation_and_lineage():
    state = initial_state()
    stable = replace(
        state,
        void_state=replace(
            state.void_state,
            phase=VoidPhase.STABLE_SUBGRID_VOID,
            mechanics_representation=MechanicsRepresentation.SUBGRID,
        ),
    )
    seeded = seed_cavity(stable, 2.0e-6)
    promoted = promote_cavity(seeded)
    before, after = seeded.void_state, promoted.void_state
    assert after.phase == VoidPhase.RESOLVED_VOID
    assert after.mechanics_representation == MechanicsRepresentation.SURROGATE_ACTIVE
    for name in (
        "radius_m",
        "area_m2",
        "birth",
        "stabilization",
        "healing",
        "ligament",
        "downstream",
        "available_inventory_area_m2",
        "consumed_inventory_area_m2",
        "rng_state",
        "length_ledgers",
    ):
        assert getattr(after, name) == getattr(before, name)


def test_06_bounded_monotonic_chain_preserves_events_and_length_ledgers():
    final, trace = run_monotonic_prototype(initial_state(), plan())
    assert [row["operation"] for row in trace] == [
        "INITIAL",
        "BIRTH_HIT",
        "EMBRYO",
        "STABILIZED",
        "INITIAL_CAVITY_SEED",
        "SUBGRID_GROWTH",
        "SURROGATE_PROMOTION",
        "RESOLVED_GROWTH",
        "LIGAMENT_FIRST_PASSAGE",
        "DOWNSTREAM_FIRST_PASSAGE",
        "CHILD_CONTINUATION",
    ]
    void = final.void_state
    assert void.phase == VoidPhase.DOWNSTREAM_FRONT_ACTIVE
    assert void.active_tip_position_m == pytest.approx(62.0e-6)
    assert void.length_ledgers.fractured_material_m == pytest.approx(12.0e-6)
    assert void.length_ledgers.free_span_m == pytest.approx(10.0e-6)
    assert void.length_ledgers.front_advance_m == pytest.approx(22.0e-6)
    assert final.fracture_state.tip_radius_m != void.radius_m
    assert void.radius_m == 5.0e-6


def test_07_connected_interval_has_no_active_tip_and_void_is_not_fractured_length():
    _, trace = run_monotonic_prototype(initial_state(), plan())
    connected = next(row for row in trace if row["operation"] == "LIGAMENT_FIRST_PASSAGE")
    assert connected["phase"] == VoidPhase.CONNECTED_VOID.value
    assert connected["active_tip_position_m"] is None
    assert connected["length_ledgers"]["fractured_material_m"] == pytest.approx(5.0e-6)
    assert connected["length_ledgers"]["free_span_m"] == 0.0
    assert connected["length_ledgers"]["connected_free_surface_extent_m"] == pytest.approx(10.0e-6)


def nodes(regime: TopologyRegime) -> tuple[ResponseMapNode, ...]:
    result = []
    for ligament in (1.0, 2.0):
        for radius in (0.1, 0.2, 0.3):
            value = 2.0 * ligament + 3.0 * radius
            result.append(
                ResponseMapNode(
                    regime=regime,
                    ligament_over_radius=ligament,
                    radius_over_scale=radius,
                    delta_G_per_load2=None if regime == TopologyRegime.CONNECTED_CAVITY else value,
                    delta_sigma_nn_per_load=value if regime != TopologyRegime.DOWNSTREAM_CHILD else None,
                    delta_sigma_tt_per_load=2.0 * value if regime != TopologyRegime.DOWNSTREAM_CHILD else None,
                    delta_tau_nt_per_load=0.0 if regime != TopologyRegime.DOWNSTREAM_CHILD else None,
                    delta_sigma_m_per_load=1.5 * value if regime != TopologyRegime.DOWNSTREAM_CHILD else None,
                    delta_compliance=0.1 * value,
                    delta_potential_energy_per_load2=0.2 * value,
                    source_identity=f"oracle:{ligament}:{radius}",
                )
            )
    return tuple(result)


def test_08_maps_are_topology_specific_additive_and_connected_has_no_G():
    pre = TensorProductResponseMap(nodes(TopologyRegime.PRECONNECTION), TopologyRegime.PRECONNECTION)
    connected = TensorProductResponseMap(
        nodes(TopologyRegime.CONNECTED_CAVITY), TopologyRegime.CONNECTED_CAVITY
    )
    backend = AlignedVoidMechanicsBackend(
        {TopologyRegime.PRECONNECTION: pre, TopologyRegime.CONNECTED_CAVITY: connected}
    )
    response = backend.evaluate(
        TopologyRegime.PRECONNECTION,
        ligament_over_radius=1.5,
        radius_over_scale=0.2,
        load=2.0,
        no_void_G_kinetic=10.0,
        no_void_compliance=1.0,
        no_void_potential_energy=2.0,
        isolated_cavity_tensor_per_load=(0.0, 0.0, 0.0, 0.0),
    )
    assert response.G_kinetic == pytest.approx(24.4)
    assert response.cavity_tensor == pytest.approx((7.2, 14.4, 0.0, 10.8))
    dormant = backend.evaluate(
        TopologyRegime.CONNECTED_CAVITY,
        ligament_over_radius=1.5,
        radius_over_scale=0.2,
        load=2.0,
        no_void_G_kinetic=None,
        no_void_compliance=1.0,
        no_void_potential_energy=2.0,
        isolated_cavity_tensor_per_load=(0.0, 0.0, 0.0, 0.0),
    )
    assert dormant.G_kinetic is None


def test_09_maps_fail_closed_outside_qualified_geometry_domain():
    model = TensorProductResponseMap(nodes(TopologyRegime.PRECONNECTION), TopologyRegime.PRECONNECTION)
    with pytest.raises(MechanicsMapDomainError, match="outside qualified domain"):
        model.evaluate(
            ligament_over_radius=0.9,
            radius_over_scale=0.2,
            load=1.0,
            no_void_G_kinetic=1.0,
            no_void_compliance=1.0,
            no_void_potential_energy=1.0,
            isolated_cavity_tensor_per_load=(0.0, 0.0, 0.0, 0.0),
        )


def test_10_no_void_limit_is_exact_and_owns_no_cavity_tensor():
    response = AlignedVoidMechanicsBackend.no_void_response(
        TopologyRegime.PRECONNECTION,
        G_kinetic=3.25,
        compliance=4.5,
        potential_energy=6.75,
    )
    assert response.G_kinetic == 3.25
    assert response.compliance == 4.5
    assert response.potential_energy == 6.75
    assert response.cavity_tensor is None
    assert response.qualification == "EXACT_NO_VOID_LIMIT"


def test_11_holdout_removes_complete_radius_family_and_scores_interpolation():
    report = held_out_geometry_family(
        nodes(TopologyRegime.PRECONNECTION),
        regime=TopologyRegime.PRECONNECTION,
        radius_over_scale=0.2,
    )
    assert report["training_nodes"] == 4
    assert report["held_out_nodes"] == 2
    for key, value in report.items():
        if key.startswith("maximum_"):
            assert value < 1.0e-12


def test_12_oracle_export_preserves_source_identity_and_fails_closed_when_incomplete(tmp_path: Path):
    static = tmp_path / "case_rows.json"
    static.write_text(json.dumps({
        "implementation_git_sha": "source-static-sha",
        "rows": [{
            "case": "aligned",
            "configuration": {
                "crack_enabled": True,
                "cavity_enabled": True,
                "crack_path_m_requested": [[0.0, 0.0], [4.0e-5, 0.0]],
                "cavity_center_m": [5.0e-5, 0.0],
                "cavity_radius_m": 5.0e-6,
                "opening_m": 4.0e-7,
            },
            "observables": {
                "reaction_top_N_per_m": 2.0,
                "compliance_m2_per_N": 3.0,
                "stored_energy_J_per_m": 4.0,
            },
        }],
    }))
    checkpoints = tmp_path / "checkpoints"
    checkpoints.mkdir()
    (checkpoints / "available.json").write_text(json.dumps({
        "state_sha256": "accepted-state",
        "topology_transaction_model_id": "v12",
        "active_tip_ids": [],
        "crack_network": {"branches": []},
        "energy_ledgers": {"latest_fem_energy_J_per_m": 1.0},
        "v12_support_state": {"complete_crack_graph_fingerprint": "topology"},
        "production_void_state_sha256": "void-state",
        "production_void_state": {
            "sites": [{"phase": "AVAILABLE_SITE", "center_m": [5.0e-5, 0.0]}],
            "cavities": [],
            "length_ledgers": {},
        },
    }))
    output = tmp_path / "oracle.json"
    result = export_v5_oracle(
        static_evidence=static,
        checkpoint_directory=checkpoints,
        output=output,
    )
    assert result["readiness"]["classification"] == "BLOCKED_ORACLE_COORDINATES_INCOMPLETE"
    assert result["rows"][0]["source_implementation_sha"] == "source-static-sha"
    assert result["rows"][1]["accepted_state_identity"] == "accepted-state"
    assert json.loads(output.read_text()) == result


def test_13_committed_retained_oracle_export_records_missing_map_coordinates():
    path = Path("analysis_outputs/oneD_v3_stateful_voiding/v5_retained_oracle_export.json")
    payload = json.loads(path.read_text())
    assert payload["readiness"]["classification"] == "BLOCKED_ORACLE_COORDINATES_INCOMPLETE"
    assert payload["readiness"]["map_ready_records"] == 0
    assert not payload["readiness"]["preconnection_complete_geometry_grid"]
    assert not payload["readiness"]["all_three_topology_regimes_present"]
    assert payload["readiness"]["missing_field_counts"] == {
            "G_kinetic_used_J_per_m2": 39,
            "cavity_sigma_m_Pa": 39,
            "cavity_sigma_nn_Pa": 39,
            "cavity_sigma_tt_Pa": 39,
            "cavity_tau_nt_Pa": 39,
            "compliance_m2_per_N": 11,
            "potential_energy_J_per_m": 1,
    }
    assert payload["readiness"]["total_records"] == 39
    assert all(row["G_kinetic_used_J_per_m2"] is None for row in payload["rows"])
    assert not any(Path(row["source_path"]).is_absolute() for row in payload["rows"])


def test_14_downstream_child_rejects_r_tip_equal_to_void_radius():
    final, _ = run_monotonic_prototype(initial_state(), plan())
    connected = replace(
        final,
        void_state=replace(
            final.void_state,
            phase=VoidPhase.CONNECTED_VOID,
            mechanics_representation=MechanicsRepresentation.CONNECTED_CAVITY,
            active_tip_position_m=None,
        ),
    )
    with pytest.raises(ValueError, match="r_tip must remain distinct"):
        nucleate_downstream_child(
            connected,
            child_length_m=1.0e-6,
            renewed_fracture_state=ReducedState(tip_radius_m=final.void_state.radius_m),
            candidate_id="far-surface",
        )


def test_15_v2_and_v3_reduced_packages_are_installable():
    packages = set(find_packages(include=["arrhenius_fracture*", "reduced_fracture_v2*", "reduced_fracture_v3*"]))
    assert "reduced_fracture_v2" in packages
    assert "reduced_fracture_v3" in packages


def test_16_far_and_small_void_limits_recover_baselines_exactly():
    model = TensorProductResponseMap(nodes(TopologyRegime.PRECONNECTION), TopologyRegime.PRECONNECTION)
    for ligament, radius in ((math.inf, 0.2), (1.5, 0.0)):
        response = model.evaluate(
            ligament_over_radius=ligament,
            radius_over_scale=radius,
            load=2.0,
            no_void_G_kinetic=10.0,
            no_void_compliance=1.0,
            no_void_potential_energy=2.0,
            isolated_cavity_tensor_per_load=(1.0, 2.0, 3.0, 4.0),
        )
        assert response.G_kinetic == 10.0
        assert response.compliance == 1.0
        assert response.potential_energy == 2.0
        assert response.cavity_tensor == (2.0, 4.0, 6.0, 8.0)
        assert response.qualification == "EXACT_FAR_OR_SMALL_VOID_REDUCTION_LIMIT"
