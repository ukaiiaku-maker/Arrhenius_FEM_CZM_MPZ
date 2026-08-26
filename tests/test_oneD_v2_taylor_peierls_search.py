from __future__ import annotations

import json

from arrhenius_fracture.emergent_gnd_contract_v913 import (
    ACTIVE_CANDIDATE_PARAMETER_FIELDS,
)
from reduced_fracture_v2.taylor_peierls_contract import (
    FROZEN_ACTIVE_FIELDS,
    SEARCH_WHITELIST,
    barrier_hashes,
    bit_identical,
    validate_candidate,
)
from scripts.run_oneD_v2_taylor_peierls_search import (
    CONTROL_IDS,
    _controls,
    sobol_population,
    transformed_bounds,
)


def test_empirical_bounds_are_transformed_and_within_observations() -> None:
    bounds = transformed_bounds()
    assert set(bounds["fields"]) == set(SEARCH_WHITELIST)
    for field, definition in bounds["fields"].items():
        assert definition["observed_min"] <= definition["lower"]
        assert definition["lower"] < definition["upper"]
        assert definition["upper"] <= definition["observed_max"]
        if "entropy" in field:
            assert definition["transform"] == "LINEAR_ADDITIVE"
        else:
            assert definition["transform"] == "LOG10_POSITIVE"


def test_small_sobol_population_is_deterministic_and_fail_closed() -> None:
    bounds = transformed_bounds()
    controls = _controls()
    for material in ("Peak", "DBTT"):
        first = sobol_population(material, 2, bounds)
        second = sobol_population(material, 2, bounds)
        assert first.full_material_sha256.tolist() == second.full_material_sha256.tolist()
        assert first.candidate_id.tolist() == second.candidate_id.tolist()
        assert first.full_material_sha256.is_unique
        assert CONTROL_IDS[material] in set(first.candidate_id)
        control = controls[material]
        expected_barriers = barrier_hashes(control)
        for _, candidate in first.iterrows():
            validate_candidate(candidate, control)
            assert barrier_hashes(candidate) == expected_barriers
            assert set(json.loads(candidate.canonical_parameter_json)) == set(
                ACTIVE_CANDIDATE_PARAMETER_FIELDS
            )
            for field in FROZEN_ACTIVE_FIELDS:
                assert bit_identical(candidate[field], control[field])


def test_material_identity_is_provider_and_lifecycle_independent() -> None:
    population = sobol_population("Peak", 2, transformed_bounds())
    assert "provider" not in population.columns
    assert "lifecycle" not in population.columns
    assert "event_length" not in population.columns
