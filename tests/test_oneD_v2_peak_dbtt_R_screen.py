from __future__ import annotations

import json

import numpy as np
import pandas as pd

from arrhenius_fracture.emergent_gnd_contract_v913 import ACTIVE_CANDIDATE_PARAMETER_FIELDS
from scripts.run_oneD_v2_peak_dbtt_R_screen import (
    CONTROL_IDS,
    OBJECTIVES,
    SEARCH_FIELDS,
    canonical_parameter_json,
    parameter_sha256,
    population,
)


def test_candidate_identity_is_full_precision_and_stable():
    frame, _ = population("Peak", 2)
    row = frame.iloc[0]
    assert json.loads(canonical_parameter_json(row)).keys() == set(ACTIVE_CANDIDATE_PARAMETER_FIELDS)
    assert parameter_sha256(row) == row.parameter_sha256
    altered = row.copy()
    altered["cleave_G00_eV"] = np.nextafter(float(row.cleave_G00_eV), np.inf)
    assert parameter_sha256(altered) != row.parameter_sha256


def test_population_contains_control_and_one_material_identity_per_hash():
    for material in ("Peak", "DBTT"):
        frame, bounds = population(material, 2)
        assert CONTROL_IDS[material] in set(frame.candidate_id)
        assert frame.parameter_sha256.is_unique
        assert frame.canonical_parameter_json.is_unique
        assert bounds["candidate_count"] == len(frame)


def test_search_does_not_put_backend_reductions_or_common_physics_in_material_schema():
    forbidden = {
        "hazard_progress_scale", "reload_gap_threshold_m", "nominal_advance_m",
        "translation_length_scale", "persistent_backstress_scale", "mpz_length_m",
        "source_zone_length_m", "shielding_orientation_factors",
    }
    assert forbidden.isdisjoint(ACTIVE_CANDIDATE_PARAMETER_FIELDS)
    assert forbidden.isdisjoint(SEARCH_FIELDS)


def test_objectives_remain_separate_pareto_dimensions():
    assert len(OBJECTIVES) == len(set(OBJECTIVES)) == 8
    assert all(name.endswith("_objective") for name in OBJECTIVES)
