from __future__ import annotations

import subprocess

import pandas as pd

from select_mpz_v9_10_4_quality_diversity import (
    SelectionConfig,
    build_response_table,
    select_quality_diverse,
)


def _details(candidate_ids):
    rows = []
    for index, candidate_id in enumerate(candidate_ids):
        for temperature in (300.0, 700.0, 900.0, 1200.0):
            rows.append(
                {
                    "candidate_id": candidate_id,
                    "T_K": temperature,
                    "K_init_proxy": 10.0 + index + 0.01 * temperature,
                    "K_plateau_proxy": 15.0 + index + 0.02 * temperature,
                    "delta_KR_proxy": 2.0 + 0.5 * index + 0.001 * temperature,
                    "early_rise_per_100um_proxy": 1.0 + index,
                    "plateau_rise_per_100um_proxy": 0.2 + 0.1 * index,
                }
            )
    return pd.DataFrame(rows)


def test_all_passers_are_hard_reserved_when_they_fit():
    candidates = pd.DataFrame(
        {
            "candidate_id": ["p0", "p1", "p2", "f0", "f1"],
            "objective": [1.0, 2.0, 3.0, 1.1, 1.2],
            "accepted_for_spatial_promotion": [True, True, True, False, False],
            "restart": [0, 0, 1, 2, 3],
            "cleave_G00_eV": [1.0, 1.1, 1.2, 8.0, 12.0],
            "emit_G00_eV": [2.0, 2.1, 2.2, 9.0, 13.0],
        }
    )
    selected, audit = select_quality_diverse(
        candidates,
        _details(candidates.candidate_id),
        SelectionConfig(count=4, preserve_restart_lineages=False),
    )
    assert {"p0", "p1", "p2"}.issubset(set(selected.candidate_id))
    assert audit["all_passers_fit_in_budget"] is True
    assert audit["all_passers_retained"] is True
    passer_reasons = selected.set_index("candidate_id").loc[
        ["p0", "p1", "p2"], "selection_reason"
    ]
    assert set(passer_reasons) == {"all_passers_reserve"}


def test_quality_diversity_avoids_objective_only_near_duplicates():
    candidates = pd.DataFrame(
        {
            "candidate_id": ["best", "duplicate", "different", "far"],
            "objective": [1.0, 1.01, 1.05, 1.2],
            "accepted_for_spatial_promotion": [True, True, True, True],
            "restart": [0, 0, 1, 2],
            "cleave_G00_eV": [1.0, 1.001, 4.0, 8.0],
            "emit_G00_eV": [2.0, 2.001, 6.0, 12.0],
        }
    )
    details = _details(candidates.candidate_id)
    # Make the first two response trajectories nearly identical as well.
    for metric in (
        "K_init_proxy",
        "K_plateau_proxy",
        "delta_KR_proxy",
        "early_rise_per_100um_proxy",
        "plateau_rise_per_100um_proxy",
    ):
        best = details.loc[details.candidate_id == "best", metric].to_numpy()
        details.loc[details.candidate_id == "duplicate", metric] = best + 1.0e-5
    selected, _ = select_quality_diverse(
        candidates,
        details,
        SelectionConfig(
            count=3,
            quality_reserve_fraction=0.34,
            quality_weight=0.10,
            preserve_restart_lineages=False,
        ),
    )
    chosen = set(selected.candidate_id)
    assert "best" in chosen
    assert {"different", "far"}.intersection(chosen)
    assert not {"best", "duplicate"}.issuperset(chosen)


def test_response_table_is_temperature_resolved():
    candidates = pd.DataFrame(
        {"candidate_id": ["a", "b"], "objective": [1.0, 2.0]}
    )
    table, columns = build_response_table(candidates, _details(["a", "b"]))
    assert len(table) == 2
    assert "response_K_init_proxy_300K" in columns
    assert "response_delta_KR_proxy_1200K" in columns


def test_selector_audit_does_not_claim_physics_changes():
    candidates = pd.DataFrame(
        {
            "candidate_id": ["a", "b"],
            "objective": [1.0, 2.0],
            "accepted_for_spatial_promotion": [True, False],
            "cleave_G00_eV": [1.0, 2.0],
        }
    )
    _, audit = select_quality_diverse(
        candidates,
        _details(candidates.candidate_id),
        SelectionConfig(count=1),
    )
    assert audit["constitutive_physics_modified"] is False
    assert audit["mechanics_closure"] == "v9_isotropic_moving_process_zone"
    assert audit["requires_spatial_promotion"] is True


def test_v9104_shell_runners_have_valid_syntax():
    for path in (
        "run_mpz_v9_10_4_dbtt_quality_diversity.sh",
        "run_mpz_v9_10_4_dbtt_spatial_promotion.sh",
    ):
        result = subprocess.run(
            ["bash", "-n", path], check=False, capture_output=True, text=True
        )
        assert result.returncode == 0, result.stderr
