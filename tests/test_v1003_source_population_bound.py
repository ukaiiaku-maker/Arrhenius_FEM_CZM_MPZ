import math

import pytest

from arrhenius_fracture.mode_i_first_passage_v10_0_3_progressive import (
    source_population_bound,
)


def test_source_bound_is_one_inventory_without_advance():
    assert source_population_bound(4.0, 0.0, 20.0e-6) == 4.0


def test_source_bound_includes_maximum_advance_refresh():
    capacity = 2.0 * 2.4387841773917582
    extension = 5.0e-6
    refresh_length = 28.24481540428955e-6
    expected = capacity * (1.0 + extension / refresh_length)
    assert source_population_bound(capacity, extension, refresh_length) == pytest.approx(
        expected, rel=1e-15
    )
    assert expected < 6.0


def test_source_bound_rejects_legacy_population_scale():
    capacity = 2.0 * 2.4387841773917582
    bound = source_population_bound(capacity, 5.0e-6, 28.24481540428955e-6)
    assert math.isfinite(bound)
    assert 93.88 > bound
