from __future__ import annotations

import numpy as np

from arrhenius_fracture.field_snapshots_v912 import map_mpz_density_to_elements


def test_front_local_mpz_inventory_maps_to_2d_elements():
    nodes = np.array([
        [0.0, -1e-6], [5e-6, -1e-6], [0.0, 1e-6],
        [5e-6, 1e-6], [10e-6, -1e-6], [10e-6, 1e-6],
    ])
    elems = np.array([[0, 1, 2], [1, 3, 2], [1, 4, 3], [4, 5, 3]])
    retained = np.zeros((1, 4))
    retained[0, 0] = 2.0
    mobile = np.zeros_like(retained)
    mobile[0, 1] = 1.0
    snap = {
        "nodes": nodes,
        "elems": elems,
        "mpz_front_states": [{
            "xy_m": [0.0, 0.0],
            "direction": [1.0, 0.0],
            "state": {
                "config": {"length_m": 10e-6, "blunting_length_m": 2e-6},
                "retained": retained.tolist(),
                "mobile": mobile.tolist(),
            },
        }],
    }
    mapped = map_mpz_density_to_elements(snap)
    assert mapped.shape == (len(elems),)
    assert np.all(mapped >= 0.0)
    assert np.max(mapped) > 0.0
    assert mapped[0] > mapped[-1]


def test_empty_mpz_inventory_maps_to_zero():
    snap = {
        "nodes": np.array([[0.0, 0.0], [1e-6, 0.0], [0.0, 1e-6]]),
        "elems": np.array([[0, 1, 2]]),
        "mpz_front_states": [],
    }
    assert np.allclose(map_mpz_density_to_elements(snap), 0.0)
