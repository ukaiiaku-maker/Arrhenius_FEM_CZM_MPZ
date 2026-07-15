from __future__ import annotations

from pathlib import Path

import numpy as np

from arrhenius_fracture.bulk_remesh_transfer_v911 import (
    install_bulk_remesh_transfer_patch,
)
from arrhenius_fracture.bulk_state_v911 import BulkPlasticityControllerV911
from arrhenius_fracture.mpz_parameterization_v911 import load_selected_row


def test_explicit_bulk_state_transfers_when_integration_point_count_changes():
    install_bulk_remesh_transfer_patch()
    root = Path(__file__).resolve().parents[1] / "mpz_v9_11_parameters"
    row = load_selected_row(
        root / "weakT" / "spatial_promotion_manifest.csv", "weakT"
    )
    ctl = BulkPlasticityControllerV911("bulk_same_pt_km", row)
    ctl.last_temperature_K = 700.0
    ctl.mobile_rho_m2 = np.array([1.0, 2.0, 3.0]) * 1.0e12
    ctl.retained_rho_m2 = np.array([10.0, 20.0, 30.0]) * 1.0e12
    transferred_retained = np.array([8.0, 12.0, 18.0, 22.0, 35.0]) * 1.0e12

    ctl._ensure_state(transferred_retained, 700.0)

    assert ctl.mobile_rho_m2.shape == transferred_retained.shape
    assert ctl.retained_rho_m2.shape == transferred_retained.shape
    assert np.array_equal(ctl.retained_rho_m2, transferred_retained)
    # Old global mobile/retained ratio was exactly 0.1.
    assert np.allclose(ctl.mobile_rho_m2, 0.1 * transferred_retained)
    summary = ctl.summary()
    assert summary["bulk_remesh_transfer_count"] == 1
    assert summary["bulk_mesh_change_rejected"] is False
    assert summary["bulk_remesh_transfer_old_size"] == 3
    assert summary["bulk_remesh_transfer_new_size"] == 5
