from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np

from arrhenius_fracture.bulk_state_v911 import (
    BulkPlasticityControllerV911,
    normalize_bulk_mode,
)
from arrhenius_fracture.mpz_parameterization_v911 import (
    apply_pt_dislocation_config,
    load_selected_row,
)


def _row(name="weakT"):
    root = Path(__file__).resolve().parents[1] / "mpz_v9_11_parameters"
    return load_selected_row(root / name / "spatial_promotion_manifest.csv", name)


def _mat():
    E = 410.0e9
    nu = 0.28
    return SimpleNamespace(
        E=E,
        nu=nu,
        G=E / (2.0 * (1.0 + nu)),
        b=2.74e-10,
        Tm=3695.0,
    )


def _disl(row):
    cfg = SimpleNamespace(
        k_store=np.sqrt(2.0),
        k_dyn=1.0,
        use_static_recovery=False,
        thermo_event_strain=1.0e-4,
        pt_equivalent_strain_factor=1.0 / np.sqrt(3.0),
    )
    apply_pt_dislocation_config(cfg, row)
    cfg.pt_mobile_fraction = 0.0
    return cfg


def test_bulk_mode_aliases_and_contract():
    assert normalize_bulk_mode("tip-source-only") == "tip_only"
    assert normalize_bulk_mode("bulk_uniform") == "bulk_same_pt_km"


def test_tip_only_disables_continuum_plasticity_without_changing_forest_field():
    row = _row("ceramic")
    ctl = BulkPlasticityControllerV911("tip_only", row)
    ep = np.zeros((3, 4))
    rho = np.full(4, 5.0e12)
    sigma = np.zeros((3, 4)); sigma[0] = 20.0e9
    out = ctl.update(ep.copy(), rho.copy(), sigma, _mat(), 700.0, 8.4, None, _disl(row), return_info=True)
    ep1, rho1, dot, info = out
    assert np.array_equal(ep1, ep)
    assert np.array_equal(rho1, rho)
    assert np.all(dot == 0.0)
    assert info["bulk_plasticity_mode"] == "tip_only"
    assert info["bulk_fixed_mobile_fraction_active"] is False


def test_exact_mobile_retained_exchange_conserves_density():
    mobile = np.array([1.0, 2.0, 3.0])
    retained = np.array([4.0, 5.0, 6.0])
    m1, r1, _, _ = BulkPlasticityControllerV911._exchange(
        mobile, retained,
        encounter_rate_s=np.array([2.0, 1.0, 0.5]),
        taylor_release_rate_s=np.array([0.5, 2.0, 1.0]),
        dt_s=3.0,
    )
    assert np.allclose(m1 + r1, mobile + retained)
    assert np.all(m1 >= 0.0)
    assert np.all(r1 >= 0.0)


def test_bulk_mode_uses_explicit_mobile_state_not_fixed_fraction():
    row = _row("weakT")
    ctl = BulkPlasticityControllerV911("bulk_same_pt_km", row)
    rho = np.full(5, 5.0e12)
    ctl._ensure_state(rho, 700.0)
    ctl.mobile_rho_m2[:] = 1.0e14
    ctl.retained_rho_m2[:] = rho

    ep = np.zeros((3, 5))
    sigma = np.zeros((3, 5)); sigma[0] = 12.0e9
    ep1, rho1, dot, info = ctl.update(
        ep, rho, sigma, _mat(), 700.0, 8.4, None, _disl(row), return_info=True
    )
    assert info["bulk_explicit_mobile_retained_state"] is True
    assert info["bulk_fixed_mobile_fraction_active"] is False
    assert np.all(rho1 >= 0.0)
    assert np.all(dot >= 0.0)
    assert ctl.mobile_rho_m2 is not None
    assert not np.allclose(ctl.mobile_rho_m2, 0.01 * ctl.retained_rho_m2)
    assert np.all(np.isfinite(ep1))


def test_existing_kocks_mecking_storage_adds_mobile_content():
    row = _row("DBTT")
    ctl = BulkPlasticityControllerV911("bulk_same_pt_km", row)
    rho = np.full(6, 5.0e12)
    ctl._ensure_state(rho, 700.0)
    before = ctl.mobile_rho_m2.copy()
    storage, _, _ = ctl._apply_storage_recovery(
        np.full(6, 1.0e-5), _mat(), _disl(row), 700.0, 1.0e-3
    )
    assert storage > 0.0
    assert np.all(ctl.mobile_rho_m2 > before)


def test_bulk_mode_keeps_manifest_specific_independent_pt_shapes():
    row = _row("weakT")
    cfg = _disl(row)
    assert np.isclose(cfg.pt_peierls_exp_a, row["peierls_exp_a"])
    assert np.isclose(cfg.pt_peierls_exp_n, row["peierls_exp_n"])
    assert np.isclose(cfg.pt_taylor_exp_a, row["taylor_exp_a"])
    assert np.isclose(cfg.pt_taylor_exp_n, row["taylor_exp_n"])
    assert cfg.pt_mobile_fraction == 0.0
