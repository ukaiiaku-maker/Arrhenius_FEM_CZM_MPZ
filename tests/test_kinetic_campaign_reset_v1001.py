from __future__ import annotations

from types import MethodType

from arrhenius_fracture.kinetic_campaign_czm import (
    CampaignKineticMPZState,
    DevelopedStateDiagnosticCZMFrontEngine,
    KineticCampaignCZMConfig,
    apply_pf_manifest_to_mpz_config,
)
from arrhenius_fracture.kinetic_campaign_czm_v1001 import (
    ResetSafeDevelopedStateDiagnosticCZMFrontEngine,
)
from arrhenius_fracture.moving_process_zone import MovingProcessZoneConfig
from arrhenius_fracture.pf_equivalent_material_manifest import load_material_manifest


def test_explicit_reset_reinstantiates_campaign_state(monkeypatch):
    manifest = load_material_manifest("weakT")
    kinetic = KineticCampaignCZMConfig()
    cfg = MovingProcessZoneConfig(length_m=100e-6, n_bins=200, n_systems=2)
    apply_pf_manifest_to_mpz_config(cfg, manifest, kinetic)

    def inherited_reset_stub(self):
        self.inherited_clock_reset_called = True
        self.mpz_state = object()  # represents the inherited v9.11 state
        self.N_em = 99.0
        self.B = 0.75
        self.a_adv = 4.0
        self.n_adv = 8
        self.W_emit = 3.0
        self.t = 2.0
        self.K_prev = 1.0
        self._lambda_c_prev = 1.0
        self._K_cleave_prev = 1.0

    monkeypatch.setattr(
        DevelopedStateDiagnosticCZMFrontEngine,
        "reset",
        inherited_reset_stub,
    )

    eng = object.__new__(ResetSafeDevelopedStateDiagnosticCZMFrontEngine)
    eng.manifest = manifest
    eng.kinetic_config = kinetic
    eng.mpz_config = cfg
    eng.b = 2.48e-10
    eng.G = 80e9
    eng._sync_compat = MethodType(
        lambda self: setattr(self, "N_em", float(self.mpz_state.retained_count)),
        eng,
    )

    eng.reset()

    assert eng.inherited_clock_reset_called
    assert isinstance(eng.mpz_state, CampaignKineticMPZState)
    assert eng.mpz_state.state_model == "kinetic_campaign_czm"
    assert eng.B == 0.0
    assert eng.a_adv == 0.0
    assert eng.n_adv == 0
    assert eng.W_emit == 0.0
    assert eng.micro_advance_total_m == 0.0
    assert eng.checkpoint_advance_total_m == 0.0


def test_campaign_state_reset_starts_from_temperature_independent_virgin_state(monkeypatch):
    manifest = load_material_manifest("DBTT")
    kinetic = KineticCampaignCZMConfig()
    cfg = MovingProcessZoneConfig(length_m=100e-6, n_bins=200, n_systems=2)
    apply_pf_manifest_to_mpz_config(cfg, manifest, kinetic)

    monkeypatch.setattr(
        DevelopedStateDiagnosticCZMFrontEngine,
        "reset",
        lambda self: setattr(self, "mpz_state", object()),
    )

    eng = object.__new__(ResetSafeDevelopedStateDiagnosticCZMFrontEngine)
    eng.manifest = manifest
    eng.kinetic_config = kinetic
    eng.mpz_config = cfg
    eng.b = 2.48e-10
    eng.G = 80e9
    eng._sync_compat = MethodType(
        lambda self: setattr(self, "N_em", float(self.mpz_state.retained_count)),
        eng,
    )

    eng.reset()
    first_capacity = eng.mpz_state.site_capacity.copy()
    first_available = eng.mpz_state.available_sites.copy()
    eng.mpz_state.available_sites[:] = 0.0
    eng.mpz_state.mobile[:] = 4.0
    eng.reset()

    assert (eng.mpz_state.site_capacity == first_capacity).all()
    assert (eng.mpz_state.available_sites == first_available).all()
    assert (eng.mpz_state.mobile == 0.0).all()
    assert (eng.mpz_state.retained == 0.0).all()
