import copy
import numpy as np

from arrhenius_fracture.moving_process_zone import MovingProcessZoneConfig, MovingProcessZoneState
from arrhenius_fracture.config import make_emergent_config
from arrhenius_fracture import sharp_front as sf


def test_finite_site_emission_is_timestep_invariant_without_recovery():
    cfg = MovingProcessZoneConfig(n_bins=10, n_systems=1, source_sites_per_system=100.0,
                                  source_recovery_rate_s=0.0,
                                  glide_nu0_s=0.0, trap_nu0_s=0.0,
                                  detrap_nu0_s=0.0, retained_recovery_nu0_s=0.0,
                                  shielding_orientation_factors=(1.0,))
    a = MovingProcessZoneState(cfg)
    b = MovingProcessZoneState(cfg)
    # Same integrated per-site hazard H=2, committed in one or ten increments.
    a.evolve(1.0, 300.0, 0.0, 2.5e-10, emission_hazard_integral=2.0)
    for _ in range(10):
        b.evolve(0.1, 300.0, 0.0, 2.5e-10, emission_hazard_integral=0.2)
    assert np.isclose(a.emitted_total, b.emitted_total, rtol=1e-12, atol=1e-12)
    assert np.isclose(a.available_sites[0], b.available_sites[0], rtol=1e-12, atol=1e-12)
    assert a.emitted_total <= 100.0 + 1e-12


def test_moving_frame_translation_conserves_active_plus_wake():
    cfg = MovingProcessZoneConfig(length_m=1.0, n_bins=10, n_systems=1,
                                  shielding_orientation_factors=(1.0,))
    s = MovingProcessZoneState(cfg)
    s.retained[0] = np.arange(1.0, 11.0)
    before = s.retained_count
    out = s.advance(0.23)
    assert np.isclose(s.retained_count + out['wake_retained'], before, rtol=1e-12, atol=1e-12)


def test_branch_split_conserves_state():
    cfg = MovingProcessZoneConfig(n_bins=8, n_systems=2)
    s = MovingProcessZoneState(cfg)
    s.mobile[:] = 2.0
    s.retained[:] = 3.0
    s.accumulated_slip[:] = 4.0
    totals = (s.mobile.sum(), s.retained.sum(), s.accumulated_slip.sum(), s.site_capacity.sum())
    c = s.split(0.37)
    assert np.allclose([s.mobile.sum()+c.mobile.sum(), s.retained.sum()+c.retained.sum(),
                        s.accumulated_slip.sum()+c.accumulated_slip.sum(),
                        s.site_capacity.sum()+c.site_capacity.sum()], totals)


def test_direct_shielding_responds_to_position():
    cfg = MovingProcessZoneConfig(length_m=1e-6, n_bins=10, n_systems=1,
                                  shielding_orientation_factors=(1.0,))
    near = MovingProcessZoneState(cfg); far = MovingProcessZoneState(cfg)
    near.retained[0, 0] = 1.0
    far.retained[0, -1] = 1.0
    assert near.shielding_K(160e9, 0.28, 2.74e-10) > far.shielding_K(160e9, 0.28, 2.74e-10) > 0.0


def _mpz_args():
    return sf._build_parser().parse_args([
        '--mode','1d','--front-state-model','moving_pz','--temperatures','300',
        '--mpz-n-bins','12','--mpz-n-systems','1','--mpz-shielding-factors','1',
        '--mpz-source-sites-per-system','20','--mpz-glide-nu0-s','0',
        '--mpz-trap-nu0-s','0','--mpz-detrap-nu0-s','0',
        '--mpz-retained-recovery-nu0-s','0','--sigma-cap-GPa','0',
    ])


def test_engine_branch_split_and_geometry_rollback():
    args = _mpz_args()
    eng = sf.build_engine(args, make_emergent_config().material)
    eng.mpz_state.retained[0, 0] = 10.0
    eng.mpz_state.accumulated_slip[0, 0] = 5.0
    eng._sync_compat()
    child = eng.clone_split(0.25)
    assert np.isclose(eng.N_em + child.N_em, 10.0)
    pre = eng.mpz_state.state_dict()
    eng.B = 1.2
    renew = eng._renew(1.0)
    assert renew['fired']
    eng.restore_geometry_veto(1)
    assert np.allclose(eng.mpz_state.retained, np.asarray(pre['retained']))
    assert eng.B >= 1.0


def test_fatigue_uses_same_front_barriers_as_monotonic():
    from arrhenius_fracture.fatigue_v1 import FatigueWaveform
    args = _mpz_args()
    eng = sf.build_engine(args, make_emergent_config().material)
    wave = FatigueWaveform(Kmax=8.0e6, R=0.1, frequency_Hz=1000.0)
    out = eng.commit_fatigue_block(wave, 300.0, cycles=1.0e-3, n_phase=32)
    # The reported peak emission barrier is evaluated by the same FrontEngine
    # barrier object used by monotonic step(), not a fatigue-only parameter set.
    sig_peak = eng.sigma_tip(wave.Kmax)
    _, _, Ge = eng.lambda_emit(sig_peak, 300.0)
    from arrhenius_fracture.config import EV_TO_J
    assert np.isclose(out['G_emit_eV'], Ge / EV_TO_J, rtol=1e-12, atol=1e-12)
    assert out['front_state_model_code'] == 1.0


def test_mixed_mode_v8_preserves_moving_pz_subclass():
    from arrhenius_fracture import mixed_mode_first_passage_v8 as mm8
    args = _mpz_args()
    mat = make_emergent_config().material
    ctx = mm8.ProductionBackendControlContext(
        1.0, 0.0, 0.0, reference_cleavage_shape=1.0
    )
    wrapped = mm8._engine_factory(sf.build_engine, ctx, sf.FrontEngine)(args, mat)
    assert getattr(wrapped, 'state_model', None) == 'moving_pz'
    assert hasattr(wrapped, 'mpz_state')
    assert hasattr(wrapped, 'step_drives')
    assert 'CalibratedMixedModeEngine' in wrapped.__class__.__name__


def test_state_serialization_round_trip_preserves_full_inventory():
    cfg = MovingProcessZoneConfig(n_bins=7, n_systems=2,
                                  shielding_orientation_factors=(1.0, 0.5))
    s = MovingProcessZoneState(cfg)
    s.available_sites[:] = [3.0, 4.0]
    s.mobile[:] = np.arange(14.0).reshape(2, 7)
    s.retained[:] = 2.0
    s.accumulated_slip[:] = 5.0
    s.emitted_total = 17.0
    s.wake_retained_total = 8.0
    r = MovingProcessZoneState.from_state_dict(s.state_dict())
    assert np.allclose(r.available_sites, s.available_sites)
    assert np.allclose(r.mobile, s.mobile)
    assert np.allclose(r.retained, s.retained)
    assert np.allclose(r.accumulated_slip, s.accumulated_slip)
    assert r.emitted_total == s.emitted_total
    assert r.wake_retained_total == s.wake_retained_total
