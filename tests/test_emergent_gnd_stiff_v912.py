"""Regression tests for near-prefactor v9.12 Arrhenius rates."""
import numpy as np

from arrhenius_fracture.emergent_gnd_dbtt_v912 import (
    CandidateParameters,
    CommonPhysics,
    EmergentGNDState,
    ExpFloorSurface,
    PTMechanism,
)


def fast_candidate() -> CandidateParameters:
    cleavage = ExpFloorSurface(2.0, 0.0, 4.0e9, 0.0, 1.0, 1.0, 0.05)
    emission = ExpFloorSurface(0.05, 0.0, 1.0e9, 0.0, 2.0, 1.0, 0.01)
    return CandidateParameters(
        candidate_id="near_prefactor",
        cleavage=cleavage,
        emission=emission,
        peierls=PTMechanism(0.05, 0.0, 2.0, 1.0, 1.0e12),
        taylor=PTMechanism(0.10, 0.0, 2.0, 1.0, 1.0e11),
        rho_source0_m2=1.0e14,
        source_refresh_length_m=5.0e-6,
        taylor_corr_rho_c_m2=1.0e14,
        taylor_corr_scale=1.0,
        recovery_nu0_s=1.0e12,
        recovery_H0_eV=0.05,
    )


def test_near_prefactor_homogeneous_segment_completes():
    state = EmergentGNDState(fast_candidate(), CommonPhysics(n_bins=1))
    totals = state.advance_time(8.4, 20.0, 1200.0)
    assert np.isclose(state.time_s, 8.4)
    assert np.all(np.isfinite(state.mobile_m2))
    assert np.all(np.isfinite(state.retained_m2))
    assert np.all(state.mobile_m2 >= 0.0)
    assert np.all(state.retained_m2 >= 0.0)
    assert totals["emitted_per_m"] >= 0.0
