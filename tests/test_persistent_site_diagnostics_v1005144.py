from __future__ import annotations

import numpy as np
import pytest

from arrhenius_fracture.persistent_site_diagnostics_v1005144 import (
    DIAGNOSTIC_MODEL,
    installed_persistent_diagnostics_v1005144,
)
from arrhenius_fracture.persistent_site_registry_v100514 import (
    select_persistent_site_row,
)
from arrhenius_fracture.persistent_site_signed_mpz_v100514 import (
    PersistentSiteSignedMPZStateV100514,
    SignedShieldingKernelV100514,
)


def state():
    candidate = select_persistent_site_row("v912_peak_0118_persistent_sites")
    active = np.zeros((2, candidate.n_bins_recommended))
    kernel = SignedShieldingKernelV100514(
        active_kernel_Pa_sqrt_m_per_signed_line=active,
        wake_kernel_Pa_sqrt_m_per_signed_line=np.zeros_like(active),
        activation_to_line_content_by_system=np.ones(2),
        metadata={
            "candidate_independent": True,
            "counts_are_signed_burgers_lines": True,
            "normalization_is_mechanically_derived": True,
        },
        source_path="synthetic_diagnostic_kernel",
    )
    return PersistentSiteSignedMPZStateV100514(
        candidate,
        kernel,
        G_Pa=160.0e9,
        nu=0.28,
        b_m=2.74e-10,
        r0_m=1.0e-6,
    )


def test_authoritative_accounting_separates_retained_from_cumulative_emission():
    model = state()
    model.mobile_positive[0, 0] = 1.0
    model.retained_negative[1, 1] = 2.0
    model.wake_mobile_positive[0, 2] = 3.0
    model.wake_retained_negative[1, 3] = 4.0
    model.emitted_total = 15.0
    model.escaped_total = 5.0
    model.recovered_total = 0.0
    model.wake_discarded_total = 0.0
    with installed_persistent_diagnostics_v1005144():
        diagnostics = model.diagnostics()
    assert diagnostics["persistent_state_diagnostic_model"] == DIAGNOSTIC_MODEL
    assert diagnostics["mpz_active_mobile_total"] == pytest.approx(1.0)
    assert diagnostics["mpz_active_retained_total"] == pytest.approx(2.0)
    assert diagnostics["mpz_wake_mobile_total"] == pytest.approx(3.0)
    assert diagnostics["mpz_wake_retained_total"] == pytest.approx(4.0)
    assert diagnostics["mpz_escaped_total"] == pytest.approx(5.0)
    assert diagnostics["mpz_emitted_total"] == pytest.approx(15.0)
    assert diagnostics["mpz_line_content_balance_error"] == pytest.approx(0.0)
    assert diagnostics["N_em_is_cumulative_emission"] is False
    assert diagnostics["legacy_N_em_semantics"] == (
        "instantaneous_active_retained_line_content"
    )
