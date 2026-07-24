from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from arrhenius_fracture import mode_i_first_passage_v10_0_5_13_barrier_only as entry
from arrhenius_fracture.persistent_site_diagnostics_v100516 import (
    installed_persistent_diagnostics_v100516,
)
from arrhenius_fracture.persistent_site_registry_v100514 import (
    select_persistent_site_row,
)
from arrhenius_fracture.persistent_site_signed_mpz_v100514 import (
    PersistentSiteSignedMPZStateV100514,
    SignedShieldingKernelV100514,
)
from arrhenius_fracture.refinement_audit_v100516 import (
    installed_refinement_audit_v100516,
)


def make_state():
    candidate = select_persistent_site_row("v912_peak_0118_persistent_sites")
    n = candidate.n_bins_recommended
    kernel = SignedShieldingKernelV100514(
        active_kernel_Pa_sqrt_m_per_signed_line=np.zeros((2, n)),
        wake_kernel_Pa_sqrt_m_per_signed_line=np.zeros((2, n)),
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
        blunting_length_m=0.5e-6,
    )


def test_discarded_slip_is_not_double_counted_as_line_content():
    state = make_state()
    state.mobile_positive[0, 0] = 2.0
    state.wake_retained_negative[1, 0] = 3.0
    state.accumulated_slip_positive[0, 0] = 6.0
    state.emitted_total = 10.0
    state.escaped_total = 1.0
    state.recovered_total = 0.0
    # Four discarded line units plus four discarded slip-history units.
    state.wake_discarded_total = 8.0
    with installed_persistent_diagnostics_v100516():
        diag = state.diagnostics()
    assert diag["mpz_wake_discarded_line_total"] == pytest.approx(4.0)
    assert diag["mpz_wake_discarded_slip_total"] == pytest.approx(4.0)
    assert diag["mpz_accounted_line_content_total"] == pytest.approx(10.0)
    assert diag["mpz_line_content_balance_relative_error"] == pytest.approx(0.0)
    assert diag["mpz_mobile_count"] == pytest.approx(2.0)
    assert diag["mpz_retained_count"] == pytest.approx(0.0)
    assert diag["mpz_local_slip_count"] > 0.0


def test_refinement_metadata_survives_missing_runtime_mesh(monkeypatch):
    radius = 330.0e-6
    mesh = SimpleNamespace(
        production_refinement_radius_m=radius,
        production_refinement_centers_m=[[0.5e-3, 0.0]],
        production_refinement_policy="synthetic_test",
    )
    monkeypatch.setattr(entry, "make_physical_refinement_mesh_v100510", lambda *a, **k: mesh)
    monkeypatch.setattr(entry, "_annotate_mesh", lambda value, *a, **k: value)
    with installed_refinement_audit_v100516():
        built = entry.make_physical_refinement_mesh_v100510(None)
        assert built is mesh
        payload = entry._mesh_payload(None, radius)
    assert payload["actual_radius_verified"] is True
    assert payload["runtime_mesh_pointer_missing"] is True
    assert payload["captured_mesh_used"] is True
