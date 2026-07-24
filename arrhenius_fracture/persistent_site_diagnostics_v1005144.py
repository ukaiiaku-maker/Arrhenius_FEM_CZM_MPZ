"""Authoritative persistent-site accounting diagnostics for v10.0.5.14.4."""
from __future__ import annotations

from contextlib import contextmanager
from typing import Iterator

import numpy as np

from .persistent_site_signed_io_v100514 import PersistentSiteSignedIOMixin

DIAGNOSTIC_MODEL = "persistent_signed_line_content_accounting_v10_0_5_14_4"


def _diagnostics_v1005144(self, *args, **kwargs):
    base = _diagnostics_v1005144._original(self, *args, **kwargs)
    active_mobile = float(
        np.sum(self.mobile_positive) + np.sum(self.mobile_negative)
    )
    active_retained = float(
        np.sum(self.retained_positive) + np.sum(self.retained_negative)
    )
    wake_mobile = float(
        np.sum(self.wake_mobile_positive) + np.sum(self.wake_mobile_negative)
    )
    wake_retained = float(
        np.sum(self.wake_retained_positive) + np.sum(self.wake_retained_negative)
    )
    wake_slip = float(
        np.sum(self.wake_slip_positive) + np.sum(self.wake_slip_negative)
    )
    active_slip = float(
        np.sum(self.accumulated_slip_positive)
        + np.sum(self.accumulated_slip_negative)
    )
    emitted = float(self.emitted_total)
    escaped = float(self.escaped_total)
    recovered = float(self.recovered_total)
    discarded = float(self.wake_discarded_total)
    accounted_line_content = (
        active_mobile
        + active_retained
        + wake_mobile
        + wake_retained
        + escaped
        + recovered
        + discarded
    )
    signed_balance = emitted - accounted_line_content
    scale = max(abs(emitted), abs(accounted_line_content), 1.0e-300)
    base.update(
        {
            "persistent_state_diagnostic_model": DIAGNOSTIC_MODEL,
            "legacy_N_em_semantics": "instantaneous_active_retained_line_content",
            "mpz_active_mobile_total": active_mobile,
            "mpz_active_retained_total": active_retained,
            "mpz_active_line_content_total": active_mobile + active_retained,
            "mpz_wake_mobile_total": wake_mobile,
            "mpz_wake_retained_total": wake_retained,
            "mpz_wake_line_content_total": wake_mobile + wake_retained,
            "mpz_active_accumulated_slip_total": active_slip,
            "mpz_wake_accumulated_slip_total": wake_slip,
            "mpz_emitted_total": emitted,
            "mpz_escaped_total": escaped,
            "mpz_recovered_total": recovered,
            "mpz_wake_discarded_total": discarded,
            "mpz_accounted_line_content_total": accounted_line_content,
            "mpz_line_content_balance_signed": signed_balance,
            "mpz_line_content_balance_error": abs(signed_balance),
            "mpz_line_content_balance_relative_error": abs(signed_balance) / scale,
            "N_em_is_cumulative_emission": False,
            "authoritative_cumulative_emission_field": "mpz_emitted_total",
        }
    )
    return base


@contextmanager
def installed_persistent_diagnostics_v1005144() -> Iterator[None]:
    old = PersistentSiteSignedIOMixin.diagnostics
    _diagnostics_v1005144._original = old
    PersistentSiteSignedIOMixin.diagnostics = _diagnostics_v1005144
    try:
        yield
    finally:
        PersistentSiteSignedIOMixin.diagnostics = old


__all__ = [
    "DIAGNOSTIC_MODEL",
    "installed_persistent_diagnostics_v1005144",
]
