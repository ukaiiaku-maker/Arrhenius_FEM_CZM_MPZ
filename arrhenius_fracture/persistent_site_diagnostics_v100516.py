"""Diagnostics-only accounting corrections for v10.0.5.16.

The legacy moving-frame counter ``wake_discarded_total`` sums losses from mobile
line content, retained line content, and accumulated slip history.  Accumulated
slip is a history variable and must not be counted a second time in the physical
line-content inventory.  This overlay separates the two without changing any
state-evolution equation or moving-frame transaction.
"""
from __future__ import annotations

from contextlib import contextmanager
from typing import Iterator

import numpy as np

from .persistent_site_signed_io_v100514 import PersistentSiteSignedIOMixin

DIAGNOSTIC_MODEL = "persistent_signed_line_and_slip_accounting_v10_0_5_16"


def _diagnostics_v100516(self, *args, **kwargs):
    base = _diagnostics_v100516._original(self, *args, **kwargs)

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
    active_slip = float(
        np.sum(self.accumulated_slip_positive)
        + np.sum(self.accumulated_slip_negative)
    )
    wake_slip = float(
        np.sum(self.wake_slip_positive) + np.sum(self.wake_slip_negative)
    )
    local_slip = float(self.local_slip_count())

    emitted = float(self.emitted_total)
    escaped = float(self.escaped_total)
    recovered = float(self.recovered_total)
    legacy_discarded = float(self.wake_discarded_total)

    present_or_escaped_line = (
        active_mobile
        + active_retained
        + wake_mobile
        + wake_retained
        + escaped
        + recovered
    )
    discarded_line = max(emitted - present_or_escaped_line, 0.0)
    discarded_slip = max(legacy_discarded - discarded_line, 0.0)
    accounted_line = present_or_escaped_line + discarded_line
    signed_balance = emitted - accounted_line
    scale = max(abs(emitted), abs(accounted_line), 1.0e-300)

    base.update(
        {
            "persistent_state_diagnostic_model": DIAGNOSTIC_MODEL,
            "legacy_N_em_semantics": "instantaneous_active_retained_line_content",
            "legacy_wake_discarded_total_semantics": (
                "discarded_mobile_plus_retained_line_content_plus_discarded_slip_history"
            ),
            "mpz_mobile_count": active_mobile,
            "mpz_retained_count": active_retained,
            "mpz_available_site_fraction": 1.0,
            "mpz_local_slip_count": local_slip,
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
            "mpz_wake_discarded_legacy_total": legacy_discarded,
            "mpz_wake_discarded_mobile_total": None,
            "mpz_wake_discarded_retained_total": None,
            "mpz_wake_discarded_line_total": discarded_line,
            "mpz_wake_discarded_slip_total": discarded_slip,
            "mpz_wake_discarded_line_total_is_balance_derived": True,
            "mpz_accounted_line_content_total": accounted_line,
            "mpz_line_content_balance_signed": signed_balance,
            "mpz_line_content_balance_error": abs(signed_balance),
            "mpz_line_content_balance_relative_error": abs(signed_balance) / scale,
            "N_em_is_cumulative_emission": False,
            "authoritative_cumulative_emission_field": "mpz_emitted_total",
            "authoritative_discarded_line_field": "mpz_wake_discarded_line_total",
        }
    )
    return base


@contextmanager
def installed_persistent_diagnostics_v100516() -> Iterator[None]:
    old = PersistentSiteSignedIOMixin.diagnostics
    _diagnostics_v100516._original = old
    PersistentSiteSignedIOMixin.diagnostics = _diagnostics_v100516
    try:
        yield
    finally:
        PersistentSiteSignedIOMixin.diagnostics = old


__all__ = [
    "DIAGNOSTIC_MODEL",
    "installed_persistent_diagnostics_v100516",
]
