from __future__ import annotations

from audit_v10_0_5_2_mpz_binding import replay_binding


def test_outer_v911_parser_factory_binding_uses_requested_200_bins():
    payload = replay_binding("DBTT", 100.0, 200)
    assert payload["binding_replay_verified"] is True
    assert payload["outer_parser_requested_mpz_n_bins"] == 200
    assert payload["active_mpz_n_bins"] == 200
    assert payload["active_mpz_length_m"] == 100.0e-6
    assert payload["active_mpz_dx_m"] == 0.5e-6
    assert payload["active_mpz_source_bin_count"] == 4
