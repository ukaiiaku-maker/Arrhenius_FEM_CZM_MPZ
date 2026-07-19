from arrhenius_fracture.field_snapshots_v9121 import authoritative_emitted_ledger


def test_authoritative_emitted_ledger_precedes_retained_compatibility_alias():
    snap = {
        "N_em": 3.0,
        "mpz_emitted_total": 17.0,
    }
    assert authoritative_emitted_ledger(snap) == 17.0


def test_authoritative_emitted_ledger_sums_front_states():
    snap = {
        "N_em": 2.0,
        "mpz_front_states": [
            {"state": {"emitted_total": 4.0}},
            {"state": {"state": {"signed_line_content_emitted_total": 6.0}}},
        ],
    }
    assert authoritative_emitted_ledger(snap) == 10.0


def test_legacy_fallback_is_only_used_when_authoritative_fields_are_absent():
    assert authoritative_emitted_ledger({"N_em": 2.5}) == 2.5
