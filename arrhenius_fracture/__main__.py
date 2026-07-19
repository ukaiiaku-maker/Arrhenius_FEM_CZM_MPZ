"""Package entry point.

Production calculations use explicit versioned runners so a command cannot
silently select an obsolete solver path.
"""
raise SystemExit(
    "Use an explicit versioned runner, for example "
    "run_mpz_v9_12_1_tip_only_material_rcurve_fullfield.py or "
    "run_v10_0_5_7_kj_audit_bracket.py."
)
