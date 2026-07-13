# Validation record

- Python compilation: passed for all Python files.
- Bash syntax: passed for both shell scripts.
- Unit/regression tests: 12 passed.
- Regression coverage includes direct loading coefficients, ratios beyond the
  old 89.9-degree cap, exact-backend calibration source checks, bracketed and
  empirical event-state control, class barrier audits, and censor-aware status.
- Plotting smoke test: passed on synthetic v6 campaign data.
- Archive layout and SHA-256 manifest: verified.

The full adaptive-CZM backend is present only in the user's active project, so the
actual production-backend calibration and fracture preflight must be executed
there. The verifier deliberately tests interfaces and control logic without
claiming that it substitutes for that physical preflight.
