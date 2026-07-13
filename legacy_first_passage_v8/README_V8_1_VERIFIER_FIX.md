# Mixed-mode FEM/CZM v8.1 verifier-only correction

This patch changes no mechanics, calibration, campaign, or plotting code.

It corrects two verifier defects in the v8 package:

1. The verifier expected the old model ID without the `_boolean_safe` suffix.
2. The verifier omitted `tests.test_mixed_mode_probe_boolean_v8`, the regression tests for the bug v8 was intended to fix.

Copy the two verifier files into the top-level `Arrhenius_FEM_CZM` directory and run the shell verifier.
